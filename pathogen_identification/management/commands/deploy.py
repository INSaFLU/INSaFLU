import os
from datetime import date

from django.core.management.base import BaseCommand
from pathogen_detection.constants_settings import ConstantsSettings
from pathogen_detection.metaruns_class import meta_orchestra


class Command(BaseCommand):
    help = "Populates the DBS"

    def add_arguments(self, parser):
        parser.add_argument(
            "--fofn",
            "-f",
            type=str,
            help="file of fastq files. one file path per line.",
        )
        parser.add_argument(
            "--fdir",
            "-d",
            type=str,
            help="dictory containing fofn files.",
        )

        parser.add_argument(
            "--nsup",
            "-n",
            type=int,
            default=1,
            help="number of parameter combinations before assembly, 0 will run all [default=0]",
        )

        parser.add_argument(
            "--nlow",
            "-m",
            type=int,
            default=1,
            help="number of parameter combinations downstrem of assembly, 0 will run all [default=0]",
        )

        parser.add_argument(
            "-o", "--odir", type=str, default="", help="Output directory."
        )

        parser.add_argument(
            "--clean",
            action="store_true",
            default=False,
            help="move output reports to final output directory, intermediate files and config files to run output directories",
        )

        parser.add_argument(
            "--fdel",
            action="store_true",
            default=False,
            help="clean output repositories, keep only report files and assembly file. Recommend for large benchmarking runs.",
        )

        parser.add_argument(
            "--ref",
            type=str,
            required=False,
            default="",
            help="reference genome for host depletion",
        )

        parser.add_argument(
            "--estimate_runs",
            action="store_true",
            default=False,
            help="estimate number of runs based on number of files and number of parameter combinations",
        )

    def handle(self, *args, **options):
        ###
        #
        os.makedirs(ConstantsSettings.project_directory, exist_ok=True)

        if options["estimate_runs"]:
            print("estimate")
            event = meta_orchestra(options["fofn"], estimate_only=True)

        if not options["odir"]:
            options["odir"] = "run_" + str(date.today())

        options["odir"] = os.path.join(
            ConstantsSettings.project_directory, options["odir"]
        )

        if options["odir"][-1] != "/":
            options["odir"] += "/"

        options["odir"] = os.path.join(os.getcwd(), options["odir"])

        if options["fofn"]:

            event = meta_orchestra(
                options["fofn"],
                sup=options["nsup"],
                down=options["nlow"],
                odir=options["odir"],
            )
            event.reference = options["ref"]

            event.data_qc()
            event.sup_deploy(options["fofn"])
            event.low_deploy()
            event.record_runs()

            # event.clean(delete=args.clean)

        elif options["fdir"]:
            if options["fdir"][-1] != "/":
                options["fdir"] += "/"

            flist = [
                options["fdir"] + x
                for x in os.listdir(options["fdir"])
                if os.path.splitext(x)[1] == ".fofn"
            ]

            for fofn in flist:
                if os.path.exists(
                    os.path.join(options["odir"], os.path.basename(fofn))
                ):
                    print("skipping {}".format(fofn))
                    # continue

                event = meta_orchestra(
                    fofn,
                    sup=options["nsup"],
                    down=options["nlow"],
                    odir=options["odir"],
                )

                event.reference = options["ref"]
                event.data_qc()
                event.sup_deploy(fofn)
                event.low_deploy()
                event.record_runs()

                if options["clean"] or options["fdel"]:
                    event.clean(delete=options["fdel"])
