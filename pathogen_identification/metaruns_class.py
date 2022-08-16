import copy
import itertools as it
import logging
import os
import pickle
import shutil
import subprocess
import sys
import time
from re import T
from threading import Thread
from typing import Type

import numpy as np
import pandas as pd
from metagen_view.settings import STATICFILES_DIRS

from pathogen_detection.constants_settings import ConstantsSettings
from pathogen_detection.object_classes import Sample_runClass
from pathogen_detection.params import (
    ACTIONS,
    ARGS_ASS,
    ARGS_CLASS,
    ARGS_ENRICH,
    ARGS_QC,
    ARGS_REMAP,
    BINARIES,
    CONSTANTS,
    DATA_TYPE,
    DIRS,
    METADATA,
    SOFTWARE,
    SOURCE,
)
from pathogen_detection.run_main import RunMain_class
from pathogen_detection.update_DBs import (
    RunIndex_Update_Retrieve_Key,
    Update_project,
    Update_QC_report,
    Update_Sample,
    Update_Sample_Runs,
    retrieve_number_of_runs,
)


class metaclass_run:
    """
    class to setup metagenomics classification run.
    """

    qcrun: Type[RunMain_class]

    def __init__(
        self,
        project_name: str,
        sample_name: str,
        id="1",
        rdir="",
        child="",
        static_dir: str = "",
    ):

        self.id = id
        self.actions = {
            "CLEANING": ACTIONS["CLEANING"],
            "SIFT": ACTIONS["PHAGE_DEPL"],
            "VIRSORT": ACTIONS["VIRSORT"],
            "QCONTROL": False,
            "ASSEMBLY": ACTIONS["ASSEMBLE"],
            "DEPLETE": ACTIONS["DEPLETE"],
            "ENRICH": ACTIONS["ENRICH"],
            "CLASSIFY": False,
            "REMAPPING": False,
        }
        self.sample_name = sample_name
        self.project_name = project_name

        self.r1 = ""
        self.r2 = ""
        self.paired = False
        self.reference = ""
        self.begin_time = time.perf_counter()
        self.exec_time = 0
        self.static_dir = static_dir

    def update_runtime(self):
        """
        update runtime in DB
        :return:
        """
        self.exec_time = time.perf_counter() - self.begin_time

    def parse_config(self, rdir, child=""):
        """
        function to return class instance by parsing config in given directory.
        :param rdir: previously spawned directory with config file.
        :param child: ID of child config file present in rdir.
        :return: self
        """
        if len(child):
            confgf = rdir + "confch_{}.sh".format(child)
            self.configf = confgf
            self.main = rdir + "main_{}.sh".format(child)
        else:
            confgf = rdir + "config.sh"
            self.configf = confgf
            self.main = rdir + "main.sh"
        settings = {}

        with open(confgf, "r") as fp:
            ln = fp.readline().strip()
            while ln:
                if ln[0] == "#":
                    ln = fp.readline().strip()
                    continue
                ln = ln.split("=")
                settings[ln[0]] = ln[1].strip('"')

                ln = fp.readline().strip()

        self.actions = {}
        self.config = {}
        self.params = {}

        for dr in ACTIONS.keys():
            if dr in settings.keys():
                self.actions[dr] = settings[dr] == "true"

        for dr in SOFTWARE.keys():
            if dr in settings.keys():
                self.config[dr] = [settings[dr]]

        pooled_params = {**ARGS_ENRICH, **ARGS_ASS, **ARGS_CLASS, **ARGS_REMAP}

        for dr in pooled_params.keys():
            if dr in settings.keys():
                self.params = [settings[dr]]

        self.dir = settings["RDIR"]
        self.id = settings["SUFFIX"]
        self.r1 = settings["INPUT"]
        self.config = pd.DataFrame(self.config)
        self.params = pd.DataFrame(self.params)

        if "PAIR" in settings.keys():
            self.r2 = settings["PAIR"]
            self.paired = True

    def config_get(self, config_frame):
        """
        function to incorporate SOFTWARE data frame into self.config.
        :param config_frame: data FRAME. columns= modules, single row, elemets = software.
        :return: self with updated config.
        """

        self.config = config_frame

        if "DEPLETION" in self.config.columns:
            if len(self.config["DEPLETION"][0]) == 0:
                self.actions["DEPLETE"] = False

        if "ENRICHMENT" in self.config.columns:

            if len(self.config["ENRICHMENT"][0]) == 0:
                self.actions["ENRICH"] = False

    def params_get(self, param_frame):
        """
        function to incorporate parameter data frame.
        :param param_frame:
        :return: None
        """
        self.params = param_frame

    def write_config(self, config="config.sh"):
        """
        write config file based on actions, dirs, config and params dfs.
        :param config:
        :return:
        """
        #
        self.configf = config
        lines = ["#!/bin/bash"]
        #
        lines.append("## INPUT")
        lines.append('{}="{}"'.format("INPUT", self.r1))
        if self.r2:
            lines.append('{}="{}"'.format("PAIR", self.r2))
        else:
            lines.append("{}=``".format("PAIR"))

        if self.reference:
            lines.append('{}="{}"'.format("REFERENCE", self.reference))

        lines.append("##")
        #
        for dr, path in SOURCE.items():
            lines.append('{}="{}"'.format(dr, path))
        lines.append("## DIRECTORIES")
        #
        for dr, path in DIRS.items():

            lines.append('{}="{}"'.format(dr, self.dir + path))
        lines.append("## ACTIONS")
        #
        for dr, g in self.actions.items():
            lines.append('{}="{}"'.format(dr, str(g).lower()))
        lines.append("## SOFTWARE")
        #
        for dr in self.config.columns:
            lines.append('{}="{}"'.format(dr, self.config[dr][0]))
        lines.append("## PARAMS")
        #
        for dr in range(self.params.shape[0]):
            lines.append('{}="{}"'.format(self.params.param[dr], self.params.value[dr]))
        lines.append("##")
        #
        lines.append('SUFFIX="{}"'.format(self.id))
        lines.append('RDIR="{}"'.format(self.dir))

        ###
        with open(config, "w") as fn:
            fn.write("\n".join(lines))

    def prep_config_dict(self):
        """
        function to prepare config dictionary for writing to config file.
        :return: config dictionary.
        """
        config = {
            "project": os.getcwd(),
            "source": SOURCE,
            "directories": {
                "root": self.dir,
            },
            "static_dir": self.static_dir,
            "actions": {},
            "bin": {},
            "threads": 6,
            "prefix": self.id,
            "type": ["SE", "PE"][int(os.path.isfile(self.r2))],
            "sample_name": self.sample_name,
            "project_name": self.project_name,
            "metadata": {
                x: os.path.join(METADATA["ROOT"], g) for x, g in METADATA.items()
            },
            "r1": self.r1,
            "r2": self.r2,
            "technology": DATA_TYPE,
            "bin": BINARIES,
        }

        for dr, g in self.actions.items():
            config["actions"][dr] = g

        for dr, g in DIRS.items():
            config["directories"][dr] = self.dir + g

        config.update(CONSTANTS)

        self.config_dict = config

    def prep_env(self, rdir=""):
        """
        from main directory bearing scripts, params.py and main.sh, create metagenome run directory

        :return:
        """
        #
        if not rdir:
            rdir = os.getcwd()

        self.dir = rdir + "{}/".format(self.id)

        for dir in DIRS.values():
            os.system("mkdir -p " + self.dir + dir)

    def spawn(self, fofn="", config="config.sh", source="main.sh", sink="main.sh"):
        """
        write config file, copy and modify instance main.sh file.
        :param fofn: file of fastas, one per line, max two.
        :param config: name of config file
        :param source: name of main*.sh file to copy and modify.
        :param sink: destinaion main*.sh file. to become instance main.
        :return:
        """
        #
        if len(fofn):
            self.read_fofn(fofn)

        self.write_config(config=config)
        self.configf = config
        #
        self.main = sink

    def read_fofn(self, fofn):
        """
        read file of files.
        :param fofn: file path
        :return:
        """
        d = 0
        with open(fofn, "r") as fn:
            ln = fn.readline().strip()
            while ln:
                if d == 0:
                    self.r1 = ln
                if d == 1:
                    self.r2 = ln
                    self.paired = True
                d += 1
                ln = fn.readline().strip()

    def run_main(self):
        """
        deploy meta classsification run
        :return:
        """
        #
        self.prep_config_dict()
        self.RunMain = RunMain_class(self.config_dict, self.params)

        self.RunMain.Run()
        self.RunMain.Summarize()

    def continue_main_run(self):
        """
        continue main run
        :return:
        """
        self.prep_config_dict()
        self.RunMain.Update(self.config_dict, self.params)

        self.RunMain.Run()
        self.RunMain.Summarize()

    def report_run_status_save(self):
        """
        Generate run summary dataclasses, update database.
        :return:
        """
        self.RunMain.generate_output_data_classes()

        Update_Sample_Runs(self.RunMain)

    def clean(self, delete=True, outf=""):

        trash = ["remap/", "reads/", "host_depletion/", "classification/"]

        if outf:
            os.makedirs(outf + self.id, exist_ok=True)

        subprocess.Popen(
            "mv {}*tsv {}".format(self.dir, self.dir + DIRS["OUTD"]), shell=True
        )
        subprocess.Popen(
            "mv {}*sh {}".format(self.dir, self.dir + DIRS["OUTD"]), shell=True
        )

        if delete:
            for dr in trash:
                subprocess.run(["rm", "-r", self.dir + dr])


class meta_orchestra:
    def __init__(self, fofn, sup=1, down=1, smax=4, odir="", estimate_only=False):
        self.sup = sup
        self.down = down
        self.processes = {}
        self.smax = smax
        self.odir = odir
        self.fofn = fofn
        self.reference = ""

        ###
        self.modules_to_stores = {
            "PREPROCESS": ARGS_QC,
            "ENRICHMENT": ARGS_ENRICH,
            "ASSEMBLY": ARGS_ASS,
            "CONTIG_CLASSIFICATION": ARGS_CLASS,
            "READ_CLASSIFICATION": ARGS_CLASS,
            "REMAPPING": ARGS_REMAP,
        }

        if estimate_only:
            self.prep_numbers()
            sys.exit(0)

        ### create directories
        self.sample_name = os.path.basename(fofn)
        self.rdir = self.odir + self.sample_name + "/"
        project_directory_path = os.path.dirname(self.odir)
        self.project_name = os.path.basename(project_directory_path)

        os.makedirs(self.rdir, exist_ok=True)

        with open(self.rdir + "parameter_combinations.txt", "w") as f:
            f.write(f"pre-assembly:\t{self.sup}\npost-assembly:\t{self.down}")
        shutil.copy(fofn, os.path.join(self.rdir + "input.fofn"))

        with open(self.rdir + "technology.txt", "w") as f:
            f.write("technology\t{}".format(DATA_TYPE))

        self.outd = self.rdir + "output/"
        os.makedirs(self.outd, exist_ok=True)

        self.staticdir = os.path.join(
            self.project_name,
            self.sample_name,
        )
        os.makedirs(
            os.path.join(ConstantsSettings.static_directory, self.staticdir),
            exist_ok=True,
        )

        ### create meta_classifier instances
        self.projects = {}
        print("updating project")
        Update_project(self.odir)
        self.total_runs = retrieve_number_of_runs(self.project_name, self.sample_name)

    def clean(self, delete: bool = True):
        for sid, sac in self.projects.items():
            sac.clean(delete=delete, outf=self.outd)

    def sample_main(
        self,
        sample=1,
        cols=["PREPROCESS", "ENRICHMENT", "ASSEMBLY", "CONTIG_CLASSIFICATION"],
    ):
        """
        sample module / software combinations from dictionaries in params.py. random.
        :param sample: how many combinations to sample. corresponds to number of metclass run directories to be created.
        :param cols: keys to sample combinations of in SOFTWARE dict in params.py
        :return: data frame.
        """
        if len(cols) == 0:
            cols = list(SOFTWARE.keys())

        venue = [SOFTWARE[x] for x in cols]
        venues = list(it.product(*venue))

        if sample > 0 and sample < len(venues):

            vex = np.random.choice(list(range(len(venues))), sample, replace=False)
            venues = [venues[x] for x in vex]

        venues = pd.DataFrame(venues)
        venues.columns = cols
        #
        return venues

    def params_extract(self, show, pooled_params, modules=[], sample=1):
        """
        takes list of software, which might or not have entries in the argument dictionaries.
        """

        if len(modules) == 0:
            modules = list(SOFTWARE.keys())

        relate = []
        new_features = []
        nvens = []
        #
        for ix, soft in enumerate(show):
            soft_module = modules[ix]
            if soft in self.modules_to_stores[soft_module].keys():
                for c, g in self.modules_to_stores[soft_module][soft].items():
                    relate.append([soft_module, soft])
                    new_features.append(c)
                    nvens.append(g)
        #
        nvens = list(it.product(*nvens))
        if sample:
            vex = np.random.choice(list(range(len(nvens))), sample, replace=False)
            nvens = [nvens[x] for x in vex]

        relate = pd.DataFrame(relate, columns=["module", "software"])
        nvens = pd.DataFrame(nvens).reset_index(drop=True)
        nvens.columns = new_features

        return nvens, relate

    def record_runs(self, outf="summary.tsv"):
        """
        record instance software and parameter informationo and print to common csv.
        :param outf: common csv file path.
        :return: None
        """
        for sd, tb in self.processes.items():
            # params = tb["params"]
            # params.to_csv(sd + ".args.tsv", sep="\t", header=True, index=False)
            with open(sd + ".runtime", "w") as f:
                f.write(str(tb["runtime"]))

    @staticmethod
    def extract_parameters(
        parameters_df_list: list, linked_db_list: list, common_index: int
    ):
        """
        extract parameters from list of data frames.
        :param parameters_df_list: list of data frames.
        :param linked_db_list: list of data frames.
        :param common_index: index of common index in data frames.

        :return: data frame.
        """

        params = (
            parameters_df_list[common_index[0]]
            .loc[[common_index[1]]]
            .reset_index(drop=True)
        )
        params = pd.DataFrame([params.columns, params.loc[0]]).T
        params = pd.concat(
            (linked_db_list[common_index[0]], params), axis=1
        ).reset_index(drop=True)
        params.columns = ["module", "software", "param", "value"]

        return params

    def generate_combinations(self, ncomb: int = 0, modules: list = []):

        hdconf = self.sample_main(sample=ncomb, cols=modules)
        params2 = {}
        paramCombs = [
            self.params_extract(hdconf.iloc[idx], params2, modules=modules, sample=0)
            for idx in range(hdconf.shape[0])
        ]
        linked_dbs = [x[1] for x in paramCombs]
        paramCombs = [x[0] for x in paramCombs]

        return hdconf, linked_dbs, paramCombs

    def combination_possibilities(self, hdconf, paramCombs, requested_combinations=0):

        suprelay = [
            [(x, y) for y in range(paramCombs[x].shape[0])]
            for x in range(hdconf.shape[0])
        ]

        suprelay = list(it.chain(*suprelay))

        if requested_combinations > len(suprelay):
            logging.info("down sampling runs to {}".format(len(suprelay)))

        if requested_combinations > 0 and requested_combinations < len(suprelay):
            vex = np.random.choice(list(range(len(suprelay))), self.down, replace=False)
            suprelay = [suprelay[x] for x in vex]
        #

        return suprelay

    def prep_numbers(self):

        hdconf, linked_dbs, paramCombs = self.generate_combinations()
        number_of_combinations = self.combination_possibilities(hdconf, paramCombs)
        number_of_combinations = len(number_of_combinations)
        print("total number of combinations: {}".format(number_of_combinations))

    def run_proc(self, srun, sac, hdconf, paramCombs, linked_db, run_prefix=""):
        """
        prepare and deploy instance of metaclassification run with only post assembly streps.

        :param sac: metaclassification run class instance.
        :param hdconf: data frame of software combinations. 1 row.
        :param srun: row index in parameter df.
        :return: None
        """

        copy_child = lambda obj: pickle.loads(pickle.dumps(obj))

        child = copy_child(sac)

        if not run_prefix:
            run_prefix = RunIndex_Update_Retrieve_Key(
                self.project_name, self.sample_name
            )

        child.id = run_prefix

        child.actions["QCONTROL"] = False
        child.actions["DEPLETE"] = False
        child.actions["ENRICH"] = False
        child.actions["ASSEMBLY"] = False
        child.actions["SIFT"] = ACTIONS["PHAGE_DEPL"]
        child.actions["VIRSORT"] = ACTIONS["VIRSORT"]
        child.actions["CLASSIFY"] = ACTIONS["CLASSIFY"]
        child.actions["REMAPPING"] = ACTIONS["REMAP"]

        if len(paramCombs[srun[0]]) == 0:
            return self

        child.config_get(hdconf.loc[[srun[0]]].reset_index(drop=True))

        params = self.extract_parameters(paramCombs, linked_db, srun)

        child.params_get(params)
        #
        child.spawn(
            fofn=child.dir + DIRS["log_dir"] + "{}_latest.fofn".format(sac.id),
            config=child.dir + "confch_{}.sh".format(child.id),
            source=child.main,
            sink=child.dir + "main_{}.sh".format(child.id),
        )

        child_params = pd.concat([sac.params, child.params], axis=0)

        child_params.to_csv(
            child.dir + child.id + ".args.tsv", sep="\t", header=True, index=False
        )

        child.continue_main_run()
        print(f"child run {child.id} finished.")
        child.update_runtime()
        child.report_run_status_save()

        with open(child.dir + child.id + ".runtime", "w") as f:
            f.write(str(child.exec_time))

        #
        ### gather report / output, mix parent and child params.
        self.processes[child.dir + child.id] = {
            "params": child_params,
            "runtime": child.exec_time,
        }

        child.RunMain = {}

    def get_run_index_name(self):
        new_name = f"run_{self.total_runs}"
        self.total_runs += 1
        return new_name

    def low_run(self, prj):
        """
        deploy children of current metaclass run instance, controls parallel thread deployment tto not exceed self.smax.
        :param prj: metaclass_run instance.
        :return: None
        """
        modules = ["CONTIG_CLASSIFICATION", "READ_CLASSIFICATION", "REMAPPING"]

        hdconf, linked_dbs, paramCombs = self.generate_combinations(self.down, modules)

        suprelay = self.combination_possibilities(hdconf, paramCombs, self.down)
        self.down = len(suprelay)
        print("low run possibilities: {}".format(self.down))

        rlow = self.down
        #
        if rlow > self.smax:
            pack = list(np.arange(0, rlow, self.smax)) + [rlow]
            lset = [list(range(pack[x], pack[x + 1])) for x in range(len(pack) - 1)]
            #
            for ix, a_args in enumerate(lset):

                threads = [
                    Thread(
                        target=self.run_proc,
                        args=(
                            suprelay[x],
                            prj,
                            hdconf.copy(),
                            paramCombs.copy(),
                            linked_dbs.copy(),
                            self.get_run_index_name(),
                        ),
                    )
                    for x in a_args
                ]
                for th in threads:
                    th.start()
                for th in threads:
                    th.join()

        else:
            a_args = list(range(rlow))
            threads = [
                Thread(
                    target=self.run_proc,
                    args=(
                        suprelay[x],
                        prj,
                        hdconf.copy(),
                        paramCombs.copy(),
                        linked_dbs.copy(),
                        self.get_run_index_name(),
                    ),
                )
                for x in a_args
            ]
            for th in threads:
                th.start()
            for th in threads:
                th.join()

    def low_deploy(self, sac=""):
        """
        deploy post assembly steps across metaclass runs.
        :param sac: metaclass_run instance.
        :return: self
        """
        if sac:
            sel = [x for x, g in self.projects.items() if g.dir == sac][0]
            self.low_run(sac)

        else:
            for idr, prj in self.projects.items():
                self.low_run(prj)

        return self

    def sup_deploy(self, fofn):
        """
        deploy metaclass runs.
        :return: self
        """
        modules = ["PREPROCESS", "ENRICHMENT", "ASSEMBLY"]

        conf, linked_dbs, paramCombs = self.generate_combinations(self.sup, modules)

        suprelay = self.combination_possibilities(conf, paramCombs, self.sup)

        self.sup = len(suprelay)
        rsup = self.sup
        print("sup run possibilities: {}".format(self.sup))

        def find_run_name(run_number):
            run_idx = run_number
            new_name = "suprun_{}".format(run_idx)

            while os.path.isdir(self.rdir + new_name):
                run_idx += 1
                new_name = "suprun_{}".format(run_idx)
            os.makedirs(self.rdir + new_name, exist_ok=True)
            return new_name

        if rsup > self.smax:
            pack = list(np.arange(0, rsup, self.smax)) + [rsup]
            lset = [list(range(pack[x], pack[x + 1])) for x in range(len(pack) - 1)]
            run_names = [find_run_name(x) for x in range(len(suprelay))]
            for a_args in lset:
                threads = [
                    Thread(
                        target=self.sup_run,
                        args=(
                            suprelay[x],
                            run_names[x],
                            conf.copy(),
                            paramCombs.copy(),
                            linked_dbs.copy(),
                        ),
                    )
                    for x in a_args
                ]
                for th in threads:
                    th.start()
                for th in threads:
                    th.join()

        else:
            a_args = list(range(rsup))
            run_names = [find_run_name(x) for x in a_args]
            threads = [
                Thread(
                    target=self.sup_run,
                    args=(
                        suprelay[x],
                        run_names[x],
                        conf.copy(),
                        paramCombs.copy(),
                        linked_dbs.copy(),
                    ),
                )
                for x in a_args
            ]
            for th in threads:
                th.start()
            for th in threads:
                th.join()

        return self

    def data_qc(self):
        """
        run quality control on input data, so only have to do it once.
        :return: self
        """
        modules = ["PREPROCESS"]
        paramsqc = {**ARGS_QC}
        conf = self.sample_main(sample=1, cols=modules)

        params = self.params_extract(conf.iloc[0], paramsqc, modules=modules, sample=1)
        softdb = params[1]
        params = params[0]

        qcrun = metaclass_run(
            self.project_name,
            self.sample_name,
            id="sampleQC",
            static_dir=self.staticdir,
        )
        qcrun.actions = {x: False for x in qcrun.actions.keys()}
        qcrun.actions["QCONTROL"] = ACTIONS["QCONTROL"]

        qcrun.config_get(conf)

        params = pd.DataFrame([params.columns, params.loc[0]]).T
        params = pd.concat((softdb, params), axis=1).reset_index(drop=True)
        params.columns = ["module", "software", "param", "value"]

        qcrun.params_get(params)
        qcrun.prep_env(rdir=self.rdir)
        qcrun.spawn(
            fofn=self.fofn, config=qcrun.dir + "config.sh", sink=qcrun.dir + "main.sh"
        )
        qcrun.run_main()
        qcrun.update_runtime()

        self.qcrun = qcrun

        Update_Sample(self.qcrun.RunMain.sample)

        Update_QC_report(self.qcrun.RunMain.sample)

        # qcrun.RunMain.sample = self.qcrun.sample

        dirlist = [
            "remap/",
            "host_depletion/",
            "classification/",
            "assembly/",
            "classification/",
            "enrichment/",
            "output/",
        ]

        for dirl in dirlist:
            if os.path.isdir(qcrun.dir + dirl):
                shutil.rmtree(qcrun.dir + dirl)

    def sup_run(self, idx, run_id, conf, paramcomb, linked_db):
        """
        prepare metaclass_run for deployment. create directories, config and main files, set up software and modules.
        :param idx: row index in software db.
        :param conf: software combination df. one line.
        :return: self
        """

        nrun = metaclass_run(
            self.project_name,
            self.sample_name,
            id=run_id,
            static_dir=self.staticdir,
        )
        nrun.reference = self.reference
        nrun.qcrun = self.qcrun.RunMain

        #
        nrun.config_get(conf.loc[[idx[0]]].reset_index(drop=True))
        #
        params = self.extract_parameters(paramcomb, linked_db, idx)

        nrun.params_get(params)
        #
        nrun.prep_env(rdir=self.rdir)
        nrun.spawn(
            fofn=self.qcrun.dir + DIRS["log_dir"] + "{}_latest.fofn".format("sampleQC"),
            config=nrun.dir + "config.sh",
            sink=nrun.dir + "main.sh",
        )  ### fofn currently flobal variable.

        shutil.copy(
            self.qcrun.dir + DIRS["log_dir"] + "reads_latest.stats",
            nrun.dir + DIRS["log_dir"],
        )

        supstart = time.perf_counter()
        copy_child = lambda obj: pickle.loads(pickle.dumps(obj))

        nrun.RunMain = copy_child(self.qcrun.RunMain)
        nrun.continue_main_run()
        supend = time.perf_counter()
        nrun.exect = supend - supstart

        self.projects[nrun.id] = nrun

        return self
