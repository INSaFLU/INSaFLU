class Host:
    remote_host: str
    remote_path: str
    remote_filename: str
    host_name: str
    common_name: str


class HomoSapiens(Host):
    def __init__(self):
        self.remote_host = "hgdownload.soe.ucsc.edu"
        self.remote_path = "/goldenPath/hg38/bigZips/"
        self.remote_filename = "hg38.fa.gz"
        self.host_name = "hg38"
        self.common_name = "human"


class MarmotaMarmota(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "/genomes/refseq/vertebrate_mammalian/Marmota_marmota/latest_assembly_versions/GCF_001458135.2_marMar/"
        self.remote_filename = "GCF_001458135.2_marMar_genomic.fna.gz"
        self.host_name = "marmota_marmota"
        self.common_name = "marmot"


class CulexPipiens(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "/genomes/refseq/invertebrate/Culex_pipiens/latest_assembly_versions/GCF_016801865.2_TS_CPP_V2/"
        self.remote_filename = "GCF_016801865.2_TS_CPP_V2_genomic.fna.gz"
        self.host_name = "culex_pipiens"
        self.common_name = "mosquito"


class AnasPlatyrhynchos(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "/genomes/refseq/vertebrate_other/Anas_platyrhynchos/latest_assembly_versions/GCF_015476345.1_ZJU1.0"
        self.remote_filename = "GCF_015476345.1_ZJU1.0_genomic.fna.gz"
        self.host_name = "anas_platyrhynchos"
        self.common_name = "duck"


class PipistrellusKuhlii(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "/genomes/refseq/invertebrate/Phlebotomus_papatasi/latest_assembly_versions/GCF_024763615.1_Ppap_2.1"
        self.remote_filename = "GCF_024763615.1_Ppap_2.1_genomic.fna.gz"
        self.host_name = "pipistrellus_kuhlii"
        self.common_name = "bat"


class PhlebotomusPapatasi(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "/genomes/refseq/invertebrate/Phlebotomus_papatasi/latest_assembly_versions/GCF_000439695.1_Ppap_1.0/"
        self.remote_filename = "GCF_000439695.1_Ppap_1.0_genomic.fna.gz"
        self.host_name = "phlebotomus_papatasi"
        self.common_name = "sandfly"


class FelisCatus(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "/genomes/refseq/vertebrate_mammalian/Felis_catus/latest_assembly_versions/GCF_018350175.1_F.catus_Fca126_mat1.0/"
        self.remote_filename = "GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna.gz"
        self.host_name = "felis_catus"
        self.common_name = "cat"


class CanisLupusFamiliaris(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "/genomes/refseq/vertebrate_mammalian/Canis_lupus_familiaris/latest_assembly_versions/GCF_011100685.1_UU_Cfam_GSD_1.0/"
        self.remote_filename = "GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna.gz"
        self.host_name = "canis_lupus_familiaris"
        self.common_name = "dog"


class CyprinusCarpio(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "/genomes/refseq/vertebrate_other/Cyprinus_carpio/latest_assembly_versions/GCF_018340385.1_ASM1834038v1/"
        self.remote_filename = "GCF_018340385.1_ASM1834038v1_genomic.fna.gz"
        self.host_name = "cyprinus_carpio"
        self.common_name = "carp"


class SusScrofa(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/vertebrate_mammalian/Sus_scrofa/latest_assembly_versions/GCF_000003025.6_Sscrofa11.1/"
        self.remote_filename = "GCF_000003025.6_Sscrofa11.1_genomic.fna.gz"
        self.host_name = "sus_scrofa"
        self.common_name = "pig"


class AedesAlbopictus(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/invertebrate/Aedes_albopictus/latest_assembly_versions/GCF_006496715.2_Aalbo_primary.1/"
        self.remote_filename = "GCF_006496715.2_Aalbo_primary.1_genomic.fna.gz"
        self.host_name = "aedes_albopictus"
        self.common_name = "mosquito"


class GallusGallus(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/vertebrate_other/Gallus_gallus/latest_assembly_versions/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/"
        self.remote_filename = (
            "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz"
        )
        self.host_name = "gallus_gallus"
        self.common_name = "chicken"


class OncorhynchusMykiss(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/vertebrate_other/Oncorhynchus_mykiss/latest_assembly_versions/GCF_013265735.2_USDA_OmykA_1.1/"
        self.remote_filename = "GCF_013265735.2_USDA_OmykA_1.1_genomic.fna.gz"
        self.host_name = "oncorhynchus_mykiss"
        self.common_name = "rainbow_trout"


class SalmoSalar(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "/genomes/refseq/vertebrate_other/Salmo_salar/latest_assembly_versions/GCF_905237065.1_Ssal_v3.1/"
        self.remote_filename = "GCF_905237065.1_Ssal_v3.1_genomic.fna.gz "
        self.host_name = "salmo_salar"
        self.common_name = "atlantic_salmon"


class BosTaurus(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.3_ARS-UCD2.0/"
        self.remote_filename = "GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz"
        self.host_name = "bos_taurus"
        self.common_name = "cow"


class NeogaleVison(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/vertebrate_mammalian/Neogale_vison/latest_assembly_versions/GCF_020171115.1_ASM_NN_V1/"
        self.remote_filename = "GCF_020171115.1_ASM_NN_V1_genomic.fna.gz"
        self.host_name = "neogale_vison"
        self.common_name = "mink"
