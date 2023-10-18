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
        self.remote_path = "genomes/refseq/vertebrate_other/Acipenser_ruthenus/latest_assembly_versions/GCF_902713425.1_fAciRut3.2_maternal_haplotype/"
        self.remote_filename = (
            "GCF_902713425.1_fAciRut3.2_maternal_haplotype_genomic.fna.gz"
        )
        self.host_name = "salmo_salar"
        self.common_name = "atlantic_salmon"


class BosTaurus(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "/genomes/refseq/vertebrate_mammalian/Bos_taurus/latest_assembly_versions/GCF_002263795.3_ARS-UCD2.0/"
        self.remote_filename = "GCF_002263795.3_ARS-UCD2.0_cds_from_genomic.fna.gz"
        self.host_name = "bos_taurus"
        self.common_name = "cow"


class NeogaleVison(Host):
    def __init__(self):
        self.remote_host = "ftp.ncbi.nlm.nih.gov"
        self.remote_path = "genomes/refseq/vertebrate_mammalian/Neogale_vison/latest_assembly_versions/GCF_020171115.1_ASM_NN_V1/"
        self.remote_filename = "GCF_020171115.1_ASM_NN_V1_cds_from_genomic.fna.gz"
        self.host_name = "neogale_vison"
        self.common_name = "mink"
