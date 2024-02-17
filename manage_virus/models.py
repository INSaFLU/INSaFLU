from django.db import models
from constants.constants import Constants
from manage_virus.constants_virus import ConstantsVirus
# Create your models here.


class Tags(models.Model):
    name = models.CharField(max_length=50, blank=True, null=True)
    def __str__(self):
        return self.name
    
    class Meta:
        ordering = ['name', ]

class UploadFile(models.Model):
    name = models.CharField(max_length=100, blank=True, null=True)
    path = models.CharField(max_length=500, blank=True, null=True)
    version = models.IntegerField(blank=True, null=True)
    creation_date = models.DateTimeField('uploaded date', auto_now_add=True)
    abricate_name = models.CharField(max_length=100, blank=True, null=True)
    def __str__(self):
        return self.name
    
    class Meta:
        ordering = ['version', ]

class SeqVirus(models.Model):
    """
    Has the virus identification
    The file must have this name <name>_<V1>.fasta, to get the version of the file
    """
    accession = models.CharField(max_length=50, blank=True, null=True)
    name = models.CharField(max_length=100, blank=True, null=True)
    kind_type = models.ForeignKey(Tags, related_name='seq_virus', blank=True, null=True, on_delete=models.CASCADE)
    file = models.ForeignKey(UploadFile, related_name='seq_virus', blank=True, null=True, on_delete=models.CASCADE)
    description = models.CharField(max_length=500, blank=True, null=True)
    def __str__(self):
        return self.name
    
    class Meta:
        ordering = ['accession', ]

class IdentifyVirus(models.Model):
    """
    identify the virus
    """
    seq_virus = models.ForeignKey(SeqVirus, related_name='identify_virus', blank=True, null=True, on_delete=models.CASCADE)
    rank = models.IntegerField(default=0)	## has the rank, first the type, then the subType
    coverage = models.CharField(max_length=10, blank=True, null=True)
    identity = models.CharField(max_length=10, blank=True, null=True)
    
    def __str__(self):
        return "name:{}     type:{}    rank:{}".format(self.seq_virus.name, self.seq_virus.kind_type, self.rank)
    
    def __eq__(self, other):
        return self.seq_virus.name == other.seq_virus.name and self.seq_virus.kind_type == other.seq_virus.kind_type
    
    def classify(self, identifyvirus_list):
        
        ### clean some possible repeated
        vect_identify_virus = []
        for identify_virus in identifyvirus_list:
            if not identify_virus in vect_identify_virus:
                vect_identify_virus.append(identify_virus)

        if len(vect_identify_virus) > 0:
            ### Corona/monkeypox
            sz_return_c = self.get_type(
                vect_identify_virus,
                ConstantsVirus.SEQ_VIRUS_GENUS,
                [ConstantsVirus.TYPE_BetaCoV],
            )
            sz_return_c += self.get_type(
                vect_identify_virus, ConstantsVirus.SEQ_VIRUS_GENOTYPE, []
            )

            ## get several species
            sz_subtype = ""
            for sub_type_available in ConstantsVirus.VECT_SUB_TYPE_other_than_influenza:
                sz_temp = self.get_type(vect_identify_virus, sub_type_available, [])
                if len(sz_temp) == 0:
                    continue
                if len(sz_subtype) > 0:
                    sz_subtype += sz_subtype + "-" + sz_temp
                else:
                    sz_subtype = sz_temp
            if len(sz_subtype) > 0:
                sz_return_c += sz_subtype if len(sz_return_c) == 0 else "-" + sz_subtype

            ### Type A
            ## the second flag is for corona
            sz_return_a = self.get_type(
                vect_identify_virus,
                ConstantsVirus.SEQ_VIRUS_TYPE,
                [ConstantsVirus.TYPE_A, ConstantsVirus.TYPE_BetaCoV],
            )
            sz_subtype = self.get_type(
                vect_identify_virus,
                ConstantsVirus.SEQ_VIRUS_SUB_TYPE,
                [ConstantsVirus.TYPE_A],
            )
            if len(sz_subtype) > 0:
                sz_return_a += sz_subtype if len(sz_return_a) == 0 else "-" + sz_subtype

            ### Type B
            sz_return_b = self.get_type(
                vect_identify_virus,
                ConstantsVirus.SEQ_VIRUS_TYPE,
                [ConstantsVirus.TYPE_B],
            )
            sz_subtype = self.get_type(
                vect_identify_virus,
                ConstantsVirus.SEQ_VIRUS_LINEAGE,
                [ConstantsVirus.TYPE_B],
            )
            if len(sz_subtype) > 0:
                sz_return_b += sz_subtype if len(sz_return_b) == 0 else "-" + sz_subtype

            sz_return = sz_return_c
            if len(sz_return_a) > 0:
                sz_return = (
                    sz_return_a
                    if len(sz_return) == 0
                    else "{}; {}".format(sz_return, sz_return_a)
                )
            if len(sz_return_b) > 0:
                sz_return = (
                    sz_return_b
                    if len(sz_return) == 0
                    else "{}; {}".format(sz_return, sz_return_b)
                )
            return sz_return
        
        return Constants.EMPTY_VALUE_TYPE_SUBTYPE


    def get_type(self, vect_identify_virus, type_to_test, vect_virus_name):
        vect_return = []
        for identify_virus in vect_identify_virus:
            if identify_virus.seq_virus.kind_type.name == type_to_test and (
                identify_virus.seq_virus.name in vect_virus_name
                or len(vect_virus_name) == 0
            ):
                vect_return.append(identify_virus.seq_virus.name)
            elif (
                type_to_test == ConstantsVirus.SEQ_VIRUS_SUB_TYPE
                and identify_virus.seq_virus.kind_type.name == type_to_test
                and ConstantsVirus.TYPE_A in vect_virus_name
            ):
                vect_return.append(identify_virus.seq_virus.name)
            elif (
                type_to_test == ConstantsVirus.SEQ_VIRUS_LINEAGE
                and identify_virus.seq_virus.kind_type.name == type_to_test
                and ConstantsVirus.TYPE_B in vect_virus_name
            ):
                vect_return.append(identify_virus.seq_virus.name)
        if type_to_test == ConstantsVirus.SEQ_VIRUS_SUB_TYPE and len(vect_return) > 2:
            return "|".join(sorted(vect_return))
        if type_to_test == ConstantsVirus.SEQ_VIRUS_HUMAN and len(vect_return) > 1:
            return "|".join(sorted(vect_return))
        return "".join(sorted(vect_return))

    def exists_type(
        self, vect_identify_virus, virus_name, type_virus=ConstantsVirus.SEQ_VIRUS_TYPE
    ):
        """
        test if exist some specific type
        """
        for identify_virus in vect_identify_virus:
            if identify_virus.seq_virus.kind_type.name == type_virus and (
                identify_virus.seq_virus.name == virus_name or len(virus_name) == 0
            ):
                return True
        return False

    def get_number_type(self, vect_identify_virus, type_to_test):
        """
        get a number for a specific type
        """
        n_return = 0
        for identify_virus in vect_identify_virus:
            if identify_virus.seq_virus.kind_type.name == type_to_test:
                n_return += 1
        return n_return

    def get_number_type_and_start_sub_type(
        self, vect_identify_virus, type_to_test, starts_sub_type
    ):
        """
        get a number for a specific type and name starts with a specific sub_type
        if (starts_sub_type == None) get all of a specific
        """
        n_return = 0
        for identify_virus in vect_identify_virus:
            if identify_virus.seq_virus.kind_type.name == type_to_test and (
                starts_sub_type == None
                or identify_virus.seq_virus.name.startswith(starts_sub_type)
            ):
                n_return += 1
        return n_return

