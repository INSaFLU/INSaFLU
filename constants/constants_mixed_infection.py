'''
Created on Dec 14, 2017

@author: mmp
'''
from manage_virus.models import IdentifyVirus
from manage_virus.constants_virus import ConstantsVirus

class ConstantsMixedInfection(object):
    '''
    classdocs
    '''
    ## [<50, 50<value<90], 
    vect_start_compare = [[78, 56],[77, 56],[76, 57],[76, 53],[55, 51],[50, 36],[49, 35],[26,11]]

    TAGS_MIXED_INFECTION_YES = 'Yes'
    TAGS_MIXED_INFECTION_MAY_BE = 'May be'
    TAGS_MIXED_INFECTION_NO = 'No'
    
    ## values to upload to database
    vect_upload_to_database = [TAGS_MIXED_INFECTION_YES, TAGS_MIXED_INFECTION_MAY_BE, TAGS_MIXED_INFECTION_NO]
    
    ### all other values are NO
    threshold_yes = 0.98	##	>=
    threshold_may_be = 0.97	##	>=
    
    def get_tag_by_value(self, value):
        """
        get tab by value
        """
        if (value >= self.threshold_yes): return self.TAGS_MIXED_INFECTION_YES
        if (value >= self.threshold_may_be): return self.TAGS_MIXED_INFECTION_MAY_BE
        return self.TAGS_MIXED_INFECTION_NO
    
    def is_alert(self, tag):
        """
        return true if is necessary generate an alert
        """
        if (tag == self.TAGS_MIXED_INFECTION_YES or tag == self.TAGS_MIXED_INFECTION_MAY_BE): return True
        return False


    def get_mixed_infection(self, vect_identify_virus_temp):
        """
        mixed infection based on the table static/mixed_infections/mixed_infections.xls
        return tuble (tag_mixed_infection, alert, message)
        tag_mixed_infection: ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO or ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES
        alert: positive number, or zero
        message: None, empty or the string with the error message
        """
        identifyvirus = IdentifyVirus()
        if vect_identify_virus_temp.count() == 0:
            return (
                ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO,
                1,
                "Warning: no typing data was obtained (possible reason: low number of influenza reads).",
            )

        vect_identify_virus = []
        for identify_virus in vect_identify_virus_temp:
            if not identify_virus in vect_identify_virus:
                vect_identify_virus.append(identify_virus)

        ## Only A not B
        if identifyvirus.exists_type(
            vect_identify_virus, ConstantsVirus.TYPE_A
        ) and not identifyvirus.exists_type(vect_identify_virus, ConstantsVirus.TYPE_B):
            #  A; #any subtype; > 0 lineage
            if (
                identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE
                )
                > 0
            ):
                return (
                    ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES,
                    1,
                    "Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.",
                )

            ## A; = 2 subtype
            if (
                identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE
                )
                == 2
            ):
                return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 0, None)
            ## A; > 2 subtype
            if (
                identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE
                )
                > 2
            ):
                return (
                    ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES,
                    1,
                    "Warning: more than two subtypes were detected for this sample, suggesting that may represent a 'mixed infection'.",
                )
            ## A < 2 subtype
            if (
                identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE
                )
                < 2
            ):
                return (
                    ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO,
                    1,
                    "Warning: an incomplete subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).",
                )

        ## Only B not A
        if not identifyvirus.exists_type(
            vect_identify_virus, ConstantsVirus.TYPE_A
        ) and identifyvirus.exists_type(vect_identify_virus, ConstantsVirus.TYPE_B):
            #  B; #any subtype; > 0 lineage
            if (
                identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE
                )
                > 0
            ):
                return (
                    ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES,
                    1,
                    "Warning: more than one type/subtypes were detected for this sample, suggesting that may represent a 'mixed infection'.",
                )

            ## B; == 1 lineage
            if (
                identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE
                )
                == 1
            ):
                return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 0, None)
            ## B; > 1 lineage
            if (
                identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE
                )
                > 1
            ):
                return (
                    ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES,
                    1,
                    "Warning: more than one lineage were detected for this sample, suggesting that may represent a 'mixed infection'.",
                )
            ## B; < 1 lineage
            if (
                identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE
                )
                < 1
            ):
                return (
                    ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO,
                    1,
                    "Warning: an incomplete lineage has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).",
                )

        ## A and B
        if identifyvirus.exists_type(
            vect_identify_virus, ConstantsVirus.TYPE_A
        ) and identifyvirus.__exists_type(vect_identify_virus, ConstantsVirus.TYPE_B):
            if (
                identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE
                )
                == 0
                and identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE
                )
                == 0
            ):
                return (
                    ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO,
                    1,
                    "Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.",
                )
            return (
                ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES,
                1,
                "Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.",
            )

        ## BEGIN corona
        if identifyvirus.exists_type(
            vect_identify_virus, "", ConstantsVirus.SEQ_VIRUS_GENUS
        ) or identifyvirus.exists_type(
            vect_identify_virus, "", ConstantsVirus.SEQ_VIRUS_SPECIES
        ):
            if (
                identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_GENUS
                )
                == 1
            ):
                if (
                    identifyvirus.get_number_type(
                        vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SPECIES
                    )
                    == 1
                ):
                    return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 0, None)
                if (
                    identifyvirus.get_number_type(
                        vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SPECIES
                    )
                    > 1
                ):
                    return (
                        ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES,
                        1,
                        "Warning: more than one human BetaCoV virus are likely present in this sample, suggesting that may represent a 'mixed infection'",
                    )
                if (
                    identifyvirus.get_number_type(
                        vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SPECIES
                    )
                    == 0
                ):
                    return (
                        ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO,
                        1,
                        "Warning: an incomplete human BetaCoV identification has been obtained (possible reasons: low number of  human BetaCoV reads, mixed infection, etc)",
                    )
            else:
                if (
                    identifyvirus.get_number_type(
                        vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SPECIES
                    )
                    == 1
                ):
                    return (
                        ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO,
                        1,
                        "Warning: an incomplete human BetaCoV identification has been obtained (possible reasons: low number of  human BetaCoV reads, mixed infection, etc)",
                    )
                else:
                    return (
                        ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES,
                        1,
                        "Warning: more than one human BetaCoV virus are likely present in this sample, suggesting that may represent a 'mixed infection'",
                    )
        ### END corona

        ## not A and not B
        if not identifyvirus.exists_type(
            vect_identify_virus, ConstantsVirus.TYPE_A
        ) and not identifyvirus.exists_type(vect_identify_virus, ConstantsVirus.TYPE_B):
            if (
                identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE
                )
                == 1
                and identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE
                )
                == 0
            ) or (
                identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_SUB_TYPE
                )
                == 1
                and identifyvirus.get_number_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE
                )
                == 0
            ):
                return (
                    ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO,
                    1,
                    "Warning: an incomplete type/subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).",
                )
            else:
                count_type_N = identifyvirus.get_number_type_and_start_sub_type(
                    vect_identify_virus,
                    ConstantsVirus.SEQ_VIRUS_SUB_TYPE,
                    ConstantsVirus.SUB_TYPE_STARTS_N,
                )
                count_type_H = identifyvirus.get_number_type_and_start_sub_type(
                    vect_identify_virus,
                    ConstantsVirus.SEQ_VIRUS_SUB_TYPE,
                    ConstantsVirus.SUB_TYPE_STARTS_H,
                )
                count_type_other = identifyvirus.get_number_type_and_start_sub_type(
                    vect_identify_virus, ConstantsVirus.SEQ_VIRUS_LINEAGE, None
                )

                if count_type_N == 1 and count_type_H == 1 and count_type_other == 0:
                    return (
                        ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO,
                        1,
                        "Warning: an incomplete type/subtype has been assigned (possible reasons: low number of influenza reads, same-subtype mixed infection, etc.).",
                    )
                elif count_type_N > 1 or count_type_H > 1 or count_type_other > 0:
                    return (
                        ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES,
                        1,
                        "Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.",
                    )
            return (
                ConstantsMixedInfection.TAGS_MIXED_INFECTION_YES,
                1,
                "Warning: more than one type/subtype were detected for this sample, suggesting that may represent a 'mixed infection'.",
            )

        ## default
        return (ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, 0, None)