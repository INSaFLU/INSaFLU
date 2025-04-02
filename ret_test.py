from pathogen_identification.utilities.utilities_views import SampleReadsRetrieve
from managing_files.models import Sample as INSaFLU_Sample
sample = INSaFLU_Sample.objects.get(pk= 25)
ret = SampleReadsRetrieve(sample)
ret.sample_televir_paths()
