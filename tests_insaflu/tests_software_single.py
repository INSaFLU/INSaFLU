'''
Created on Oct 28, 2017

@author: mmp
'''

from django.test import TestCase
from django.conf import settings 
from constants.constantsTestsCase import ConstantsTestsCase
from constants.constants_mixed_infection import ConstantsMixedInfection
from constants.constants import Constants, TypePath, FileType, FileExtensions
from constants.meta_key_and_values import MetaKeyAndValue
from constants.tag_names_constants import TagNamesConstants
from utils.software import Software, Contigs2Sequences
from utils.software_pangolin import SoftwarePangolin
from constants.software_names import SoftwareNames
from utils.utils import Utils
from utils.parse_out_files import ParseOutFiles
from utils.result import DecodeObjects, Coverage, CountHits, GeneticElement, MaskingConsensus
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference, Statistics
from manage_virus.uploadFiles import UploadFiles
from manage_virus.models import UploadFile
from managing_files.manage_database import ManageDatabase
from django.test.utils import override_settings
from utils.tree import CreateTree
import os, filecmp, csv
from utils.parse_in_files import ParseInFiles
from utils.result import DecodeObjects, MixedInfectionMainVector
from managing_files.models import CountVariations, MixedInfections
from utils.mixed_infections_management import MixedInfectionsManagement
from manage_virus.constants_virus import ConstantsVirus
from manage_virus.models import Tags, SeqVirus, IdentifyVirus
from settings.default_software import DefaultSoftware
from settings.models import Software as Software2, Parameter
from utils.parse_coverage_file import GetCoverage
from plotly.figure_factory._dendrogram import scs
from managing_files.models import Software as SoftwareModel
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class Test(TestCase):

	### static
	software = Software()
	software_names = SoftwareNames()
	software_pangolin = SoftwarePangolin()
	utils = Utils()
	constants = Constants()
	constants_tests_case = ConstantsTestsCase()

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass


	def tearDown(self):
		pass
	
	def test_mask_sequences_by_position(self):
		
		seq_ref =       SeqRecord(Seq("AACA-AC--AAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AA-AAACAC--C"), id="xpto")
		mask_sites = "5,6"
		mask_from_beginning = "20"
		mask_from_end = "20"
		mask_range = "1-20"
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(12, len(str(sequence_consensus.seq)))
		self.assertEqual("AAANNNNNNNNN", str(sequence_consensus.seq))
		
		seq_ref =       SeqRecord(Seq("AACA-AC--AAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AA-AAACAC--C"), id="xpto")
		mask_sites = "1,2"
		mask_from_beginning = "-1"
		mask_from_end = "-1"
		mask_range = ""
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(9, len(str(sequence_consensus.seq)))
		self.assertEqual("AAANNCACC", str(sequence_consensus.seq))

		seq_ref =       SeqRecord(Seq("AACAACAAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AAAAACACC"), id="xpto")
		mask_sites = ""
		mask_from_beginning = "2"
		mask_from_end = "2"
		mask_range = ""
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(11, len(str(sequence_consensus.seq)))
		self.assertEqual("AAANNCACCNN", str(sequence_consensus.seq))
		
		seq_ref =       SeqRecord(Seq("AACAACAAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AACAACAAA"), id="xpto")
		mask_sites = ""
		mask_from_beginning = "2"
		mask_from_end = "2"
		mask_range = ""
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(9, len(str(sequence_consensus.seq)))
		self.assertEqual("NNCAACANN", str(sequence_consensus.seq))
			
		seq_ref =       SeqRecord(Seq("AACAACAAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AACAACAAA"), id="xpto")
		mask_sites = ""
		mask_from_beginning = "2"
		mask_from_end = "2"
		mask_range = "0-40"
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(9, len(str(sequence_consensus.seq)))
		self.assertEqual("NNNNNNNNN", str(sequence_consensus.seq))
		
		seq_ref =       SeqRecord(Seq("AACAACCAAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AACAACAAA"), id="xpto")
		mask_sites = "5,6"
		mask_from_beginning = "0"
		mask_from_end = "0"
		mask_range = "0"
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(9, len(str(sequence_consensus.seq)))
		self.assertEqual("AACANNAAA", str(sequence_consensus.seq))
		
		seq_ref =       SeqRecord(Seq("AACAAAAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AACAACAAA"), id="xpto")
		mask_sites = "5,6"
		mask_from_beginning = "0"
		mask_from_end = "0"
		mask_range = "0"
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(9, len(str(sequence_consensus.seq)))
		self.assertEqual("AACANNAAA", str(sequence_consensus.seq))
		
		seq_ref =       SeqRecord(Seq("AACAAAAAAAAAAAAAAA"), id="xpto")
		seq_consensus = SeqRecord(Seq("AACAACAAAAAAAAAAAAA"), id="xpto")
		mask_sites = "5,6"
		mask_from_beginning = "0"
		mask_from_end = "0"
		mask_range = "0"
		sequence_consensus = self.software.mask_sequence(seq_ref, seq_consensus, mask_sites, mask_from_beginning, mask_from_end, mask_range)
		self.assertEqual(19, len(str(sequence_consensus.seq)))
		self.assertEqual("AACANNAAAAAAAAAAAAA", str(sequence_consensus.seq))
		
		temp_ref_file = self.utils.get_temp_file("ref_test", ".fasta")
		temp_consensus_file = self.utils.get_temp_file("consensus_test", ".fasta")
		vect_record_out = [SeqRecord(Seq("AACAAAAAAAAAAAAAAA"), id="xpto")]
		with open(temp_ref_file, "w") as handle_fasta_out:
			SeqIO.write(vect_record_out, handle_fasta_out, "fasta")
		vect_record_out = [SeqRecord(Seq("AACAACAAAAAAAAAAAAA"), id="xpto")]
		with open(temp_consensus_file, "w") as handle_fasta_out:
			SeqIO.write(vect_record_out, handle_fasta_out, "fasta")
		
		masking_consensus = MaskingConsensus()
		masking_consensus.set_mask_sites("7")
		masking_consensus.set_mask_from_beginning("2")
		masking_consensus.set_mask_from_ends("2")
		masking_consensus.set_mask_regions("10-11")
		genetic_element = GeneticElement()
		genetic_element.set_mask_consensus_element("xpto", masking_consensus)
		self.software.mask_sequence_by_sites(temp_ref_file, temp_consensus_file, genetic_element)
		dt_records = SeqIO.to_dict(SeqIO.parse(temp_consensus_file, "fasta"))
		self.assertEqual("NNCAACNAANNAAAAANNA", str(dt_records["xpto"].seq))
		
		vect_record_out = [SeqRecord(Seq("AACAACAAAAAAAAAAAA"), id="xpto")]
		with open(temp_consensus_file, "w") as handle_fasta_out:
			SeqIO.write(vect_record_out, handle_fasta_out, "fasta")
		masking_consensus = MaskingConsensus()
		masking_consensus.set_mask_sites("7")
		masking_consensus.set_mask_from_beginning("2")
		masking_consensus.set_mask_from_ends("3")
		masking_consensus.set_mask_regions("10-14")
		genetic_element = GeneticElement()
		genetic_element.set_mask_consensus_element("xpto", masking_consensus)
		self.software.mask_sequence_by_sites(temp_ref_file, temp_consensus_file, genetic_element)
		dt_records = SeqIO.to_dict(SeqIO.parse(temp_consensus_file, "fasta"))
		self.assertEqual("NNCAACNAANNNNNANNN", str(dt_records["xpto"].seq))
		
		os.unlink(temp_consensus_file)
		os.unlink(temp_ref_file)
