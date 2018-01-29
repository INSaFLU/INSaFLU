'''
Created on Jan 23, 2018

@author: mmp
'''
import unittest, os, time, filecmp
from utils.utils import Utils
from utils.process_SGE import ProcessSGE
from constants.constants import Constants, TypePath, FileType, FileExtensions
from django.conf import settings 
from django.test.utils import override_settings
from constants.constantsTestsCase import ConstantsTestsCase
from django.contrib.auth.models import User
from managing_files.models import Sample, Project, ProjectSample, Reference, TagName, TagNames
from utils.collect_extra_data import CollectExtraData
from managing_files.manage_database import ManageDatabase
from constants.software_names import SoftwareNames
from constants.meta_key_and_values import MetaKeyAndValue
from utils.result import DecodeObjects, Coverage
from manage_virus.uploadFiles import UploadFiles
from constants.constants_mixed_infection import ConstantsMixedInfection

class Test(unittest.TestCase):

	software_names = SoftwareNames()
	constants_tests_case = ConstantsTestsCase()
	utils = Utils()

	def setUp(self):
		self.baseDirectory = os.path.join(getattr(settings, "STATIC_ROOT", None), ConstantsTestsCase.MANAGING_TESTS)
		pass

	def test_set_script_run_sge(self):
		"""
		
		"""
		process_SGE = ProcessSGE()
		out_dir = self.utils.get_temp_dir()
		self.assertEquals(None, process_SGE.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_GLOBAL, []))
		self.utils.remove_dir(out_dir)
		
	def test_set_script_run_sge_1(self):
		"""
		
		"""
		process_SGE = ProcessSGE()
		out_dir = self.utils.get_temp_dir()
		vect_command = ['vai.te e vira-te']
		path_file = process_SGE.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_GLOBAL, vect_command)
		self.assertEquals(os.path.join(out_dir, ProcessSGE.FILE_NAME_SCRIPT_SGE), path_file)
		self.assertTrue(os.path.getsize(path_file) > 70)
		self.utils.remove_dir(out_dir)
		
	def test_submitte_job(self):
		"""
		submit a job
		"""
		process_SGE = ProcessSGE()
		out_dir = self.utils.get_temp_dir()
		temp_file = self.utils.get_temp_file_from_dir(out_dir, 'test_sge', '.txt')
		if os.path.exists(temp_file): os.unlink(temp_file)

		vect_command = ['echo $HOSTNAME > ' + temp_file]
		vect_command.append('echo "start waiting" >> ' + temp_file)
		vect_command.append('date >> ' + temp_file)
		vect_command.append('sleep 2s')
		vect_command.append('date >> ' + temp_file)
		vect_command.append('echo "end" >> ' + temp_file)
		path_file = process_SGE.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_GLOBAL, vect_command)
		try:
			sge_id = process_SGE.submitte_job(path_file)
		except:
			self.fail('Fail to submit the task')
		self.assertTrue(sge_id != None)
		
		n_count = 0
		while True:
			if (process_SGE.is_finished(sge_id)): break
			time.sleep(1)
			n_count += 1
			if (n_count == 100):
				self.fail('Wait to much time until end of the SGE process')

		self.assertTrue(os.path.getsize(temp_file) > 80 and os.path.getsize(temp_file) < 100)
		self.utils.remove_dir(out_dir)
		
	def test_submitte_job_1(self):
		"""
		submit a job
		"""
		process_SGE = ProcessSGE()
		out_dir = self.utils.get_temp_dir()
		temp_file = self.utils.get_temp_file_from_dir(out_dir, 'test_sge', '.txt')
		if os.path.exists(temp_file): os.unlink(temp_file)

		vect_command = ['echo $HOSTNAME > ' + temp_file]
		vect_command.append('echo "start waiting" >> ' + temp_file)
		vect_command.append('date >> ' + temp_file)
		vect_command.append('sleep 2s')
		vect_command.append('date >> ' + temp_file)
		vect_command.append('echo "end" >> ' + temp_file)
		path_file = process_SGE.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_GLOBAL, vect_command, True)
		try:
			sge_id = process_SGE.submitte_job(path_file)
		except:
			self.fail('Fail to submit the task')
		self.assertTrue(sge_id != None)
		
		n_count = 0
		while True:
			if (process_SGE.is_finished(sge_id)): break
			time.sleep(1)
			n_count += 1
			if (n_count == 100):
				self.fail('Wait to much time until end of the SGE process')

		self.assertFalse(os.path.exists(temp_file))


	def test_submitte_job_2(self):
		"""
		submit a job
		"""
		process_SGE = ProcessSGE()
		out_dir = self.utils.get_temp_dir()
		temp_file = self.utils.get_temp_file_from_dir(out_dir, 'test_sge', '.txt')
		if os.path.exists(temp_file): os.unlink(temp_file)

		vect_command = ['echo $HOSTNAME > ' + temp_file]
		vect_command.append('echo "start waiting" >> ' + temp_file)
		vect_command.append('date >> ' + temp_file)
		vect_command.append('sleep 2s')
		vect_command.append('date >> ' + temp_file)
		vect_command.append('echo "end" >> ' + temp_file)
		vect_command.append('rm xpto_xpt')
		path_file = process_SGE.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_GLOBAL, vect_command, True)
		try:
			sge_id = process_SGE.submitte_job(path_file)
		except:
			self.fail('Fail to submit the task')
		self.assertTrue(sge_id != None)
		
		n_count = 0
		while True:
			if (process_SGE.is_finished(sge_id)): break
			time.sleep(1)
			n_count += 1
			if (n_count == 100):
				self.fail('Wait to much time until end of the SGE process')

		self.assertTrue(os.path.exists(temp_file))
		self.assertTrue(os.path.getsize(temp_file) > 80 and os.path.getsize(temp_file) < 100)
		self.utils.remove_dir(out_dir)

	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_create_tree_and_alignments(self):
		"""
 		test global method
 		"""
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		path_destination = os.path.join(getattr(settings, "MEDIA_ROOT_TEST", None), ConstantsTestsCase.DIR_PROJECTS)
		self.utils.make_path(os.path.join(getattr(settings, "MEDIA_ROOT_TEST", None), ConstantsTestsCase.DIR_PROJECTS))
		cmd = 'cp -r {}/{}/* {}'.format(self.baseDirectory, ConstantsTestsCase.DIR_PROJECTS, path_destination)
		os.system(cmd)
		
		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME + '5000')
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME + '5000'
			user.id = 5000
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stage_2_ test_create_tree"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.display_name = ref_name
			reference.reference_fasta.name = fasta_file
			reference.reference_fasta_name = os.path.basename(fasta_file)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
		
		project_name = "several_names_test_create_tree"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.id= 5000
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()
		
		## get all fastq files
		vect_files = self.constants_tests_case.get_all_fastq_files(self.baseDirectory)
		
		tag_name_name = 'xpto'
		try:
			tag_name = TagName.objects.get(name='xpto')
		except TagName.DoesNotExist as e:
			tag_name = TagName()
			tag_name.owner = user
			tag_name.name = tag_name_name
			tag_name.is_meta_data = False
			tag_name.save()
		
		ProjectSample.objects.all().delete()
		Sample.objects.all().delete()
		temp_dir = self.utils.get_temp_dir()
		count = 0
		for vect_file in vect_files:
			self.utils.copy_file(vect_file[0], os.path.join(temp_dir, os.path.basename(vect_file[0])))
			self.utils.copy_file(vect_file[1], os.path.join(temp_dir, os.path.basename(vect_file[1])))
			
			n_id = 5000 + count + 1
			sample_name = "_".join(os.path.basename(vect_file[0]).split('_')[0:2])
			try:
				
				sample = Sample.objects.get(pk = n_id)
			except Sample.DoesNotExist as e:
				sample = Sample()
				sample.id = n_id
				sample.name = sample_name
				sample.is_valid_1 = True
				sample.file_name_1 = os.path.basename(vect_file[0])
				sample.path_name_1.name = os.path.join(temp_dir, os.path.basename(vect_file[0]))
				sample.is_valid_2 = False
				sample.file_name_2 = os.path.basename(vect_file[1])
				sample.path_name_2.name = os.path.join(temp_dir, os.path.basename(vect_file[1]))
				sample.owner = user
				sample.is_ready_for_projects = True
				sample.is_obsolete = False
				sample.type_subtype = 'xpto, zpto'
				sample.save()

				## add tag names to sample
				tag_names = TagNames()
				tag_names.value = tag_name_name + " " + tag_name_name
				tag_names.tag_name = tag_name
				tag_names.sample = sample
				tag_names.save()
			
			n_id = 5000 + count + 1
			try:
				project_sample = ProjectSample.objects.get(pk = n_id)
			except ProjectSample.DoesNotExist as e:
				## create project_sample
				project_sample = ProjectSample()
				project_sample.id = n_id
				project_sample.sample = sample
				project_sample.project = project
				project_sample.is_finished = True
				project_sample.is_deleted = False
				project_sample.is_error = False
				project_sample.save()
			count += 1
		
		### test group set tag names
		self.assertEquals(1, project_sample.sample.get_tag_names().count()) 
		
		## launch the process
		process_SGE = ProcessSGE()
		out_dir = self.utils.get_temp_dir()
		temp_file = self.utils.get_temp_file_from_dir(out_dir, 'test_sge', '.txt')
		if os.path.exists(temp_file): os.unlink(temp_file)

		vect_command = ['python3 {} collect_global_files --project_id {} --user_id {} --settings fluwebvirus.settings_test'.format(\
			os.path.join(settings.BASE_DIR, 'manage.py'), project.pk, user.pk)]
		path_file = process_SGE.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_GLOBAL, vect_command)
		try:
			sge_id = process_SGE.submitte_job(path_file)
		except:
			self.fail('Fail to submit the task')
		self.assertTrue(sge_id != None)
		
		n_count = 0
		while True:
			if (process_SGE.is_finished(sge_id)): break
			time.sleep(5)
			n_count += 1
			if (n_count == 100):
				self.fail('Wait to much time until end of the SGE process')

		### instance 
		collect_extra_data = CollectExtraData();
		
		### samples test
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output.csv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_COMMA)
		
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_sample_output.tsv")
		out_file = collect_extra_data.collect_sample_table(project, Constants.SEPARATOR_TAB)
		
		self.assertTrue(os.path.exists(out_file))
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### collect variations from snippy
		out_file = collect_extra_data.collect_variations_snippy(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_snippy.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### collect variations from freebayes
		out_file = collect_extra_data.collect_variations_freebayes(project)
		expected_file_samples = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_GLOBAL_PROJECT, "insa_flu_variations_freebayes.tsv")
		self.assertTrue(filecmp.cmp(out_file, expected_file_samples))
		if (os.path.exists(out_file)): os.unlink(out_file)
		
		### get sample result file
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(out_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_process_second_stage_analysis_single_file(self):
		"""
 		test global method
 		"""
		self.assertEquals(getattr(settings, "MEDIA_ROOT_TEST", None), getattr(settings, "MEDIA_ROOT", None))
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT_TEST", None))
		self.utils.make_path(getattr(settings, "MEDIA_ROOT_TEST", None))

		gb_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_GBK)
		fasta_file = os.path.join(self.baseDirectory, ConstantsTestsCase.MANAGING_DIR, ConstantsTestsCase.MANAGING_FILES_FASTA)
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		ref_name = "second_stagis_single_file2_2"
		try:
			reference = Reference.objects.get(name=ref_name)
		except Reference.DoesNotExist:
			reference = Reference()
			reference.name = ref_name
			reference.reference_fasta.name = fasta_file
			reference.reference_fasta_name = os.path.basename(fasta_file)
			reference.reference_genbank.name = gb_file
			reference.reference_genbank_name = os.path.basename(gb_file)
			reference.owner = user
			reference.save()
			
		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
			
		sample_name = "run_snippyis_single_file1_3"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = False
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()

		project_name = "file_naais_single_filee_3_4"
		try:
			project = Project.objects.get(name=project_name)
		except Project.DoesNotExist:
			project = Project()
			project.name = project_name
			project.reference = reference
			project.owner = user
			project.save()
		
		## create project_sample
		project_sample = ProjectSample()
		project_sample.sample = sample
		project_sample.project = project
		project_sample.save()
		
		### copy to trimmomatic files
		self.utils.copy_file(file_1, project_sample.sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True))
		self.utils.copy_file(file_2, project_sample.sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False))
		
		manageDatabase = ManageDatabase()
		metaKeyAndValue = MetaKeyAndValue()
		meta_key_project_sample = metaKeyAndValue.get_meta_key_queue_by_project_sample_id(project_sample.id)
		manageDatabase.set_project_sample_metakey(project_sample, user, meta_key_project_sample, MetaKeyAndValue.META_VALUE_Queue, "meta_sample.description")
		
		## launch the process
		process_SGE = ProcessSGE()
		out_dir = self.utils.get_temp_dir()
		temp_file = self.utils.get_temp_file_from_dir(out_dir, 'test_sge', '.txt')
		if os.path.exists(temp_file): os.unlink(temp_file)

		vect_command = ['python3 {} second_stage_snippy --project_sample_id {} --user_id {} --settings fluwebvirus.settings_test'.format(\
			os.path.join(settings.BASE_DIR, 'manage.py'), project_sample.pk, user.pk)]
		path_file = process_SGE.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_GLOBAL, vect_command)
		try:
			sge_id = process_SGE.submitte_job(path_file)
		except:
			self.fail('Fail to submit the task')
		self.assertTrue(sge_id != None)
		
		n_count = 0
		while True:
			if (process_SGE.is_finished(sge_id)): break
			time.sleep(5)
			n_count += 1
			if (n_count == 200):
				self.fail('Wait to much time until end of the SGE process')

		### start testing				
		try:
			project_sample = ProjectSample.objects.get(pk=project_sample.id)
		except ProjectSample.DoesNotExist:
			self.fail("Must exist")
		self.assertTrue(project_sample.is_finished)
		
		### test the files
		print(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FA, SoftwareNames.SOFTWARE_SNIPPY_name))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FA, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CONSENSUS_FASTA, SoftwareNames.SOFTWARE_SNIPPY_name).index(SoftwareNames.SOFTWARE_SNIPPY_name.lower()) != -1)
		
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_BAM_BAI, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertFalse(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_DEPTH_GZ_TBI, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_CSV, SoftwareNames.SOFTWARE_SNIPPY_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_SNIPPY_name)))
		## freebayes
		self.assertTrue(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name).index(SoftwareNames.SOFTWARE_FREEBAYES_name.lower()) != -1)
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_VCF_GZ, SoftwareNames.SOFTWARE_FREEBAYES_name)))

		## test freebayes clean
		self.assertTrue(os.path.exists(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)))
		self.assertTrue(os.path.getsize(project_sample.get_file_output(TypePath.MEDIA_ROOT, FileType.FILE_TAB, SoftwareNames.SOFTWARE_FREEBAYES_name)) > 10)

		## test consensus file
		self.assertTrue(os.path.exists(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)))
		self.assertTrue(os.path.getsize(project_sample.get_consensus_file(TypePath.MEDIA_ROOT)) > 10)

		### human file name, snippy tab
		self.assertTrue(os.path.exists(project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB,  self.software_names.get_snippy_name())))
		self.assertTrue(os.path.getsize(project_sample.get_file_output_human(TypePath.MEDIA_ROOT, FileType.FILE_TAB,  self.software_names.get_snippy_name())) > 10)
		
		### set flag that is finished
		list_meta = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Count_Hits, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Count_Hits, list_meta[0].meta_tag.name)
		
		### get the hits value
		decode_coverage = DecodeObjects()
		count_hits = decode_coverage.decode_result(list_meta[0].description)
		self.assertEquals(0, count_hits.get_hits_50_90())
		self.assertEquals(3, count_hits.get_hits_less_50())
		self.assertEquals(3, count_hits.get_total_50_50_90())
		self.assertEquals(125, count_hits.get_total())
		
		self.assertEquals(3, project_sample.count_variations.var_less_50)
		self.assertEquals(0, project_sample.count_variations.var_bigger_50_90)
		self.assertEquals(122, project_sample.count_variations.var_bigger_90)
		self.assertEquals(0, project_sample.alert_first_level)
		self.assertEquals(0, project_sample.alert_second_level)
		
		### test mixed infections
		self.assertEquals('0.815100969586423', '{}'.format(project_sample.mixed_infections.average_value))
		self.assertEquals('No', project_sample.mixed_infections.tag.name)
		self.assertFalse(project_sample.mixed_infections.has_master_vector)
		
		list_meta = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Snippy_Freebayes, None)
		self.assertEquals(1, len(list_meta))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Snippy_Freebayes, list_meta[0].meta_tag.name)
		
		decode_result = DecodeObjects()
		result = decode_result.decode_result(list_meta[0].description)
		self.assertTrue(result is not None)
		self.assertEquals("Snippy-3.2-dev; (--mapqual 20 --mincov 10 --minfrac 0.51)", result.get_software(self.software_names.get_snippy_name()))
		self.assertEquals("Freebayes-v1.1.0-54-g49413aa; (-p 2 -q 20 -m 20 --min-coverage 100 --min-alternate-fraction 0.01 --min-alternate-count 10 -V)",\
 						result.get_software(self.software_names.get_freebayes_name()))

		meta_value = manageDatabase.get_project_sample_metakey(project_sample, MetaKeyAndValue.META_KEY_Coverage, MetaKeyAndValue.META_VALUE_Success)
		decode_coverage = DecodeObjects()
		coverage = decode_coverage.decode_result(meta_value.description)
		self.assertEqual(coverage.get_coverage('PA', Coverage.COVERAGE_ALL), "527.4")
		self.assertEqual(coverage.get_coverage('MP', Coverage.COVERAGE_ALL), "2198.8")
		self.assertEqual(coverage.get_coverage('HA', Coverage.COVERAGE_ALL), "1449.3")
		self.assertEqual(coverage.get_coverage('PA', Coverage.COVERAGE_MORE_9), "100.0")
		self.assertEqual(coverage.get_coverage('MP', Coverage.COVERAGE_MORE_9), "100.0")
		self.assertEqual(coverage.get_coverage('HA', Coverage.COVERAGE_MORE_9), "100.0")
		self.assertEqual(coverage.get_coverage('PA', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('MP', Coverage.COVERAGE_MORE_0), "100.0")
		self.assertEqual(coverage.get_coverage('HA', Coverage.COVERAGE_MORE_0), "100.0")
		
		lst_meta_sample = manageDatabase.get_project_sample_metakey(project_sample, meta_key_project_sample, None)
		self.assertEquals(2, len(lst_meta_sample))
		self.assertEquals(MetaKeyAndValue.META_VALUE_Queue, lst_meta_sample[1].value)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, lst_meta_sample[0].value)
		
		### check if the images exist
		for gene in self.utils.get_elements_from_db(reference, user):
			output_image = project_sample.get_global_file_by_element(TypePath.MEDIA_ROOT, ProjectSample.PREFIX_FILE_COVERAGE, gene, FileExtensions.FILE_PNG)
			self.assertTrue(os.path.exists(output_image))

		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(out_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))


	@override_settings(MEDIA_ROOT=getattr(settings, "MEDIA_ROOT_TEST", None))
	def test_run_fastq_and_trimmomatic_and_identify_species(self):
		"""
		Test run fastq and trimmomatic all together
		"""

		uploadFiles = UploadFiles()
		to_test = True
		(version, file) = uploadFiles.get_file_to_upload(to_test)
		self.assertEqual("2", version)
		self.assertEqual(os.path.join(self.baseDirectory, "db/type_identification/test_db_influenza_typing_v2.fasta"), file)
		uploadFiles.upload_file(version, file)	## upload file
		
		file_1 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_1)
		file_2 = os.path.join(self.baseDirectory, ConstantsTestsCase.DIR_FASTQ, ConstantsTestsCase.FASTQ1_2)

		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = False
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()

		temp_dir = self.utils.get_temp_dir()
		self.utils.copy_file(file_1, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1))
		self.utils.copy_file(file_2, os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2))
			
		sample_name = "run_fastq_and_trimmomatic"
		try:
			sample = Sample.objects.get(name=sample_name)
		except Sample.DoesNotExist:
			sample = Sample()
			sample.name = sample_name
			sample.is_valid_1 = True
			sample.file_name_1 = ConstantsTestsCase.FASTQ1_1
			sample.path_name_1.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_1)
			sample.is_valid_2 = True
			sample.file_name_2 = ConstantsTestsCase.FASTQ1_2
			sample.path_name_2.name = os.path.join(temp_dir, ConstantsTestsCase.FASTQ1_2)
			sample.owner = user
			sample.save()
		
		self.assertFalse(sample.is_ready_for_projects)
		
		### set the job
		taskID = "xpto_task" 
		manageDatabase = ManageDatabase()
		manageDatabase.set_sample_metakey(sample, user, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Queue, taskID)

		## launch the process
		process_SGE = ProcessSGE()
		out_dir = self.utils.get_temp_dir()
		temp_file = self.utils.get_temp_file_from_dir(out_dir, 'test_sge', '.txt')
		if os.path.exists(temp_file): os.unlink(temp_file)

		vect_command = ['python3 {} run_trimmomatic_species --sample_id {} --user_id {} --settings fluwebvirus.settings_test'.format(\
			os.path.join(settings.BASE_DIR, 'manage.py'), sample.pk, user.pk)]
		path_file = process_SGE.set_script_run_sge(out_dir, Constants.QUEUE_SGE_NAME_GLOBAL, vect_command)
		try:
			sge_id = process_SGE.submitte_job(path_file)
		except:
			self.fail('Fail to submit the task')
		self.assertTrue(sge_id != None)
		
		n_count = 0
		while True:
			if (process_SGE.is_finished(sge_id)): break
			time.sleep(5)
			n_count += 1
			if (n_count == 200):
				self.fail('Wait to much time until end of the SGE process')

		### start testing				
		meta_sample = manageDatabase.get_sample_metakey_last(sample, MetaKeyAndValue.META_KEY_Queue_TaskID, MetaKeyAndValue.META_VALUE_Success)
		self.assertTrue(meta_sample != None)
		self.assertEquals(taskID, meta_sample.description)

		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_trimmomatic_file(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, False)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, os.path.basename(sample.get_fastq_output(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, True)))))
		self.assertTrue(os.path.exists(os.path.join(temp_dir, Constants.DIR_PROCESSED_PROCESSED, os.path.basename(sample.get_fastq_trimmomatic(TypePath.MEDIA_ROOT, False)))))
		
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)
		
		### number of sequences
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Number_And_Average_Reads, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Number_And_Average_Reads, list_meta[0].meta_tag.name)

		decodeResultAverageAndNumberReads = DecodeObjects()
		result_average = decodeResultAverageAndNumberReads.decode_result(list_meta[0].description)
		self.assertEqual('41254', result_average.number_file_1)
		self.assertEqual('142.0', result_average.average_file_1)
		self.assertEqual('41254', result_average.number_file_2)
		self.assertEqual('139.1', result_average.average_file_2)

		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Fastq_Trimmomatic, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Fastq(0.11.5), Trimmomatic(0.27)", list_meta[0].description)
		
		sample = Sample.objects.get(pk=sample.id)
		self.assertTrue(sample.is_ready_for_projects)
		self.assertTrue(sample.get_is_ready_for_projects())
		
		file_abricate = sample.get_abricate_output(TypePath.MEDIA_ROOT)
		self.assertTrue(os.path.exists(sample.get_abricate_output(TypePath.MEDIA_ROOT)))
		manageDatabase = ManageDatabase()
		list_meta = manageDatabase.get_sample_metakey(sample, MetaKeyAndValue.META_KEY_Identify_Sample, None)
		self.assertTrue(len(list_meta) == 1)
		self.assertEquals(MetaKeyAndValue.META_VALUE_Success, list_meta[0].value)
		self.assertEquals(MetaKeyAndValue.META_KEY_Identify_Sample, list_meta[0].meta_tag.name)
		self.assertEquals("Success, Spades(3.11.1), Abricate(0.8-dev)", list_meta[0].description)
		self.assertEquals("A-H3N2", sample.type_subtype)
		self.assertEquals("A-H3N2", sample.get_type_sub_type())
		if (os.path.exists(file_abricate)): os.unlink(file_abricate)
		
		### mixed infections
		self.assertEquals(0, sample.number_alerts)
		self.assertEquals(ConstantsMixedInfection.TAGS_MIXED_INFECTION_NO, sample.mixed_infections_tag.name)
		
		## remove all files
		self.utils.remove_dir(temp_dir)
		self.utils.remove_dir(out_dir)
		self.utils.remove_dir(getattr(settings, "MEDIA_ROOT", None))
