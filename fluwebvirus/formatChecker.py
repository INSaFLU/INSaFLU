from django.db.models import FileField
from django.forms import forms
from django.template.defaultfilters import filesizeformat
from django.utils.translation import ugettext_lazy as _
import logging

class ContentTypeRestrictedFileField(FileField):
	"""
	Same as FileField, but you can specify:
		* content_types - list containing allowed content_types. Example: ['application/pdf', 'image/jpeg', 'video/x-msvideo', 'video/mp4', 'audio/mpeg', 'txt/css', 
						'application/octet-stream']
		* max_upload_size - a number indicating the maximum file size allowed for upload.
			#	2.5MB -   2621440
			#	5MB   -   5242880
			#	10MB  -  10485760
			#	20MB  -  20971520
			#	50MB  -  52428800
			#	100MB - 104857600
			#	250MB - 214958080
			#	500MB - 429916160
	"""
	
	## logging
	logger_debug = logging.getLogger("fluWebVirus.debug")
	logger_production = logging.getLogger("fluWebVirus.production")
	
	def __init__(self, *args, **kwargs):
		if ("max_upload_size" in kwargs): self.max_upload_size = kwargs.pop("max_upload_size")
		else: self.max_upload_size = 10
		if ("content_types" in kwargs): self.content_types = kwargs.pop("content_types")
		else: self.content_types = ['application/pdf', 'txt/css']
		
		super(ContentTypeRestrictedFileField, self).__init__(*args, **kwargs)

	def clean(self, *args, **kwargs):
		data = super(ContentTypeRestrictedFileField, self).clean(*args, **kwargs)
		try:
			file = data.file
			content_type = file.content_type
			
			### Important to catch the content_type 
			self.logger_production.info("Read '{}' size '{}' content type: {}".format(file.name, file._size, content_type))
			self.logger_debug.info("Read '{}' size '{}' content type: {}".format(file.name, file._size, content_type))
			if content_type in self.content_types:
				if file._size > self.max_upload_size:
					self.logger_debug.warning(_('Please, keep file size under %s. Current file size %s') % (filesizeformat(self.max_upload_size), filesizeformat(file._size)))
					self.logger_production.warning(_('Please, keep file size under %s. Current file size %s') % (filesizeformat(self.max_upload_size), filesizeformat(file._size)))
					raise forms.ValidationError(_('Please, keep file size under %s. Current file size %s') % (filesizeformat(self.max_upload_size), filesizeformat(file._size)))
			else:
				self.logger_debug.warning(_("File type '{}' not supported.".format(content_type)))
				self.logger_production.warning(_("File type '{}' not supported.".format(content_type)))
				raise forms.ValidationError(_("File type '{}' not supported.".format(content_type)))
		except AttributeError:
			raise forms.ValidationError(_("File not supported."))
		return data