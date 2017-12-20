from django.contrib import admin

# Register your models here.

from .models import TagName, DataSet, VaccineStatus

class DataSetAdmin(admin.ModelAdmin):
	fields = ('name', )
	class Meta:
		model = DataSet
	
	def has_add_permission(self, request):
		return True

	def save_model(self, request, obj, form, change):
		obj.owner = request.user
		super(DataSetAdmin, self).save_model(request, obj, form, change)


class VacineStatusAdmin(admin.ModelAdmin):
	fields = ('name', )
	class Meta:
		model = VaccineStatus
	
	def has_add_permission(self, request):
		return True

	def save_model(self, request, obj, form, change):
		obj.owner = request.user
		super(VacineStatusAdmin, self).save_model(request, obj, form, change)

admin.site.register(TagName)
admin.site.register(DataSet, DataSetAdmin)
admin.site.register(VaccineStatus)