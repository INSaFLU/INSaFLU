from settings.models import Technology

from pathogen_identification.models import Projects
from django.utils.html import format_html


class ReferenceProjectForm(forms.ModelForm):
    """
    Create a new project
    """

    error_css_class = "error"

    class Meta:
        model = Projects
        exclude = ()

        fields = ("technology", "name", "description")

    def __init__(self, *args, **kwargs):
        super(ReferenceProjectForm, self).__init__(*args, **kwargs)

        field_text = [
            ("name", "Name", "Unique identifier for this project", True),
        ]
        for x in field_text:
            self.fields[x[0]].label = x[1]
            self.fields[x[0]].help_text = x[2]
            self.fields[x[0]].required = x[3]

    def clean(self):
        """
        clean names provide warning if spaces are used
        """

        cleaned_data = super(ReferenceProjectForm, self).clean()
        from django.core.exceptions import ValidationError

        if self.name == "":
            raise ValidationError("Empty error message")

        return cleaned_data
