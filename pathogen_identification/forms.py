from pathogen_identification.models import Projects


class ReferenceProjectForm(forms.ModelForm):
    """
    Create a new project
    """

    error_css_class = "error"

    class Meta:
        model = Projects
        exclude = ()

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
        Clean all together because it's necessary to compare the genbank and fasta files
        """
        cleaned_data = super(ReferenceProjectForm, self).clean()
        return cleaned_data
