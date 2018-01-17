from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.models import User
from django.utils.translation import ugettext_lazy as _
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Div, Layout, ButtonHolder, Submit, Button, HTML
from django.urls import reverse
from django import forms
from django.core.exceptions import ValidationError
from django.conf import settings

## https://simpleisbetterthancomplex.com/tutorial/2016/06/13/how-to-send-email.html
class RegistrationForm(UserCreationForm):

	first_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
	last_name = forms.CharField(max_length=30, required=False, help_text='Optional.')
	email = forms.EmailField(max_length=254, required=True, help_text='Required. Inform a valid email address. It is going to be used to valid your account.')
	institution = forms.CharField(max_length=100, required=False, help_text='Optional. Affiliation institutional ')

	class Meta:
		model = User
		fields = ('username', 'first_name', 'last_name', 'email', 'password1', 'password2', )
		
	def __init__(self, *args, **kwargs):
		super(RegistrationForm, self).__init__(*args, **kwargs)

		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		self.helper.layout = Layout(
			Div(
				Div('username', css_class="col-sm-4"),
				Div('first_name', css_class="col-sm-4"),
				Div('last_name', css_class="col-sm-4"),
				css_class='row'
			), 
			Div(
				Div('email', css_class="col-sm-4"),
				Div('password1', css_class="col-sm-4"),
				Div('password2', css_class="col-sm-4"),
				css_class='row'
			),
			Div(
				Div('institution', css_class="col-sm-4"),
				css_class='row'
			),
			Div(
				HTML('<div class="g-recaptcha insa-recaptcha" data-sitekey="{}"></div>'.format(settings.SITE_KEY)),
				css_class='row'
			),
			ButtonHolder(
				Submit('register', 'Register', css_class='btn-primary'),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('dashboard'))),
			)
		)
	
	def clean(self):
		"""
		Test the account name
		"""
		cleaned_data = super(RegistrationForm, self).clean()
		username = cleaned_data['username']
		email = cleaned_data['email']
		try:
			User.objects.get(username__iexact=username)
			self.add_error('username', ValidationError(_("There's an account with this username."), code='invalid'))
		except User.DoesNotExist:
			pass
		
		try:
			User.objects.get(email__iexact=email)
			self.add_error('email', ValidationError(_("There's an account with this email."), code='invalid'))
		except User.DoesNotExist:
			pass
			
		return cleaned_data


class ResetPasswordForm(forms.ModelForm):
	"""
	reset the password
	"""
	email = forms.EmailField(label=_("E-mail"), max_length=254, required=True, 
			widget=forms.TextInput(attrs={
				"type": "email",
		 	   	"size": "30",
			 	"placeholder": _("E-mail address"), }),
			help_text='Required. Inform a valid email address. It is going to be used to reset your password.')
	class Meta:
		model = User
		fields = ('email', )
		
	def __init__(self, *args, **kwargs):
		super(ResetPasswordForm, self).__init__(*args, **kwargs)
		
		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		self.helper.layout = Layout(
			Div(
				Div('email', css_class="col-sm-4"),
				css_class='row'
			),
			Div(
				HTML('<div class="g-recaptcha insa-recaptcha" data-sitekey="{}"></div>'.format(settings.SITE_KEY)),
				css_class='row'
			),
			ButtonHolder(
				Submit('send', 'Send', css_class='btn-primary'),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('dashboard'))),
			)
		)
	
class ChangePasswordForm(forms.ModelForm):

	password1 = forms.CharField(label=_("Password"), required=True, widget=forms.PasswordInput)
	password2 = forms.CharField(label=_("Retype Password"), required=True, widget=forms.PasswordInput)
	
	class Meta:
		model = User
		fields = ('password1', 'password2', )
		
	def __init__(self, *args, **kwargs):
		super(ChangePasswordForm, self).__init__(*args, **kwargs)

		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		self.helper.layout = Layout(
			'password1',
			'password2',
			ButtonHolder(
				Submit('register', 'Change Password', css_class='btn-primary'),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('dashboard'))),
			)
		)

	def clean(self):
		"""
		Test the account name
		"""
		cleaned_data = super(ChangePasswordForm, self).clean()
		password1 = cleaned_data.get('password1')
		password2 = cleaned_data.get('password2')
		if (password1 and password2) and password1 != password2:
			self.add_error(
				'password2', _("You must type the same password each time.")
			)
		return cleaned_data


	
class LoginForm(AuthenticationForm):
	def __init__(self, *args, **kwargs):
		super(LoginForm, self).__init__(*args, **kwargs)

		field_text= [
			('username', 'User name/Email', 'A valid user name or email', True),
			('password', 'Password', '', True),
		]
		for x in field_text:
			self.fields[x[0]].label = x[1]
			self.fields[x[0]].help_text = x[2]
			self.fields[x[0]].required = x[3]
			
		self.helper = FormHelper()
		self.helper.form_method = 'POST'
		self.helper.layout = Layout(
			'username',
			'password',
			ButtonHolder(
				Submit('login', 'Login', css_class='btn-primary'),
				Button('register', 'Register', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('register'))),
				Button('reset_password', 'Reset password', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('reset_password'))),
				Button('cancel', 'Cancel', css_class='btn-secondary', onclick='window.location.href="{}"'.format(reverse('dashboard'))),
			)
		)
