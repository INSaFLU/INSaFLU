from __future__ import absolute_import

from django.views import generic
from django.contrib.auth import authenticate, login, logout
from django.core.urlresolvers import reverse_lazy
from braces.views import AnonymousRequiredMixin, FormValidMessageMixin, LoginRequiredMixin, MessageMixin
from ipware.ip import get_ip
from log_login.models import LoginHistory
from managing_files.models import DataSet
from constants.constants import Constants
from django.contrib import messages
from fluwebvirus.forms import RegistrationForm, LoginForm, ResetPasswordForm, ChangePasswordForm
from django.contrib.auth.models import User
from django.contrib.sites.shortcuts import get_current_site
from django.shortcuts import redirect
from django.utils.encoding import force_bytes
from django.utils.http import urlsafe_base64_encode
from django.template.loader import render_to_string
from fluwebvirus.tokens import account_activation_token
from django.utils.encoding import force_text
from django.utils.http import urlsafe_base64_decode
from django.conf import settings
import urllib, json

class HomePageView(generic.TemplateView):
	"""
	Home page
	"""
	template_name = 'home.html'
	
	def get_context_data(self, **kwargs):
		context = super(HomePageView, self).get_context_data(**kwargs)
		context['nav_dashboard'] = True			
		context['add_google_analytis'] = settings.ADD_GOOGLE_ANALYTICS			
		context['not_show_breadcrumbs'] = True	## to not show breadcrumbs
		context['is_authenticated'] = self.request.user.is_authenticated
		return context

class SignUpView(AnonymousRequiredMixin, FormValidMessageMixin, generic.CreateView):
	"""
	SignUpView
	"""
	form_class = RegistrationForm
	model = User
	template_name = 'accounts/signup.html'

	def get_context_data(self, **kwargs):
		context = super(SignUpView, self).get_context_data(**kwargs)
		context['nav_modal'] = True	## short the size of modal window
		context['not_show_breadcrumbs'] = True	## to not show breadcrumbs
		return context
	
	def form_valid(self, form):
		
		if form.is_valid():
			### Begin reCAPTCHA validation '''
			recaptcha_response = self.request.POST.get('g-recaptcha-response')
			url = 'https://www.google.com/recaptcha/api/siteverify'
			values = {
				'secret': settings.GOOGLE_RECAPTCHA_SECRET_KEY,
				'response': recaptcha_response
			}
			data = urllib.parse.urlencode(values).encode("utf-8")
			req = urllib.request.Request(url, data)
			response = urllib.request.urlopen(req).read()
			result = json.loads(response.decode('utf-8'))
			### End reCAPTCHA validation
			
			if (not result['success']):
				messages.warning(self.request, "Wrong reCAPTCHA. Please, try again.", fail_silently=True)
			else:
				user = form.save(commit=False)
				user.is_active = False
				user.save()
				user.profile.institution = form.cleaned_data['institution']
				user.profile.save()
				
				current_site = get_current_site(self.request)
				subject = 'Activate your InsaFlu Account'
				message = render_to_string('accounts/account_activation_email.html', {
					'user': user,
					'domain': current_site.domain,
					'uid': urlsafe_base64_encode(force_bytes(user.pk)),
					'token': account_activation_token.make_token(user),
				})
				user.email_user(subject, message)
				messages.success(self.request, "An email was sent to validate your account. Please, follow the link in the e-mail.", fail_silently=True)
			return redirect('dashboard')
			
	## static method
	form_valid_message = ""


class ResetPasswordView(AnonymousRequiredMixin, FormValidMessageMixin, generic.CreateView):
	"""
	Reset password
	"""
	form_class = ResetPasswordForm
	model = User
	template_name = 'accounts/reset_password.html'

	def get_context_data(self, **kwargs):
		context = super(ResetPasswordView, self).get_context_data(**kwargs)
		context['nav_modal'] = True	## short the size of modal window
		context['not_show_breadcrumbs'] = True	## to not show breadcrumbs
		return context
	
	def form_valid(self, form):
		
		if form.is_valid():
			email_name = form.cleaned_data['email']
			if (len(email_name.strip()) == 0 or email_name.lower() == Constants.DEFAULT_USER_EMAIL.lower() or\
					email_name.lower() == Constants.USER_ANONYMOUS_EMAIL.lower()):
				messages.warning(self.request, "The account '{}' does not exist in database.".format(email_name), fail_silently=True)
				return redirect('dashboard')
			
			### 
			try:
				### Begin reCAPTCHA validation '''
				recaptcha_response = self.request.POST.get('g-recaptcha-response')
				url = 'https://www.google.com/recaptcha/api/siteverify'
				values = {
					'secret': settings.GOOGLE_RECAPTCHA_SECRET_KEY,
					'response': recaptcha_response
				}
				data = urllib.parse.urlencode(values).encode("utf-8")
				req = urllib.request.Request(url, data)
				response = urllib.request.urlopen(req).read()
				result = json.loads(response.decode('utf-8'))
				### End reCAPTCHA validation
				
				if (result['success']):
					user = User.objects.get(email__iexact=email_name, is_active=True)
					## invalidate current account
					user.is_active = False
					user.save()
					current_site = get_current_site(self.request)
					subject = 'Reseting password  in your INSaFlu Account'
					message = render_to_string('accounts/account_reset_pass_email.html', {
						'user': user,
						'domain': current_site.domain,
						'uid': urlsafe_base64_encode(force_bytes(user.pk)),
						'token': account_activation_token.make_token(user),
					})
					user.email_user(subject, message, from_email=settings.DEFAULT_FROM_EMAIL, auth_user=settings.EMAIL_HOST_USER,\
							auth_password=settings.EMAIL_HOST_PASSWORD, html_message=message)
					messages.success(self.request, "An email was sent to change your account. Please, follow the link.", fail_silently=True)
				else:
					messages.warning(self.request, "Wrong reCAPTCHA. Please, try again.", fail_silently=True)
					
			except User.DoesNotExist as e:
				messages.warning(self.request, "The account '{}' does not exist in database or is disabled.".format(email_name), fail_silently=True)
			
			return redirect('dashboard')
			
	## static method
	form_valid_message = ""
	
class ChangePasswordView(AnonymousRequiredMixin, FormValidMessageMixin, generic.CreateView):
	"""
	Change the passwords
	"""
	form_class = ChangePasswordForm
	model = User
	template_name = 'accounts/change_password.html'

	def get_context_data(self, **kwargs):
		context = super(ChangePasswordView, self).get_context_data(**kwargs)
		context['nav_modal'] = True	## short the size of modal window
		context['not_show_breadcrumbs'] = True	## to not show breadcrumbs
		return context
	
	def form_valid(self, form):
		
		if form.is_valid() and Constants.SESSION_KEY_USER_ID in self.request.session:
			try:
				user = User.objects.get(pk=self.request.session[Constants.SESSION_KEY_USER_ID])
				user.set_password(form.cleaned_data['password1'])
				user.is_active = True
				user.profile.email_confirmed = True
				user.save()
				
				## clean value in session
				del self.request.session[Constants.SESSION_KEY_USER_ID]
				
				messages.success(self.request, "Congratulations, your account has a new password.", fail_silently=True)
				return redirect('dashboard')
			except User.DoesNotExist as e:
				pass

		messages.error(self.request, "Error, fail to change passwords.", fail_silently=True)
		return redirect('dashboard')
			
	## static method
	form_valid_message = ""

class LoginView(AnonymousRequiredMixin, FormValidMessageMixin, generic.FormView):
	"""
	Login
	"""
	form_class = LoginForm
	success_url = reverse_lazy('home')
	template_name = 'accounts/login.html'

	def get_context_data(self, **kwargs):
		context = super(LoginView, self).get_context_data(**kwargs)
		context['nav_modal'] = True	## short the size of modal window
		context['not_show_breadcrumbs'] = True	## to not show breadcrumbs
		return context
	
	def form_valid(self, form):
		username = form.cleaned_data['username']
		password = form.cleaned_data['password']
		
		###
		if ('login_anonymous' in form.cleaned_data and form.cleaned_data['login_anonymous']):
			username = Constants.USER_ANONYMOUS
			password = Constants.USER_ANONYMOUS_PASS
		
		user = authenticate(username=username, password=password)
		### if none try it with email
		if (user is None): user = authenticate(email=username, password=password)
		if user is not None and user.is_active and user.profile.email_confirmed:
			login(self.request, user)
			
			## set login history
			login_history = LoginHistory()
			login_history.ip = get_ip(self.request)
			login_history.owner = self.request.user
			login_history.operation = LoginHistory.LOGIN_IN
			login_history.description = get_all_info(self.request)
			login_history.save()
			
			if (DataSet.objects.filter(owner__id=user.id).count() == 0):
				### need to create it generic
				dataSet = DataSet()
				dataSet.name = Constants.DATA_SET_GENERIC
				dataSet.owner = self.request.user
				dataSet.save()
				
			return super(LoginView, self).form_valid(form)
		else:
			return self.form_invalid(form)

	## static method
	form_valid_message = "You've been logged in. Welcome back!"
	
	## dinamic method instead
#	def get_form_valid_message(self):
#		return u"{0} created!".format(self.object.title)


class LogOutView(LoginRequiredMixin, MessageMixin, generic.RedirectView):
	"""
	Logout
	"""
	url = reverse_lazy('home')
	def get(self, request, *args, **kwargs):
		## set login history
		login_history = LoginHistory()
		## need to set a proxy if they use one...
		## get_trusted_ip(request, trusted_proxies=['23.91.45.15'])
		login_history.ip = get_ip(self.request)
		login_history.owner = self.request.user
		login_history.operation = LoginHistory.LOGIN_OUT
		login_history.save()
		logout(request)
		
		self.messages.success("You've been logged out. Come back soon!")
		return super(LogOutView, self).get(request, *args, **kwargs)

def activate(request, uidb64, token):
	try:
		uid = force_text(urlsafe_base64_decode(uidb64))
		user = User.objects.get(pk=uid)
	except (TypeError, ValueError, OverflowError, User.DoesNotExist):
		user = None

	if user is not None and account_activation_token.check_token(user, token):
		user.is_active = True
		user.profile.email_confirmed = True
		user.save()
		login(request, user)
		
		## set login history
		login_history = LoginHistory()
		login_history.ip = get_ip(request)
		login_history.owner = user
		login_history.operation = LoginHistory.LOGIN_IN
		login_history.description = get_all_info(request)
		login_history.save()
		
		if (DataSet.objects.filter(owner__id=user.id).count() == 0):
			### need to create it generic
			dataSet = DataSet()
			dataSet.name = Constants.DATA_SET_GENERIC
			dataSet.owner = user
			dataSet.save()
		
		messages.success(request, "Congratulations, your account was confirmed. Enjoy the experience.", fail_silently=True)
		return redirect('dashboard')
	else:
		messages.error(request, "Error, your account activation token is invalid.", fail_silently=True)
		return redirect('dashboard')
	
def reset_password_key(request, uidb64, token):
	try:
		uid = force_text(urlsafe_base64_decode(uidb64))
		user = User.objects.get(pk=uid)
	except (TypeError, ValueError, OverflowError, User.DoesNotExist):
		user = None

	if user is not None and account_activation_token.check_token(user, token):
		request.session[Constants.SESSION_KEY_USER_ID] = uid
		return redirect('change_password')
	else:
		messages.error(request, "Error, your reset password token is invalid.", fail_silently=True)
		return redirect('dashboard')


def get_all_info(request):
	"""
	return all info about user
	"""
	sz_return = "is_mobile: {}".format(request.user_agent.is_mobile) # returns True
	sz_return += ";   is_tablet: {}".format(request.user_agent.is_tablet) # returns False
	sz_return += ";   is_touch_capable: {}".format(request.user_agent.is_touch_capable) # returns True
	sz_return += ";   is_pc: {}".format(request.user_agent.is_pc) # returns False
	sz_return += ";   is_bot: {}".format(request.user_agent.is_bot) # returns False

	# Accessing user agent's browser attributes
	sz_return += ";   browser_family: {}".format(request.user_agent.browser.family)  # returns 'Mobile Safari'
	sz_return += ";   browser_version: {}".format(request.user_agent.browser.version_string)   # returns '5.1'

	# Operating System properties
	sz_return += ";   os_family: {}".format(request.user_agent.os.family)  # returns 'iOS'
	sz_return += ";   os_version: {}".format(request.user_agent.os.version_string)  # returns '5.1'

	# Device properties
	sz_return += ";   device_family: {}".format(request.user_agent.device.family)  # returns 'iPhone'
	return sz_return
