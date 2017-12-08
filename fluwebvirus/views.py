from __future__ import absolute_import

from django.views import generic
from django.contrib.auth import authenticate, login, logout
from django.core.urlresolvers import reverse_lazy
from braces.views import AnonymousRequiredMixin, FormValidMessageMixin, LoginRequiredMixin, MessageMixin
from ipware.ip import get_ip
from log_login.models import LoginHistory
from .forms import RegistrationForm
from django.contrib.auth.models import User
from .forms import LoginForm


class HomePageView(generic.TemplateView):
	"""
	Home page
	"""
	template_name = 'home.html'
	
	def get_context_data(self, **kwargs):
		context = super(HomePageView, self).get_context_data(**kwargs)
		context['nav_dashboard'] = True	## short the size of modal window
		return context

class SignUpView(AnonymousRequiredMixin, FormValidMessageMixin, generic.CreateView):
	"""
	SignUpView
	"""
	form_class = RegistrationForm
	model = User
	template_name = 'accounts/signup.html'


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
		return context
	
	def form_valid(self, form):
		username = form.cleaned_data['username']
		password = form.cleaned_data['password']
		user = authenticate(username=username, password=password)
		if user is not None and user.is_active:
			login(self.request, user)
			
			## set login history
			login_history = LoginHistory()
			login_history.ip = get_ip(self.request)
			login_history.owner = self.request.user
			login_history.operation = LoginHistory.LOGIN_IN
			login_history.description = get_all_info(self.request)
			login_history.save()
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


