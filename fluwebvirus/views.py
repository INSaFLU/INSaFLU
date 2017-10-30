from __future__ import absolute_import

from django.views import generic
from django.contrib.auth import authenticate, login, logout
from django.core.urlresolvers import reverse_lazy
from braces.views import AnonymousRequiredMixin, FormValidMessageMixin, LoginRequiredMixin, MessageMixin

from .forms import RegistrationForm
from django.contrib.auth.models import User
from .forms import LoginForm


class HomePageView(generic.TemplateView):
	"""
	Home page
	"""
	template_name = 'home.html'


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
		logout(request)
		self.messages.success("You've been logged out. Come back soon!")
		return super(LogOutView, self).get(request, *args, **kwargs)
