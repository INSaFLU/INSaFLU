"""fluwebvirus URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.11/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import include, url
from django.contrib import admin
from fluwebvirus.views import LoginView, HomePageView, LogOutView, SignUpView
from fluwebvirus.views import ResetPasswordView, ChangePasswordView, GetMessageConfirmEmailView
from django.conf import settings
from django.conf.urls.static import static
from fluwebvirus.views import activate, reset_password_key

urlpatterns = [
    url('^$', HomePageView.as_view(), name='home'),
    url(r'^accounts/register/$', SignUpView.as_view(), name='register'),
    url(r'^accounts/reset_password/$', ResetPasswordView.as_view(), name='reset_password'),
    url(r'^accounts/get_message_confirm_email/$', GetMessageConfirmEmailView.as_view(), name='get_message_confirm_email'),
    url(r'^accounts/change_password/$', ChangePasswordView.as_view(), name='change_password'),
    url(r'^accounts/login/$', LoginView.as_view(), name='login'),
    url(r'^accounts/logout/$', LogOutView.as_view(), name='logout'),
    url(r'^accounts/activate/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$', activate, name='activate'),
    url(r'^accounts/reset_password_key/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$', reset_password_key, name='reset_password_key'),
    url(r'^accounts/change_password/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$', ChangePasswordView.as_view(), name='change_password'),
    url(r'^dashboard/$', HomePageView.as_view(), name='dashboard'),

    url(r'^admin/', admin.site.urls),
    url(r'^managing_files/', include('managing_files.urls')),
    url(r'^phylogeny/', include('phylogeny.urls')),
]

urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)