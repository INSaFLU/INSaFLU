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

from django.conf import settings
from django.urls import include
from django.conf.urls import re_path
from django.conf.urls.static import static
from django.contrib import admin

from fluwebvirus.views import (
    ChangePasswordView,
    GetMessageConfirmEmailView,
    HomePageView,
    LoginView,
    LogOutView,
    ResetPasswordView,
    SignUpView,
    activate,
    reset_password_key,
)

urlpatterns = []
if settings.ADMIN_ENABLED:
    urlpatterns += [
        url(r"^admin/", admin.site.urls),
    ]

urlpatterns += [
    re_path("^$", HomePageView.as_view(), name="home"),
    re_path(r"^accounts/register/$", SignUpView.as_view(), name="register"),
    re_path(
        r"^accounts/reset_password/$",
        ResetPasswordView.as_view(),
        name="reset_password",
    ),
    re_path(
        r"^accounts/get_message_confirm_email/$",
        GetMessageConfirmEmailView.as_view(),
        name="get_message_confirm_email",
    ),
    re_path(
        r"^accounts/change_password/$",
        ChangePasswordView.as_view(),
        name="change_password",
    ),
    re_path(r"^accounts/login/$", LoginView.as_view(), name="login"),
    re_path(r"^accounts/logout/$", LogOutView.as_view(), name="logout"),
    re_path(
        r"^accounts/activate/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$",
        activate,
        name="activate",
    ),
    re_path(
        r"^accounts/reset_password_key/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$",
        reset_password_key,
        name="reset_password_key",
    ),
    re_path(
        r"^accounts/change_password/(?P<uidb64>[0-9A-Za-z_\-]+)/(?P<token>[0-9A-Za-z]{1,13}-[0-9A-Za-z]{1,20})/$",
        ChangePasswordView.as_view(),
        name="change_password",
    ),
    re_path(r"^dashboard/$", HomePageView.as_view(), name="dashboard"),
    re_path(r"^managing_files/", include("managing_files.urls")),
    re_path(r"^pathogen_identification/", include("pathogen_identification.urls")),
    re_path(r"^phylogeny/", include("phylogeny.urls")),
    re_path(r"^settings/", include("settings.urls")),
    #    url(r"^settings_pf/", include("settings_pf.urls")),
    re_path(r"^datasets/", include("datasets.urls")),
]

if settings.DEBUG is True:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
