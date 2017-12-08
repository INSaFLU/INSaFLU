"""
WSGI config for fluwebvirus project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.11/howto/deployment/wsgi/
"""

import os, sys

from django.core.wsgi import get_wsgi_application

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fluwebvirus.settings")

path = '/usr/local/web_site'
if path not in sys.path:
	sys.path.append(path)
	
application = get_wsgi_application()
