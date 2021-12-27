'''
Created on 01/03/2020

@author: mmp
'''

from django.conf import settings
from constants.constantsTestsCase import ConstantsTestsCase
from django.test import TestCase
from django.contrib.auth.models import User
from django.utils.http import urlsafe_base64_encode
from django.template.loader import render_to_string
from fluwebvirus.tokens import account_activation_token
from django.utils.encoding import force_bytes
from django.test.utils import override_settings

class Test(TestCase):


	@override_settings(EMAIL_BACKEND='django.core.mail.backends.smtp.EmailBackend')
	def test_send_test_email(self):
		"""
		Only run this to a real configuration test
		"""
		try:
			user = User.objects.get(username=ConstantsTestsCase.TEST_USER_NAME)
		except User.DoesNotExist:
			user = User()
			user.first_name = "Manel"
			user.last_name = "Maria"
			user.username = ConstantsTestsCase.TEST_USER_NAME
			user.is_active = True
			user.email = settings.EMAIL_DESTINATION_TO_SEND_A_TEST
			user.password = ConstantsTestsCase.TEST_USER_NAME
			user.save()
			
		message = render_to_string('accounts/account_activation_email.txt', {
				'user': user,
				'domain': "insa.min-saude.pt",
				'uid': urlsafe_base64_encode(force_bytes(user.pk)),
				'token': account_activation_token.make_token(user),
			})
		
		messageHTML = render_to_string('accounts/account_activation_email.html', {
				'user': user,
				'domain': "insa.min-saude.pt",
				'uid': urlsafe_base64_encode(force_bytes(user.pk)),
				'token': account_activation_token.make_token(user),
			})
		
		print("Email Backen: {}".format(settings.EMAIL_BACKEND))
		print("Try to send email to: {}".format(user.email))
		subject = "Test INSaFLU email, it's the same of activation account. It is only used for test performs."
		user.email_user(subject, message, from_email=settings.DEFAULT_FROM_EMAIL, html_message=messageHTML)
		print("Please, check the account '{}'. If everything worked fine you must and email in previous mail box.".format(settings.EMAIL_DESTINATION_TO_SEND_A_TEST))

