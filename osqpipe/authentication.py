'''
Authentication classes for Django Rest Framework.

Classes:
  ExpiringTokenAuthentication: Authentication using extended authtoken model.
'''

from rest_framework import exceptions
from rest_framework.authentication import TokenAuthentication
from rest_framework.authtoken.models import Token

from django.conf import settings
from django.utils import timezone
from datetime import timedelta

TOKEN_EXPIRE_HOURS = getattr(settings, 'REST_FRAMEWORK', {}).get('TOKEN_EXPIRE_HOURS', 24)

class ExpiringTokenAuthentication(TokenAuthentication):
  '''
  Extends default token auth to have time-based expiration.
  Based on http://stackoverflow.com/questions/14567586/
  '''
  def authenticate_credentials(self, key):
    '''
    Attempt token authentication using the provided key.
    '''
    try:
      token = self.model.objects.get(key=key)
    except self.model.DoesNotExist:
      raise exceptions.AuthenticationFailed('Invalid token')

    if not token.user.is_active:
      raise exceptions.AuthenticationFailed('User inactive or deleted')

    if token.created < timezone.now() - timedelta(hours=TOKEN_EXPIRE_HOURS):
      raise exceptions.AuthenticationFailed('Token has expired')

    return (token.user, token)
