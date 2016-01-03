'''
Custom project-level permissions.
'''

from .models import Project, Library
from rest_framework import permissions

class IsProjectMember(permissions.BasePermission):
  '''
  Custom permission to only allow project members to view it.
  '''
  def has_object_permission(self, request, view, obj):
    '''
    Method assumes that the passed object is a Project instance.
    '''
    # At the moment we only ever want read-only access.
    if request.method not in permissions.SAFE_METHODS:
      return False

    # Access permissions are only allowed to project members.
    if type(obj) is Project:
      return request.user in obj.people.all()
    elif type(obj) is Library:
      return request.user in [ person
                               for prj in obj.projects.all()
                               for person in prj.people.all() ]
    else:
      return False