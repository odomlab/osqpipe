#
# Copyright 2018 Odom Lab, CRUK-CI, University of Cambridge
#
# This file is part of the osqpipe python package.
#
# The osqpipe python package is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# The osqpipe python package is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the osqpipe python package.  If not, see
# <http://www.gnu.org/licenses/>.

'''
Custom project-level permissions.
'''

from .models import Project, Library, Lane, Alignment, Sample, MergedAlignment
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
    elif type(obj) is Lane:
      return request.user in [ person
                               for prj in obj.library.projects.all()
                               for person in prj.people.all() ]
    elif type(obj) is Alignment:
      return request.user in [ person
                               for prj in obj.lane.library.projects.all()
                               for person in prj.people.all() ]
    elif type(obj) is MergedAlignment:
      return request.user in [ person
                               for aln in obj.alignments.all()
                               for prj in aln.lane.library.projects.all()
                               for person in prj.people.all() ]
    elif type(obj) is Sample:
      return request.user in [ person
                               for lib in obj.library_set.all()
                               for prj in lib.projects.all()
                               for person in prj.people.all() ]
    else:
      return False
