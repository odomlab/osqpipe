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

from django.conf.urls import include, url
from django.views.generic import TemplateView
from django.contrib.auth.views import login as authlogin
from django.contrib.auth.views import logout as authlogout
from . import views
from django.core.urlresolvers import reverse_lazy

# REST API router config.
from rest_framework.routers import DefaultRouter
router = DefaultRouter()
router.register(r'projects',  views.ProjectViewSet, base_name='project')
router.register(r'libraries', views.LibraryViewSet, base_name='library')
router.register(r'lanes',     views.LaneViewSet,    base_name='lane')
router.register(r'alignments', views.AlignmentViewSet, base_name='alignment')
router.register(r'mergedalignments', views.MergedAlignmentViewSet, base_name='mergedalignment')

# REST API file download links. This can be extended as necessary to
# Alnfile, QCfile etc.
restdownloadurls = [
  url(r'^download/lanefile/(?P<pk>\d+)$',
      views.RESTFileDownloadView.as_view(),
      name='lanefile-download',
      kwargs={'cls': 'lanefile'}
    ),
  url(r'^download/alnfile/(?P<pk>\d+)$',
      views.RESTFileDownloadView.as_view(),
      name='alnfile-download',
      kwargs={'cls': 'alnfile'}
    ),
  url(r'^download/mergedalnfile/(?P<pk>\d+)$',
      views.RESTFileDownloadView.as_view(),
      name='mergedalnfile-download',
      kwargs={'cls': 'mergedalnfile'}
    ),
]

# Repository app URLs.
urlpatterns = [

  url(r'^$',
      TemplateView.as_view(
      template_name='repository/home.html'),
      name='repo-home',
      ),

  url(r'^project/$',
      views.ProjectListView.as_view(),
      name='project-list',
      ),

  url(r'^library/(?P<project>\w+)/$',
      views.LibraryListView.as_view(),
      name='library-list',
      ),
  url(r'^library/(?P<project>\w+)/search$',
      views.LibrarySearchView.as_view(),
      name='library-search',
      ),
  url(r'^library/(?P<slug>\w+)$',
      views.LibraryDetailView.as_view(),
      name='library-detail',
      ),
  url(r'^library/(?P<slug>\w+)/edit$',
      views.LibraryEditView.as_view(),
      name='library-edit',
      ),
  url(r'^lane/(?P<pk>\d+)$',
      views.LaneDetailView.as_view(),
      name='lane-detail',
      ),
  url(r'^lane/qualplot/(?P<pk>\d+)$',
      views.QualplotDetailView.as_view(),
      name='lane-qualplot',
      ),
  url(r'^sample/(?P<pk>\d+)$',
      views.SampleDetailView.as_view(),
      name='sample-detail',
      ),
  url(r'^download/(?P<cls>alnfile|lanefile|qcfile|alnqcfile|peakfile|mergedalnfile|histologyimagefile)/(?P<pk>\d+)$',
      views.FileDownloadView.as_view(),
      name='file-download',
      ),

  url(r'^tempfile/(?P<file>[A-Za-z0-9_\.-]+)$',
      views.TempfileDownloadView.as_view(),
      name='tempfile-link',
      ),

  url(r'^genome/$',
      views.GenomeListView.as_view(),
      name='genome-list',
      ),

  url(r'^default_genome/$',
      views.DefaultGenomeListView.as_view(),
      name='default-genome-list',
      ),

  # Our login and logout urls are managed within this application but
  # use the django.contrib.auth backend.
  url(r'^login$',  authlogin, name='auth_login'),
  url(r'^logout$', authlogout, {'next_page': reverse_lazy('repo-home')}, name='auth_logout'),
  url(r'^denied/?$',
      TemplateView.as_view(
      template_name='repository/denied.html'),
      name='denied'),

  # The following define our REST API URLs. Beware that the API namespace
  # given here is also reused in serializers.py.
  url(r'^api/', include(router.urls + restdownloadurls, namespace='api')),
  url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
  url(r'^api-token-auth/', views.obtain_expiring_auth_token),
]
