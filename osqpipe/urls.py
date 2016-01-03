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
  url(r'^download/(?P<cls>alnfile|lanefile|qcfile|peakfile|mergedalnfile|histologyimagefile)/(?P<pk>\d+)$',
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

  # The following define our REST API URLs.
  url(r'^api/', include(router.urls, namespace='api')),
  url(r'^api/download/(?P<cls>lanefile)/(?P<pk>\d+)$',
      views.FileDownloadView.as_view(),
      name='api:file-download'), # FIXME is this appropriate?
  url(r'^api-auth/', include('rest_framework.urls', namespace='rest_framework')),
]
