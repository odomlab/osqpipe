from django.conf.urls import patterns, include, url
from django.views.generic import DetailView, ListView, UpdateView, TemplateView
from models import Project
from views import ProjectListView, LibraryListView, GenomeListView, LibraryDetailView, LaneDetailView, QualplotDetailView, FileDownloadView, TempfileDownloadView, LibrarySearchView
from django.core.urlresolvers import reverse_lazy

# Examples:
# url(r'^$', 'cs_pipeline.views.home', name='home'),
# url(r'^cs_pipeline/', include('cs_pipeline.foo.urls')),

# Repository app
urlpatterns = patterns(
  '',

  url(r'^$',
      TemplateView.as_view(
      template_name='repository/home.html'),
      name='repo-home',
      ),

  url(r'^project/$',
      ProjectListView.as_view(),
      name='project-list',
      ),

  url(r'^library/(?P<project>\w+)/$',
      LibraryListView.as_view(),
      name='library-list',
      ),
  url(r'^library/(?P<project>\w+)/search$',
      LibrarySearchView.as_view(),
      name='library-search',
      ),
  url(r'^library/(?P<slug>\w+)$',
      LibraryDetailView.as_view(),
      name='library-detail',
      ),
  url(r'^lane/(?P<pk>\d+)$',
      LaneDetailView.as_view(),
      name='lane-detail',
      ),
  url(r'^lane/qualplot/(?P<pk>\d+)$',
      QualplotDetailView.as_view(),
      name='lane-qualplot',
      ),
  url(r'^download/(?P<cls>alnfile|lanefile|qcfile|peakfile)/(?P<pk>\d+)$',
      FileDownloadView.as_view(),
      name='file-download',
      ),

  url(r'^tempfile/(?P<file>[A-Za-z0-9_\.-]+)$',
      TempfileDownloadView.as_view(),
      name='tempfile-link',
      ),

  url(r'^genome/$',
      GenomeListView.as_view(),
      name='genome-list',
      ),

  # Our login and logout urls are managed within this application but
  # use the django.contrib.auth backend.
  url(r'^login$',  'django.contrib.auth.views.login', name='auth_login'),
  url(r'^logout$', 'django.contrib.auth.views.logout', {'next_page': reverse_lazy('repo-home')}, name='auth_logout'),
  url(r'^denied/?$',
      TemplateView.as_view(
      template_name='repository/denied.html'),
      name='denied'),
)
