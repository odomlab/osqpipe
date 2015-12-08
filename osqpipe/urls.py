from django.conf.urls import include, url
from django.views.generic import DetailView, ListView, UpdateView, TemplateView
from django.contrib.auth.views import login as authlogin
from django.contrib.auth.views import logout as authlogout
from models import Project
from views import ProjectListView, LibraryListView, \
    GenomeListView, DefaultGenomeListView, \
    LibraryDetailView, LaneDetailView, QualplotDetailView, \
    FileDownloadView, TempfileDownloadView, LibrarySearchView, \
    LibraryEditView
from django.core.urlresolvers import reverse_lazy

# Examples:
# url(r'^$', 'cs_pipeline.views.home', name='home'),
# url(r'^cs_pipeline/', include('cs_pipeline.foo.urls')),

# Repository app
app_name    = 'repository'

urlpatterns = [

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
  url(r'^library/(?P<slug>\w+)/edit$',
      LibraryEditView.as_view(),
      name='library-edit',
      ),
  url(r'^lane/(?P<pk>\d+)$',
      LaneDetailView.as_view(),
      name='lane-detail',
      ),
  url(r'^lane/qualplot/(?P<pk>\d+)$',
      QualplotDetailView.as_view(),
      name='lane-qualplot',
      ),
  url(r'^download/(?P<cls>alnfile|lanefile|qcfile|peakfile|mergedalnfile|histologyimagefile)/(?P<pk>\d+)$',
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

  url(r'^default_genome/$',
      DefaultGenomeListView.as_view(),
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
]
