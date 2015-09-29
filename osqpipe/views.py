import os
import re

from urllib import urlencode

from django.http import Http404, HttpResponse
from django.utils.translation import ugettext as _
from django.shortcuts import get_object_or_404, redirect
from django.views.generic.edit import FormMixin
from collections import OrderedDict

from models import Library, Project, Genome, Lane, Alnfile, Lanefile, QCfile,\
    Peakfile, MergedAlignment, MergedAlnfile
from forms import SimpleSearchForm, LibrarySearchForm, LibraryEditForm,\
    LibraryProjectPicker

from django.utils.encoding import smart_str
from django.contrib.auth.views import login as django_login, logout as django_logout
from django.core.urlresolvers import reverse

from django.contrib import messages

from urlparse import urlparse # FIXME workaround for qualplot output
from qualplot import plot_pfqual_values, plot_all_qual_values
from pipeline.config import Config
from mimetypes import guess_type

from viewclasses import MyListView, MyDetailView, MyFormView, FilterMixin,\
    FormListView, RestrictedFileDownloadView

CONFIG = Config()

class ProjectListView(MyListView):

  context_object_name = 'projects'
  template_name       = 'repository/project/list.html'
  paginate_by         = 25

  def get_queryset(self):
    '''Restrict the queryset to just those projects assigned to the
    logged-in user.'''
    
    self.queryset = Project.objects.filter(people=self.request.user)
    return super(ProjectListView, self).get_queryset()    

# Note that the order of mixins/superclasses here is important.
class LibraryListView(FilterMixin, FormListView):
  context_object_name = 'libraries'
  template_name       = 'repository/library/list.html'
  paginate_by         = 25
  form_class          = SimpleSearchForm
  allow_empty         = True

  _project = None

  # Note the use of __icontains to provide something like wildcard
  # support. To add a field to the general search initiated by
  # e.g. SimpleSearchForm, just add it here. By convention I think
  # it's best if the search field matches the display field. If the
  # display field is not an intuitive search term it should perhaps be
  # fixed in the database.
  allowed_filters = {
    'code'     : 'code__icontains',
    'libtype'  : 'libtype__name__icontains',
    'strain'   : 'strain__name__icontains',
    'genome'   : 'genome__code__icontains',
    'tissue'   : 'tissue__name__icontains',
    'sex'      : 'sex__name__iexact', # otherwise 'male' is not a useful search term.
    'factor'   : 'factor__name__icontains',
    'antibody' : 'antibody__name__icontains',
    'experiment' : 'chipsample__icontains',
    'individual' : 'individual__icontains',
    'genomicsid' : 'lane__genomicssampleid__icontains',
    'flowcell'   : 'lane__flowcell__icontains',
    'accession'  : 'lane__external_records__accession__icontains',
    }

  def get(self, request, *args, **kwargs):
    self._project = get_object_or_404(Project, code=self.kwargs['project'])

    # Per-project user authorization.
    if self.request.user not in self._project.people.all():
      return redirect('denied')

    return super(LibraryListView, self).get(request, *args, **kwargs)

  def get_queryset(self):
    '''Restrict the queryset to the requested project, and call the
    superclass custom filtering code.'''
    if self._project is None:
      self._project = get_object_or_404(Project, code=self.kwargs['project'])

    self.request.session['session_project'] = self._project.code

    # We'd like to order such that do34 comes before do332. The
    # extra__*_prefix fields here do just that. The select_related
    # clause allows us to greatly reduce the number of SQL queries a
    # given page view triggers.
    self.queryset = Library.objects\
        .select_related('genome', 'tissue', 'libtype', 'strain', 'factor', 'antibody')\
        .prefetch_related('lane_set', 'lane_set__facility', 'lane_set__lanefile_set',
                          'lane_set__alignment_set', 'lane_set__alignment_set__genome',
                          'lane_set__alignment_set__provenance__program', 'lane_set__alignment_set__lane',
                          'lane_set__lanefile_set__filetype',
                          'lane_set__alignment_set__alnfile_set',
                          'lane_set__alignment_set__alnfile_set__filetype')\
        .filter(projects=self._project)\
        .order_by('extra__code_text_prefix',
                  'extra__code_numeric_suffix')
    return super(LibraryListView, self).get_queryset()

  def get_context_data(self, *args, **kwargs):
    '''Add extra context to control the view.'''
    context = super(LibraryListView, self).get_context_data(*args, **kwargs)
    context['linked_filetypes'] = ('fasta', 'fastq', 'bam')
    context['project_code']     = self.kwargs['project']
    return context

class GenomeListView(FilterMixin, FormListView):
  model               = Genome
  context_object_name = 'genomes'
  template_name       = 'repository/genome/list.html'
  paginate_by         = 25
  form_class          = SimpleSearchForm
  allow_empty         = True

  allowed_filters = {
    'code'        : 'code__icontains',
    'commonname'  : 'common_name__icontains',
    'sciname'     : 'scientific_name__icontains',
    'version'     : 'version__icontains',
    }

class DefaultGenomeListView(MyListView):
  context_object_name = 'genomes'
  template_name       = 'repository/genome/list_defaults.html'
  paginate_by         = 9999 # Short list should never really need pagination.
  allow_empty         = True

  def get_queryset(self):
    self.queryset = Genome.objects.filter(code__in=CONFIG.genome_synonyms.values())
    return super(DefaultGenomeListView, self).get_queryset()
  
  def get_context_data(self, *args, **kwargs):
    '''Map genome synonyms to actual genome objects.'''
    context = super(DefaultGenomeListView, self).get_context_data(*args, **kwargs)
    syn     = CONFIG.genome_synonyms
    gen     = dict( (x.code, x) for x in self.get_queryset() )
    syn     = OrderedDict( (k, gen.setdefault(syn[k], '')) for k in sorted(syn) ) 
    context['synonyms'] = syn
    return context

class LibrarySearchView(MyFormView):
  template_name = 'repository/library/search.html'
  form_class    = LibrarySearchForm

  def form_valid(self, form):

    # This is quite ugly code. FIXME if possible?
    search   = { key: form.cleaned_data[key] for key in form.cleaned_data.keys()
                 if form.cleaned_data[key] != '' }
    paramstr = urlencode(search)
    return redirect("%s?%s" % (self.get_success_url(), paramstr))

  def get_success_url(self):
    return reverse('library-list', kwargs={'project':self.kwargs['project']})

  def get(self, request, *args, **kwargs):

    project = get_object_or_404(Project, code=self.kwargs['project'])

    # Per-project user authorization. This might be overkill for a search form.
    if self.request.user not in project.people.all():
      return redirect('denied')

    # Probably not actually necessary, this is just defensive. It does
    # mean that hard links to search pages will set the parent breadcrumbs
    # correctly.
    self.request.session['session_project'] = project.code

    return super(LibrarySearchView, self).get(request, *args, **kwargs)

class LibraryEditView(MyFormView):

  model         = Library
  slug_field    = 'code'
  template_name = 'repository/library/edit.html'
  form_class    = LibraryEditForm

  def post(self, request, *args, **kwargs):

    # From ProcessFormMixin
    form_class = self.get_form_class()
    self.form = self.get_form(form_class)

    # Process the form data here, set an appropriate message.
    if self.form.is_valid():
      library  = get_object_or_404(self.model, code=self.kwargs['slug'])
      library.comment = self.form.cleaned_data['comment']
      library.bad     = self.form.cleaned_data['bad']
      library.save()
    else:
      messages.error(request, "Error in library form submission.")

    return redirect(self.get_success_url())

  def get_success_url(self):
    return reverse('library-detail', kwargs={'slug':self.kwargs['slug']})

  def get(self, request, *args, **kwargs):

    self.object = get_object_or_404(self.model, code=self.kwargs['slug'])

    # Per-project user authorization.
    allowed_users = [ person for project in self.object.projects.all()
                             for person  in project.people.all() ]
    if self.request.user not in allowed_users:
      return redirect('denied')

    # Hard links to edit pages will set the parent breadcrumbs
    # correctly.
    self.request.session['session_library'] = self.object.code

    return super(LibraryEditView, self).get(request, *args, **kwargs)

  def get_form_kwargs(self):

    kwargs   = super(LibraryEditView, self).get_form_kwargs()

    library = get_object_or_404(self.model, code=self.kwargs['slug']) # FIXME dupe

    # Prior data.
    kwargs['initial'] = dict( (key, vars(library)[key]) for key in ('comment','bad') )
    
    return kwargs
  
class LibraryDetailView(FormMixin, MyDetailView):
  context_object_name = 'library'
  template_name       = 'repository/library/detail.html'
  model               = Library
  slug_field          = 'code'
  form_class          = LibraryProjectPicker

  def get(self, request, *args, **kwargs):
    self.object = get_object_or_404(self.model, code=self.kwargs['slug'])

    # Per-project user authorization.
    allowed_users = [ person for project in self.object.projects.all()
                             for person  in project.people.all() ]
    if self.request.user not in allowed_users:
      return redirect('denied')

    self.request.session['session_library'] = self.object.code
    
    # From ProcessFormMixin
    form_class = self.get_form_class()
    self.form  = self.get_form(form_class)
    context    = self.get_context_data(form=self.form)

    # Include any MergedAlignment objects linked to this library in the output.
    context['merged_alignments'] = MergedAlignment.objects.filter(alignments__lane__library=self.object).distinct()
    
    return self.render_to_response(context)

  def assign_projects(self, request):

    object  = get_object_or_404(self.model, code=self.kwargs['slug'])
    default = Project.objects.get(code=CONFIG.defaultproject)
    wanted  = self.form.cleaned_data['projects']
    allowed = request.user.project_set.all()
    current = object.projects.all()
    changecount = 0

    # Test that the logged-in user is a member of the default project
    # (i.e., is a lab member) before allowing any changes.
    if default not in allowed:
      messages.error(request, "You are not authorised to change library project assignments.")
      return

    try:

      # Add new project assignments. Test that the user is a member of
      # all the projects for which we make changes.
      for proj in wanted:
        if proj not in current and proj in allowed:
          object.projects.add(proj)
          changecount += 1

      # Remove old, unwanted assignments (but never the default project!).
      for proj in current:
        if proj not in wanted and proj in allowed:
          if proj == default:
            messages.warning(request, "Cannot remove library from default project.")
          else:
            object.projects.remove(proj)
            changecount += 1

      if changecount > 0:
        messages.info(request, "Successfully reassigned library to projects.")
      else:
        messages.info(request, "No changes made to library project assignment.")
          
    except:
      messages.error(request, "Unable to assign library to projects.")

  def post(self, request, *args, **kwargs):

    # From ProcessFormMixin
    form_class = self.get_form_class()
    self.form = self.get_form(form_class)

    # Process the form data here, set an appropriate message.
    if self.form.is_valid():
      self.assign_projects(request)
    else:
      messages.error(request, "Error in project form submission.")

    # We just redirect to the same page.
    return redirect(request.path)
    
  def get_form_kwargs(self):

    kwargs   = super(LibraryDetailView, self).get_form_kwargs()
    projects = self.request.user.project_set.all().order_by('name')
    kwargs['projects'] = projects

    object = get_object_or_404(self.model, code=self.kwargs['slug']) # FIXME dupe

    # Currently selected projects
    kwargs['initial'] = { 'projects': object.projects.all() }
    
    return kwargs

class LaneDetailView(MyDetailView):
  context_object_name = 'lane'
  template_name       = 'repository/lane/detail.html'
  model               = Lane

  def get(self, request, *args, **kwargs):
    object = get_object_or_404(self.model, id=self.kwargs['pk'])

    # Per-project user authorization.
    allowed_users = [ person for project in object.library.projects.all()
                             for person  in project.people.all() ]
    if self.request.user not in allowed_users:
      return redirect('denied')

    self.request.session['session_lane'] = object.pk
    self.request.session['session_library'] = object.library.code
    return super(LaneDetailView, self).get(request, *args, **kwargs)

class QualplotDetailView(LaneDetailView):
  template_name       = 'repository/lane/qualplot.html'

  def get_context_data(self, *args, **kwargs):
    '''Add extra context to control the view.'''
    context = super(QualplotDetailView, self).get_context_data(*args, **kwargs)
    lane    = get_object_or_404(self.model, id=self.kwargs['pk'])
    url = urlparse(plot_all_qual_values(lane))
    context['qualplot_file_all'] = url.path.split("/")[-1]
    url = urlparse(plot_pfqual_values(lane))
    context['qualplot_file_pf']  = url.path.split("/")[-1]
    return context

# See also http://djangosnippets.org/snippets/2549/
class FileDownloadView(RestrictedFileDownloadView):

  def get(self, request, *args, **kwargs):
    cls = self.kwargs['cls']

    if cls == 'alnfile':
      model = Alnfile
      lanerel = lambda x: x.alignment.lane
    elif cls == 'lanefile':
      model = Lanefile
      lanerel = lambda x: x.lane
    elif cls == 'qcfile':
      model = QCfile
      lanerel = lambda x: x.laneqc.lane
    elif cls == 'peakfile':
      model = Peakfile
      lanerel = lambda x: x.peakcalling.factor_align.lane
    elif cls == 'mergedalnfile':
      model = MergedAlnfile
      lanerel = lambda x: x.alignment.alignments.all()[0].lane # First alignment determines access FIXME.
    else:
      raise ValueError("Unrecognised file class for download: %s" % (cls))

    fobj = get_object_or_404(model, id=self.kwargs['pk'])
    filepath = fobj.repository_file_path

    # Per-project user authorization.
    lane = lanerel(fobj)
    allowed_users = [ person for project in lane.library.projects.all()
                             for person  in project.people.all() ]
    if self.request.user not in allowed_users:
      return redirect('denied')

    # File not found is a standard 404 situation. We could conceivably
    # use 503 for archive unavailable, although it's not a perfect fit
    # either.
    if not os.path.exists(filepath):
      if fobj.archive:
        raise Http404(_(u"Requested file not found: %s. File archive may be unavailable (%s)"
                        % (filepath, fobj.archive.name)))
      else:
        raise Http404(_(u"Requested file not found: %s" % (filepath,)))

    # N.B. no need to set Content-Length as mod_xsendfile does this for us (tested).
    fname = fobj.filename
    if fobj.filetype.gzip:
      fname += ".gz"
    if fobj.filetype.code == 'pdf' and not fobj.filetype.gzip:
      mtype    = guess_type(fname)
      response = HttpResponse(mimetype=mtype[0])
    else:
      response = HttpResponse(mimetype='application/force-download')
      response['Content-Disposition'] = 'attachment; filename=%s' % smart_str(fname)
    response['X-Sendfile']          = smart_str(filepath)
    # You can also set any other required headers: Cache-Control, etc.

    return response

# Requires a login, but we can't limit access per project for
# temporary files.
class TempfileDownloadView(RestrictedFileDownloadView):
  
  slug_field = 'file'

  def get(self, request, *args, **kwargs):

    fname = self.kwargs['file']

    filepath = os.path.join(CONFIG.httptmpdir, fname)
    if not os.path.exists(filepath):
      raise Http404(_(u"Requested temporary file not found: %s" % (fname,)))

    mtype    = guess_type(fname)
    response = HttpResponse(mimetype=mtype[0])
    response['Content-Encoding'] = mtype[1]
    response['X-Sendfile']       = smart_str(filepath)

    return response

# Wrapper methods for the standard django.contrib.auth.views login and
# logout methods which at the moment just add success messages.
def login(request, *args, **kwargs):
  response = django_login(request, *args, **kwargs)
  if request.user.is_authenticated():
    messages.info(request, "Logged in successfully.")
  elif request.method == 'POST':  # Not quite ideal; FIXME
    messages.error(request, "Login Failed.")
  return response
  
def logout(request, *args, **kwargs):
  response = django_logout(request, *args, **kwargs)
  if not request.user.is_authenticated():
    messages.info(request, "Successfully logged out.")
  else:
    messages.error(request, "Failed to log out.")
  return response
  
  
