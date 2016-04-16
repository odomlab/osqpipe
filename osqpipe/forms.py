from django import forms
from .models import Project, Sex, Libtype

def model_names_choices_list(model):
  choices  = [('', '<select>')]
  choices += [ (x.name, x.name.replace("_", " "))
               for x in model.objects.all() ]
  return choices

class SimpleSearchForm(forms.Form):
  query = forms.CharField()

class LibrarySearchForm(forms.Form):

  # N.B. currently these are being kept in sync with the
  # views.LibraryListView.allowed_filters dict.

  code       = forms.CharField(required=False, label='Library code')
  experiment = forms.CharField(required=False)

  # Note that the current search implementation will become confused
  # if we have libtypes which are substrings of other libtypes.
  libtype    = forms.ChoiceField(required=False, label='Type, e.g. ChIP-Seq',
                                 choices=model_names_choices_list(Libtype))
  
  strain     = forms.CharField(required=False)
  individual = forms.CharField(required=False)
  tissue     = forms.CharField(required=False)

  sex        = forms.ChoiceField(required=False,
                                 choices=model_names_choices_list(Sex))

  genome     = forms.CharField(required=False, label='Requested genome')
  factor     = forms.CharField(required=False)
  antibody   = forms.CharField(required=False)
  genomicsid = forms.CharField(required=False, label='Genomics SLX ID')
  flowcell   = forms.CharField(required=False, label='Flowcell ID')
  accession  = forms.CharField(required=False, label='Public Accession')

class LibraryEditForm(forms.Form):

  bad          = forms.BooleanField(required=False, label='Library Failed')
  comment      = forms.CharField(required=False, label='Comments', widget=forms.Textarea)

class LibraryProjectPicker(forms.Form):

  projects = forms.ModelMultipleChoiceField(queryset=Project.objects.none(),
                                            label='Project membership')

  # We need to be able to set the project listing based on the
  # logged-in user. See http://stackoverflow.com/a/4880869
  def __init__(self, *args, **kwargs):
    qs = kwargs.pop('projects')
    super(LibraryProjectPicker, self).__init__(*args, **kwargs)
    self.fields['projects'].queryset = qs
    self.fields['projects'].help_text = '' # see https://code.djangoproject.com/ticket/9321
