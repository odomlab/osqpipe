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
  condition  = forms.CharField(required=False)
  antibody   = forms.CharField(required=False)
  genomicsid = forms.CharField(required=False, label='Genomics SLX ID')
  flowcell   = forms.CharField(required=False, label='Flowcell ID')
  accession  = forms.CharField(required=False, label='Public Accession')

class LibraryEditForm(forms.Form):

  bad          = forms.BooleanField(required=False, label='Library Failed')
  comment      = forms.CharField(required=False, label='Comments', widget=forms.Textarea)

class LibraryProjectChoiceField(forms.ModelMultipleChoiceField):

  def label_from_instance(self, obj):
    if hasattr(obj, 'is_frozen') and obj.is_frozen:
      slug = "%s *"
    else:
      slug = "%s"
    return slug % obj

class LibraryProjectPicker(forms.Form):

  projects = LibraryProjectChoiceField(queryset=Project.objects.none(),
                                       label='Project membership')

  # We need to be able to set the project listing based on the
  # logged-in user. See http://stackoverflow.com/a/4880869
  def __init__(self, *args, **kwargs):
    qs = kwargs.pop('projects')
    super(LibraryProjectPicker, self).__init__(*args, **kwargs)
    self.fields['projects'].queryset = qs
    self.fields['projects'].help_text = '' # see https://code.djangoproject.com/ticket/9321
