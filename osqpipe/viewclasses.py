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

import re
from urllib import urlencode
from copy import copy

from django.shortcuts import render_to_response, redirect
from django.views.generic import View, ListView, DetailView, FormView
from django.views.generic.edit import FormMixin
from django.utils.decorators import method_decorator
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse_lazy
from django.utils.translation import ugettext as _
from django.http import Http404

'''Base classes and mixins used to create the actual view classes in
views.py.'''

# These three classes are used to populate our view context with some
# details from the session; this is used to provide e.g. breadcrumbs
# via sitetree.
class SessionContextMixin(object):

  def copy_session_to_context(self, context):
    session_re = re.compile('^session_')
    for item in self.request.session.keys():
      if session_re.match(item):
        context[item] = self.request.session[item]
    return context

class MyListView(ListView, SessionContextMixin):

  def get_context_data(self, *args, **kwargs):
    '''Add session context items for navigation.'''
    context = super(MyListView, self).get_context_data(*args, **kwargs)
    context = self.copy_session_to_context(context)
    return context

  @method_decorator(login_required(login_url=reverse_lazy('auth_login'), redirect_field_name=None))
  def dispatch(self, *args, **kwargs):
    return super(MyListView, self).dispatch(*args, **kwargs)

class MyDetailView(DetailView, SessionContextMixin):

  def get_context_data(self, *args, **kwargs):
    '''Add session context items for navigation.'''
    context = super(MyDetailView, self).get_context_data(*args, **kwargs)
    context = self.copy_session_to_context(context)
    return context
  
  @method_decorator(login_required(login_url=reverse_lazy('auth_login'), redirect_field_name=None))
  def dispatch(self, *args, **kwargs):
    return super(MyDetailView, self).dispatch(*args, **kwargs)

class MyFormView(FormView, SessionContextMixin):

  def get_context_data(self, *args, **kwargs):
    '''Add session context items for navigation.'''
    context = super(MyFormView, self).get_context_data(*args, **kwargs)
    context = self.copy_session_to_context(context)
    return context
  
  @method_decorator(login_required(login_url=reverse_lazy('auth_login'), redirect_field_name=None))
  def dispatch(self, *args, **kwargs):
    return super(MyFormView, self).dispatch(*args, **kwargs)

# See http://stackoverflow.com/a/7012420 for the originator of this
# rather cool mixin idea.
class FilterMixin(object):

  def get_queryset_filters(self):
    '''Pull out allowed_filters terms into a dict.'''
    filters = {}
    for item in self.allowed_filters:
      if item in self.request.GET:
        filters[self.allowed_filters[item]] = self.request.GET[item]
    return filters

  def get_queryset(self):
    '''Apply any terms defined in the subclass allowed_filters attribute.'''

    if self.queryset is None:
      self.queryset = self.model.objects

    # Support a broad-based "query" term from a SimpleSearchForm; use
    # OR to combine queries across all allowed_filters.
    if 'query' in self.request.GET:
      qs = self.queryset.none()
      for term in self.allowed_filters.values():
        qs = qs | self.queryset.filter(**{ term:self.request.GET['query'] })
      return qs.distinct()
    else:

      # The alternative is to combine specific queries using AND. This
      # may be used in e.g. an advanced search interface.
      return super(FilterMixin, self).get_queryset().filter(**self.get_queryset_filters()).distinct()

# Note that the order of mixins/superclasses here is important.
class FormListView(MyListView, FormMixin):

  # See http://stackoverflow.com/a/9423819
  def get(self, request, *args, **kwargs):

    # From ProcessFormMixin
    form_class = self.get_form_class()
    self.form = self.get_form(form_class)

    # From BaseListView
    self.object_list = self.get_queryset()
    allow_empty = self.get_allow_empty()
    if not allow_empty and len(self.object_list) == 0:
      raise Http404(_(u"Empty list and '%(class_name)s.allow_empty' is False.")
                    % {'class_name': self.__class__.__name__})

    args = request.path_info

    # Sort out our paramstr variable which is used to construct pagination links.
    params = copy(self.request.GET)
    if 'page' in params:
      del(params['page'])
    paramstr = urlencode(params)
    if len(paramstr):
      paramstr += '&'  # Will be appended by "page=2" etc.

    context = self.get_context_data(object_list=self.object_list,
                                    form=self.form,
                                    paramstr=paramstr)
    return self.render_to_response(context)

  def post(self, request, *args, **kwargs):

    # From ProcessFormMixin
    form_class = self.get_form_class()
    self.form = self.get_form(form_class)

    self.success_url = request.path
    
    if self.form.is_valid():
      self.success_url += "?%s" % (urlencode(self.form.cleaned_data),)
      return self.form_valid(self.form)
    else:

      # Empty search returns a list of all objects.
      return redirect(self.success_url)

class RestrictedFileDownloadView(View):
  @method_decorator(login_required(login_url=reverse_lazy('auth_login'), redirect_field_name=None))
  def dispatch(self, *args, **kwargs):
    return super(RestrictedFileDownloadView, self).dispatch(*args, **kwargs)

