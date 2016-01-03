'''
Serializer classes used by the REST API.
'''

from rest_framework import serializers
from .models import Project, Library

class ProjectSerializer(serializers.ModelSerializer):
  libraries = serializers.HyperlinkedRelatedField(many=True,
                                                  view_name='api:library-detail',
                                                  read_only=True)
  class Meta:
    model     = Project
    fields    = ('code', 'description', 'libraries')

class LibrarySerializer(serializers.ModelSerializer):
  libtype = serializers.StringRelatedField()
  class Meta:
    model     = Library
    fields    = ('code','libtype')