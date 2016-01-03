'''
Serializer classes used by the REST API.
'''

from rest_framework import serializers
from .models import Project, Library, Source, Sample, Lane, Lanefile

class SourceSerializer(serializers.ModelSerializer):
  '''
  General Source-related metadata.
  '''
  sex     = serializers.StringRelatedField()
  strain  = serializers.StringRelatedField()
  species = serializers.StringRelatedField()
  class Meta:
    model     = Source
    fields    = ('name', 'sex', 'strain', 'species')

class SampleSerializer(serializers.ModelSerializer):
  '''
  Samples bridge Library to Source.
  '''
  source  = SourceSerializer(read_only=True)
  tissue  = serializers.StringRelatedField()
  class Meta:
    model     = Sample
    fields    = ('name', 'source', 'tissue')

class ProjectSerializer(serializers.ModelSerializer):
  '''
  The top-level Project serializer.
  '''
  libraries = serializers.HyperlinkedRelatedField(many=True,
                                                  view_name='api:library-detail',
                                                  read_only=True)
  class Meta:
    model     = Project
    fields    = ('code', 'description', 'libraries')

class LibrarySerializer(serializers.ModelSerializer):
  '''
  Library serialization.
  '''
  libtype  = serializers.StringRelatedField()
  factor   = serializers.StringRelatedField()
  sample   = SampleSerializer(read_only=True)
  lane_set = serializers.HyperlinkedRelatedField(many=True,
                                                 view_name='api:lane-detail',
                                                 read_only=True)
  
  class Meta:
    model     = Library
    fields    = ('code', 'libtype', 'factor', 'sample', 'lane_set')

class LanefileSerializer(serializers.ModelSerializer):
  '''
  This Lanefile serializer links back to the main UI file download view.
  '''
  filetype = serializers.StringRelatedField()
  download = serializers.HyperlinkedIdentityField(view_name='api:lanefile-download')
  
  class Meta:
    model  = Lanefile
    fields = ('filename','checksum','filetype','download')

class LaneSerializer(serializers.ModelSerializer):
  '''
  A bridging serializer between Library and Lanefile.
  '''
  machine  = serializers.StringRelatedField()
  facility = serializers.StringRelatedField()
  lanefile_set = LanefileSerializer(read_only=True, many=True)
  
  class Meta:
    model = Lane
    fields = ('machine', 'facility', 'flowcell', 'flowlane', 'lanefile_set')

