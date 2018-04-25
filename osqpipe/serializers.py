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

'''
Serializer classes used by the REST API.
'''

from rest_framework import serializers
from .models import Project, Library, Source, Sample, Lane, Lanefile, \
  Alignment, Alnfile, MergedAlignment, MergedAlnfile

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
  url       = serializers.HyperlinkedIdentityField(view_name='api:project-detail')

  class Meta:
    model     = Project
    fields    = ('code', 'url', 'description', 'libraries')

class LibrarySerializer(serializers.ModelSerializer):
  '''
  Library serialization.
  '''
  libtype  = serializers.StringRelatedField()
  factor   = serializers.StringRelatedField()
  sample   = SampleSerializer(read_only=True)
  adapter  = serializers.StringRelatedField()
  adapter2 = serializers.StringRelatedField()
  lane_set = serializers.HyperlinkedRelatedField(many=True,
                                                 view_name='api:lane-detail',
                                                 read_only=True)
  
  class Meta:
    model     = Library
    fields    = ('code', 'libtype', 'factor', 'adapter', 'adapter2',
                 'sample', 'lane_set')

class LanefileSerializer(serializers.ModelSerializer):
  '''
  This Lanefile serializer links back to the main UI file download view.
  '''
  filetype = serializers.StringRelatedField()
  download = serializers.HyperlinkedIdentityField(view_name='api:lanefile-download')
  
  class Meta:
    model  = Lanefile
    fields = ('filename_on_disk','checksum','filetype','download')

class LaneSerializer(serializers.ModelSerializer):
  '''
  A bridging serializer between Library and Lanefile.
  '''
  library  = serializers.StringRelatedField()
  machine  = serializers.StringRelatedField()
  facility = serializers.StringRelatedField()
  status   = serializers.StringRelatedField()
  lanefile_set = LanefileSerializer(read_only=True, many=True)
  alignment_set = serializers.HyperlinkedRelatedField(many=True,
                                                      view_name='api:alignment-detail',
                                                      read_only=True)
  
  class Meta:
    model = Lane
    fields = ('library', 'machine', 'facility', 'flowcell', 'flowlane',
              'rundate', 'genomicssampleid', 'paired', 'readlength',
              'total_passedpf', 'lanefile_set', 'alignment_set', 'status')

class MergedAlnfileSerializer(serializers.ModelSerializer):
  '''
  This MergedAlnfile serializer links back to the main UI file download view.
  '''
  filetype = serializers.StringRelatedField()
  download = serializers.HyperlinkedIdentityField(view_name='api:mergedalnfile-download')

  class Meta:
    model  = MergedAlnfile
    fields = ('filename_on_disk','checksum','filetype','download')

class MergedAlignmentSerializer(serializers.ModelSerializer):
  '''
  A bridging serializer between Sample and MergedAlnfile.
  '''
  genome   = serializers.StringRelatedField()
  mergedalnfile_set = MergedAlnfileSerializer(read_only=True, many=True)

  class Meta:
    model = MergedAlignment
    fields = ('genome','mergedalnfile_set',)

class AlnfileSerializer(serializers.ModelSerializer):
  '''
  This Alnfile serializer links back to the main UI file download view.
  '''
  filetype = serializers.StringRelatedField()
  download = serializers.HyperlinkedIdentityField(view_name='api:alnfile-download')
  
  class Meta:
    model  = Alnfile
    fields = ('filename_on_disk','checksum','filetype','download')

class AlignmentSerializer(serializers.ModelSerializer):
  '''
  A bridging serializer between Lane and Alnfile.
  '''
  genome   = serializers.StringRelatedField()
  alnfile_set = AlnfileSerializer(read_only=True, many=True)
  mergedalignment_set = serializers.HyperlinkedRelatedField(many=True,
                                                            view_name='api:mergedalignment-detail',
                                                            read_only=True)
  
  class Meta:
    model = Alignment
    fields = ('genome', 'total_reads', 'mapped', 'munique', 'headtrim', 'tailtrim',
              'alnfile_set', 'mergedalignment_set')

