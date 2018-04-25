# -*- coding: utf-8 -*-
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

from __future__ import unicode_literals

from django.db import models, migrations

def populate_species(apps, schema_editor):

    Species = apps.get_model('osqpipe', 'Species')

    speciesdata = """\
Ambystoma mexicanum,Axolotl,8296
Anolis carolinensis,Green anole,28377
Arabidopsis thaliana,Thale cress,3702
Bos taurus,Cow,9913
Caenorhabditis elegans,Nematode,6239
Callithrix jacchus,Marmoset,9483
Canis familiaris,Dog,9615
Cavia porcellus,Guinea Pig,10141
Chlorocebus sabaeus,Vervet,60711
Daubentonia madagascariensis,Aye-aye,31869
Equus caballus,Horse,9796
Felis catus,Cat,9685
Gallus gallus,Chicken,9031
Gorilla gorilla,Gorilla,9593
Heter glaber,Naked mole-rat,10181
Homo sapiens,Human,9606
Macaca mulatta,Rhesus macaque,9544
Microcebus murinus,Gray mouse lemur,30608
Monodelphis domestica,Opossum,13616
Mus caroli,Ryukyu mouse,10089
Mus musculus,Mouse,10090
Mus musculus AJ,AJ mouse,10090_AJ
Mus musculus castaneus,Southeastern Asian house mouse,10091
Mus musculus TC1,Tc1 Mouse,10090_TC1
Mus pahari,Shrew mouse,10093
Mus spretus,Algerian mouse,10096
Mustela putorius furo,Ferret,9669
Myotis lucifugus,Microbat,59463
Nomascus leucogenys,Gibbon,61853
Oryctolagus cuniculus,Rabbit,9986
Otolemur garnettii,Small-eared galago,30611
Ovis aries,Sheep,9940
Pan troglodytes,Chimpanzee,9598
Papio hamadryas,Baboon,9557
Pongo abelii,Sumatran orangutan,9601
Pongo pygmaeus,Bornean orangutan,9600
Rattus norvegicus,Rat,10116
Saccharomyces cerevisiae,Yeast,4932
Sarcophilus harrisii,Tasmanian devil,9305
Sorex araneus,Common shrew,42254
Suncus murinus,Asian house shrew,9378
Sus scrofa,Pig,9823
Tetraodon nigroviridis,Spotted green pufferfish,99883
Trachypithecus francoisi,Francois Langur,54180
Trichuris muris,Whipworm (Mouse Nematode Parasite),70415
Tupaia belangeri,Tree Shrew,37347
Tursiops truncatus,Bottlenosed dolphin,9739
Xenopus tropicalis,Western clawed frog,8364\
"""
    
    for row in speciesdata.split("\n"):
        (sci, com, tax) = [ x.strip() for x in row.split(',') ]
        Species.objects.get_or_create(scientific_name=sci,
                                      common_name=com,
                                      accession=tax)

def relink_genomes(apps, schema_editor):

    Species = apps.get_model('osqpipe', 'Species')
    Genome = apps.get_model('osqpipe', 'Genome')

    for gen in Genome.objects.all():
        sp = Species.objects.get(scientific_name=gen.scientific_name)
        gen.species = sp
        gen.save()

class Migration(migrations.Migration):

    dependencies = [
        ('osqpipe', '0002_auto_20151218_1338'),
    ]

    operations = [
        migrations.RunPython(populate_species),
        migrations.RunPython(relink_genomes)
    ]
