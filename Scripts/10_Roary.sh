#!/bin/bash

# Producing core-genome alignment with Roary.
# As input it takes annotated gff3 files produced by Prokka.

cd Data/fna/gff

#For genus Actinoplanes and Ensifer
roary -f ./roary-i75 -e -n -v -i 75 -p 80 *.gff 

#For genus Azospirillum and Cupriavidus
#roary -f ./roary-i75 -e -n -v -i 75 -g 70000 -p 80 *.gff 

#For genus Massilia
#roary -f ./roary-i75 -e -n -v -i 75 -g 58000 -p 80 *.gff