#!/bin/bash

# Make phylogenetic trees for each genus separatly,
# using multiple sequence alignmnet of 16S rRNA sequence and its variable regions.

cd Data

for file in *.fasta
do
iqtree -s $file -B 1000 -T AUTO -ntmax 15
done

mkdir Trees
mv *.treefile Trees