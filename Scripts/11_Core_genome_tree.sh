#!/bin/bash

# Make core-genome tree using iqtree tool.
# As input it takes core-genome alignment produced by Roary.

cd Data/fna/gff/roary-i75
iqtree -s core_gene_alignment.aln -B 1000 -T AUTO -ntmax 20