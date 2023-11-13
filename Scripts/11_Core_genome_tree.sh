#!/bin/bash

cd Data/fna/gff/roary-i75
iqtree -s core_gene_alignment.aln -B 1000 -T AUTO -ntmax 20