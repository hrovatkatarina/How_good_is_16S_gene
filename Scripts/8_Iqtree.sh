#!/bin/bash

cd Data

for file in *.fasta
do
iqtree -s $file -B 1000 -T AUTO -ntmax 15
done

mkdir Trees
mv *.treefile Trees