#!/bin/bash
cd Data/fna
for file in *.fna; do barrnap "$file" > "${file%.fna}.gff" --outseq "${file%.fna}.gff"; done
cd ..
mkdir gff
mv fna/*.gff gff
rm fna/*.fna.fai
cd gff
for file in *.gff; do sed -i "s/^>/>${file%.*} /" "$file"; done 
#this command appends genome_id to each 16S sequence