#!/bin/bash

cd Data/fna

for file in *.fna
do
name=$(basename "$file" .fna)
prokka --kingdom Bacteria --outdir "$name" --prefix "$name" --locustag "$name" "$file" --cpus 20
done
mkdir gff
find -name "*.gff" -exec mv -t gff {} +