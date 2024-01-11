#!/bin/bash

# Download quality complete and WGS genomes from BV-BRC database and move fna files to fna directory

cd Data

for i in `cat genome_list.txt`; do wget -qN "ftp://ftp.bvbrc.org/genomes/$i/$i.fna";
done

for i in `cat genome_list_wgs.txt`; do wget -qN "ftp://ftp.bvbrc.org/genomes/$i/$i.fna";
done

mkdir fna
mv *.fna fna
