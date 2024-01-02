#!/bin/bash

for i in `cat genome_list.txt`; do wget -qN "ftp://ftp.bvbrc.org/genomes/$i/$i.fna";
done

