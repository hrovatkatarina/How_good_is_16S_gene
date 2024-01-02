# Taxonomic resolution of different 16S rRNA variable regions varies strongly across plant-associated bacteria



## Description

This repository contains code used to perform the analysis reported in Hrovat et al. 'Taxonomic resolution of different 16S rRNA variable regions varies strongly across plant-associated bacteria'.


## Data availability

All genome IDs of genomes used in analysis are provided in [Genomes](./Genomes) and can be downloaded for each genus separatly from BV-BRC database using command:

```bash
for i in `cat genus.txt`; do wget -qN "ftp://ftp.bvbrc.org/genomes/$i/$i.fna";
done
```

## Analysis

Code and data to reproduce Figure 1 a-f is in [Figure1](./Figure1)
