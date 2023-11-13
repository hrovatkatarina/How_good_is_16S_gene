# How good is the 16S rRNA gene marker for plant related microbes?



## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Figures](#figures)

## Installation

Clone the repository: 

```bash
git clone https://github.com/hrovatkatarina/
```

List of Dependencies:


## Usage

1. Download excel file from BV-BRC database, where in search you type the name of the genus of interest. Apply a filter on Genus column, where you select genomes with your genus name. Move excel file to Data directory.

2. Selection of quality complete and WGS genomes. As <your-excel-file.xslx> type the name of your excel file.

```bash
Rscript 1_Filtering.R <Data/your-excel-file.xslx>
```

Outputs: 
- Data/genome_list.txt 
- Data/genome_list_wgs.txt

3. Download fna files from BV-BRC database. Files will be downloaded in new fna directory.

```bash
bash 2_Get_genomes.sh
```
Output:
-Data/fna/*.fna

4. Predict 16S rRNA sequences with barrnap. Retrieve the longest 16S sequence predicted for each genome, rermove outliers - too long or too short sequences and write fasta file with 16S sequences of all genomes in genus. fna files of the genomes, where there was no 16S rRNA sequence found or sequence was too short or too long, will be moved to a fna_out directory. These genomes are excluded from further analysis.

```bash
bash 3_Predict_16S.sh
```
Output:
-gff directory with .gff files

```bash
python 4_16S_seq.py
```
Outputs:
- Data/missing_16S.txt (list of fna files with missing 16S sequences)
- Data/short_16S.txt (list of fna files with too short or too long 16S sequences)
- Data/fna_out (.fna files where 16S sequences couldn't be found or they were too long/short)
- barrnap_longest_16S.fasta (all 16S sequences combined in one file)

Note: if more genera are analyzed, combine all sequences in one fasta file.

5. Multiple sequence alignmnet using ssu-align. 

```bash
bash 5_ssu-align.sh
```
Output
- Data/myseqs/myseqs.bacteria.mask.afa (multiple sequence alignmnet og 16S rRNA sequences)

6. Entropy across 16S rRNA sequences.
Note: for here on MSA of our study is used. File name: full16S_MSAall.fasta

```bash
python 6_Entropy.py
```
Outputs:
- Data/{genus}_new.fasta (MSA for each genus in separate fasta file) - you don't need this files - DELETE!
- Data/Smooth_std_entropy.csv (Smoothen standardized entropy for all genera) - Use this for entropy plots
- Data/Std_entropy.csv (Standardized entropy for all genera)

7. Phylogenetic trees of 16S and its V-regions.

Get alignments of each regions of 16S and save them in separate fasta files for each genus.

```bash
Rscript 7_Regions.R
```
Outputs: 
- Data/genus_region.fasta

Calculate phylogenetic trees using iqtree.

```bash
bash 8_Iqtree.sh
```
Ouput:
- Data/Trees/*.treefile

8. Core genome phylogeny

Gene annotations with prokka

```bash
bash 9_Prokka.sh
```
Note: Run this for each genus separatly.
 
Ouput:
- Data/fna/gff/*.gff (annotated genes)

Core-genome alignment with Roary

```bash
bash 10_Roary.sh
```
Ouput:
- Data/fna/gff/roary-i75/core_gene_alignmnet.aln

Tree from core-genome alignment using iqtree

```bash
bash 11_Core_genome_tree.sh
```
Output:
- Data/fna/gff/roary-i75/core_gene_alignmnet.treefile  

9. ANI: Calculate ANI from fna files and Calculate bootstrap supported dendrogram of ANI matrix using pvclust.
Todo: Organize input data and modify scripts.

```bash
bash 12_ANIclustermap.txt
```
Output:
- Data/ANIclustermap_result/ANIclustermap_matrix.tsv

```bash
Rscript pvclust.R
```
Outputs:
- Data/ANIclustermap_result/{genus}/ANI_tree_BT.newick

Note: This script calculates dendrograms for all genera - each genus has its own directory with fan files in fna dir. 
Todo: Update!

Calculating ANI groups from ANI matrix.

```bash
python One_per_ANI.py
```
Output:
- cluster_table_95.tsv

Calculating ANI for each new  ANI group.

```bash
bash ANIclustermap.txt
```
10. Distance calculation.
Todo: MCI calculation is missing, distance caluclation for Figure 1f is not finished yet.

```bash
Rscript Distance_calculations.R
```
Outputs:
- Coregenome_distance.csv (ANI, 16S and V-regions vs core-genome tree. For figure 1b)
- JRF_distance_ANI.csv (Full 16S, V-regions for all genera vs. ANI. For figure 1e)
- Distance_ANI_groups.csv (16S and one per ANI group vs ANI (for genera with a lot of genomes): For figure 1f)

## Figures

Data and code for figures are in Figure1 directory.

