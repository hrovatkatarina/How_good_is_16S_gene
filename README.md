# Taxonomic resolution of different 16S rRNA variable regions varies strongly across plant-associated bacteria

## 🖋️ Citation

If you use this resource, please cite:

 > Hrovat, K., Dutilh, B. E., Medema, M. H., & Melkonian, C. (2024). Taxonomic resolution of different 16S rRNA variable regions varies strongly across plant-associated bacteria. In ISME Communications. Oxford University Press (OUP). https://doi.org/10.1093/ismeco/ycae034

## 🌿🦠 Abstract

Plant-microbiome research plays a pivotal role in understanding the relationships between plants and their associated microbial communities, with implications for agriculture and ecosystem dynamics. Metabarcoding analysis on variable regions of the 16S ribosomal RNA (rRNA) gene remains the dominant technology to study microbiome diversity in this field. However, the choice of the targeted variable region might affect the outcome of the microbiome studies. In our in-silico analysis, we have evaluated whether the targeted variable region has an impact on taxonomic resolution in 16 plant-related microbial genera. Through a comparison of 16S rRNA gene variable regions with whole-genome data, our findings suggest that the V1-V3 region is generally a more suitable option than the widely used V3-V4 region for targeting microbiome analysis in plant-related genera. However, sole reliance on one region could introduce detection biases for specific genera. Thus, we are suggesting that while transitioning to full-length 16S rRNA gene and whole-genome sequencing for plant-microbiome analysis, the usage of genus-specific variable regions can achieve more precise taxonomic assignments. More broadly, our approach provides a blueprint to identify the most discriminating variable regions of the 16S rRNA gene for genus of interest.


## 🧬 Data availability

All genome IDs of genomes used in analysis are provided in [Genomes](./Genomes) and can be downloaded for each genus separatly from BV-BRC database using command:

```bash
for i in `cat genus.txt`; do wget -qN "ftp://ftp.bvbrc.org/genomes/$i/$i.fna";
done
```

## 🔎 Analysis

The code and data necessary to reproduce the analysis and generate each panel in Figure 1 (a-f) can be accessed in the Figure1 directory.

Clone repo:
```
$ git clone git@github.com:hrovatkatarina/How_good_is_16S_gene.git
```
See repo structure:
```
$ tree -L 3
.
├── CITATION.bib
├── Figure1
│   ├── Figure1a
│   │   ├── full16S_MSAall.fasta
│   │   ├── PCA_full16S.png
│   │   └── PCA_full16S.R
│   ├── Figure1b
│   │   ├── Boxplot_1b.png
│   │   ├── Boxplot_1b.R
│   │   ├── Coregenome_distance.csv
│   │   └── Distance_calculation_core.R
│   ├── Figure1c
│   │   ├── 6_Entropy.py
│   │   ├── Cupriavidus_Entropy_plot.png
│   │   ├── Data
│   │   ├── entropy_newgenomes.csv
│   │   ├── Entropy_plot.R
│   │   ├── full16S_MSAall.fasta
│   │   ├── Rhizobium_Entropy_plot.png
│   │   ├── Smooth_std_entropy.csv
│   │   ├── Std_entropy.csv
│   │   ├── Supp_fig1.png
│   │   └── Supplementary_figure1.R
│   ├── Figure1d
│   │   ├── ANIclustermap_matrix.tsv
│   │   ├── ANI_tree_BT.newick
│   │   ├── Comparing_trees_1d.R
│   │   ├── Cupriavidus_MSA_V1V3_bt.fasta.treefile
│   │   ├── Cupriavidus_MSA_V6V8_bt.fasta.treefile
│   │   ├── CupriavidusV1V3_1d.png
│   │   ├── CupriavidusV6V8_1d.png
│   │   └── Rplots.pdf
│   ├── Figure1e
│   │   ├── Distance_calculation_heatmap.R
│   │   ├── Heatmap_1e.png
│   │   ├── Heatmap_1e.R
│   │   └── JRF_distance_ANI.csv
│   └── Figure1f
│       ├── Boxplot_1f.png
│       ├── Boxplot_1f.R
│       ├── Distance_ANI_groups.csv
│       ├── JRF_distance_ANI.csv
│       └── Pearson_correlation.png
├── Genomes
│   ├── 2_Get_genomes.sh
│   ├── Actinoplanes.txt
│   ├── Agrobacterium.txt
│   ├── Aminobacter.txt
│   ├── Azospirillum.txt
│   ├── Azotobacter.txt
│   ├── Bacillus.txt
│   ├── Bradyrhizobium.txt
│   ├── Burkholderia.txt
│   ├── Citrobacter.txt
│   ├── Cupriavidus.txt
│   ├── Duganella.txt
│   ├── Ensifer.txt
│   ├── Enterobacter.txt
│   ├── Kitasatospora.txt
│   ├── Massilia.txt
│   ├── Mesorhizobium.txt
│   ├── Micromonospora.txt
│   ├── Nitrospirillum.txt
│   ├── Paenibacillus.txt
│   ├── Paraburkholderia.txt
│   ├── Pseudomonas.txt
│   ├── Pseudoxanthomonas.txt
│   ├── Rahnella.txt
│   ├── Ralstonia.txt
│   ├── Rhizobium.txt
│   ├── Rhodopseudomonas.txt
│   ├── Serratia.txt
│   ├── Sinorhizobium.txt
│   ├── Streptomyces.txt
│   ├── Xanthomonas.txt
│   └── Xylella.txt
├── README.md
└── Scripts
    ├── 10_Roary.sh
    ├── 11_Core_genome_tree.sh
    ├── 12_ANIclustermap.txt
    ├── 13_ANI_clusters.py
    ├── 14_One_per_ANI.py
    ├── 15_pvclust.R
    ├── 16_Distance_calculation.R
    ├── 1_Filtering.R
    ├── 2_Get_genomes.sh
    ├── 3_Predict_16S.sh
    ├── 4_16S_seq.py
    ├── 5_ssu-align.sh
    ├── 6_Entropy.py
    ├── 7_Regions.R
    ├── 8_Iqtree.sh
    ├── 9_Prokka.sh
    └── Data
        ├── Actinoplanes
        ├── Agrobacterium
        ├── Aminobacter
        ├── Azospirillum
        ├── Azotobacter
        ├── Bacillus
        ├── Bradyrhizobium
        ├── Burkholderia
        ├── Citrobacter
        ├── Coregenome_distance.csv
        ├── Cupriavidus
        ├── Duganella
        ├── Ensifer
        ├── Enterobacter
        ├── full16S_MSAall.fasta
        ├── JRF_disatnce_ANI.csv
        ├── Kitasatospora
        ├── Massilia
        ├── Mesorhizobium
        ├── Micromonospora
        ├── Nitrospirillum
        ├── Paenibacillus
        ├── Paraburkholderia
        ├── Pseudomonas
        ├── Pseudoxanthomonas
        ├── Rahnella
        ├── Ralstonia
        ├── Rhizobium
        ├── Rhodopseudomonas
        ├── Serratia
        ├── Sinorhizobium
        ├── Smooth_std_entropy.csv
        ├── Std_entropy.csv
        ├── Streptomyces
        ├── Xanthomonas
        └── Xylella

42 directories, 89 files

```
