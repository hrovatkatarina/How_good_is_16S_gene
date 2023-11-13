import pandas as pd
import os

# list of directories
genera = ["Actinoplanes", "Azospirillum", "Bacillus", "Bradyrhizobium", "Burkholderia", "Cupriavidus", "Ensifer",
          "Enterobacter", "Massilia", "Mesorhizobium", "Pseudomonas", "Rhizobium", "Serratia", "Streptomyces",
          "Xanthomonas", "Xylella"]

for genus in genera:
    directory = f"{genus}\\ANIclustermap_result"

    matrix_tsv_file = os.path.join(directory, "ANIclustermap_matrix.tsv") # ANIclustermap_matrix.tsv
    cluster_ani_thr = 95.0
    cluster_tsv_file = os.path.join(directory, "cluster_table_95.tsv")  # Output cluster table file

    # Parse cluster ANI matrix
    df = pd.read_table(matrix_tsv_file)
    cluster_id = 1
    cluster_base_idx = 0
    cluster_size_record = 1
    genome_name2cluster_id = {}
    for i, genome_name in enumerate(df.columns):
        cluster_candidate_df = df.iloc[cluster_base_idx : i + 1, cluster_base_idx : i + 1]
        ani_thr_match_count = (cluster_candidate_df > cluster_ani_thr).sum().sum()
        if ani_thr_match_count != cluster_size_record**2:
            cluster_id += 1
            cluster_base_idx = i
            cluster_size_record = 1

        genome_name2cluster_id[genome_name] = cluster_id
        cluster_size_record += 1

    # Output cluster table
    cluster_table_dict = {
        "genome": genome_name2cluster_id.keys(),
        "cluster_id": genome_name2cluster_id.values(),
    }
    cluster_table_df = pd.DataFrame(cluster_table_dict)
    #print("Table created", genus)
    cluster_table_df.to_csv(cluster_tsv_file, sep="\t", index=False)