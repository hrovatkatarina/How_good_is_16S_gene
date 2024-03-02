from Bio import SeqIO, AlignIO
import numpy as np
import pandas as pd

# Save MSA for each genus in separate fasta file
infile = "Data/full16S_MSAall.fasta"

genera_set = set()

for record in SeqIO.parse(infile, "fasta"):
    genus = record.description.split()[1]
    genera_set.add(genus)
    output_file = f"Data/{genus}_new.fasta"
    with open(output_file, "a") as f:
        SeqIO.write(record, f, "fasta")

entropy_matrix = np.empty((0, 1494))
entropy_matrix_s = np.empty((0, 1494))
genera = list(genera_set)
print(genera)

for genus in genera:
    #infile = "Cupriavidus_new.fasta"
    infile = f"Data/{genus}_new.fasta"
    alignment_infile = AlignIO.read(infile, 'fasta')

    # convert the alignment to a numpy array
    alignment = np.array([list(rec.seq) for rec in alignment_infile], dtype='S1')

    # Calculate the entropy
    entropy = []
    for i in range(alignment.shape[1]):
    # loops over the columns of the alignment, alignment_array.shape[1] gives the length of the alignment,
    # which is the number of columns in the array
        column = alignment[:, i]  # extracts the column
        unique, counts = np.unique(column, return_counts=True)
    # unique: unique elements in the column counts: the number of times each unique item appears
        freq = counts / len(column)
        entropy.append(-np.sum(freq * np.log2(freq)))  # sum entropy for each unique nc in the column

    # standardize the entropy values
    entropy_std = entropy / np.max(entropy)

    # apply a moving average filter to the entropy values
    window_size = 70
    entropy_smooth = np.convolve(entropy, np.ones(window_size) / window_size, mode='same')

    # standardize the smooth entropy values
    entropy_smooth_std = entropy_smooth / np.max(entropy_smooth)

    entropy_matrix = np.vstack((entropy_matrix, entropy_std))
    entropy_matrix_s = np.vstack((entropy_matrix_s, entropy_smooth_std))
#print(entropy_matrix)

# Save std entropy to csv
entropy_data = {}
for i in range(1, 1495):
    entropy_data[i] = entropy_matrix[:, i-1].tolist()

entropy_df = pd.DataFrame(entropy_data, index=genera)
entropy_df.to_csv('Data/Std_entropy.csv', index=True, header=False)
#print(entropy_df.head())

# Save smooth std entropy to csv
entropy_data_s = {}
for i in range(1, 1495):
    entropy_data_s[i] = entropy_matrix_s[:, i-1].tolist()

entropy_df_s = pd.DataFrame(entropy_data_s, index=genera)
entropy_df_s.to_csv('Data/Smooth_std_entropy.csv', index=True, header=False)
print("CSV created")



