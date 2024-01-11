import os
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as pylab
import numpy as np
import shutil

# This script retrieves the longest 16S sequence predicted with barrnap (in gff directory) for each genome, 
# rermove outliers - too long or too short sequences and 
# write fasta file with 16S sequences of all genomes in each genus (barrnap_longest_16S.fasta). 
# fna files of the genomes, where there was no 16S rRNA sequence found or sequence was too short or too long, 
# will be moved to a fna_out directory. These genomes are excluded from further analysis.

directory_path = "Data/gff"
missing_16s_fna = "Data/missing_16S.txt"

with open(missing_16s_fna, 'a'):
   pass

# Dictionary to store the longest gene from each file
longest_genes = {}

# Iterate over the files in the directory ending in ".gff"
for file_path in glob.glob(os.path.join(directory_path, "*.gff")):
    with open(file_path, "r") as handle:
        # Initialize variables to keep track of the longest 16S rRNA gene in the file
        longest_gene_length = 0
        longest_gene_header = ""
        longest_gene_seq = ""

        # Parse the file using Biopython
        for record in SeqIO.parse(handle, "fasta"):
            # Check if the description line contains the 16S rRNA gene (check the description)
            if "16S_rRNA" in record.description:
                # Check if the gene is longer than the current longest gene
                if len(record.seq) > longest_gene_length:
                    longest_gene_length = len(record.seq)
                    longest_gene_header = record.description
                    longest_gene_seq = record.seq

        # Check if 16S gene was found in the file
        if longest_gene_length > 0:
            # Add the longest gene from the file to the dictionary
            longest_genes[file_path] = {"header": longest_gene_header, "seq": longest_gene_seq}
        else:
            print("Couldn't find 16S gene:", file_path)
            genome_id = os.path.basename(file_path)
            genome_id_fna = genome_id.replace(".gff", ".fna")
            with open(missing_16s_fna, "a") as missing_file:
                missing_file.write(genome_id_fna + "\n")

# Initialize a list to store the SeqRecord objects
seq_records = []

# Create a new SeqRecord object for each longest gene and add it to the list
for file_path, gene in longest_genes.items():
    file_name_with_ext = os.path.basename(file_path)
    file_name, file_ext = os.path.splitext(file_name_with_ext)
    seq_record = SeqRecord(gene["seq"], id=file_name, description=gene["header"]) #With SeqRecord return an object whereas ID = file path, Decription: header, Seq: sequence
    seq_records.append(seq_record)

# lengths of 16s genes in list
lengths = []
for record in seq_records:
    seq_length = len(record.seq)
    lengths.append(seq_length)

print(len(lengths))
set_order = (set(lengths))
print(sorted(set_order))

pylab.boxplot(lengths)
pylab.ylabel("Length (bp)")
pylab.title("Length distribution of 16S gene sequences")
pylab.show()

# write new fasta file - set the threshold based on the one from boxplot
output_file = "Data/barrnap_longest_16S.fasta"
short_16s_fna = "Data/short_16s.txt"

with open(short_16s_fna, 'a'):
   pass

# Calculate the IQR of the lengths
Q1 = np.percentile(lengths, 25)
Q3 = np.percentile(lengths, 75)
IQR = Q3 - Q1

# Define the lower and upper thresholds
lower_threshold = (Q1 - 1.5 * IQR) - 10
print("Lower threshold: ", lower_threshold)
upper_threshold = (Q3 + 1.5 * IQR) + 10
print("Upper threshold: ", upper_threshold)

with open(output_file, "w") as out_file, open(short_16s_fna, "w") as new_file:
    for record in seq_records:
        if lower_threshold <= len(record.seq) <= upper_threshold:
            SeqIO.write(record, out_file, "fasta")
        else:
            print("Sequence was shorter than threshold:", record.id, len(record.seq))
            genome_id = os.path.basename(record.id)
            genome_id_fna = genome_id + ".fna"
            new_file.write(genome_id_fna + "\n")

# Define the source and destination directories
source_dir = "Data/fna"
dest_dir = "Data/fna_out"
if not os.path.exists(dest_dir):
    os.mkdir(dest_dir)

# Open the files and read the contents
with open("Data/missing_16s.txt", "r") as missing_file, open("Data/short_16s.txt", "r") as short_file:
    missing_files = missing_file.readlines()
    short_files = short_file.readlines()

# Move the files
for file in missing_files + short_files:
    file = file.strip()
    source_path = os.path.join(source_dir, file)
    dest_path = os.path.join(dest_dir, file)
    shutil.move(source_path, dest_path)
