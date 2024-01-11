from Bio import SeqIO
import random

# Randomly select sequence from MSA
# For Bacillus, Streptomyces, Rhizobium, Pseudomonas, Bradyrhizobium

def select_one_per_ANI_group(genus, iteration):
    """
    This function randomly selects one genome from each ANI group,
    repeats this process 10 times
    and writes the 16S sequences of selected genomes to a fasta file.
    """

    dir_ANI_table = f"C:Data/{genus}/ANIclustermap_result_new/cluster_table_95.tsv"

    # Read the tab file and store genome IDs and group numbers in a dictionary
    group_dict = {}
    with open(dir_ANI_table) as tab_file:
        next(tab_file)
        for line in tab_file:
            genome_id, group_number = line.strip().split("\t")
            group_dict[genome_id] = group_number

    # Create a list of unique group numbers
    group_numbers = list(set(group_dict.values()))

    # Randomly select one genome ID from each group
    selected_genome_ids = []
    for group_number in group_numbers:
        genomes_in_group = [genome_id for genome_id, group in group_dict.items() if group == group_number]
        random_genome_id = random.choice(genomes_in_group)
        selected_genome_ids.append(random_genome_id)

    # Iterate over the fasta files and filter records based on selected genome IDs
    fasta_file = f"Data/{genus}_full16S.fasta"

    filtered_records = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in selected_genome_ids:
            filtered_records.append(record)

    output_file = f"Data/{genus}/{genus}_one_per_ANI_group_{iteration}.fasta"
    SeqIO.write(filtered_records, output_file, "fasta")
    #print(output_file, "file created, number of seq=", len(filtered_records))

    filenames_output = f"Data/{genus}/{genus}_filenames_oneperANI_{iteration}.txt"
    with open(filenames_output, 'w') as outfile:
        for number in selected_genome_ids:
            outfile.write(str(number) + ".fna\n")
        #print(outfile, "created")


genera = ["Bacillus", "Pseudomonas", "Streptomyces", "Rhizobium", "Bradyrhizobium"]

for i in range(1, 11):
    for genus in genera:
        select_one_per_ANI_group(genus, i)
