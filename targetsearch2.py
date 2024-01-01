from Bio import SeqIO
from collections import Counter

# Step 1: Load and Process CFTR Gene Data
file_path = "E:\\Sharanya\\bioppr\\gene1.fna"
gene_sequence = SeqIO.read(file_path, 'fasta').seq

records = list(SeqIO.parse(file_path, "fasta"))
gene_symbol = "CFTR"
record = records[0]

gene_length = len(record.seq)
print(f"The length of {gene_symbol} gene is: {gene_length} nucleotides")

# Count individual nucleotides
nucleotide_counts = Counter(record.seq)

# Step 2: Simulate CRISPR Model
def simulate_CRISPR(gene_sequence):
    pam_possibilities = ['AGG', 'CGG', 'TGG', 'GGG']
    all_target_sites = []

    for pam in pam_possibilities:
        target_sites = [i for i in range(len(gene_sequence) - len(pam) + 1) if str(gene_sequence[i:i + len(pam)]) == pam]
        all_target_sites.append(target_sites)

    return all_target_sites

# Simulate CRISPR and find potential target sites for each possibility
potential_target_sites = simulate_CRISPR(gene_sequence)

# Step 3: CRISPR Model Logic
def crispr_model_logic(potential_target_sites):
    model_output = [(len(target_sites), target_sites) for target_sites in potential_target_sites]
    return model_output

# Simulate the CRISPR model based on the logic
crispr_output = crispr_model_logic(potential_target_sites)

# Display the count and locations of target sites
for i, (count, locations) in enumerate(crispr_output):
    for location in locations:
        nucleotide_countsbef = Counter(str(gene_sequence[location-20:location]))
        gc_contentbef = (nucleotide_countsbef["G"] + nucleotide_countsbef["C"]) / 20 * 100
        nucleotide_countsaft = Counter(str(gene_sequence[location+3:location+23]))
        gc_contentaft = (nucleotide_countsaft["G"] + nucleotide_countsaft["C"]) / 20 * 100
        if gc_contentbef < 40 or gc_contentaft < 40:
            print(f"PAM Sequence {i + 1}: Location = {location}")