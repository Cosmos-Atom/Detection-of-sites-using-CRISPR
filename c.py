from Bio import SeqIO
import random

# Step 1: Load and Process CFTR Gene Data
file_path = 'gene1.fna'
gene_sequence = str(SeqIO.read(file_path, 'fasta').seq)

# Step 2: Simulate CRISPR Model
def simulate_CRISPR(gene_sequence):
    pam_possibilities = ['AGG', 'CGG', 'TGG', 'GGG'] #NGG format
    all_target_sites = []

    for pam in pam_possibilities:
        target_sites = [i for i in range(len(gene_sequence) - len(pam) + 1) if gene_sequence[i:i+len(pam)] == pam]
        all_target_sites.append(target_sites)

    return all_target_sites

# Simulate CRISPR and find potential target sites for each possibility
potential_target_sites = simulate_CRISPR(gene_sequence)

# Step 3: CRISPR Model Logic
def crispr_model_logic(potential_target_sites):
    model_output = [len(target_sites) for target_sites in potential_target_sites]
    return model_output

# Simulate the CRISPR model based on the logic
crispr_output = crispr_model_logic(potential_target_sites)
print("CRISPR Model Output for Each Possibility:", crispr_output)
