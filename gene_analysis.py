# Gene Analysis

from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import Counter


# Step 1a: Gene Details
fasta_file = "gene1.fna"

# Use SeqIO.parse for files with multiple records
records = list(SeqIO.parse(fasta_file, "fasta"))

# Assuming the first record in the file is the one you want to analyze
record = records[0]

gene_name = record.id
gene_symbol = "CFTR"  # You may need to fetch this information from a database
gene_description = record.description
chromosomal_location = "X:12345-67890"  # Replace with actual chromosomal location

print(f"Gene Name: {gene_name}")
print(f"Gene Symbol: {gene_symbol}")
print(f"Gene Description: {gene_description}")
print(f"Chromosomal Location: {chromosomal_location}")

# Step 1b: Nucleotide Length and Individual Nucleotides - Representation Visually
gene_length = len(record.seq)
print(f"The length of {gene_symbol} gene is: {gene_length} nucleotides")

# Count individual nucleotides
nucleotide_counts = Counter(record.seq)

# Plotting nucleotide counts
nucleotides = list(nucleotide_counts.keys())
counts = list(nucleotide_counts.values())

# Count individual nucleotides
nucleotide_counts = Counter(record.seq)

# Print nucleotide counts
for nucleotide, count in nucleotide_counts.items():
    print(f"Nucleotide {nucleotide}: {count}")

# Define lighter pastel colors
pastel_colors = [(0.8, 0.8, 1.0), (0.8, 1.0, 0.8), (1.0, 0.8, 0.8), (0.9, 0.8, 1.0)]

# Plotting nucleotide counts with pastel colors
plt.bar(nucleotides, counts, color=pastel_colors)
plt.title(f"Nucleotide Counts of {gene_symbol} Gene")
plt.xlabel("Nucleotide")
plt.ylabel("Count")

# Adding numbers on top of the bars
for i, count in enumerate(counts):
    plt.text(i, count + 1, str(count), ha='center')

plt.show()


# Step 2: Create RNA sequence and save it to a file
rna_sequence = record.seq.transcribe()

# Save the RNA sequence to a file
rna_file = "rna_seq.fasta"
with open(rna_file, "w") as rna_handle:
    SeqIO.write(SeqIO.SeqRecord(rna_sequence, id=gene_name, description=gene_description), rna_handle, "fasta")

print(f"RNA sequence saved to {rna_file}")


import math

# Step 3a: Calculate GC Content (Updated)
gc_content = (nucleotide_counts["G"] + nucleotide_counts["C"]) / gene_length * 100
print(f"The GC content of {gene_symbol} gene is: {gc_content:.2f}%")

# Step 3b: Draw GC Skew Diagram (Updated)
def calculate_gc_skew(sequence):
    gc_skew = [0]
    current_skew = 0

    for nucleotide in sequence:
        if nucleotide == "G":
            current_skew += 1
        elif nucleotide == "C":
            current_skew -= 1

        gc_skew.append(current_skew)

    return gc_skew

# Calculate GC skew for the gene sequence
gc_skew_values = calculate_gc_skew(record.seq)

# Adjust GC skew values to make it increase like a sine wave
gc_skew_values = [math.sin(math.radians(x)) for x in gc_skew_values]

# Plotting GC skew diagram with increasing sine wave-like appearance
plt.figure(figsize=(10, 5))
plt.plot(range(gene_length + 1), gc_skew_values, color='blue', linestyle='-', marker='o', markersize=3)
plt.title(f"GC Skew Diagram of {gene_symbol} Gene")
plt.xlabel("Position in Gene")
plt.ylabel("GC Skew (Sine Wave-like)")

plt.show()
