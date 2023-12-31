from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Blast import NCBIWWW

# Load your gene DNA sequence
gene_sequence = SeqIO.read("gene1.fna", "fasta")

# Break down the sequence into smaller fragments
fragment_size = 500  # Adjust the size as needed
gene_fragments = [gene_sequence.seq[i:i+fragment_size] for i in range(0, len(gene_sequence.seq), fragment_size)]

# Perform BLAST search for each fragment
for i, fragment in enumerate(gene_fragments):
    result_handle = NCBIWWW.qblast("blastn", "nt", fragment)
    blast_results = result_handle.read()
    with open(f"blast_results_fragment_{i+1}.xml", "w") as result_file:
        result_file.write(blast_results)
    result_handle.close()
