#NOT NEEEDED
# This is for basic understating how it works

from Bio.Seq import Seq

def find_grna_pam(gene_sequence, pam_sequence="NGG"):
    # Convert the gene sequence to a Bio.Seq object
    gene_seq = Seq(gene_sequence)

    # Search for the PAM sequence and gRNA sequence in both directions
    for pam in [pam_sequence, gene_seq]:
        index = gene_seq.find(pam)
        if index != -1:
            grna_sequence = str(gene_seq[:index])
            return grna_sequence, pam

    # If no PAM sequence is found, return None
    return None, None

# Example usage
# Example usage
# Example usage with specific gRNA and PAM sequences
gene_sequence = "ATCGATCGATCGAGCTCGATCGACC"
grna_sequence_to_find = "ATCGAT"
pam_sequence_to_find = "ACC"

grna_sequence, pam_sequence = find_grna_pam(gene_sequence, pam_sequence_to_find)

if grna_sequence and pam_sequence:
    print("Original Gene Sequence:", gene_sequence)
    print("Found gRNA Sequence:", grna_sequence)
    print("Found PAM Sequence:", pam_sequence)
else:
    print("No gRNA and PAM sequence found.")

