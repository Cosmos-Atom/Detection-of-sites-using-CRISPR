import random

def read_fasta(file_path):
    """
    Read a FASTA file and return a list of sequences.
    """
    records = [record.seq for record in SeqIO.parse(file_path, "fasta")]
    return records

def get_random_motifs(dna, k):
    """
    Randomly select k-mers as motifs from each sequence in the DNA.
    """
    motifs = [dna[i:i+k] for i in range(len(dna) - k + 1)]
    return motifs

def score_motifs(motifs):
    """
    Calculate the score of a set of motifs based on consensus.
    """
    consensus = ""
    for i in range(len(motifs[0])):
        column = [motif[i] for motif in motifs]
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for base in column:
            counts[base] += 1
        consensus += max(counts, key=counts.get)
    score = sum([1 for i in range(len(consensus)) if consensus[i] != motifs[0][i]])
    return score

def randomized_motif_search(dna, k, t, max_iter=1000):
    """
    Perform Randomized Motif Search to find motifs in DNA sequences.
    """
    best_motifs = get_random_motifs(random.choice(dna), k)
    best_score = score_motifs(best_motifs)

    for _ in range(max_iter):
        motifs = get_random_motifs(random.choice(dna), k)
        while True:
            profile = create_profile(motifs)
            motifs = [profile_most_probable(dna_seq, k, profile) for dna_seq in dna]
            current_score = score_motifs(motifs)
            if current_score < best_score:
                best_motifs = motifs
                best_score = current_score
            else:
                break

    return best_motifs

def create_profile(motifs):
    """
    Create a profile matrix from a set of motifs.
    """
    profile = {base: [0] * len(motifs[0]) for base in 'ACGT'}
    for i in range(len(motifs[0])):
        column = [motif[i] for motif in motifs]
        for base in 'ACGT':
            profile[base][i] = (column.count(base) + 1) / (len(column) + 4)
    return profile

def profile_most_probable(dna_seq, k, profile):
    """
    Find the most probable k-mer in a DNA sequence based on a profile matrix.
    """
    max_prob = -1
    most_probable = ""
    for i in range(len(dna_seq) - k + 1):
        kmer = dna_seq[i:i+k]
        prob = 1
        for j in range(k):
            prob *= profile[kmer[j]][j]
        if prob > max_prob:
            max_prob = prob
            most_probable = kmer
    return most_probable

if __name__ == "__main__":
    from Bio import SeqIO

    # File path
    sequence_file = "gene1.fna"

    # Read sequences
    sequences = read_fasta(sequence_file)

    # Set parameters
    k = 8  # length of motif
    t = len(sequences)  # number of sequences

    # Perform Randomized Motif Search
    best_motifs = randomized_motif_search(sequences, k, t)

    # Print the best motifs
    print("Best motifs found:")
    for motif in best_motifs:
        print(motif)
