import random
import time

def symbolToNumber(symbol):
    if symbol == "A":
        return 0
    elif symbol == "C":
        return 1
    elif symbol == "G":
        return 2
    elif symbol == "T":
        return 3

def numberToSymbol(x):
    if x == 0:
        return "A"
    elif x == 1:
        return "C"
    elif x == 2:
        return "G"
    elif x == 3:
        return "T"

# making probability profile matrix
def profileForm(motifs):
    k = len(motifs[0])
    profile = [[1 for i in range(k)] for j in range(4)]
    for x in motifs:
        for i in range(len(x)):
            j = symbolToNumber(x[i])
            profile[j][i] += 1
    for x in profile:
        for i in range(len(x)):
            x[i] = x[i] / len(motifs)  # normalize the counts by dividing each count by the total number of motifs
    return profile

def consensus(profile):
    result = ""
    for i in range(len(profile[0])):
        max_val = 0
        loc = 0
        for j in range(4):
            if profile[j][i] > max_val:  # choosing the nucleotide with the best probability in the profile matrix
                loc = j
                max_val = profile[j][i]
        result += numberToSymbol(loc)
    return result

def score(motifs):
    profile = profileForm(motifs)
    cons = consensus(profile)
    score = 0
    for x in motifs:
        for i in range(len(x)):
            if cons[i] != x[i]:
                score += 1
    return score

def greedyMotifSearch(dna, k, t):
    bestMotifs = [string[:k] for string in dna]  # Initialize bestMotifs with the first k-mers from each sequence
    bestPositions = [0] * t  # Initialize bestPositions with zeros

    for i in range(len(dna[0]) - k + 1):
        motifs = [dna[0][i:i+k]]  # Initialize motifs with the current k-mer from the first sequence
        positions = [i]  # Initialize positions with the starting position of the current k-mer

        for j in range(1, t):
            profile = profileForm(motifs)  # Create a profile matrix from the current motifs
            motif = profileMostProbable(dna[j], k, profile)  # Find the most probable k-mer in the j-th sequence based on the profile
            motifs.append(motif)  # Add the most probable k-mer to motifs
            positions.append(dna[j].find(motif))  # Record the starting position of the most probable k-mer in the j-th sequence

        if score(motifs) < score(bestMotifs):  # If the score of the current motifs is lower than the score of the best motifs
            bestMotifs = motifs  # Update bestMotifs with the current motifs
            bestPositions = positions  # Update bestPositions with the current starting positions

    return bestMotifs, bestPositions

def profileMostProbable(text, k, profile):
    max_prob = -1
    most_probable = text[:k]
    # to keep track of the curr max prob, and corres probable k-mer
    for i in range(len(text) - k + 1):
        # iterate over all possible k-mers
        pattern = text[i:i+k]
        prob = 1.0
        for j in range(k):
            l = symbolToNumber(pattern[j])
            prob *= profile[l][j]
            # prod of the probabilities of individual symbols at each pos in the k-mer
        if prob > max_prob:
            max_prob = prob
            most_probable = pattern
    return most_probable

# Calculate consensus motif directly from motifs
def calculateConsensusFromMotifs(motifs):
    profile = profileForm(motifs)
    return consensus(profile)

# Validate motifs using consensus motif
def validateMotifs(bestMotifs, profile):
    consensusFromMotifs = calculateConsensusFromMotifs(bestMotifs)
    consensusDirect = consensus(profile)

    print("Consensus Motif from Motifs:", consensusFromMotifs)
    print("Consensus Motif Directly:", consensusDirect)

    if consensusFromMotifs == consensusDirect:
        print("Validation: Motifs are consistent with the consensus.")
    else:
        print("Validation: Motifs are inconsistent with the consensus.")

# Function to measure runtime
def measure_runtime(func, *args):
    start_time = time.time()
    result = func(*args)
    end_time = time.time()
    runtime = end_time - start_time
    return result[0], result[1], runtime


# Example usage
with open("gene3.fna") as file:
    k, t, _ = [int(x) for x in file.readline().split()] 
    dna = [line.rstrip() for line in file]

# Measure runtime of the greedyMotifSearch function
bestMotifs, bestPositions, runtime = measure_runtime(greedyMotifSearch, dna, k, t)

# Print runtime
print(f"Runtime: {runtime} seconds")

# Validate motifs using consensus motif
validateMotifs(bestMotifs, profileForm(bestMotifs))

# Print motifs and positions
for motif, position in zip(bestMotifs, bestPositions):
    print(f"Motif: {motif}, Starting Position: {position}")
