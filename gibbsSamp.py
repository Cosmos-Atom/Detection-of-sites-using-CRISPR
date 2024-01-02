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

def profileRandom(k, profile, text):
    # probabilistically choose a position for a new motif in a given sequence.
    probs = []
    for i in range(0, len(text) - k + 1):
        prob = 1.0  # for each pos i calc the prblity of the k-mer motif starting at i
        pattern = text[i:i + k]
        for j in range(k):
            l = symbolToNumber(pattern[j])
            prob *= profile[l][j]  # probability is calculated by multiplying the probabilities of each nucleotide in the k-mer according to the given profile
        probs.append(prob)  # probs list contains the problities of each possible motif starting posititon in the text
    r = myRandom(probs)  # randomly samples an index r from the list of prblities
    return r

def profileForm(motifs):
    k = len(motifs[0])
    profile = [[1 for i in range(k)] for j in range(4)]  # Create a 2D list (profile) with 4 rows (for nucleotides A, C, G, and T) and k columns. Initialize all elements to 1.
    for x in motifs:
        for i in range(len(x)):
            j = symbolToNumber(x[i])
            profile[j][i] += 1
            # Count Nucleotide Occurrences in Motifs
    for x in profile:
        for i in range(len(x)):
            x[i] = x[i] / len(motifs)  # normalize the counts to probabilities by dividing each count by the total number of motifs
    return profile

def consensus(profile):
    # profile is a prblity profile matrix representing the likelihood of each pos in the motif
    result = ""
    for i in range(len(profile[0])):
        # iterate over each pos in the motif and for each position iterate over the rows of the profile matrix(A,T,C,G)
        max_val = 0  # for finding the nucleotide with the highest probability
        loc = 0
        for j in range(4):
            if profile[j][i] > max_val:
                loc = j  # Update loc to the index of the nucleotide with the maximum probability.
                max_val = profile[j][i]
        result += numberToSymbol(loc)  # appending the consensus nucleotide to the res
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

def myRandom(dist):
    # dist is a list of floats representing prblity distribution
    s = 0.0
    for x in dist:
        s += x
    i = random.random()  # random no. btw 0 and 1
    partial = 0.0  # var to store the cumulative prblity
    for x in range(len(dist)):
        partial += dist[x]
        if partial / s >= i:  # checking if the cumulative prblity normalised by the total prblity is greater than or eql to the random no i
            return x

def calculateConsensusFromMotifs(motifs):
    profile = profileForm(motifs)
    return consensus(profile)

def validateMotifsGibbs(bestMotifs, dna):
    profile = profileForm(bestMotifs)

    consensusFromMotifs = calculateConsensusFromMotifs(bestMotifs)
    consensusDirect = consensus(profile)

    print("Consensus Motif from Motifs:", consensusFromMotifs)
    print("Consensus Motif Directly:", consensusDirect)

    if consensusFromMotifs == consensusDirect:
        print("Validation: Motifs are consistent with the consensus.")
    else:
        print("Validation: Motifs are inconsistent with the consensus.")

def gibbsSampler(dna, k, t, n):
    bestMotifs = []
    bestPositions = []
    motifs = []
    for x in range(t):
        # for each seq in teh input DNA randomly select starting motifs
        i = random.randint(0, len(dna[x]) - k)
        motifs.append(dna[x][i:i + k])
    bestMotifs = motifs[:]
    for i in range(n):
        # randomly choose a seq index
        j = random.randint(0, t - 1)
        profile = profileForm(motifs[:j] + motifs[j + 1:])  # build a profile matrix excluding the jth seq
        r = profileRandom(k, profile, dna[j])
        motifs[j] = dna[j][r:r + k]  # using the profile sample a new motif for seq j
        if score(motifs) < score(bestMotifs):
            bestMotifs = motifs[:]  # update the best motif if the current motif has a lower score
    for i in range(t):
        pos = dna[i].find(bestMotifs[i])
        bestPositions.append(pos)  # return the best motifs and corresponding pos
    return bestMotifs, bestPositions

# Function to measure runtime
def measure_runtime(func, *args):
    start_time = time.time()
    result = func(*args)
    end_time = time.time()
    runtime = end_time - start_time
    return result[0], result[1], runtime

# Example usage
with open("gene3.fna") as file:
    k, t, n = [int(x) for x in file.readline().split()]
    dna = [line.rstrip() for line in file]

# Measure runtime of the gibbsSampler function
bestMotifs, bestPositions, runtime = measure_runtime(gibbsSampler, dna, k, t, n)

# Print runtime
print(f"Runtime: {runtime} seconds")

# Print motifs and positions
for motif, position in zip(bestMotifs, bestPositions):
    print(f"Motif: {motif}, Starting Position: {position}")

# Validate motifs using consensus motif
validateMotifsGibbs(bestMotifs, dna)
