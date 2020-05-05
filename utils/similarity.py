
# Function to get residues close to the mutation site
def slice_sequence(true_i, seq_original, seq_comparison, L=5):
    lower_bound = max(true_i-L,0)
    upper_bound = min(true_i+L+1, len(seq_original))
    return seq_original[lower_bound:upper_bound], seq_comparison[lower_bound:upper_bound]

# Function to compare two sequences by proportion sequence identity (higher is more similar)
def sequence_identity(seq_0, seq_1):
    assert len(seq_0) == len(seq_1)
    identities = 0
    for i, _ in enumerate(seq_0):
        if seq_0[i] == seq_1[i]:
            identities += 1
    return identities/len(seq_0)

# Function to parse similarity matrix into python-friendly object
def parse_similarity_matrix(similarity_matrix_copied):
    return [[float(x) for x in line.split(" ")] for line in similarity_matrix_copied.split("\n")]

# Function to determine residue similarity (lower is more similar)
def residue_similarity(r_0, r_1, similarity_matrix, residue_lookups): 
    indices = [residue_lookups.index(r) for r in [r_0, r_1]]
    indices = sorted(indices)
    indices.reverse()
    lookup = similarity_matrix[indices[0]][indices[1]]
    return lookup

# Function to determine average residue similarity in an area (lower is more similar)
def sequence_similarity_score(seq_0, seq_1, similarity_matrix, residue_lookups): 
    assert len(seq_0) == len(seq_1)
    identities = []
    for i, _ in enumerate(seq_0):
        identities.append(residue_similarity(seq_0[i], seq_1[i], similarity_matrix, residue_lookups))
    return sum(identities)/len(seq_0)

