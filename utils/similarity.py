
def slice_sequence(true_i, seq_original, seq_comparison, L=5):
    """
    Function to get residues near the mutation site. 
    Args:
        true_i (int): index of the mutation site on the original sequence, after correcting for spaces
        seq_original, seq_comparison (Biopython sequences): EcFtn and PfFtn sequences respectively 
    Returns:
        tuple (EcFtn, PfFtn) of Biopython sequences 
    """
    lower_bound = max(true_i-L,0)
    upper_bound = min(true_i+L+1, len(seq_original))
    return seq_original[lower_bound:upper_bound], seq_comparison[lower_bound:upper_bound]

def sequence_identity(seq_0, seq_1):
    """
    Function to compare two sequences by proportion sequence identity (higher is more similar)
    Args:
        seq_0, seq_1: Biopython sequences
    Returns:
        float: proportion sequence identity
    """
    assert len(seq_0) == len(seq_1)
    identities = 0
    for i, _ in enumerate(seq_0):
        if seq_0[i] == seq_1[i]:
            identities += 1
    return identities/len(seq_0)

def parse_similarity_matrix(similarity_matrix_copied):
    """
    Function to parse similarity matrix into python-friendly object. 
    Args:
        similarity_matrix_copied: str (read from similarity_matrix.txt)
    Returns:
        list of list of floats: parsed similarity matrix
    """
    return [[float(x) for x in line.split(" ")] for line in similarity_matrix_copied.split("\n")]

def residue_similarity(r_0, r_1, similarity_matrix): 
    """
    Function to determine similarity between two amino acids (lower is more similar)
    Args:
        r_0, r_1: str of length 1 (amino acid residues)
        similarity_matrix: returned from parse_similarity_matrix
    Returns:
        float: residue similarity score
    """
    residue_lookups = "A R N D C Q E G H I K L M P F S T W Y V -".split(" ")

    indices = [residue_lookups.index(r) for r in [r_0, r_1]]
    indices = sorted(indices)
    indices.reverse()
    lookup = similarity_matrix[indices[0]][indices[1]]
    return lookup

def sequence_similarity_score(seq_0, seq_1, similarity_matrix): 
    """
    Function to determine mean residue similarity in an area (lower is more similar)
    Args:
        seq_0, seq_1: Biopython sequence
        similarity_matrix: returned from parse_similarity_matrix
    Returns:
        float: mean residue similarity score 
    """
    assert len(seq_0) == len(seq_1)
    identities = []
    for i, _ in enumerate(seq_0):
        identities.append(residue_similarity(seq_0[i], seq_1[i], similarity_matrix))
    return sum(identities)/len(seq_0)

