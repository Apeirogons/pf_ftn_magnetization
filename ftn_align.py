# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# We have a list of mutations done in E. coli ferritin done to increase magnetizability, which we want to translate to P. furiosus.
# Here's the mutations: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5126674/figure/f2/ 

# %%
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from os import listdir
from utils.sequence import remove_duplicates, padding, seq_print
from utils.similarity import parse_similarity_matrix, slice_sequence, sequence_identity, sequence_similarity_score
from utils.utils import parse_other_known_sites, parse_mutations
import numpy as np


# %%
# reading sequences. [0] allows us to get each sequence in the Seq format as opposed to the SingleLetterAlphabet format. 
fastas = listdir("manual_sequences")
sequences = [AlignIO.read("manual_sequences/"+pdb, "fasta")[0] for pdb in fastas] 
unique_sequences = remove_duplicates(sequences)
unique_sequences = padding(unique_sequences)

print(" ")
for seq in unique_sequences: 
    seq_print(seq)


# %%
#### It's also interesting to note that 1EUM and 4XGS have really similar sequences.
seq_names = [seq.id[:4] for seq in sequences]
int_seq_0 = sequences[seq_names.index("1EUM")]
int_seq_1 = sequences[seq_names.index("4XGS")]
for i, residue in enumerate(int_seq_0):
    if int_seq_0[i] != int_seq_1[i]:
        print("Difference at point " + residue + str(i+1)) 
print('')
seq_print(int_seq_0)
seq_print(int_seq_1)


# %%
# While the E. coli paper doesn't give the full sequence of the ferritin used, it does mention what the mutations were at certain points. There seems to be only one sequence which matches all the sites given (therefore, this is probably the sequence the researchers used). This is the sequence we need to compare to the P. furiosus sequence.

with open("sites/mutations.txt", "r") as m:
    mutations = parse_mutations(m.read())

with open("sites/other_known_sites.txt", "r") as m:
    other_known_sites = parse_other_known_sites(m.read())

all_known_sites = other_known_sites
all_known_mutations = []
for mutation_list in mutations:
    for mutation in mutation_list:
        site = mutation[:-1]
        if site not in all_known_sites:
            all_known_sites.append(site)
        if mutation not in all_known_mutations:
            all_known_mutations.append(mutation)
all_known_sites.remove(("X", 0)) 
print(all_known_sites)

for seq in unique_sequences:
    if all([seq.seq[ks[1]-1] == ks[0] for ks in all_known_sites]):
        seq_print(seq)


# %%
padded = padding([sequences[seq_names.index("1EUM")], sequences[seq_names.index("2JD7")]])
print("Unaligned: ")
for seq in padded: 
    seq_print(seq)


# %%
# In order to align the sequences with this package, we have to add gaps (-). It seems wasteful to write, then align from file, then read again, but this process is still pretty fast. 
    
unaligned = MultipleSeqAlignment(padded)
AlignIO.write(unaligned, "alignment.fasta", "fasta")
muscle_cline = MuscleCommandline("muscle.exe", input="alignment.fasta", out = "alignment_results.fasta")
stdout, stderr = muscle_cline()
aligned = AlignIO.read("alignment_results.fasta", "fasta")

print("")
print("Aligned: ")
for seq in aligned: 
    print(seq.id[:4] + "|" +seq.seq)


# %%
# Obtaining all of the sequences near the mutation sites. Can change the number to the desired search radius (i.e. 7 residues to the left and right of the mutation).
slices = []
for site in all_known_mutations:
    slices.append([site,slice_sequence(site[1]-1, aligned[0].seq, aligned[1].seq, 7)])


# %%
# Comparing the sequences by sequence identity, it seems that in the region of L18 (on the E. coli), the structures are relatively similar.
identity_scores = [("".join([str(x) for x in sl[0]]), sequence_identity(sl[1][0], sl[1][1])) for sl in slices]
sorted_identity = sorted(identity_scores,key=lambda x: x[1])
sorted_identity.reverse()
for x in sorted_identity:
    print(x)


# %%
# However, comparing sequences by identity loses quite a bit of information (certain residues are similar). So, it might make sense to use this information.
# I copy-pasted this similarity matrix, and wrote some functions to parse it.  https://sci-hub.tw/10.1016/S1093-3263(98)80002-8 

with open("utils/similarity_matrix.txt", "r") as s:
    similarity_matrix_copied = s.read()

similarity_matrix = parse_similarity_matrix(similarity_matrix_copied)


# %%
# This is the comparison of areas near the mutation site using the sequence similarity scores, using the knowledge that certain residues are similar.
# Here, the lower the score, the more similar the area.

similarity_scores = [("".join([str(x) for x in sl[0]]), sequence_similarity_score(sl[1][0], sl[1][1], similarity_matrix)) for sl in slices]
sorted_similarity = sorted(similarity_scores,key=lambda x: x[1])
for x in sorted_similarity:
    print(x)


# %%
# Print the sequences of the most similar near-mutation regions

arg_sorted = np.argsort([x[1] for x in similarity_scores])
print(slices[arg_sorted[0]])

