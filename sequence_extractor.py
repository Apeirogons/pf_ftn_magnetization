# Extracts sequences from pdb_structures to sequence_from_structure. Unfortunately, these sequences seem incomplete for some reason (maybe some residues at beginning+end got cut off)
# so for the purpose of ftn_align, we will use the fasta files manually downloaded from PDB.

from os import listdir
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PPBuilder
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

structures = []
pdb_ids = []
structures_dir = "pdb_structures"
parser = MMCIFParser()

# Read structures from IO
for item in listdir(structures_dir):
    if item.find('.') == -1:
        for subitem in listdir(structures_dir + "/" + item):
            print("Parsing " + subitem)
            structures.append(parser.get_structure(subitem[:4], structures_dir + "/" + item + "/" + subitem))
            pdb_ids.append(subitem[:4])

# Extract peptide sequences and write to sequence_from_structure
ppb = PPBuilder()
for i, structure in enumerate(structures):
    pdb_id = pdb_ids[i]
    print(pdb_id)
    peptides = ppb.build_peptides(structure)
    seqs = []
    for peptide in peptides:
        seqs.append(peptide.get_sequence())
    sorted_seqs = sorted(seqs, key=len)

    AlignIO.write(MultipleSeqAlignment([SeqRecord(sorted_seqs[-1], id=pdb_id)]), "sequence_from_structure/" + pdb_id +".fasta", "fasta")