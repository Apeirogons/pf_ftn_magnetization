We want to increase the magnetizability of ferritin derived from P. furiosus, however, we only have mutations increasing magnetizability in E. coli. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5126674/figure/f2/

To generalize this, we want to find areas of high similarity between P. furiosus and E. coli ferritin, which would be good starting points to induce mutations to increase the magnetizability of the former.

Sequences: 2JD7 and 2JD8 are PfFtn, all else are EcFtn.

Maintained by Apeirogons (Matthew So)

Folders:
- manual_sequences contains fasta sequences manually downloaded from PDB. 
- pdb_sequences contains PDB structure+sequences downloaded from PDB.
- sequence_from_structure contains fasta sequences automatically extracted from pdb_sequences by running sequence_extractor.py
- sites contains a list of mutations found in Supplementary Information of the linked paper; other_known_sites contains sites known to be in the EcFtn found in the paper (found in Figure 2A)
- utils contains commonly used functions as well as the amino acid similarity table.
- unused contains unused files, do not expect these to be readable or even usable. 

Others:
- sequence_extractor.py used to extract sequences from pdb_structures
- muscle.exe downloaded from online, used for sequence alignment
- ftn_align.ipynb is the main development notebook
- ftn_align.py is a script converted from ftn_align.ipynb
- alignment.fasta and alignment_results.fasta should be temp files used for sequence alignment, but it does make some sense to leave them in the main directory to read their outputs without running ftn_align.py.