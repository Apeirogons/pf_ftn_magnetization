# Detect and remove duplicate sequences.
def remove_duplicates(sequences):
    unique_sequences = []
    for seq in sequences:
      in_unique = False
      for unique_seq in [unique for unique in unique_sequences]:
          if seq.seq == unique_seq.seq:
              in_unique = True
              print("The sequence of " + seq.id[:4] + " is the same as " + unique_seq.id[:4] + ", removing "+ seq.id[:4])
              break
      if not in_unique:
         unique_sequences.append(seq)
    return unique_sequences

# Pad sequences with - to make them the same length for processing with MUSCLE.
def padding(sequences):
    padded_sequences = []
    max_length = max([len(seq) for seq in sequences])
    for seq in sequences:
        padded_sequences.append(seq + (max_length-len(seq))*'-')
    return padded_sequences

# Function to print a Sequence object in a more accessible way.
def seq_print(seq):
    print(seq.id[:4] + "|" + seq.seq)