import re

import util

protein_seq = util.parse_sequence_from_file("../data/01-sars2_spike_vaccine.fasta", 6)
seq_list = util.split_sequence(protein_seq[1:], 6)
log_file = "../experiments/resulti_6.log"  # update path if needed

sequences = []

with open(log_file, "r") as f:
    for line in f:
        match = re.match(r"Sequence:\s+(\S+)", line.strip())
        if match:
            seq = match.group(1)
            if seq not in sequences:  # deduplicate
                sequences.append(seq)

for i, seq in enumerate(sequences, 1):
    print(f"{i}: {seq}", seq in seq_list)

print(f"\nTotal: {len(sequences)} sequences")
print(seq_list)