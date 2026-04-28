import data_process_util

uid = "P03138"

seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))

print("Full length:", len(seq))


# --------------------------
# 1. OPTIONAL: signal peptide (better than fixed 20)
# --------------------------

def detect_signal_peptide(seq):
    # simple hydrophobic N-terminal rule
    window = 20
    segment = seq[:window]
    hydrophobic = sum(a in "AILMFWYV" for a in segment)
    return hydrophobic > 12  # heuristic


if detect_signal_peptide(seq):
    signal_removed = seq[20:]  # approximation OK
else:
    signal_removed = seq

# --------------------------
# 2. KEEP full antigen (recommended)
# --------------------------

vaccine_antigen = signal_removed

# --------------------------
# 3. save FASTA
# --------------------------

head_line = ">HBV_HBsAg_mRNA_vaccine"
data_process_util.save("../data/15-HBV_HBsAg_vaccine.fasta", head_line, vaccine_antigen)
print("Saved HBV mRNA vaccine sequence")
