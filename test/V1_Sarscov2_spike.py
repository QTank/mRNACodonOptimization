import data_process_util

# --------------------------
# 1. fetch Spike protein
# --------------------------

uid = "P0DTC2"
seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))

# --------------------------
# 2. optional: remove transmembrane region (soluble ectodomain)
# --------------------------

spike_ectodomain = seq[:-120]  # approximate TM + cytoplasmic tail removal

# --------------------------
# 3. optional: RBD region extraction
# --------------------------

rbd = seq[320:550]  # approximate RBD region

# --------------------------
# 4. choose vaccine antigen
# --------------------------

vaccine_antigen = spike_ectodomain  # or rbd

# --------------------------
# 5. save FASTA
# --------------------------

head_line = ">SARSCoV2_Spike_vaccine_antigen"

data_process_util.save("../data/01-sars2_spike_vaccine.fasta.bk", head_line,
                       vaccine_antigen)
print("Saved SARS-CoV-2 Spike mRNA vaccine antigen sequence")
