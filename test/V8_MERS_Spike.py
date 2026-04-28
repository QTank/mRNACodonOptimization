import data_process_util

# --------------------------
# 1. fetch MERS Spike
# --------------------------

uid = "K9N5Q8"
seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))

# --------------------------
# 2. remove transmembrane region
# --------------------------

ectodomain = seq[:-60]  # spike TM + tail approx removal

# --------------------------
# 3. optional: keep only RBD region (simplified heuristic)
# --------------------------

rbd = ectodomain[int(len(ectodomain) * 0.3):int(len(ectodomain) * 0.6)]

# --------------------------
# 4. choose vaccine antigen
# --------------------------

vaccine_antigen = rbd

# --------------------------
# 5. save FASTA
# --------------------------

head_line = ">MERS_Spike_vaccine_antigen"
data_process_util.save("../data/08-mers_vaccine_antigen.fasta", head_line, vaccine_antigen)
print("Saved MERS mRNA vaccine antigen sequence")
