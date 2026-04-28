import data_process_util

# --------------------------
# 1. fetch Rabies G protein
# --------------------------

uid = "P08667"

seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))

# --------------------------
# 2. remove transmembrane region (C-terminal anchor)
# --------------------------

# Rabies G TM ~ last 40 aa (typical lyssavirus GP)
ectodomain = seq[:-40]

# --------------------------
# 3. optional stabilization trimming (remove flexible tail)
# --------------------------
vaccine_antigen = ectodomain[:500]

# --------------------------
# 4. save vaccine-ready sequence
# --------------------------

head_line = '>Rabies_G_vaccine_antigen'
data_process_util.save("../data/06-rabies_vaccine_antigen.fasta", head_line, vaccine_antigen)
print("Saved rabies mRNA vaccine antigen sequence")
