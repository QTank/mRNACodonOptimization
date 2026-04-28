import data_process_util

# --------------------------
# 1. fetch Lassa GP
# --------------------------

uid = "P08669"
seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))

# --------------------------
# 2. remove transmembrane region
# --------------------------

ectodomain = seq[:-55]  # GP2 TM region approx

# --------------------------
# 3. extract GP1-like region
# --------------------------

gp1_region = ectodomain[:500]

# --------------------------
# 4. choose vaccine antigen
# --------------------------

vaccine_antigen = gp1_region

# --------------------------
# 5. save FASTA
# --------------------------

head_line = ">Lassa_GP_vaccine_antigen"
data_process_util.save("../data/12-lassa_vaccine_antigen.fasta", head_line, vaccine_antigen)
print("Saved Lassa mRNA vaccine antigen sequence")
