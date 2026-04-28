import data_process_util

# --------------------------
# 1. fetch Ebola GP sequence
# --------------------------

uid = "Q05320"

seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))

# --------------------------
# 2. remove TM region (GP2 anchor)
# --------------------------

ectodomain = seq[:-50]  # approximate TM + tail removal

# --------------------------
# 3. remove flexible C-terminal tail
# --------------------------

vaccine_antigen = ectodomain[:500]

# --------------------------
# 4. save vaccine-ready sequence
# --------------------------

head_line = ">Ebola_GP_vaccine_antigen"
data_process_util.save("../data/07-ebola_vaccine_antigen.fasta", head_line, vaccine_antigen)
print("Saved Ebola mRNA vaccine antigen sequence")
