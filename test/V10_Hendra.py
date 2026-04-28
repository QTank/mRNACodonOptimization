import data_process_util

# --------------------------
# 1. fetch Hendra G protein
# --------------------------

uid = "Q8V4Q9"  # Hendra glycoprotein G (example UniProt)
seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))

# --------------------------
# 2. remove TM region (approx)
# --------------------------

ectodomain = seq[:-50]  # TM + tail removal heuristic

# --------------------------
# 3. extract head domain (optional)
# --------------------------

head_domain = ectodomain[:500]

# --------------------------
# 4. choose vaccine antigen
# --------------------------

vaccine_antigen = ectodomain  # or head_domain


# --------------------------
# 5. save FASTA
# --------------------------

head_line = ">Hendra_G_vaccine_antigen"
data_process_util.save("../data/10-hendra_vaccine_antigen.fasta", head_line, vaccine_antigen)
print("Saved Hendra mRNA vaccine antigen sequence")