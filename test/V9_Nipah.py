import data_process_util

# --------------------------
# 1. fetch Nipah G protein
# --------------------------

uid = "Q9IH63"
seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))

# --------------------------
# 2. remove transmembrane region
# --------------------------

ectodomain = seq[:-50]  # approximate TM + tail removal

# --------------------------
# 3. extract receptor-binding region (heuristic)
# --------------------------

rbd_like = ectodomain[:500]  # G head domain approx

# --------------------------
# 4. choose vaccine antigen
# --------------------------

vaccine_antigen = rbd_like  # or rbd_like

# --------------------------
# 5. save FASTA
# --------------------------

head_line = ">Nipah_G_vaccine_antigen"
data_process_util.save("../data/09-nipah_vaccine_antigen.fasta", head_line, vaccine_antigen)
print("Saved Nipah mRNA vaccine antigen sequence")
