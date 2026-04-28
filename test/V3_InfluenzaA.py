import data_process_util

# --------------------------
# 1. fetch Influenza HA protein
# --------------------------

uid = "P03437"  # example H3 HA
seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))

# --------------------------
# 2. remove transmembrane region (ectodomain)
# --------------------------

ha_ectodomain = seq[:-60]  # TM + cytoplasmic tail approx

# --------------------------
# 3. extract HA head domain (optional)
# --------------------------

ha_head = seq[:330]  # antigenic head region approx

# --------------------------
# 4. choose vaccine antigen
# --------------------------

vaccine_antigen = ha_ectodomain  # or ha_head

# --------------------------
# 5. save FASTA
# --------------------------

head_line = ">InfluenzaA_HA_vaccine_antigen"
data_process_util.save("../data/03-influenza_ha_vaccine.fasta", head_line, vaccine_antigen)
print("Saved Influenza HA mRNA vaccine antigen sequence")
