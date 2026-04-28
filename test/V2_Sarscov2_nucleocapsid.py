import data_process_util

# --------------------------
# 1. fetch SARS-CoV-2 N protein
# --------------------------

uid = "P0DTC9"
seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))
# --------------------------
# 2. N protein is internal → no TM removal needed
# --------------------------

n_full = seq

# --------------------------
# 3. optional: immunodominant region (T-cell rich)
# --------------------------

tcell_region = n_full[50:350]

# --------------------------
# 4. choose vaccine antigen
# --------------------------

vaccine_antigen = tcell_region

# --------------------------
# 5. save FASTA
# --------------------------

head_line = ">SARS2_Nucleocapsid_vaccine_antigen"
data_process_util.save("../data/02-sars2_n_vaccine.fasta", head_line, vaccine_antigen)
print("Saved SARS-CoV-2 Nucleoprotein protein mRNA vaccine antigen sequence")
