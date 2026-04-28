import requests
from textwrap import wrap
import data_process_util
# --------------------------
# 1. fetch Marburg GP
# --------------------------

uid = "Q05320"
seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))

# --------------------------
# 2. remove transmembrane region
# --------------------------

ectodomain = seq[:-60]  # GP2 TM region removal (approx)

# --------------------------
# 3. extract GP1-like region (optional)
# --------------------------

gp1_region = ectodomain[:550]

# --------------------------
# 4. choose vaccine antigen
# --------------------------

vaccine_antigen = gp1_region

# --------------------------
# 5. save FASTA
# --------------------------

head_line = ">Marburg_GP_vaccine_antigen"
data_process_util.save("../data/11-marburg_vaccine_antigen.fasta", head_line, vaccine_antigen)
print("Saved Marburg mRNA vaccine antigen sequence")