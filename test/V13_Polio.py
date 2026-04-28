import data_process_util

uid = "P03300"  # Poliovirus VP1
seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))
# --------------------------
# Extract VP1 region
# --------------------------

VP1_START = 520
VP1_END = 817

vp1 = seq[VP1_START:VP1_END]

print("VP1 length:", len(vp1))

# sanity check
if 250 <= len(vp1) <= 350:
    print("VP1 extraction looks correct ✔")
else:
    print("Check region — something may be wrong ❌")

# --------------------------
# Save FASTA
# --------------------------

head_line = ">Poliovirus_VP1_mRNA_vaccine"
data_process_util.save("../data/13-Poliovirus_VP1.fasta", head_line, vp1)
print("Saved VP1 sequence")
