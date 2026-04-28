import data_process_util

# --------------------------
# 1. fetch VP1 sequence
# --------------------------


# Norovirus VP1 (GI.1 example)
uid = "P33471"
seq = data_process_util.get_seq(data_process_util.fetch_fasta(uid))

print("VP1 length:", len(seq))


# --------------------------
# 2. choose vaccine antigen
# --------------------------

def build_vaccine_sequence(seq, mode="full"):
    """
    mode:
    - full: full VP1 (recommended)
    - p_domain: P domain only
    """

    if mode == "full":
        return seq

    elif mode == "p_domain":
        # structural VP1 P domain boundary (~225 aa)
        return seq[225:]

    else:
        raise ValueError("mode must be full or p_domain")


# --------------------------
# 3. select antigen
# --------------------------

vaccine_antigen = build_vaccine_sequence(seq, mode="p_domain")

print("Vaccine protein length:", len(vaccine_antigen))

# --------------------------
# 4. save FASTA
# --------------------------

head_line = "Norovirus_VP1_vaccine_antigen"
data_process_util.save("../data/14-Norovirus_VP1_mRNA_vaccine.fasta", head_line,
                       vaccine_antigen)
print("Saved Norovirus_VP1 vaccine sequence")
