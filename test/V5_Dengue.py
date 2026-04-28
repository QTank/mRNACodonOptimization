import requests
import data_process_util

# --------------------------
# detect TM
# --------------------------
uid = "P17763"


def detect_tm(seq, window=20):
    hydrophobic = set("AILMFWYV")

    for i in range(len(seq) - window):
        segment = seq[i:i + window]

        score = sum(
            aa in hydrophobic
            for aa in segment
        ) / window

        if score >= 0.65:
            return i

    return None


# --------------------------
# extract Envelope region
# --------------------------

def extract_E_region(uid):
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
    data = requests.get(url).json()
    seq = data["sequence"]["value"]
    features = data["features"]

    for f in features:
        if "Envelope protein" in f.get("description", ""):
            start = f["location"]["start"]["value"]
            end = f["location"]["end"]["value"]

            E_seq = seq[start - 1:end]

            return E_seq

    return None


# --------------------------
# save fasta
# --------------------------

print("\nProcessing:", "DENV1")
E_seq = extract_E_region(uid)
print("E length:", len(E_seq))
tm_start = detect_tm(E_seq)
ectodomain = E_seq[:tm_start]
print("Ectodomain length:", len(ectodomain))

head_line = ">DENV1_E_ectodomain_vaccine_antigen"
data_process_util.save(f"../data/05-DENV1_E_ectodomain.fasta", head_line, ectodomain)
print("Save DENV1 E ectodomain mRNA vaccine antigen sequence")
