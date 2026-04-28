import requests
import data_process_util
from textwrap import wrap

uid = "A0A024B7W1"

# --------------------------
# 1. fetch JSON annotation
# --------------------------

url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"

data = requests.get(url).json()


# --------------------------
# 2. find Envelope region
# --------------------------

def find_E_region(data):
    features = data["features"]

    for f in features:

        if "Envelope" in f.get("description", ""):
            start = f["location"]["start"]["value"]
            end = f["location"]["end"]["value"]

            return start, end

    return None, None


E_START, E_END = find_E_region(data)

print("E region:", E_START, E_END)

# --------------------------
# 3. get sequence
# --------------------------

seq = data["sequence"]["value"]

E_seq = seq[E_START - 1:E_END]

print("E length:", len(E_seq))


# --------------------------
# 4. detect TM
# --------------------------

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


tm_start = detect_tm(E_seq)

ectodomain = E_seq[:tm_start]

print("Ectodomain length:", len(ectodomain))

# --------------------------
# 5. save
# --------------------------

head_line = ">Zika_E_ectodomain"
data_process_util.save("../data/04-Zika_E_ectodomain.fasta", head_line, ectodomain)

print("Save Zika E ectodomain mRNA vaccine antigen sequence")
