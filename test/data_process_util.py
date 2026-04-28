import requests
from textwrap import wrap
import mRNACodonOptimization.src.util as util

import re


def fetch_fasta(uid):
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    return requests.get(url).text


def get_seq(fasta):
    return "".join([l for l in fasta.split("\n") if not l.startswith(">")])


def save(name, head_line, s):
    with open(name, "w") as f:
        f.write(f"{head_line}\n")
        f.write("\n".join(wrap(s, 60)))


def parse_file(name):
    return util.parse_sequence_from_file(name)


def extract_dna_by_method(name, method_name):
    """
    method_name examples:
    - VQE, optimized DNA
    - SA, optimized DNA
    - BRUTE, optimized DNA
    """
    text = util.parse_sequence_from_file(name)
    pattern = rf"{method_name}:\s*([ATCG\n ]+)"
    match = re.search(pattern, text)

    if match:
        dna = match.group(1)
        return dna.replace("\n", "").replace(" ", "")
    return None

