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
    with open(name, 'r') as f:
        text = f.read()
    if method_name in ['VQE', 'QAOA', 'SA']:
        opt_dna = re.search(rf'{method_name},\s+optimized\s+DNA:\s*\n\s*([ATCG]+)', text)
        if opt_dna:
            opt_seq = opt_dna.group(1)
            print(f"{method_name} Success! Length: {len(opt_seq)} bp")
            return opt_seq
        else:
            print("No match found")
    return None

