import re
import util
import numpy as np

def extract_mrna_lists(text):
    """
    Extract mRNA sequences into three lists:
    - vqe_list
    - sa_list
    - brute_list
    """

    vqe_list = []
    qaoa_list = []
    sa_list = []
    brute_list = []

    # Find all VQE mRNA sequences
    vqe_matches = re.findall(
        r"\[VQE\][\s\S]*?mRNA\s*:\s*([AUGC]+)",
        text
    )
    qaoa_matches = re.findall(
        r"\[QAOA\][\s\S]*?mRNA\s*:\s*([AUGC]+)",
        text
    )

    # Find all SA mRNA sequences
    sa_matches = re.findall(
        r"\[SA\][\s\S]*?mRNA\s*:\s*([AUGC]+)",
        text
    )

    # Find all BRUTE mRNA sequences
    brute_matches = re.findall(
        r"\[BRUTE\][\s\S]*?mRNA\s*:\s*([AUGC]+)",
        text
    )

    vqe_list.extend(vqe_matches)
    qaoa_list.extend(qaoa_matches)
    sa_list.extend(sa_matches)
    brute_list.extend(brute_matches)

    return vqe_list, qaoa_list, sa_list, brute_list


def parse_time(text):

    # Pattern to capture method name and time
    pattern = r"\[(VQE|SA|QAOA|BRUTE)\][\s\S]*?- Time\s+:\s+([\d.]+)s"

    matches = re.findall(pattern, text)

    # Store results
    times = {"QAOA": [], "BRUTE": [], "VQE": [], "SA": []}

    for method, time in matches:
        times[method].append(float(time))

    # Print results
    print("Parsed Times:")
    for method in times:
        if times[method]:
            print(f"{method}: {np.mean(times[method]):.2f}, {times[method]}")

# Example usage
if __name__ == "__main__":
    with open("../experiments/03-influenza_ha_vaccine_qaoa.log", "r") as f:
        text = f.read()
    vqe_list, qaoa_list, sa_list, brute_list = extract_mrna_lists(text)
    print("VQE list:", len(vqe_list))
    vqe_rna = "".join(vqe_list)
    print("VQE, optimized RNA:")
    print(vqe_rna)
    print("VQE, optimized DNA:")
    print(vqe_rna.replace("U", "T"))

    print("QAOA list:", len(qaoa_list))
    qaoa_rna = "".join(qaoa_list)
    print("VQE, optimized RNA:")
    print(qaoa_rna)
    print("VQE, optimized DNA:")
    print(qaoa_rna.replace('U', "T"))

    print("\nSA list:")
    sa_rna = "".join(sa_list)
    print("SA, optimized RNA:")
    print(sa_rna)
    print("SA, optimized DNA:")
    print(util.convert_rna_to_dna(sa_rna))
    print("\nBRUTE list:")
    brute_rna = "".join(brute_list)
    print("BRUTE, optimized RNA:")
    print(brute_rna)
    print("BRUTE, optimized DNA:")
    print(util.convert_rna_to_dna(brute_rna))

    parse_time(text)