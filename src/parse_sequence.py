import re
import util
def extract_mrna_lists(text):
    """
    Extract mRNA sequences into three lists:
    - vqe_list
    - sa_list
    - brute_list
    """

    vqe_list = []
    sa_list = []
    brute_list = []

    # Find all VQE mRNA sequences
    vqe_matches = re.findall(
        r"\[VQE\][\s\S]*?mRNA\s*:\s*([AUGC]+)",
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
    sa_list.extend(sa_matches)
    brute_list.extend(brute_matches)

    return vqe_list, sa_list, brute_list


# Example usage
if __name__ == "__main__":
    with open("../experiments/resulti_5.log", "r") as f:
        text = f.read()
    vqe_list, sa_list, brute_list = extract_mrna_lists(text)
    print("VQE list:", len(vqe_list))
    vqe_rna = "".join(vqe_list)
    print("VQE, optimized RNA:")
    print(vqe_rna)
    print("VQE, optimized DNA:")
    print(vqe_rna.replace("U", "T"))
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