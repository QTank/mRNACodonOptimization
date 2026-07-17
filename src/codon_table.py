GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

CODON_USAGE_RAW = {
    "h_sapiens_9606": {  # Homo sapiens
        "TTT": 17.6, "TTC": 20.3, "TTA": 7.7, "TTG": 12.9,
        "CTT": 13.2, "CTC": 19.6, "CTA": 7.2, "CTG": 39.6,
        "ATT": 16.0, "ATC": 20.8, "ATA": 7.5, "ATG": 22.0,
        "GTT": 11.0, "GTC": 14.5, "GTA": 7.1, "GTG": 28.1,
        "TCT": 15.2, "TCC": 17.7, "TCA": 12.2, "TCG": 4.4,
        "CCT": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9,
        "ACT": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1,
        "GCT": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4,
        "TAT": 12.2, "TAC": 15.3, "TAA": 1.0, "TAG": 0.8,
        "CAT": 10.9, "CAC": 15.1, "CAA": 12.3, "CAG": 34.2,
        "AAT": 17.0, "AAC": 19.1, "AAA": 24.4, "AAG": 31.9,
        "GAT": 21.8, "GAC": 25.1, "GAA": 29.0, "GAG": 39.6,
        "TGT": 10.6, "TGC": 12.6, "TGA": 1.6, "TGG": 13.2,
        "CGT": 4.5, "CGC": 10.4, "CGA": 6.2, "CGG": 11.4,
        "AGT": 12.1, "AGC": 19.5, "AGA": 12.2, "AGG": 12.0,
        "GGT": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5,
    },
    "e_coli_316407": {  # E. coli K12
        "TTT": 22.1, "TTC": 16.6, "TTA": 13.9, "TTG": 13.4,
        "CTT": 11.0, "CTC": 10.6, "CTA": 3.9, "CTG": 52.6,
        "ATT": 30.3, "ATC": 25.1, "ATA": 4.4, "ATG": 27.9,
        "GTT": 18.3, "GTC": 15.3, "GTA": 10.9, "GTG": 26.4,
        "TCT": 8.5, "TCC": 8.6, "TCA": 7.2, "TCG": 8.9,
        "CCT": 7.0, "CCC": 5.5, "CCA": 8.6, "CCG": 23.2,
        "ACT": 9.0, "ACC": 22.8, "ACA": 7.7, "ACG": 14.5,
        "GCT": 15.3, "GCC": 25.5, "GCA": 20.5, "GCG": 32.9,
        "TAT": 16.1, "TAC": 12.1, "TAA": 2.0, "TAG": 0.2,
        "CAT": 12.5, "CAC": 9.3, "CAA": 15.3, "CAG": 28.9,
        "AAT": 17.9, "AAC": 21.6, "AAA": 33.6, "AAG": 10.3,
        "GAT": 32.1, "GAC": 19.1, "GAA": 39.1, "GAG": 18.7,
        "TGT": 5.2, "TGC": 6.1, "TGA": 1.0, "TGG": 13.9,
        "CGT": 20.9, "CGC": 22.0, "CGA": 3.6, "CGG": 5.4,
        "AGT": 8.8, "AGC": 15.2, "AGA": 2.1, "AGG": 1.2,
        "GGT": 24.7, "GGC": 29.7, "GGA": 8.1, "GGG": 11.1,
    },
}


def get_codons_table(organism="e_coli_316407", raw_usage=None):
    """
    Parameters
    ----------
    organism : str
         "h_sapiens_9606", "e_coli_316407"。
    raw_usage : dict[str, float], optional

    Returns
    -------
    dict[str, dict[str, float]]
    """
    if raw_usage is None:
        if organism not in CODON_USAGE_RAW:
            raise ValueError(
                f"内置数据没有物种 '{organism}', 可用: {list(CODON_USAGE_RAW)}, "
                f"或者通过 raw_usage 参数传入自定义数据。"
            )
        raw_usage = CODON_USAGE_RAW[organism]

    grouped = {}
    for codon, aa in GENETIC_CODE.items():
        grouped.setdefault(aa, {})[codon] = raw_usage.get(codon, 0.0)

    codons_table = {}
    for aa, codon_counts in grouped.items():
        total = sum(codon_counts.values())
        if total == 0:
            codons_table[aa] = {c: 0.0 for c in codon_counts}
        else:
            codons_table[aa] = {c: round(v / total, 2) for c, v in codon_counts.items()}

    return codons_table


def get_optimal_codons(organism="e_coli_316407", raw_usage=None):
    codons_table = get_codons_table(organism, raw_usage=raw_usage)
    optimal = {aa: max(freqs, key=freqs.get) for aa, freqs in codons_table.items() if freqs}
    return optimal, codons_table


if __name__ == "__main__":
    if False:
        for org in ("h_sapiens_9606", "e_coli_316407"):
            test_pct = pct.get_codons_table(org)
            test_py = get_codons_table(org)
