import math
from collections import defaultdict
import python_codon_tables as pct

import data_process_util


# ---------------------------------------------------------------------------
# Helpers: build codon->amino-acid map and RSCU from the library
# ---------------------------------------------------------------------------

def _build_codon_to_aa(table_name: str) -> dict:
    """
    Return a flat {codon: amino_acid} dict from a python-codon-tables table.
    Stop codons map to '*'.
    """
    raw = pct.get_codons_table(table_name)  # {aa: {codon: freq, ...}, ...}
    codon_to_aa = {}
    for aa, codons in raw.items():
        for codon in codons:
            codon_to_aa[codon] = aa
    return codon_to_aa


def _freq_table_to_rscu(table_name: str) -> dict:
    """
    Convert python-codon-tables relative-frequency values to RSCU.

    The library stores *relative synonymous frequencies* that sum to 1 per AA
    (i.e. f_i = count_i / total_counts_for_AA).  RSCU is defined as:

        RSCU_i = f_i * n_synonymous

    so the uniform RSCU for every codon of an n-fold degenerate AA equals 1.
    """
    raw = pct.get_codons_table(table_name)  # {aa: {codon: freq}}
    rscu = {}
    for aa, codons in raw.items():
        n = len(codons)
        for codon, freq in codons.items():
            rscu[codon] = freq * n
    return rscu


def _rscu_to_w(rscu: dict, codon_to_aa: dict) -> dict:
    """
    Compute relative adaptiveness w for each codon:
        w_i = RSCU_i / max(RSCU of synonymous group)
    Stop codons ('*') are excluded.
    """
    aa_codons = defaultdict(list)
    for codon, aa in codon_to_aa.items():
        if aa != '*':
            aa_codons[aa].append(codon)

    w = {}
    for aa, codons in aa_codons.items():
        max_rscu = max(rscu.get(c, 0.0) for c in codons)
        if max_rscu == 0:
            max_rscu = 1.0
        for c in codons:
            w[c] = rscu.get(c, 0.0) / max_rscu
    return w


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def list_organisms() -> list:
    """Return all organism table names available in python-codon-tables."""
    return pct.available_codon_tables_names


def calculate_cai(
        sequence: str,
        organism: str = 'e_coli_316407',
        custom_rscu: dict = None,
) -> dict:
    """
    Calculate the Codon Adaptation Index (CAI) for a DNA coding sequence.

    Parameters
    ----------
    sequence : str
        DNA coding sequence (5'->3', ATG ... stop optional).
    organism : str
        Organism table name from python-codon-tables
        (default: 'e_coli_316407').  Ignored for RSCU when custom_rscu given.
    custom_rscu : dict, optional
        User-supplied {codon: RSCU_value} table.  When provided, the
        organism argument is still used for the codon->AA mapping.

    Returns
    -------
    dict with keys:
        cai        – float, final CAI score [0, 1]
        n_codons   – int, number of sense codons scored (Met/Trp excluded)
        codon_ws   – list of (codon, w) tuples in sequence order
        w_table    – dict {codon: w} for the full table
        rscu       – dict {codon: RSCU} used
        organism   – str, resolved organism name
        warnings   – list of warning strings
    """
    warnings_list = []

    # Validate organism
    available = pct.available_codon_tables_names
    if organism not in available:
        raise ValueError(
            f"Unknown organism '{organism}'.\nAvailable: {available}"
        )

    # Build codon->AA map from the library (replaces hardcoded CODON_TABLE)
    codon_to_aa = _build_codon_to_aa(organism)

    # Derive RSCU (from library frequencies, or use custom table)
    rscu = custom_rscu if custom_rscu is not None else _freq_table_to_rscu(organism)

    # Derive w table
    w_table = _rscu_to_w(rscu, codon_to_aa)

    # Clean sequence

    seq = sequence.upper().replace(' ', '').replace('\n', '').replace('\r', '')
    bad_chars = sorted({b for b in seq if b not in 'ATCG'})
    if bad_chars:
        warnings_list.append(f"Non-ATCG characters ignored: {bad_chars}")
        seq = ''.join(b for b in seq if b in 'ATCG')

    if len(seq) % 3 != 0:
        trim = len(seq) % 3
        warnings_list.append(
            f"Sequence length {len(seq)} not divisible by 3; trimming last {trim} base(s)."
        )
        seq = seq[: len(seq) - trim]

    # Score each codon
    codon_ws = []
    for i in range(0, len(seq), 3):
        codon = seq[i: i + 3]
        if len(codon) < 3:
            break

        aa = codon_to_aa.get(codon)
        if aa is None:
            warnings_list.append(f"Unrecognised codon '{codon}' at pos {i + 1}, skipped.")
            continue
        if aa == '*':
            continue  # stop codon — excluded by convention

        # Single-codon AAs (Met=ATG, Trp=TGG) always w=1;
        # including them would inflate CAI — standard practice excludes them.
        if aa in ('M', 'W'):
            continue

        w = w_table.get(codon, 0.0)
        if w <= 0:
            warnings_list.append(
                f"Codon '{codon}' (AA={aa}) has w=0. Using epsilon to avoid log(0)."
            )
            w = 1e-10

        codon_ws.append((codon, w))

    if not codon_ws:
        raise ValueError(
            "No scoreable codons found after excluding Met, Trp, and stop codons."
        )

    # CAI = geometric mean of all per-codon w values
    log_sum = sum(math.log(w) for _, w in codon_ws)
    cai = math.exp(log_sum / len(codon_ws))

    return {
        'cai': round(cai, 6),
        'n_codons': len(codon_ws),
        'codon_ws': codon_ws,
        'w_table': w_table,
        'rscu': rscu,
        'organism': organism,
        'warnings': warnings_list,
    }


def compute_rscu_from_sequences(sequences: list, organism: str = 'e_coli_316407') -> dict:
    """
    Compute a custom RSCU table from a list of reference coding sequences
    (e.g. a set of highly-expressed genes for your organism).

    Parameters
    ----------
    sequences : list[str]
        DNA coding sequences used as the reference set.
    organism : str
        Used only to obtain the codon->AA mapping from python-codon-tables.

    Returns
    -------
    dict {codon: RSCU}
    """
    codon_to_aa = _build_codon_to_aa(organism)

    aa_codons = defaultdict(list)
    for codon, aa in codon_to_aa.items():
        if aa != '*':
            aa_codons[aa].append(codon)

    total_counts = defaultdict(int)
    for seq in sequences:
        s = seq.upper().replace(' ', '').replace('\n', '')
        for i in range(0, len(s) - 2, 3):
            codon = s[i: i + 3]
            if len(codon) == 3 and all(b in 'ATCG' for b in codon):
                if codon_to_aa.get(codon, '*') != '*':
                    total_counts[codon] += 1

    rscu = {}
    for aa, codons in aa_codons.items():
        total = sum(total_counts[c] for c in codons)
        n_syn = len(codons)
        for c in codons:
            if total == 0:
                rscu[c] = 1.0
            else:
                expected = total / n_syn
                rscu[c] = total_counts[c] / expected if expected > 0 else 0.0
    return rscu


# ---------------------------------------------------------------------------
# Pretty-print report
# ---------------------------------------------------------------------------

def print_cai_report(result: dict, seq_len: int = None):
    cai = result['cai']
    print("=" * 57)
    print("  CAI ANALYSIS REPORT")
    print("=" * 57)
    print(f"  Organism              : {result['organism']}")
    if seq_len is not None:
        print(f"  Sequence length       : {seq_len} bp")
    print(f"  Codons scored         : {result['n_codons']}  (Met/Trp/stops excluded)")
    print(f"  CAI score             : {cai:.4f}")
    print()

    if result['warnings']:
        print("  Warnings:")
        for w in result['warnings']:
            print(f"    ⚠  {w}")
        print()

    if cai >= 0.80:
        level = "HIGH  – well-optimized for expression"
    elif cai >= 0.60:
        level = "MODERATE – reasonable expression expected"
    elif cai >= 0.40:
        level = "LOW – consider codon optimization"
    else:
        level = "VERY LOW – strong codon bias mismatch"
    print(f"  Interpretation        : {level}")
    print("=" * 57)


def gc_content(dna_seq):
    if len(dna_seq) == 0: return 0
    dna_seq = dna_seq.upper()

    g = dna_seq.count('G')
    c = dna_seq.count('C')

    return (g + c) / len(dna_seq) * 100


def run_cai_compare(file_name):
    solver_list = ['SA']
    for solver in solver_list:
        seq = data_process_util.extract_dna_by_method(file_name, solver + ", optimized DNA")
        result = calculate_cai(seq, organism='e_coli_316407')
        print(f"CAI({solver}): {result['cai']:.4f}, GC Content: {gc_content(seq):.2f}")
    return result

# ---------------------------------------------------------------------------
# Calculate CAI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from pathlib import Path
    path = Path("../experiments/")
    txt_files = [f.name for f in path.glob("*.log")]
    val = 0
    for file_name in txt_files:
        print(file_name)
        r = run_cai_compare("../experiments/" + file_name)
        val += r['cai']
        print()

    print(f"{val/len(txt_files)}")