import vqe_solver
import qaoa_solver
from denseCodon import DenseCodon
from oneHotCodon import OneHotCodon
import util
import time, json
from hamiltonian import CodonOptimizer
import sa_solver
import brute_force
import numpy as np

def run_optimization(sequence, config, type_opt="vqe", encoding_type="dense"):
    start_time = time.time()

    if type_opt == "vqe":
        print("Starting VQE (Variational Quantum Eigensolver) on a Quantum Simulator...\n")
        codon_opt = CodonOptimizer(sequence, config, DenseCodon, "dense")
        qubit_op = codon_opt.create_qubit_op()
        bitstring, energy = vqe_solver.get_min(qubit_op, config['vqe_settings'])
        final_mRNA_sequence = util.decode_bitstring(sequence, bitstring,
                                                    table_name=config['metadata']['table_name'])

    elif type_opt == "qaoa":
        print("Starting QAOA on a Quantum Simulator...\n")
        codon_opt = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")
        qubit_op = codon_opt.create_qubit_op()
        bitstring = qaoa_solver.get_min(qubit_op, config['qaoa_settings'])
        final_mRNA_sequence = util.decode_one_hot_bitstring(sequence, bitstring,
                                                            table_name=config['metadata']['table_name'])

    elif type_opt == "sa":
        print("Starting Simulated Annealing optimization...\n")
        if encoding_type == "dense":
            codon_opt = CodonOptimizer(sequence, config, DenseCodon, "dense")
            qubit_op = codon_opt.create_qubit_op()
            bitstring, energy = sa_solver.get_min(qubit_op, codon_opt.qubit_len, config['sa_settings'], codon_opt.codon_list)
            final_mRNA_sequence = util.decode_bitstring(sequence, bitstring,
                                                        table_name=config['metadata']['table_name'])
        else:
            codon_opt = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")
            qubit_op = codon_opt.create_qubit_op()
            bitstring, energy = sa_solver.get_min(qubit_op, codon_opt.qubit_len, config['sa_settings'], codon_opt.codon_list)
            final_mRNA_sequence = util.decode_one_hot_bitstring(sequence, bitstring,
                                                            table_name=config['metadata']['table_name'])

    elif type_opt == "brute":
        print("Starting Brute Force optimization...\n")
        if encoding_type == "dense":
            codon_opt = CodonOptimizer(sequence, config, DenseCodon, "dense")
            qubit_op = codon_opt.create_qubit_op()
            bitstring, energy = brute_force.get_min(qubit_op, codon_opt.codon_list)
            final_mRNA_sequence = util.decode_bitstring(sequence, bitstring,
                                                        table_name=config['metadata']['table_name'])
        else:
            codon_opt = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")
            qubit_op = codon_opt.create_qubit_op()
            bitstring, energy = brute_force.get_min(qubit_op, codon_opt.codon_list)
            final_mRNA_sequence = util.decode_one_hot_bitstring(sequence, bitstring,
                                                            table_name=config['metadata']['table_name'])

    else:
        raise ValueError(f"Unknown type_opt: {type_opt}")
 
    elapsed = time.time() - start_time

    print(f"  - Sequence  : {sequence}")
    print(f"  - Qubits    : {codon_opt.qubit_len}")
    print(f"  - Bitstring : {bitstring}")
    print(f"  - mRNA      : {final_mRNA_sequence}")
    print(f"  - Energy    : {energy:.4f}")
    print(f"  - Time      : {elapsed:.2f}s\n")

    final_string_index = 0
    for index in range(len(sequence)):
        amino_acid = sequence[index]
        chunk = bitstring[final_string_index:final_string_index + codon_opt.codon_list[index].encoding_qubit_len]
        codon = final_mRNA_sequence[index * 3: (index + 1) * 3] if final_mRNA_sequence else "???"
        print(f"    {amino_acid} → {chunk} → {codon}")
        final_string_index += codon_opt.codon_list[index].encoding_qubit_len
    print()

    return final_mRNA_sequence, energy, elapsed


def mrna_to_one_hot_bitstring(mrna_sequence, protein_sequence, table_name='e_coli_316407'):
    """
    Convert a decoded mRNA sequence back to a one-hot bitstring.
    This allows fair energy comparison across solvers that use different encodings.

    Example:
        protein = "FM", mrna = "UUUAUG"
        F has codons [TTT, TTC] → TTT is index 0 → one-hot "10"
        M has codons [ATG]      → ATG is index 0 → one-hot "1"
        result: "101"
    """
    import python_codon_tables as pct
    codon_table = pct.get_codons_table(table_name)
    bitstring = ""

    # mRNA uses U, codon table uses T
    dna_sequence = mrna_sequence.replace("U", "T")

    for i, amino in enumerate(protein_sequence):
        codon = dna_sequence[i * 3: (i + 1) * 3]
        codon_list = list(codon_table[amino].keys())
        encoding_len = len(codon_list)

        if codon not in codon_list:
            print(f"Warning: codon {codon} not found for amino {amino}")
            return None

        idx = codon_list.index(codon)
        one_hot = ["0"] * encoding_len
        one_hot[idx] = "1"
        bitstring += "".join(one_hot)

    return bitstring

def mrna_to_dense_bitstring(mrna_sequence, protein_sequence, table_name='e_coli_316407'):
    """Convert mRNA back to dense (binary) bitstring for energy evaluation."""
    import python_codon_tables as pct
    codon_table = pct.get_codons_table(table_name)
    bitstring = ""
    dna_sequence = mrna_sequence.replace("U", "T")

    for i, amino in enumerate(protein_sequence):
        codon = dna_sequence[i * 3: (i + 1) * 3]
        codon_list = list(codon_table[amino].keys())
        encoding_len = round(np.log2(len(codon_list)))

        if encoding_len == 0:
            continue  # Met/Trp — no bits needed

        idx = codon_list.index(codon)
        bitstring += f"{idx:0{encoding_len}b}"

    return bitstring

def compare_solvers(sequence, config):
    """Run VQE, SA, and Brute Force on the same sequence and compare energy."""
    print("=" * 60)
    print(f"Sequence: {sequence}")
    print("=" * 60)

    results = {}
    for solver in ["vqe", "sa", "brute"]:
        print(f"[{solver.upper()}]")
        mRNA, energy, elapsed = run_optimization(sequence, config, type_opt=solver)
        results[solver] = {"mRNA": mRNA, "energy": energy, "time": elapsed}

    brute_energy = results["brute"]["energy"]

    print("=" * 60)
    print("SUMMARY")
    print(f"{'Solver':<10} {'mRNA':<25} {'Energy':>12} {'Gap':>10} {'Time':>8}")
    print("-" * 70)
    for solver, res in results.items():
        gap = res["energy"] - brute_energy
        match = "✓" if abs(gap) < 1e-6 else f"+{gap:.2f}"
        print(f"{solver.upper():<10} {str(res['mRNA']):<25} {res['energy']:>12.4f} {match:>10} {res['time']:>6.2f}s")
    print("=" * 60)

    return results


def codon_optimization_experiment():
    with open("config.json", "r") as f:
        config = json.load(f)

    data_file_name = config['metadata']['data_file']
    chunk_size = config['metadata']['chunk_size']
    sequences = util.split_sequence_from_file(data_file_name, chunk_size)

    all_results = []
    for i, seq in enumerate(sequences):
        print(f"\n{'#' * 60}")
        print(f"Sequence {i + 1}/{len(sequences)}")
        results = compare_solvers(seq, config)
        all_results.append({"sequence": seq, **results})

    # Aggregate
    n = len(all_results)
    print("\n" + "=" * 60)
    print("AGGREGATE RESULTS")
    print("=" * 60)
    for solver in ["dense", "sa"]:
        matches = sum(1 for r in all_results
                      if abs(r[solver]["energy"] - r["brute"]["energy"]) < 1e-6)
        avg_gap = sum(r[solver]["energy"] - r["brute"]["energy"] for r in all_results) / n
        avg_time = sum(r[solver]["time"] for r in all_results) / n
        print(f"{solver.upper():<8} optimal: {matches}/{n} | avg energy gap: {avg_gap:+.4f} | avg time: {avg_time:.2f}s")

    bf_avg = sum(r["brute"]["time"] for r in all_results) / n
    print(f"{'BRUTE':<8} (ground truth)                            | avg time: {bf_avg:.2f}s")
    print("=" * 60)


if __name__ == '__main__':
    codon_optimization_experiment()