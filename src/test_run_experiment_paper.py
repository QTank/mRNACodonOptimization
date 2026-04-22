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
import qiskit_util


def run_optimization(sequence, config, type_opt="vqe", encoding_type="one-hot"):
    start_time = time.time()

    if type_opt == "vqe":
        print("Starting VQE (Variational Quantum Eigensolver) on a Quantum Simulator...\n")
        codon_opt = CodonOptimizer(sequence, config, DenseCodon, "dense")
        qubit_op = codon_opt.create_qubit_op()
        sampler = qiskit_util.get_backend(inject_noise=True)
        bitstring, energy = vqe_solver.get_min(qubit_op, config['vqe_settings'], sampler)
        final_mRNA_sequence = util.decode_bitstring(sequence, bitstring,
                                                    table_name=config['metadata']['table_name'])

    elif type_opt == "qaoa":
        print("Starting QAOA on a Quantum Simulator...\n")
        codon_opt = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")
        qubit_op = codon_opt.create_qubit_op()
        sampler = qiskit_util.get_backend(inject_noise=True)
        bitstring = qaoa_solver.get_min(qubit_op, config['qaoa_settings'], sampler)
        energy = util.evaluate_energy(qubit_op, bitstring)
        final_mRNA_sequence = util.decode_one_hot_bitstring(sequence, bitstring,
                                                            table_name=config['metadata']['table_name'])

    elif type_opt == "sa":
        print("Starting Simulated Annealing optimization...\n")
        if encoding_type == "dense":
            codon_opt = CodonOptimizer(sequence, config, DenseCodon, "dense")
            qubit_op = codon_opt.create_qubit_op()
            bitstring, energy = sa_solver.get_min(qubit_op, codon_opt.qubit_len, config['sa_settings'],
                                                  codon_opt.codon_list)
            final_mRNA_sequence = util.decode_bitstring(sequence, bitstring,
                                                        table_name=config['metadata']['table_name'])
        else:
            codon_opt = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")
            qubit_op = codon_opt.create_qubit_op()
            bitstring, energy = sa_solver.get_min(qubit_op, codon_opt.qubit_len, config['sa_settings'],
                                                  codon_opt.codon_list)
            final_mRNA_sequence = util.decode_one_hot_bitstring(sequence, bitstring,
                                                                table_name=config['metadata']['table_name'])

    elif type_opt == "brute":
        print("Starting Brute Force optimization...\n")
        if encoding_type == "dense":
            codon_opt = CodonOptimizer(sequence, config, DenseCodon, "dense")
            qubit_op = codon_opt.create_qubit_op()
            bitstring, energy = brute_force.brute_force_search(qubit_op, codon_opt.qubit_len)
            final_mRNA_sequence = util.decode_bitstring(sequence, bitstring,
                                                        table_name=config['metadata']['table_name'])
        else:
            codon_opt = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")
            qubit_op = codon_opt.create_qubit_op()
            bitstring, energy = brute_force.brute_force_search(qubit_op, codon_opt.qubit_len)
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

    if False:
        final_string_index = 0
        for index in range(len(sequence)):
            amino_acid = sequence[index]
            chunk = bitstring[final_string_index:final_string_index + codon_opt.codon_list[index].encoding_qubit_len]
            codon = final_mRNA_sequence[index * 3: (index + 1) * 3] if final_mRNA_sequence else "???"
            print(f"    {amino_acid} → {chunk} → {codon}")
            final_string_index += codon_opt.codon_list[index].encoding_qubit_len

    return final_mRNA_sequence, energy, elapsed


def compare_solvers(sequence, config, solver_list):
    """Run VQE, SA, and Brute Force on the same sequence and compare energy."""
    print("=" * 60)
    print(f"Sequence: {sequence}")
    print("=" * 60)

    results = {}
    for solver in solver_list:
        print(f"[{solver.upper()}]")
        if solver == 'brute':
            mRNA, energy, elapsed = run_optimization(sequence, config, type_opt=solver, encoding_type='dense')
        else:
            mRNA, energy, elapsed = run_optimization(sequence, config, type_opt=solver, encoding_type='one-hot')
        results[solver] = {"mRNA": mRNA, "energy": energy, "time": elapsed}

    reference_energy = results["brute"]["energy"]

    print("=" * 60)
    print("SUMMARY")
    print(f"{'Solver':<10} {'mRNA':<25} {'Energy':>12} {'Gap':>10} {'Time':>8}")
    print("-" * 70)
    for solver, res in results.items():
        gap = res["energy"] - reference_energy
        match = "✓" if abs(gap) < 1e-6 else f"+{gap:.2f}"
        print(f"{solver.upper():<10} {str(res['mRNA']):<25} {res['energy']:>12.4f} {match:>10} {res['time']:>6.2f}s")
    print("=" * 60)

    return results


def codon_optimization_experiment():
    with open("../config.json", "r") as f:
        config = json.load(f)

    data_file_name = config['metadata']['dataset']
    chunk_size = config['metadata']['chunk_size']
    dataset_seq = util.parse_fasta(data_file_name)
    for k, seq in dataset_seq.items():
        if k != "V17_HepatitisB_Surface": continue
        print(f"start optimizing {k}, splite into fragments with the length of {chunk_size}")
        sequences = [seq[i:i + chunk_size] for i in range(0, len(seq), chunk_size)]
        output_name = seq + ".txt"
        codon_optimization_sequences(sequences, config, output_name)


def codon_optimization_sequences(sequences, config, output_name):
    all_results = []
    solver_list = ['vqe', 'sa', 'brute']
    final_rna_strings = {solver: "" for solver in solver_list}

    for i, seq in enumerate(sequences):
        print(f"\n{'#' * 60}")
        print(f"Sequence {i + 1}/{len(sequences)}")
        results = compare_solvers(seq, config, solver_list)
        all_results.append({"sequence": seq, **results})

        for solver in solver_list:
            final_rna_strings[solver] += results[solver]["mRNA"]

    with open(output_name, "w") as f:
        # Aggregate
        n = len(all_results)
        write_text("\n" + "=" * 60, f)
        write_text("AGGREGATE RESULTS", f)
        write_text("=" * 60, f)
        for solver in solver_list:
            matches = sum(1 for r in all_results if abs(r[solver]["energy"] - r["brute"]["energy"]) < 1e-6)
            avg_gap = sum(r[solver]["energy"] - r["brute"]["energy"] for r in all_results) / n
            avg_time = sum(r[solver]["time"] for r in all_results) / n
            write_text(
                f"{solver.upper():<8} optimal: {matches}/{n} | avg energy gap: {avg_gap:+.4f} | avg time: {avg_time:.2f}s",
                f)

        bf_avg = sum(r["brute"]["time"] for r in all_results) / n
        write_text(f"{'BRUTE':<8} (ground truth)                            | avg time: {bf_avg:.2f}s", f)
        write_text("=" * 60, f)

        for solver in solver_list:
            rna = final_rna_strings[solver]
            write_text(f"{solver.upper()}, optimized RNA:\n {rna}", f)
            write_text(f"{solver.upper()}, optimized DNA:\n "f"{util.convert_rna_to_dna(rna)}", f)


def write_text(text, file):
    print(text)
    file.write(text + "\n")


if __name__ == '__main__':
    codon_optimization_experiment()
