import vqe_solver
import qaoa_solver
from denseCodon import DenseCodon
from oneHotCodon import OneHotCodon
import util
import time, json
from hamiltonian import CodonOptimizer
import sa_solver
import brute_force
import qiskit_util
from datetime import datetime
import numpy as np


def run_optimization(sequence, config, type_opt="vqe", encoding_type="one-hot"):
    start_time = time.time()

    if type_opt == "vqe":
        write_text("Starting VQE (Variational Quantum Eigensolver) on a Quantum Simulator...\n")
        codon_opt = CodonOptimizer(sequence, config, DenseCodon, "dense")
        qubit_op = codon_opt.create_qubit_op()
        sampler = qiskit_util.get_backend(inject_noise=True)
        bitstring, energy = vqe_solver.get_min(qubit_op, config['vqe_settings'], sampler)
        final_mRNA_sequence = util.decode_bitstring(sequence, bitstring,
                                                    table_name=config['metadata']['table_name'])

    elif type_opt == "qaoa":
        write_text("Starting QAOA on a Quantum Simulator...\n")
        codon_opt = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")
        qubit_op = codon_opt.create_qubit_op()
        sampler = qiskit_util.get_backend(inject_noise=True)
        bitstring = qaoa_solver.get_min(qubit_op, config['qaoa_settings'], sampler)
        energy = util.evaluate_energy(qubit_op, bitstring)
        final_mRNA_sequence = util.decode_one_hot_bitstring(sequence, bitstring,
                                                            table_name=config['metadata']['table_name'])

    elif type_opt == "sa":
        write_text("Starting Simulated Annealing optimization...\n")
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
        write_text("Starting Brute Force optimization...\n")
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

    write_text(f"  - Sequence  : {sequence}")
    if type_opt in ['vqe', 'qaoa']:
        write_text(f"  - Qubits    : {codon_opt.qubit_len}")
    write_text(f"  - Bitstring : {bitstring}")
    write_text(f"  - mRNA      : {final_mRNA_sequence}")
    write_text(f"  - Energy    : {energy:.4f}")
    write_text(f"  - Time      : {elapsed:.2f}s\n")

    return final_mRNA_sequence, energy, elapsed


def compare_solvers(sequence, config, solver_list, f):
    """Run VQE, SA, and Brute Force on the same sequence and compare energy."""
    write_text("=" * 60, f)
    write_text(f"Sequence: {sequence}", f)
    write_text("=" * 60, f)

    results = {}
    for solver in solver_list:
        write_text(f"[{solver.upper()}]", f)
        if solver == 'brute':
            mRNA, energy, elapsed = run_optimization(sequence, config, type_opt=solver, encoding_type='dense')
        else:
            mRNA, energy, elapsed = run_optimization(sequence, config, type_opt=solver, encoding_type='one-hot')
        results[solver] = {"mRNA": mRNA, "energy": energy, "time": elapsed}

    reference_energy = results["brute"]["energy"]

    write_text("=" * 60, f)
    write_text("SUMMARY", f)
    write_text(f"{'Solver':<10} {'mRNA':<25} {'Energy':>12} {'Gap':>10} {'Time':>8}", f)
    write_text("-" * 70, f)
    for solver, res in results.items():
        gap = res["energy"] - reference_energy
        match = "✓" if abs(gap) < 1e-6 else f"+{gap:.2f}"
        write_text(
            f"{solver.upper():<10} {str(res['mRNA']):<25} {res['energy']:>12.4f} {match:>10} {res['time']:>6.2f}s")
    write_text("=" * 60)

    return results


def codon_optimization_experiment(data_file_name):
    with open("../config.json", "r") as f:
        config = json.load(f)

    output_path = config['metadata']['output_path']
    chunk_size = config['metadata']['chunk_size']
    data = util.parse_sequence_from_file(data_file_name)
    sequences = util.split_sequence(data, chunk_size)
    for i, seq in enumerate(sequences):
        print(i+1, seq)
        if False:
            sequences = [seq[i:i + chunk_size] for i in range(0, len(seq), chunk_size)]
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_name = f"{output_path}test_{data_file_name}_{timestamp}.txt"
            codon_optimization_sequences(sequences, config, output_name)


def codon_optimization_sequences(sequences, config, output_name):
    all_results = []
    solver_list = ['vqe', 'sa', 'brute']
    final_rna_strings = {solver: "" for solver in solver_list}
    with open(output_name, "w") as f:
        for i, seq in enumerate(sequences):
            write_text(f"\n{'#' * 60}")
            write_text(f"Sequence {i + 1}/{len(sequences)}")
            results = compare_solvers(seq, config, solver_list, f)
            all_results.append({"sequence": seq, **results})

            for solver in solver_list:
                final_rna_strings[solver] += results[solver]["mRNA"]

        # Aggregate
        n = len(all_results)
        write_text(f"The length of fragment: {len(sequences[0])}", f)
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
    codon_optimization_experiment("../data/02-sars2_n_vaccine.fasta")
