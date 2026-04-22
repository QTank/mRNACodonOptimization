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


def run_optimization(sequence, config, type_opt="vqe", encoding_type="one-hot"):
    start_time = time.time()

    codon_opt, bitstring, energy = _solve(sequence, config, type_opt, encoding_type)
    final_mRNA = _decode(sequence, codon_opt, bitstring, encoding_type, config)

    elapsed = time.time() - start_time
    _print_result(sequence, codon_opt, bitstring, final_mRNA, energy, elapsed)

    return final_mRNA, energy, elapsed


def _solve(sequence, config, type_opt, encoding_type):
    """Build the qubit operator and run the requested solver."""
    encoding_class, enc_label = _get_encoding(type_opt, encoding_type)
    codon_opt = CodonOptimizer(sequence, config, encoding_class, enc_label)
    qubit_op = codon_opt.create_qubit_op()

    match type_opt:
        case "vqe":
            print("Starting VQE (Variational Quantum Eigensolver) on a Quantum Simulator...\n")
            sampler = qiskit_util.get_backend(inject_noise=True)
            bitstring, energy = vqe_solver.get_min(qubit_op, config["vqe_settings"], sampler)

        case "qaoa":
            print("Starting QAOA on a Quantum Simulator...\n")
            sampler = qiskit_util.get_backend(inject_noise=True)
            bitstring = qaoa_solver.get_min(qubit_op, config["qaoa_settings"], sampler)
            energy = util.evaluate_energy(qubit_op, bitstring)

        case "sa":
            print("Starting Simulated Annealing optimization...\n")
            bitstring, energy = sa_solver.get_min(
                qubit_op, codon_opt.qubit_len, config["sa_settings"], codon_opt.codon_list
            )

        case "brute":
            print("Starting Brute Force optimization...\n")
            bitstring, energy = brute_force.brute_force_search(qubit_op, codon_opt.qubit_len)

        case _:
            raise ValueError(f"Unknown type_opt: {type_opt}")

    return codon_opt, bitstring, energy


def _get_encoding(type_opt, encoding_type):
    """Return (encoding_class, label). VQE always uses dense; others respect encoding_type."""
    if type_opt == "vqe":
        return DenseCodon, "dense"
    if encoding_type == "dense":
        return DenseCodon, "dense"
    return OneHotCodon, "one-hot"


def _decode(sequence, codon_opt, bitstring, encoding_type, config):
    """Decode a bitstring back to an mRNA sequence."""
    table = config["metadata"]["table_name"]
    if isinstance(codon_opt.codon_list[0], DenseCodon):
        return util.decode_bitstring(sequence, bitstring, table_name=table)
    return util.decode_one_hot_bitstring(sequence, bitstring, table_name=table)


def _print_result(sequence, codon_opt, bitstring, final_mRNA, energy, elapsed):
    print(f"  - Sequence  : {sequence}")
    print(f"  - Qubits    : {codon_opt.qubit_len}")
    print(f"  - Bitstring : {bitstring}")
    print(f"  - mRNA      : {final_mRNA}")
    print(f"  - Energy    : {energy:.4f}")
    print(f"  - Time      : {elapsed:.2f}s\n")


def compare_solvers(sequence, config, solver_list):
    """Run all solvers on the same sequence and compare energies."""
    print("=" * 60)
    print(f"Sequence: {sequence}")
    print("=" * 60)

    results = {
        solver: dict(zip(("mRNA", "energy", "time"),
                         run_optimization(sequence, config, type_opt=solver, encoding_type="one-hot")))
        for solver in solver_list
    }

    reference_energy = results.get("brute", next(iter(results.values())))["energy"]

    print("=" * 60)
    print("SUMMARY")
    print(f"{'Solver':<10} {'mRNA':<25} {'Energy':>12} {'Gap':>10} {'Time':>8}")
    print("-" * 70)
    for solver, res in results.items():
        gap = res["energy"] - reference_energy
        match_str = "✓" if abs(gap) < 1e-6 else f"+{gap:.2f}"
        print(f"{solver.upper():<10} {str(res['mRNA']):<25} {res['energy']:>12.4f} {match_str:>10} {res['time']:>6.2f}s")
    print("=" * 60)

    return results


def codon_optimization_experiment():
    with open("../config.json") as f:
        config = json.load(f)

    data_file = config["metadata"]["data_file"]
    chunk_size = config["metadata"]["chunk_size"]
    sequences = util.split_sequence_from_file(data_file, chunk_size)

    solver_list = ["vqe", "sa", 'brute']
    all_results = []
    final_rna = {solver: "" for solver in solver_list}

    for i, seq in enumerate(sequences):
        print(f"\n{'#' * 60}")
        print(f"Sequence {i + 1}/{len(sequences)}")
        results = compare_solvers(seq, config, solver_list)
        all_results.append({"sequence": seq, **results})
        for solver in solver_list:
            final_rna[solver] += results[solver]["mRNA"]

    _print_aggregate(all_results, solver_list)

    for solver in solver_list:
        rna = final_rna[solver]
        print(f"{solver.upper()}, optimized RNA:\n  {rna}")
        print(f"{solver.upper()}, optimized DNA:\n  {util.convert_rna_to_dna(rna)}")


def _print_aggregate(all_results, solver_list):
    n = len(all_results)
    reference_key = "brute" if "brute" in solver_list else solver_list[0]

    print("\n" + "=" * 60)
    print("AGGREGATE RESULTS")
    print("=" * 60)

    for solver in solver_list:
        matches = sum(1 for r in all_results if abs(r[solver]["energy"] - r[reference_key]["energy"]) < 1e-6)
        avg_gap = sum(r[solver]["energy"] - r[reference_key]["energy"] for r in all_results) / n
        avg_time = sum(r[solver]["time"] for r in all_results) / n
        print(f"{solver.upper():<8} optimal: {matches}/{n} | avg energy gap: {avg_gap:+.4f} | avg time: {avg_time:.2f}s")

    bf_avg = sum(r[reference_key]["time"] for r in all_results) / n
    print(f"Reference {reference_key:<8} (ground truth)                        | avg time: {bf_avg:.2f}s")
    print("=" * 60)

if __name__ == '__main__':
    codon_optimization_experiment()