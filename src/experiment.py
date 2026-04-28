import time
import json
import logging
from datetime import datetime

import util
import vqe_solver
import qaoa_solver
import sa_solver
import brute_force
import qiskit_util

from denseCodon import DenseCodon
from oneHotCodon import OneHotCodon
from hamiltonian import CodonOptimizer


# =========================================
# Logging Setup
# =========================================

def setup_logging():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = f"experiment_{timestamp}.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )


logger = logging.getLogger(__name__)


# =========================================
# Core Optimization
# =========================================

def build_optimizer(sequence, config, encoding_type):
    if encoding_type == "dense":
        encoding_cls = DenseCodon
    else:
        encoding_cls = OneHotCodon

    codon_opt = CodonOptimizer(
        sequence,
        config,
        encoding_cls,
        encoding_type
    )

    qubit_op = codon_opt.create_qubit_op()

    return codon_opt, qubit_op


def decode_result(sequence, bitstring, config, encoding_type):
    table_name = config['metadata']['table_name']

    if encoding_type == "dense":
        return util.decode_bitstring(
            sequence,
            bitstring,
            table_name=table_name
        )

    return util.decode_one_hot_bitstring(
        sequence,
        bitstring,
        table_name=table_name
    )


def run_optimization(sequence, config, type_opt="vqe", encoding_type="one-hot"):
    start_time = time.time()

    logger.info(f"Starting {type_opt.upper()} optimization")

    codon_opt, qubit_op = build_optimizer(
        sequence,
        config,
        encoding_type
    )

    try:
        if type_opt == "vqe":
            sampler = qiskit_util.get_backend(inject_noise=True)

            bitstring, energy = vqe_solver.get_min(
                qubit_op,
                config['vqe_settings'],
                sampler
            )

        elif type_opt == "qaoa":
            sampler = qiskit_util.get_backend(inject_noise=True)

            bitstring = qaoa_solver.get_min(
                qubit_op,
                config['qaoa_settings'],
                sampler
            )

            energy = util.evaluate_energy(qubit_op, bitstring)

        elif type_opt == "sa":
            bitstring, energy = sa_solver.get_min(
                qubit_op,
                codon_opt.qubit_len,
                config['sa_settings'],
                codon_opt.codon_list
            )

        elif type_opt == "brute":
            bitstring, energy = brute_force.brute_force_search(
                qubit_op,
                codon_opt.qubit_len
            )

        else:
            raise ValueError(f"Unknown type_opt: {type_opt}")

    except Exception:
        logger.exception(f"{type_opt.upper()} solver failed")
        raise

    final_mRNA_sequence = decode_result(
        sequence,
        bitstring,
        config,
        encoding_type
    )

    elapsed = time.time() - start_time

    logger.info(f"Sequence  : {sequence}")
    logger.info(f"Qubits    : {codon_opt.qubit_len}")
    logger.info(f"Bitstring : {bitstring}")
    logger.info(f"mRNA      : {final_mRNA_sequence}")
    logger.info(f"Energy    : {energy:.4f}")
    logger.info(f"Time      : {elapsed:.2f}s")

    return final_mRNA_sequence, energy, elapsed


# =========================================
# Solver Comparison
# =========================================

def compare_solvers(sequence, config, solver_list):
    logger.info("=" * 60)
    logger.info(f"Sequence: {sequence}")
    logger.info("=" * 60)

    results = {}

    for solver in solver_list:
        logger.info(f"[{solver.upper()}]")

        if solver == 'brute':
            encoding = "dense"
        else:
            encoding = "one-hot"

        mRNA, energy, elapsed = run_optimization(
            sequence,
            config,
            type_opt=solver,
            encoding_type=encoding
        )

        results[solver] = {
            "mRNA": mRNA,
            "energy": energy,
            "time": elapsed
        }

    reference_energy = results["brute"]["energy"]

    logger.info("=" * 60)
    logger.info("SUMMARY")

    for solver, res in results.items():
        gap = res["energy"] - reference_energy

        if abs(gap) < 1e-6:
            match = "✓"
        else:
            match = f"{gap:+.2f}"

        logger.info(
            f"{solver.upper():<10}"
            f"{res['energy']:>12.4f}"
            f"{match:>10}"
            f"{res['time']:>6.2f}s"
        )

    logger.info("=" * 60)

    return results


# =========================================
# Experiment Runner
# =========================================

def codon_optimization_experiment():
    with open("../config.json", "r") as f:
        config = json.load(f)

    dataset = util.parse_fasta(
        config['metadata']['dataset']
    )

    chunk_size = config['metadata']['chunk_size']
    output_path = config['metadata']['output_path']

    solver_list = ['vqe', 'sa', 'brute']

    for name, seq in dataset.items():

        logger.info(f"Start optimizing {name}")

        sequences = [
            seq[i:i + chunk_size]
            for i in range(
                0,
                len(seq),
                chunk_size
            )
        ]

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        output_name = (
            f"{output_path}"
            f"{name}_{timestamp}.txt"
        )

        codon_optimization_sequences(
            sequences,
            config,
            output_name,
            solver_list
        )


def codon_optimization_sequences(sequences, config, output_name, solver_list):
    all_results = []

    final_rna_strings = {
        solver: ""
        for solver in solver_list
    }

    for i, seq in enumerate(sequences):

        logger.info("#" * 60)
        logger.info(
            f"Sequence {i + 1}/"
            f"{len(sequences)}"
        )

        results = compare_solvers(seq, config, solver_list)

        all_results.append({
            "sequence": seq,
            **results
        })

        for solver in solver_list:
            final_rna_strings[solver] += results[solver]["mRNA"]

    write_output_file(
        output_name,
        all_results,
        final_rna_strings,
        solver_list
    )


# =========================================
# Output
# =========================================

def write_output_file(output_name, all_results, final_rna_strings, solver_list):
    with open(output_name, "w") as f:
        n = len(all_results)

        write_text("\n" + "=" * 60, f)
        write_text("AGGREGATE RESULTS", f)
        write_text("=" * 60, f)

        for solver in solver_list:
            matches = sum(
                1
                for r in all_results
                if abs(
                    r[solver]["energy"]
                    - r["brute"]["energy"]
                ) < 1e-6
            )

            avg_gap = sum(
                r[solver]["energy"]
                - r["brute"]["energy"]
                for r in all_results
            ) / n

            avg_time = sum(
                r[solver]["time"]
                for r in all_results
            ) / n

            write_text(
                f"{solver.upper():<8}"
                f" optimal: {matches}/{n}"
                f" | avg energy gap: {avg_gap:+.4f}"
                f" | avg time: {avg_time:.2f}s",
                f
            )


def write_text(text, file):
    logger.info(text)

    file.write(text + "\n")


# =========================================

if __name__ == '__main__':
    setup_logging()

    codon_optimization_experiment()
