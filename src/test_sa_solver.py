import vqe_solver
import qaoa_solver
from denseCodon import DenseCodon

from oneHotCodon import OneHotCodon
import util
import time, json

from hamiltonian import CodonOptimizer
import sa_solver

def test_sa_solver():
    with open("config.json", "r") as f:
        config = json.load(f)

    # Test sequence HAIHV
    test_sequence = "HAIHV"
    print("Testing Simulated Annealing Solver with sequence:", test_sequence)
    codon_opt = CodonOptimizer(test_sequence, config, DenseCodon, "dense")
    qubit_op = codon_opt.create_qubit_op()
    bitstring, energy = sa_solver.get_min(qubit_op, codon_opt.qubit_len, config['vqe_settings'])
    final_mRNA_sequence = util.decode_bitstring(test_sequence, bitstring,
                                               table_name=config['metadata']['table_name'])
    print("Results:")
    print(f"  - Input Amino Acid Sequence: {test_sequence}")
    print(f"  - The required number of qubits: {codon_opt.qubit_len}")
    print(f"  - Final Binary String Representation: {bitstring}")
    print(f"  - Optimized mRNA Codon Sequence: {final_mRNA_sequence}")
    print(f"  - Energy of the solution: {energy:.4f}\n")


if __name__ == '__main__':
    test_sa_solver()