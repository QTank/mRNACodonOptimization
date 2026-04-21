from src import vqe_solver
from src import qaoa_solver
from src.denseCodon import DenseCodon

from src.oneHotCodon import OneHotCodon
from src import util
import time, json

from src.hamiltonian import CodonOptimizer


def run_optimization(sequence, config, type_opt="dense"):
    start_time = time.time()
    if type_opt == "dense":
        print("Starting the VQE (Variational Quantum Eigensolver) on a Quantum Simulator...\n")
        codon_opt = CodonOptimizer(sequence, config, DenseCodon, "dense")
        qubit_op = codon_opt.create_qubit_op()
        bitstring, _ = vqe_solver.get_min(qubit_op, config['vqe_settings'])
        final_mRNA_sequence = util.decode_bitstring(sequence, bitstring,
                                               table_name=config['metadata']['table_name'])
    else:
        print("Starting QAOA on a Quantum Simulator...\n")
        codon_opt = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")
        qubit_op = codon_opt.create_qubit_op()
        bitstring = qaoa_solver.get_min(qubit_op, config['qaoa_settings'])
        final_mRNA_sequence = util.decode_one_hot_bitstring(sequence, bitstring,
                                                            table_name=config['metadata']['table_name'])
    print("Results:")
    print(f"  - Input Amino Acid Sequence: {sequence}")
    print(f"  - The required number of qubit: {codon_opt.qubit_len}")
    print(f"  - Final Binary String Representation: {bitstring}")
    print(f"  - Optimized mRNA Codon Sequence: {final_mRNA_sequence}")
    print(f"  - Running time: {time.time() - start_time:.2f}s \n")
    print("Optimized Codon Mapping:")
    print("  Each amino acid is mapped to its optimized codon:")
    final_string_index = 0
    for index in range(len(sequence)):
        amino_acid = sequence[index]
        print(f"    {amino_acid:} → {bitstring[final_string_index:final_string_index + codon_opt.codon_list[index].encoding_qubit_len]} → {final_mRNA_sequence[index * 3: (index+1) * 3]}")
        final_string_index += codon_opt.codon_list[index].encoding_qubit_len

    return final_mRNA_sequence



if __name__ == '__main__':
    with open("config.json", "r") as f:
        config = json.load(f)

    # Test sequence HAIHV
    test_sequence = "HAIHV"
    final_sequence = run_optimization(test_sequence, config, "dense")
    final_sequence = run_optimization(test_sequence, config, "one-hot")



