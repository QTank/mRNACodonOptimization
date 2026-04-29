import vqe_solver
import qaoa_solver
from denseCodon import DenseCodon

from oneHotCodon import OneHotCodon
import util
import time, json

from hamiltonian import CodonOptimizer
import sa_solver
import brute_force

def run_optimization(sequence, config, type_opt="dense"):
    start_time = time.time()
    if type_opt == "dense":
        print("Starting the VQE (Variational Quantum Eigensolver) on a Quantum Simulator...\n")
        codon_opt = CodonOptimizer(sequence, config, DenseCodon, "dense")
        qubit_op = codon_opt.create_qubit_op()
        bitstring, _ = vqe_solver.get_min(qubit_op, config['vqe_settings'])
        final_mRNA_sequence = util.decode_bitstring(sequence, bitstring,
                                               table_name=config['metadata']['table_name'])
    if type_opt == "one-hot":
        print("Starting QAOA on a Quantum Simulator...\n")
        codon_opt = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")
        qubit_op = codon_opt.create_qubit_op()
        bitstring = qaoa_solver.get_min(qubit_op, config['qaoa_settings'])
        final_mRNA_sequence = util.decode_one_hot_bitstring(sequence, bitstring,
                                                            table_name=config['metadata']['table_name'])
    
    if type_opt == "sa":
        if False:  # Set to True to print detailed results for SA
            print("Starting Simulated Annealing optimization...\n")
            codon_opt = CodonOptimizer(sequence, config, DenseCodon, "dense")
            qubit_op = codon_opt.create_qubit_op()
            bitstring, _ = sa_solver.get_min(qubit_op, codon_opt.qubit_len, config['sa_settings'])
            final_mRNA_sequence = util.decode_bitstring(sequence, bitstring,
                                               table_name=config['metadata']['table_name'])
        else:
            print("Starting Simulated Annealing optimization...\n")
            codon_opt = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")
            qubit_op = codon_opt.create_qubit_op()
            bitstring, _ = sa_solver.get_min(qubit_op, codon_opt.qubit_len, config['sa_settings'], codon_opt.codon_list)
            final_mRNA_sequence = util.decode_one_hot_bitstring(sequence, bitstring,
                                                                table_name=config['metadata']['table_name'])

    if type_opt == "brute-force":
        print("Starting brute-force optimization...\n")
        codon_opt = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")
        qubit_op = codon_opt.create_qubit_op()
        bitstring, _ = brute_force.get_min(qubit_op, codon_opt.codon_list)
        final_mRNA_sequence = util.decode_one_hot_bitstring(sequence, bitstring,
                                                            table_name=config['metadata']['table_name'])
        
    print_results = False
    if print_results:
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


def codon_optimization_experiment():
    with open("config.json", "r") as f:
        config = json.load(f)

    data_file_name = config['metadata']['data_file']
    chunk_size = config['metadata']['chunk_size']
    sequences = util.split_sequence_from_file(data_file_name, chunk_size)
    vqe_start_time = time.time()
    vqe_results = []
    qaoa_results = []
    sa_results = []
    for i, seq in enumerate(sequences):
        print(f"{i+1}/{len(sequences)} Optimizing sequence: {seq}")
        vqe_results.append(run_optimization(seq, config, type_opt="dense"))
    print(f"Total running time for VQE: {time.time() - vqe_start_time:.2f}s\n")
    
    qaoa_start_time = time.time()
    if False:
        for i, seq in enumerate(sequences):
            print(f"{i+1}/{len(sequences)} Optimizing sequence: {seq}")
            qaoa_results.append(run_optimization(seq, config, type_opt="one-hot"))
        print(f"Total running time for QAOA: {time.time() - qaoa_start_time:.2f}s\n")

    sa_start_time = time.time()
    for i, seq in enumerate(sequences):
        print(f"{i+1}/{len(sequences)} Optimizing sequence: {seq}")
        sa_results.append(run_optimization(seq, config, type_opt="sa"))
    print(f"Total running time for Simulated Annealing: {time.time() - sa_start_time:.2f}s\n")   

    print("Optimizations completed. Summary of results:")
    print(f"  - VQE Optimized mRNA Sequences: {''.join(vqe_results)}, time taken: {qaoa_start_time - vqe_start_time:.2f}s")
    print(f"  - QAOA Optimized mRNA Sequences: {''.join(qaoa_results)}, time taken: {sa_start_time - qaoa_start_time:.2f}s")
    print(f"  - Simulated Annealing Optimized mRNA Sequences: {''.join(sa_results)}, time taken: {time.time() - sa_start_time:.2f}s")


if __name__ == '__main__':
    codon_optimization_experiment()
