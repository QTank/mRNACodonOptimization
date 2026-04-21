import numpy as np
from typing import Set
from qiskit.opflow import PauliOp, I, Z
import python_codon_tables as pct
from qiskit.quantum_info import SparsePauliOp


def build_full_identity(num_qubits: int) -> PauliOp:
    full_identity = I
    for _ in range(1, num_qubits):
        full_identity = I ^ full_identity
    return full_identity


def build_pauli_z_op(num_qubits: int, pauli_z_indices: Set[int]) -> PauliOp:
    if 0 in pauli_z_indices:
        operator = Z
    else:
        operator = I
    for i in range(1, num_qubits):
        if i in pauli_z_indices:
            operator = operator ^ Z
        else:
            operator = operator ^ I

    return operator


def build_indicator_qubit(qubit_len: int, pauli_z_index: int) -> PauliOp:
    return 0.5 * build_full_identity(qubit_len) - 0.5 * build_pauli_z_op(qubit_len, {pauli_z_index})


def int_to_binary(index: int, length: int) -> str:
    return f"{index:0{length}b}"


def get_gc_count(codon: str) -> int:
    return codon.count('G') + codon.count('C')


def get_neighbour_repetition(codon1, codon2):
    if not codon1 or not codon2:
        return 0

    last_base = codon1[-1]
    if codon2[0] != last_base:
        return 0

    left_count = 1
    for i in range(len(codon1) - 2, -1, -1):
        if codon1[i] == last_base:
            left_count += 1
        else:
            break

    right_count = 1
    for base in codon2[1:]:
        if base == last_base:
            right_count += 1
        else:
            break

    total_repeat = left_count + right_count
    return total_repeat


def decode_bitstring(protein_sequence, bitstring, table_name='e_coli_316407'):
    mRNA_sequence = ""
    start_index = 0
    codon_table = pct.get_codons_table(table_name)
    try:
        for amino in protein_sequence:
            codon_list = list(codon_table[amino].keys())
            encoding_len = round(np.log2(len(codon_list)))

            if encoding_len == 0:
                mRNA_sequence += codon_list[0]
            else:
                loc = int(bitstring[start_index:start_index + encoding_len], 2)
                start_index += encoding_len
                mRNA_sequence += codon_list[loc]
    except:
        print(f"the bitstring {bitstring} is invalid!")

    return mRNA_sequence.replace("T", "U")


def decode_one_hot_bitstring(protein_sequence, bitstring, table_name='e_coli_316407'):
    mRNA_sequence = ""
    start_index = 0
    codon_table = pct.get_codons_table(table_name)
    try:
        for amino in protein_sequence:
            codon_list = list(codon_table[amino].keys())
            encoding_len = len(codon_list)
            valid_index = [2 ** i for i in range(encoding_len)]

            loc = int(bitstring[start_index:start_index + encoding_len], 2)
            if loc not in valid_index:
                print(f"the bitstring {bitstring} is invalid!")
                return mRNA_sequence

            codon_index = bitstring[start_index:start_index + encoding_len].find("1")
            mRNA_sequence += codon_list[codon_index]
            start_index += encoding_len
    except:
        print(f"the bitstring {bitstring} is invalid!")

    return mRNA_sequence.replace("T", "U")




def split_sequence_from_file(filename, chunk_size):
    with open(filename, "r") as f:
        sequence = f.read().strip()

    sequence = sequence.replace("\n", "")

    chunks = [
        sequence[i:i+chunk_size]
        for i in range(0, len(sequence), chunk_size)
    ]

    return chunks


def parse_fasta(file_path):
    sequences = {}
    current_name = None
    current_seq = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                # save previous record
                if current_name is not None:
                    sequences[current_name] = "".join(current_seq)

                # start new record
                current_name = line[1:]  # remove ">"
                current_seq = []
            else:
                current_seq.append(line)

        # save last record
        if current_name is not None:
            sequences[current_name] = "".join(current_seq)

    return sequences


def convert_rna_to_dna(rna_sequence):
    return rna_sequence.replace("U", "T")


def evaluate_energy(qubit_op, bitstring):
    energy = 0
    n = len(bitstring)
    qubit_op = normalize_operator(qubit_op)

    for pauli, coef in zip(qubit_op.paulis, qubit_op.coeffs):
        z_index = np.where(pauli.z)[0]
        parity = 1.0

        for idx in z_index:
            val = -1 if bitstring[n - 1 - idx] == '1' else 1
            parity *= val

        energy += coef.real * parity
    return energy


def normalize_operator(qubit_op):
    if isinstance(qubit_op, SparsePauliOp):
        return qubit_op

    try:
        return qubit_op.primitive
    except Exception:
        raise RuntimeError("Failed to convert qubit_op to SparsePauliOp.")