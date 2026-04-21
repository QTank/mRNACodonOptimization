import itertools
import numpy as np

def evaluate_energy(qubit_op, bitstring: str) -> float:
    """
    Evaluate energy of a computational basis state.

    Assumes Hamiltonian is diagonal in Z-basis.
    Works with SparsePauliOp.

    Parameters
    ----------
    qubit_op : SparsePauliOp
        Hamiltonian operator.
    bitstring : str
        Binary string state.

    Returns
    -------
    float
        Energy value.
    """
    n = len(bitstring)

    energy = 0.0

    for pauli, coef in zip(
        qubit_op.paulis,
        qubit_op.coeffs
    ):

        # Indices where Z acts
        z_indices = np.where(pauli.z)[0]
        parity = 1.0

        for idx in z_indices:
            val = 1 if bitstring[n - 1 - idx] == '1' else -1
            parity *= val

        energy += coef.real * parity

    return energy


def normalize_operator(qubit_op):
    """
    Ensure operator is SparsePauliOp.
    """

    from qiskit.quantum_info import SparsePauliOp

    if isinstance(qubit_op, SparsePauliOp):
        return qubit_op

    try:
        return qubit_op.primitive
    except Exception as exc:
        raise RuntimeError(
            "Failed to convert qubit_op to SparsePauliOp."
        ) from exc
    
def get_min(qubit_op, codon_list):
    """
    Brute force solver: exhaustively checks all valid one-hot states.
    
    Parameters
    ----------
    qubit_op : SparsePauliOp
        Hamiltonian.
    codon_list : list
        List of OneHotCodon objects defining group boundaries.
    
    Returns
    -------
    tuple[str, float]
        (best_bitstring, best_energy)
    """
    sparse_op = normalize_operator(qubit_op)
    qubit_len = sum(c.encoding_qubit_len for c in codon_list)

    # Build group slices
    groups = []
    idx = 0
    for codon in codon_list:
        length = codon.encoding_qubit_len
        groups.append((idx, idx + length))  # one-hot position choices
        idx += length

    best_bitstring = None
    best_energy = float('inf')

    # Iterate over all combinations: for each group, pick which bit is active
    group_sizes = [end - start for start, end in groups]
    
    for combo in itertools.product(*[range(s) for s in group_sizes]):
        # Build bitstring from this combination
        bits = [0] * qubit_len
        for group_idx, active_offset in enumerate(combo):
            start, _ = groups[group_idx]
            bits[start + active_offset] = 1
        
        bitstring = "".join(map(str, bits))
        energy = evaluate_energy(sparse_op, bitstring)
        
        if energy < best_energy:
            best_energy = energy
            best_bitstring = bitstring

    return best_bitstring, float(best_energy)