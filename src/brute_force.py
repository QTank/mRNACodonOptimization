import util
import numpy as np


def prepare_operator_bitmask(qubit_op):
    """
    Convert Pauli Z operators into integer bitmasks.
    """

    qubit_op = util.normalize_operator(qubit_op)

    masks = []
    coeffs = []

    for pauli, coef in zip(qubit_op.paulis, qubit_op.coeffs):
        mask = 0
        z_index = np.where(pauli.z)[0]

        for idx in z_index:
            mask |= (1 << idx)

        masks.append(mask)                 # Python int
        coeffs.append(float(coef.real))   # Python float

    return masks, coeffs


def evaluate_energy_int(masks, coeffs, bitint):
    energy = 0.0

    for mask, coef in zip(masks, coeffs):
        parity_bits = bitint & mask   # Now always valid
        ones = parity_bits.bit_count()
        parity = -1 if (ones % 2) else 1
        energy += coef * parity

    return energy


# ---------- BRUTE FORCE ----------
def brute_force_search(qubit_op, qubit_len):
    masks, coeffs = prepare_operator_bitmask(qubit_op)
    min_energy = float("inf")
    best_int = 0

    total_states = 1 << qubit_len

    for i in range(total_states):
        energy = evaluate_energy_int(masks, coeffs, i)
        if energy < min_energy:
            min_energy = energy
            best_int = i

    best_string = format(best_int, f"0{qubit_len}b")
    return best_string, min_energy
