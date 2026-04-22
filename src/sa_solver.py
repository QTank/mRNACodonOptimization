"""
sa_solver.py
------------
Classical solver for codon optimization Hamiltonians using
Simulated Annealing (scikit-opt).

Acts as a drop-in classical replacement for VQE / QAOA solvers.
"""

from __future__ import annotations
import util
import numpy as np

# ---------------------------------------------------------------------------
# Bitstring Utilities
# ---------------------------------------------------------------------------

def array_to_bitstring(x: np.ndarray) -> str:
    """
    Convert continuous vector to binary bitstring.

    Example:
        [0.6, 0.2, 0.9] -> '101'
    """
    return "".join("1" if v >= 0.5 else "0" for v in x)


def bitstring_to_array(bitstring: str) -> np.ndarray:
    """
    Convert bitstring to numpy array.

    Example:
        '0110' -> [0,1,1,0]
    """
    return np.array([float(b) for b in bitstring], dtype=float)

# ---------------------------------------------------------------------------
# Simulated Annealing Solver
# ---------------------------------------------------------------------------

DEFAULT_SA_CONFIG = {
    "T_max": 100,
    "T_min": 1e-7,
    "L": 300,
    "max_stay_counter": 150,
}


def get_min(
    qubit_op,
    qubit_len: int,
    sa_config: dict | None = None, codon_list=None
):
    """
    Minimise Hamiltonian using Simulated Annealing.

    Parameters
    ----------
    qubit_op : PauliSumOp | SparsePauliOp
        Hamiltonian.
    qubit_len : int
        Number of qubits.
    sa_config : dict | None
        Optional SA parameters.

    Returns
    -------
    tuple[str, float]
        (best_bitstring, best_energy)
    """

    # Merge configs
    config = DEFAULT_SA_CONFIG.copy()

    if sa_config:
        config.update(sa_config)

    # Objective function
    groups = []
    idx = 0
    for codon in codon_list:
        length = codon.encoding_qubit_len
        groups.append((idx, idx + length))
        idx += length

    def random_valid_state():
        """Initialize with a valid one-hot state per group."""
        bits = [0] * qubit_len
        for start, end in groups:

            # Fix invalid group bounds
            start = max(0, start)
            end = min(qubit_len, end)

            if start >= qubit_len:
                continue  # skip invalid group

            if end - start <= 1:
                chosen = start
            else:
                chosen = np.random.randint(start, end)

            bits[chosen] = 1
        return bits

    def neighbour(current):
        """Pick a random group and swap the active bit to another position."""
        bits = current.copy()
        swappable = [(start, end) for start, end in groups if end - start > 1]
        if not swappable:
            return bits  # nothing to swap
        
        start, end = swappable[np.random.randint(len(swappable))]
        active = next(i for i in range(start, end) if bits[i] == 1)
        choices = [i for i in range(start, end) if i != active]
        new_active = np.random.choice(choices)
        bits[active] = 0
        bits[new_active] = 1
        return bits

    current = random_valid_state()
    current_bs = "".join(map(str, current))
    current_energy = util.evaluate_energy(qubit_op, current_bs)
    best, best_energy = list(current), current_energy

    T = config["T_max"]
    T_min = config["T_min"]
    L = config["L"]
    alpha = (T_min / T) ** (1.0 / L)

    for _ in range(L):
        for _ in range(qubit_len):
            nbr = neighbour(current)
            nbr_bs = "".join(map(str, nbr))
            e = util.evaluate_energy(qubit_op, nbr_bs)
            delta = e - current_energy
            if delta < 0 or np.random.rand() < np.exp(-delta / T):
                current, current_energy = nbr, e
                if e < best_energy:
                    best, best_energy = list(nbr), e
        T *= alpha

    return "".join(map(str, best)), float(best_energy)