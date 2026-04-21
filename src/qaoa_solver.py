from qiskit_algorithms import QAOA
from qiskit_algorithms.optimizers import COBYLA
from qiskit.primitives import Sampler
from qiskit.circuit.library import EfficientSU2
import numpy as np


def get_min(qubit_op, qaoa_config):
    optimizer = COBYLA(maxiter=qaoa_config.get('maxiter', 50))
    ansatz = EfficientSU2(
        su2_gates=['rx'],
        entanglement=qaoa_config['ansatz'].get('entanglement', 'circular'),
        reps=qaoa_config['ansatz'].get('reps', 2)
    )

    qaoa = QAOA(sampler=Sampler(), optimizer=optimizer, reps=qaoa_config.get('reps', 2))
    qaoa.ansatz = ansatz
    result = qaoa.compute_minimum_eigenvalue(qubit_op)

    if hasattr(result, 'best_measurement'):
        best_bitstring = result.best_measurement['bitstring']
    else:
        if hasattr(result, 'eigenstate'):
            state = result.eigenstate
            probabilities = state.probabilities()
            best_index = np.argmax(probabilities)
            best_bitstring = format(best_index, f'0{qubit_op.len}b')
        else:
            print("No solution found in result")
            return None

    return best_bitstring