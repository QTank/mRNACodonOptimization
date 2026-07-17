from qiskit_algorithms import QAOA
from qiskit_algorithms.optimizers import COBYLA
from qiskit.primitives import StatevectorSampler as Sampler
from qiskit.circuit.library import EfficientSU2
import numpy as np
from qiskit_aer import AerSimulator
from qiskit import transpile



class QAOACustomAnsatz(QAOA):
    """QAOA subclass that allows a fixed, user-supplied ansatz instead of
    the auto-generated QAOAAnsatz."""
    def _check_operator_ansatz(self, operator):
        # Skip QAOA's forced rebuild; keep whatever ansatz was assigned.
        pass

def get_min(qubit_op, qaoa_config, sampler=None):
    if sampler is None:
        sampler = Sampler()

    optimizer = COBYLA(maxiter=qaoa_config.get('maxiter', 50))
    ansatz = EfficientSU2(
        num_qubits=qubit_op.num_qubits,
        su2_gates=['rx'],
        entanglement=qaoa_config['ansatz'].get('entanglement', 'circular'),
        reps=qaoa_config['ansatz'].get('reps', 2)
    ).decompose(reps=3)


    qaoa = QAOACustomAnsatz(sampler=sampler, optimizer=optimizer, reps=qaoa_config.get('reps', 2))
    qaoa.ansatz = ansatz
    result = qaoa.compute_minimum_eigenvalue(qubit_op)

    optimal_circuit = qaoa.ansatz.assign_parameters(result.optimal_parameters)
    backend = AerSimulator()

    # Transpile for backend
    transpiled = transpile(optimal_circuit, backend=backend, optimization_level=2)

    # Gate statistics
    statistics = [transpiled.count_ops(), sum(transpiled.count_ops().values()), transpiled.depth()]

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
    print(f"parameters len: {len(result.optimal_parameters)}")
    return best_bitstring, statistics
