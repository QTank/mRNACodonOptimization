from qiskit_algorithms import QAOA
from qiskit_algorithms.optimizers import COBYLA
from qiskit.primitives import StatevectorSampler as Sampler
import numpy as np
from qiskit import transpile
from qiskit.circuit.library import QAOAAnsatz

class QAOACustomAnsatz(QAOA):
    """QAOA subclass that allows a fixed, user-supplied ansatz instead of
    the auto-generated QAOAAnsatz."""
    def _check_operator_ansatz(self, operator):
        # Skip QAOA's forced rebuild; keep whatever ansatz was assigned.
        pass

def get_min(qubit_op, qaoa_config, sampler=None, fake_backend=None):
    if sampler is None:
        sampler = Sampler()
    reps = 2
    optimizer = COBYLA(maxiter=qaoa_config.get('maxiter', 50))
    ansatz = QAOAAnsatz(
        cost_operator=qubit_op,
        reps=reps
    ).decompose(reps=reps)


    qaoa = QAOACustomAnsatz(sampler=sampler, optimizer=optimizer, reps=qaoa_config.get('reps', 2))
    qaoa.ansatz = ansatz
    result = qaoa.compute_minimum_eigenvalue(qubit_op)

    optimal_circuit = ansatz.assign_parameters(result.optimal_parameters)
    # Transpile for backend
    if fake_backend is not None and qubit_op.num_qubits > fake_backend.configuration().n_qubits:
        print(f"Warning: {qubit_op.num_qubits} qubits exceeds backend capacity "
              f"({fake_backend.configuration().n_qubits}); transpiling without backend constraints.")
        transpiled = transpile(optimal_circuit, optimization_level=1)
    else:
        transpiled = transpile(optimal_circuit, backend=fake_backend, optimization_level=1)
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
