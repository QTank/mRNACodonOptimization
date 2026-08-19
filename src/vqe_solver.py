from qiskit_algorithms.optimizers import COBYLA
from qiskit_algorithms.minimum_eigensolvers import SamplingVQE
from qiskit.primitives import StatevectorSampler as Sampler
from qiskit.circuit.library import EfficientSU2
import qiskit_util
from qiskit import transpile


def get_min(qubit_op, vqe_config, sampler=None, fake_backend=None):
    if sampler is None:
        sampler = Sampler()

    optimizer = COBYLA(maxiter=vqe_config.get('maxiter', 50))

    ansatz = EfficientSU2(
        num_qubits=qubit_op.num_qubits,
        su2_gates=['ry'],
        entanglement=vqe_config['ansatz'].get('entanglement', 'circular'),
        reps=vqe_config['ansatz'].get('reps', 2)
    ).decompose(reps=2)

    counts = []
    values = []

    def store_intermediate_result(eval_count, parameters, mean, std):
        counts.append(eval_count)
        values.append(mean)

    # initialize VQE using CVaR with alpha = 0.1 or 0.05 which can be set by user,
    # this method set alpha = 0.1
    vqe = SamplingVQE(
        sampler,
        ansatz=ansatz,
        optimizer=optimizer,
        aggregation=0.05,
        callback=store_intermediate_result,
    )
    raw_result = vqe.compute_minimum_eigenvalue(qubit_op)

    optimal_circuit = ansatz.assign_parameters(raw_result.optimal_parameters)

    transpiled = transpile(optimal_circuit, backend=fake_backend, optimization_level=1)
    print(f"parameters len: {len(raw_result.optimal_parameters)}")
    # Gate statistics
    statistics = [transpiled.count_ops(), sum(transpiled.count_ops().values()), transpiled.depth()]
    return raw_result.best_measurement['bitstring'], raw_result.best_measurement['value'].real, statistics
