from qiskit.algorithms.optimizers import COBYLA
from qiskit.algorithms.minimum_eigensolvers import SamplingVQE
from qiskit.primitives import Sampler
from qiskit.circuit.library import EfficientSU2


def get_min(qubit_op, vqe_config):
    optimizer = COBYLA(maxiter=vqe_config.get('maxiter', 50))
    ansatz = EfficientSU2(
        su2_gates=['rx'],
        entanglement=vqe_config.get('entanglement', 'circular'),
        reps=vqe_config.get('reps', 2)
    )

    counts = []
    values = []

    def store_intermediate_result(eval_count, parameters, mean, std):
        counts.append(eval_count)
        values.append(mean)

    # initialize VQE using CVaR with alpha = 0.1 or 0.05 which can be set by user,
    # this method set alpha = 0.1
    vqe = SamplingVQE(
        Sampler(),
        ansatz=ansatz,
        optimizer=optimizer,
        aggregation=vqe_config.get('alpha', 0.1),
        callback=store_intermediate_result,
    )
    raw_result = vqe.compute_minimum_eigenvalue(qubit_op)
    # We only need the final qubit and the lowest value
    return raw_result.best_measurement['bitstring'], raw_result.best_measurement['value']
