# Use this if you want to run on a fake/real backend
from logging import config

from qiskit_ibm_runtime import SamplerV2 as Sampler
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke
from qiskit.circuit.library import EfficientSU2
backend = FakeSherbrooke()
sampler = Sampler(mode=backend)

from qiskit import transpile
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA, SLSQP

def run_on_fake_device(hamiltonian, ansatz, sampler, backend):
    # 1. Transpile the ansatz for the fake backend's gate set
    # This is critical for accurate noise simulation
    transpiled_ansatz = transpile(ansatz, backend=backend)
    
    # 2. Setup VQE with the noisy sampler
    optimizer = SLSQP(maxiter=50)
    vqe = VQE(sampler=sampler, ansatz=transpiled_ansatz, optimizer=optimizer)
    
    # 3. Run
    result = vqe.compute_minimum_eigenvalue(hamiltonian)
    return result

from src.denseCodon import DenseCodon

from src.oneHotCodon import OneHotCodon
from src import util
import time, json

from src.hamiltonian import CodonOptimizer



from qiskit_algorithms import QAOA
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke
from qiskit_ibm_runtime import SamplerV2

def solve_with_qaoa_fake(qubit_op):
    # 1. Initialize Fake Backend
    backend = FakeSherbrooke()
    
    # 2. Setup the Sampler
    # In Qiskit 2.x, SamplerV2 is the optimized choice for backends
    sampler = SamplerV2(mode=backend)
    
    # 3. Hamiltonian & Optimizer

    
    optimizer = COBYLA(maxiter=50)
    
    # 4. Run QAOA
    # Note: reps=2 matches your original JSON 'qaoa_settings'
    qaoa = QAOA(sampler=sampler, optimizer=optimizer, reps=2)
    result = qaoa.compute_minimum_eigenvalue(qubit_op)
    
    return result

with open("config.json", "r") as f:
    config = json.load(f)
# Test sequence HAIHV
test_sequence = "HAIHV"
ansatz = EfficientSU2(
    su2_gates=['rx'],
    entanglement=config['vqe_settings'].get('entanglement', 'circular'),
    reps=config['vqe_settings'].get('reps', 2)
)
codon_opt = CodonOptimizer(test_sequence, config, DenseCodon, "dense")
qubit_op = codon_opt.create_qubit_op()
solve_with_qaoa_fake(qubit_op)
print("Running on a fake quantum device (simulating noise)...\n")
