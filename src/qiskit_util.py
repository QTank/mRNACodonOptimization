import warnings
from qiskit_aer.noise import NoiseModel
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke, FakeManilaV2, FakeSydneyV2, FakeOslo, FakeWashingtonV2
from qiskit_aer.primitives import SamplerV2 as Sampler
from qiskit import transpile
from qiskit_aer import AerSimulator

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

from denseCodon import DenseCodon
from oneHotCodon import OneHotCodon
from hamiltonian import CodonOptimizer
import util

_FAKE_BACKENDS = {
    "FakeManila": FakeManilaV2,
    "FakeSydney": FakeSydneyV2,
    "FakeSherbrooke": FakeSherbrooke,
    "FakeOslo": FakeOslo,
    "FakeWashington": FakeWashingtonV2,
}


def get_backend(backend_name='FakeSherbrooke', inject_noise=False):
    if backend_name not in _FAKE_BACKENDS:
        raise ValueError("Please input the correct backend name.")

    fake_backend = _FAKE_BACKENDS[backend_name]()

    if inject_noise:
        noise_model = NoiseModel.from_backend(fake_backend)
        aer_backend = AerSimulator(noise_model=noise_model)
        sampler = Sampler.from_backend(aer_backend)
        sampler.options.default_shots = 1024
    else:
        sampler = Sampler()
    return sampler, fake_backend


def analyze_circuit(circuit, backend_name):
    if backend_name not in _FAKE_BACKENDS:
        raise ValueError("Please input the correct backend name.")

    fake_backend = _FAKE_BACKENDS[backend_name]()

    transpiled_circuit = transpile(
        circuit,
        backend=fake_backend,
        optimization_level=3
    )

    ops = transpiled_circuit.count_ops()
    depth = transpiled_circuit.depth()

    return ops.values(), depth


def estimate_qubits(sequence, config, encode_type):
    if encode_type == 'dense':
        codon = CodonOptimizer(sequence, config, DenseCodon, encode_type)
    else:
        codon = CodonOptimizer(sequence, config, OneHotCodon, "one-hot")

    return codon.qubit_len


def estimate_qubits_from_protein(protein, chunk_size, config, encode_type):
    sequences = util.split_sequence(protein, chunk_size)

    qubits_list = []
    for seq in sequences:
        qubit_len = estimate_qubits(seq, config, encode_type)
        qubits_list.append(qubit_len)

    return qubits_list