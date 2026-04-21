import warnings
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel
from qiskit.providers.fake_provider import FakeManila, FakeSydney, FakeSherbrooke, FakeOslo, FakeTokyo
from qiskit.primitives import Sampler
from qiskit.primitives import BackendSampler
from qiskit import transpile
warnings.filterwarnings("ignore", category=UserWarning)


_FAKE_BACKENDS = {
    "FakeManila": FakeManila,
    "FakeSydney": FakeSydney,
    "FakeSherbrooke": FakeSherbrooke,
    "FakeOslo": FakeOslo,
    "FakeTokyo": FakeTokyo
}


def get_backend(backend_name='FakeOslo', inject_noise=False):
    if backend_name not in _FAKE_BACKENDS:
        raise ValueError("Please input the correct backend name.")

    fake_backend = _FAKE_BACKENDS[backend_name]()

    if inject_noise:
        noise_model = NoiseModel.from_backend(fake_backend)
        backend = AerSimulator(noise_model=noise_model)
        sampler = BackendSampler(backend)
    else:
        sampler = Sampler()
    return sampler


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
