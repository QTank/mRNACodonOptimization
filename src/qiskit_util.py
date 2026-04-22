import warnings
from qiskit_aer.noise import NoiseModel
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke, FakeManila, FakeSydney, FakeOslo, FakeTokyo
from qiskit_aer.primitives import Sampler
from qiskit import transpile
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

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
        sampler = Sampler(
            backend_options={
                "noise_model": noise_model
            },
            run_options={
                "shots": 1024
            }
        )
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
