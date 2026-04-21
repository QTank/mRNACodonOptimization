from hamiltonian import CodonOptimizer
from denseCodon import DenseCodon
from oneHotCodon import OneHotCodon
from sa_solver import get_min
import util

protein_sequence = "HAIHV"     # replace with your sequence

config = {
    "metadata": {
        "table_name": "e_coli_316407",
        "gc_content": 0.5,
    },
    "weights": {
        "usage": 0.3,
        "target_gc": 0.9,
        "repeated_nucleotides": 1.0,
        "redundant_encoding": 13000,
    },
}

sa_config = {
    "T_max": 50,        # starting temperature
    "T_min": 1e-2,       # stopping temperature
    "L": 50,            # Markov chain length per temperature step
    "max_stay_counter": 100,
    # "q": 0.95          # cooling rate (scikit-opt uses L to control cooling)
}

# ── 2. Choose encoding: DenseCodon (fewer qubits) or OneHotCodon ─────────────
codon_class = DenseCodon    # swap to OneHotCodon if preferred
codon_type  = "dense"       # or "one-hot"

# ── 3. Build the Hamiltonian ─────────────────────────────────────────────────
optimizer = CodonOptimizer(protein_sequence, config, codon_class, codon_type)
qubit_op  = optimizer.create_qubit_op()
qubit_len = optimizer.qubit_len
H_spare = qubit_op.primitive

print(f"Total qubits : {qubit_len}")
print(f"Qubit operator type:{type(H_spare)}")

# ── 4. Run SA ────────────────────────────────────────────────────────────────
best_bitstring, best_value = get_min(qubit_op, qubit_len, sa_config)

print(f"\nBest bitstring : {best_bitstring}")
print(f"Best energy    : {best_value:.4f}")

# ── 5. Decode to mRNA ────────────────────────────────────────────────────────
if codon_type == "dense":
    mrna = util.decode_bitstring(protein_sequence, best_bitstring)
else:
    mrna = util.decode_one_hot_bitstring(protein_sequence, best_bitstring)

print(f"mRNA sequence  : {mrna}")