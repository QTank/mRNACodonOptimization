import re
import numpy as np


def run(text):
    results = []

    # Split each sequence block
    sequence_blocks = re.split(r"#{10,}", text)

    for seq_block in sequence_blocks:

        # Get sequence name
        seq_match = re.search(r"Sequence:\s*(\w+)", seq_block)

        if not seq_match:
            continue

        sequence = seq_match.group(1)

        # Find each algorithm block separately
        algo_blocks = re.findall(r"\[(VQE|QAOA)\](.*?)(?=\n\[|\Z)",
                                 seq_block, re.DOTALL)

        for algorithm, block in algo_blocks:

            depth_match = re.search(r"Circuit depth:\s*(\d+)", block)
            gates_match = re.search(r"Total gates\s*:\s*(\d+)", block)
            qubit_match = re.search(r"Qubits\s*:\s*(\d+)", block)

            if qubit_match and depth_match and gates_match:
                qubits = int(qubit_match.group(1))
                depth = int(depth_match.group(1))
                total_gates = int(gates_match.group(1))

                results.append((sequence, algorithm, qubits, depth, total_gates))

    if False:
        # Print parsed table
        print(f"{'Sequence':<10} {'Algorithm':<10} {'Depth':<10} {'Total Gates':<12}")
        print("-" * 50)

        for seq, algo, depth, gates in results:
            print(f"{seq:<10} {algo:<10} {depth:<10} {gates:<12}")

    # Mean calculation
    depths = {"VQE": [], "QAOA": []}
    gates = {"VQE": [], "QAOA": []}
    qubits = {"VQE": [], "QAOA": []}

    for seq, algo, q, depth, total_gates in results:
        qubits[algo].append(q)
        depths[algo].append(depth)
        gates[algo].append(total_gates)

    print("\nMean Values")
    print("-" * 50)

    for algo in ["VQE", "QAOA"]:

        mean_depth = sum(depths[algo]) / len(depths[algo])
        mean_gates = sum(gates[algo]) / len(gates[algo])
        mean_qubit = sum(qubits[algo]) / len(qubits[algo])
        print(f"{algo}")
        print(f"  Max  Qubits        : {max(qubits[algo])}")
        print(f"  Mean Circuit Depth : {np.ceil(mean_depth)}")
        print(f"  Mean Total Gates   : {np.ceil(mean_gates)}")
        print()


def get_optimal_time(text):
    pattern = re.compile(
        r"(VQE|QAOA)\s+optimal:\s*(\d+)\s*/\s*(\d+)\s*\|.*?avg time:\s*([\d.]+)s"
    )

    results = []

    for m in pattern.finditer(text):
        algo = m.group(1)
        optimal = int(m.group(2))
        total = int(m.group(3))
        avg_time = float(m.group(4))

        results.append((algo, optimal, total, avg_time))

    # Print results
    print(f"{'Algo':<6} {'Optimal':<10} {'Total':<8} {'Rate':<8} {'Avg Time(s)':<12}")
    print("-" * 40)

    for algo, opt, total, t in results:
        print(f"{algo:<6} {opt:<10} {total:<8} {opt / total:<8.2f} {t * total:<12.2f}")
    print()


from pathlib import Path

path = Path("../figures/")
log_files = [f.name for f in path.glob("*.log")]
for file_name in log_files:
    print(f"{file_name}")
    with open("../figures/" + file_name, "r") as f:
        text = f.read()
    run(text)
    get_optimal_time(text)
