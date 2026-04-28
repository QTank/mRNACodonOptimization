import re

import util

log_file = "../experiments/resulti_6.log"  # update path if needed

results = []
current_seq = None
current_solver = None
seq_set = []
with open(log_file, "r") as f:
    for line in f:
        line = line.strip()

        # Detect sequence header
        seq_match = re.match(r"Sequence:\s+(\S+)", line)
        if seq_match:
            current_seq = seq_match.group(1)

        # Detect solver block
        if line in ("[VQE]", "[SA]", "[BRUTE]"):
            current_solver = line.strip("[]")

        # Extract energy
        energy_match = re.match(r"-\s+Energy\s+:\s+([\d.]+)", line)
        if energy_match and current_seq and current_solver:
            energy = float(energy_match.group(1))

        # Extract time (last field in block — finalize record here)
        time_match = re.match(r"-\s+Time\s+:\s+([\d.]+)s", line)
        if time_match and current_seq and current_solver:
            time_val = float(time_match.group(1))
            results.append({
                "sequence": current_seq,
                "solver": current_solver,
                "energy": energy,
                "time_s": time_val,
            })

# Print results
print(f"{'Sequence':<20} {'Solver':<8} {'Energy':>12} {'Time (s)':>10}")
print("-" * 54)
for r in results:
    print(f"{r['sequence']:<20} {r['solver']:<8} {r['energy']:>12.4f} {r['time_s']:>10.2f}")

print(f"\nTotal records: {len(results)}")


data = util.split_sequence_from_file("../data/data", 6)
for key in data:
    print(key, key in results)