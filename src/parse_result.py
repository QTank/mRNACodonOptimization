import re
from collections import defaultdict

methods = ["VQE", "SA", "BRUTE"]

stats = {
    m: {"opt": 0, "gap": 0.0, "time": 0.0}
    for m in methods
}

count = 0

with open("../experiments/resulti_5.log", 'r') as f:
    text = f.read()

blocks = text.split("############################################################")

vqe_rna = ""
sa_rna = ""
brute_rna = ""
for block in blocks:
    if "Sequence:" not in block:
        continue

    count += 1

    energies = {}
    times = {}

    for m in methods:
        e = re.search(fr"\[{m}\].*?Energy\s*:\s*([0-9.]+)", block, re.S)
        t = re.search(fr"\[{m}\].*?Time\s*:\s*([0-9.]+)s", block, re.S)
        s = re.search(rf"\[{m}\].*?mRNA\s*:\s*([AUCG]+)", block, re.S)

        if e:
            energies[m] = float(e.group(1))
        if t:
            times[m] = float(t.group(1))


        if s:
            if m == 'VQE':
                vqe_rna += s.group(1)
            if m == 'SA':
                sa_rna += s.group(1)
            if m == 'BRUTE':
                brute_rna += s.group(1)

    if "BRUTE" not in energies:
        continue

    ref = energies["BRUTE"]

    for m in methods:
        if m in energies:
            if abs(energies[m] - ref) < 1e-9:
                stats[m]["opt"] += 1
            stats[m]["gap"] += (energies[m] - ref)
            stats[m]["time"] += times.get(m, 0.0)

print("=" * 60)
print("AGGREGATE RESULTS")
print("=" * 60)

if count == 0:
    print("❌ No valid sequences parsed. Check log structure.")
else:
    for m in methods:
        print(
            f"{m:6} optimal: {stats[m]['opt']}/{count} | "
            f"avg energy gap: {stats[m]['gap'] / count:+.4f} | "
            f"avg time: {stats[m]['time'] / count:.2f}s"
        )

import util
print(vqe_rna)
print(util.convert_rna_to_dna(vqe_rna))
print(sa_rna)
print(util.convert_rna_to_dna(sa_rna))
print(brute_rna)
print(util.convert_rna_to_dna(brute_rna))
