import matplotlib.pyplot as plt

plt.style.use('seaborn-v0_8-whitegrid')

plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'legend.fontsize': 11
})

data = {
    "Dataset": [
        "SARS-CoV-2",
        "Influenza",
        "Zika",
        "DENV1",
        "Rabies",
        "Ebola",
        "Mers",
        "Nipah",
        "Hendra",
        "Marburg",
        "Lassa",
        "Poliovirus",
        "Norovirus",
        "HBV"
    ],

    "VQE_Qubits": [9, 9, 9, 8, 9, 8, 9, 9, 8, 8, 9, 9, 9, 9],
    "QAOA_Qubits": [18, 18, 18, 16, 18, 16, 18, 18, 16, 16, 18, 18, 18, 18],

    "VQE_Depth": [15, 15., 15., 15., 15, 15,
                  16, 15, 14, 15, 15, 15,
                  15, 16],
    "QAOA_Depth": [145, 133, 138, 137, 148, 152,
                   162, 147, 118, 147, 143, 146,
                   151, 160],

    "VQE_Gates": [33, 32, 33, 33, 33, 34,
                  35, 34, 30, 34, 32, 34,
                  32, 35],
    "QAOA_Gates": [297, 286, 295, 293, 316, 324,
                   349, 312, 237, 313, 291, 312,
                   295, 357],

    "VQE_Time": [88.55, 38.87, 36.48, 41.08, 38.88, 40.08,
                 31.20, 36.52, 23.69, 44.16, 35.04, 23.76,
                 29.43, 35.10],
    "QAOA_Time": [423.50, 196.04, 170.24, 180.12, 208.98, 205.41,
                  166.40, 172.64, 87.55, 200.56, 173.74, 113.85,
                  126.44, 200.20],
}

labels = data["Dataset"]

plt.figure(dpi=500)
plt.plot(labels, data["VQE_Qubits"], marker='o', label="VQE Qubits")
plt.plot(labels, data["QAOA_Qubits"], marker='s', label="QAOA Qubits")

plt.xticks(rotation=45)
plt.ylabel("Total Gates")
plt.title("VQE vs QAOA Qubits")
plt.legend()
plt.tight_layout()
plt.savefig("qubits.png", dpi=600, bbox_inches='tight')

# ---- Circuit Depth Plot ----
plt.figure(dpi=500)

plt.plot(labels, data["VQE_Depth"], marker='o', label="VQE Depth")
plt.plot(labels, data["QAOA_Depth"], marker='s', label="QAOA Depth")

plt.xticks(rotation=45)
plt.ylabel("Circuit Depth")
plt.title("VQE vs QAOA Circuit Depth")
plt.legend()
plt.tight_layout()
plt.savefig("circuit_depth.png", dpi=600, bbox_inches='tight')

# ---- Gate Count Plot ----
plt.figure(dpi=500)
plt.plot(labels, data["VQE_Gates"], marker='o', label="VQE Gates")
plt.plot(labels, data["QAOA_Gates"], marker='s', label="QAOA Gates")

plt.xticks(rotation=45)
plt.ylabel("Total Gates")
plt.title("VQE vs QAOA Gate Count")
plt.legend()
plt.tight_layout()
plt.savefig("gate_count.png", dpi=600, bbox_inches='tight')

plt.figure(dpi=500)
plt.plot(labels, data["VQE_Time"], marker='o', label="VQE Time")
plt.plot(labels, data["QAOA_Time"], marker='s', label="QAOA Time")

plt.xticks(rotation=45)
plt.ylabel("Time (s)")
plt.title("VQE vs QAOA Runtime")
plt.legend()
plt.tight_layout()
plt.savefig("runtime.png", dpi=600, bbox_inches='tight')
