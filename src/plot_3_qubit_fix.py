import json

import matplotlib.pyplot as plt

import util
from denseCodon import DenseCodon
from hamiltonian import CodonOptimizer
from oneHotCodon import OneHotCodon

CONFIG_PATH = "../config.json"
DATA_PATH = "../data/03-influenza_ha_vaccine.fasta"
OUTPUT_PATH = "qubit_with_fragment.png"

CHUNK_SIZE_RANGE = range(3, 13)


def load_config(config_path: str) -> dict:
    with open(config_path, "r") as f:
        return json.load(f)


def get_sequences_list(data_file_name: str, chunk_size: int) -> list:
    data = util.parse_sequence_from_file(data_file_name)
    return util.split_sequence(data, chunk_size)


def calculate_qubits(sequences: list, config: dict, encode_type: str) -> list[int]:
    codon_class = DenseCodon if encode_type == "dense" else OneHotCodon
    return [
        CodonOptimizer(sequence, config, codon_class, encode_type).qubit_len
        for sequence in sequences
    ]


def collect_qubit_counts(data_file_name: str, config: dict) -> tuple[list, list]:
    dense_max_qubits = []
    one_hot_max_qubits = []

    for chunk_size in CHUNK_SIZE_RANGE:
        sequences = get_sequences_list(data_file_name, chunk_size)
        dense_max_qubits.append(max(calculate_qubits(sequences, config, "dense")))
        one_hot_max_qubits.append(max(calculate_qubits(sequences, config, "one-hot")))
        print()
    return dense_max_qubits, one_hot_max_qubits

def collect_qubits_dataset(config):
    file_name_list = ["01-sars2_spike_vaccine.fasta", "03-influenza_ha_vaccine.fasta",
                 "04-Zika_E_ectodomain.fasta", "05-DENV1_E_ectodomain.fasta",
                 "06-rabies_vaccine_antigen.fasta", "07-ebola_vaccine_antigen.fasta",
                 "08-mers_vaccine_antigen.fasta", "10-hendra_vaccine_antigen.fasta",
                 "11-marburg_vaccine_antigen.fasta", "12-lassa_vaccine_antigen.fasta",
                 "13-Poliovirus_VP1.fasta", "14-Norovirus_VP1_mRNA_vaccine.fasta",
                 "15-HBV_HBsAg_vaccine.fasta"]

    for file_name in file_name_list:
        dense_max_qubits, one_hot_max_qubits = collect_qubit_counts(file_name, config)
def plot_qubits(data_file_name: str) -> None:
    config = load_config(CONFIG_PATH)

    dense_max_qubits, one_hot_max_qubits = collect_qubit_counts(data_file_name, config)

    x = list(CHUNK_SIZE_RANGE)
    Dataset = [
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
    ]

    plt.figure(dpi=500)
    plt.plot(x, dense_max_qubits, marker="o", label="dense")
    plt.plot(x, one_hot_max_qubits, marker="s", label="one-hot")
    plt.xlabel("Fragment size (codons)")
    plt.ylabel("Qubits")
    plt.title("Required number of qubits")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(OUTPUT_PATH)


if __name__ == "__main__":
    plot_qubits(DATA_PATH)