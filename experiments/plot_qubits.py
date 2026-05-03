import matplotlib.pyplot as plt
import time, csv
import numpy as np


def plot_qubits(dense_qubits, one_hot_qubits):
    x = range(3, 21)
    plt.plot(x, one_hot_qubits[0], marker='o', label='maximum on QAOA')
    plt.plot(x, dense_qubits[0], marker='*', label='maximum on our method')
    plt.plot(x, one_hot_qubits[1], marker='h', label='mean on QAOA')
    plt.plot(x, dense_qubits[1], marker="v", label="mean on our method")
    plt.xlabel("Length of each fragment")
    plt.ylabel("Number of qubits")
    # plt.title("the required number of gates with one-hot and dense encoding")
    plt.xlim(2, 21)
    plt.xticks(np.arange(3, 21, 2))
    plt.legend()
    plt.show()


def plot_qubits_split(dense_qubits, one_hot_qubits, output_file):
    x = range(3, 21)
    fig, axes = plt.subplots(1, 2)

    axes[1].plot(x, one_hot_qubits[0], marker='o', label='maximum on QAOA')
    axes[1].plot(x, dense_qubits[0], marker='*', label='maximum on our method')
    axes[0].plot(x, one_hot_qubits[1], marker='h', label='mean on QAOA')
    axes[0].plot(x, dense_qubits[1], marker="v", label="mean on our method")
    # axes[1].xlabel("Length of each fragment")
    # axes[1].ylabel("Number of qubits")
    # plt.title("the required number of gates with one-hot and dense encoding")
    # axes[1].xlim(2, 21)
    # axes[1].xticks(np.arange(3, 21, 2))

    # axes[1].legend()
    # plt.savefig(output_file)
    plt.figure(dpi=400)
    plt.show()
