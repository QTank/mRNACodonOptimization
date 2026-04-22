# Quantum mRNA Codon Optimization via VQE

## Overview
This project implements a hybrid quantum–classical approach for mRNA codon optimization using the Variational Quantum Eigensolver (VQE). It focuses on mapping biological constraints into a qubit Hamiltonian and solving the combinatorial optimization problem efficiently.

**Key biological constraints:**
- Codon usage bias  
- GC content  
- Repetition penalty  
- Valid codon mapping

**Quantum highlights:**
- Dense encoding for qubit efficiency  
- Hamiltonian formulation of constraints  
- Hybrid VQE optimization (Qiskit)

---

## Pipeline

1. Encode amino acids into qubits (dense encoding)  
2. Construct Hamiltonian with constraints  
3. Run VQE with classical optimizer  
4. Decode optimal qubit state into valid mRNA sequence  

---

## Installation

```bash
git clone https://github.com/QTank/mRNACodonOptimization.git
cd mRNACodonOptimization
pip install -r requirements.txt