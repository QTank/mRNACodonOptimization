import numpy as np
import python_codon_tables as pct
import util

import vqe_solver, qaoa_solver, sa_solver, brute_force

class CodonOptimizer:
    """
    Unified class for Codon Optimization using Qubit Operators.
    Supports both One-Hot and Binary encoding depending on the Codon class passed.
    """

    def __init__(self, protein_sequence, config, codon_class, codon_type):
        self.protein_seq = protein_sequence
        self.n_len = 3 * len(protein_sequence)

        # Metadata & Weights
        meta = config.get('metadata', {})
        weights = config.get('weights', {})

        self.table_name = meta.get('table_name', 'e_coli_316407')
        self.target_gc = meta.get('gc_content', 0.5)

        self.w_usage = weights.get('usage', 0.3)
        self.w_gc = weights.get('target_gc', 0.9)
        self.w_rep = weights.get('repeated_nucleotides', 1.0)
        self.penalty = weights.get('redundant_encoding', 13000)

        self.codon_table = pct.get_codons_table(self.table_name)
        self.qubit_len = self.calculate_qubit_len(codon_type)
        self.codon_list = self._build_codon_list(codon_class)

    def calculate_qubit_len(self, codon_type):
        count = 0
        for amino in self.protein_seq:
            if codon_type == "one-hot":
                count += len(self.codon_table[amino])
            elif codon_type == "dense":
                count += round(np.log2(len(self.codon_table[amino])))
            else:
                raise("Codon_type error!")
        return count

    def _build_codon_list(self, codon_class):
        """Instantiates codon objects for each amino acid in the sequence."""
        codons = []
        current_idx = 0
        for amino in self.protein_seq:
            # We calculate qubit length locally or pass a total placeholder
            c = codon_class(amino, current_idx, self.qubit_len)
            codons.append(c)
            current_idx += c.encoding_qubit_len

        return codons

    def create_usage_term(self, epsilon=0.05):
        """H_usage: Penalizes codons with low frequency in the target organism."""
        h_usage = 0
        for codon in self.codon_list:
            # Skip Amino Acids with only one codon (Met, Trp)
            if len(self.codon_table.get(codon.amino, {})) <= 1:
                continue

            for seq, (op, freq) in codon.indicator_dict.items():
                score = -np.log10(freq + epsilon)
                h_usage += op * round(score, 3)

        return h_usage * self.w_usage * (self.n_len ** 2)

    def create_gc_term(self):
        """H_gc: Penalizes deviation from the target GC content."""
        h_gc_sum = 0
        for codon in self.codon_list:
            for seq, (op, _) in codon.indicator_dict.items():
                h_gc_sum += op * util.get_gc_count(seq)

        identity = util.build_full_identity(self.qubit_len)
        target = identity * self.target_gc * self.n_len
        return ((h_gc_sum - target) ** 2) * self.w_gc

    def create_repetition_term(self):
        """H_rep: Penalizes repeated nucleotides at codon junctions."""
        h_rep = 0
        for i in range(1, len(self.codon_list)):
            prev, curr = self.codon_list[i - 1], self.codon_list[i]

            for p_seq, (p_op, _) in prev.indicator_dict.items():
                for c_seq, (c_op, _) in curr.indicator_dict.items():
                    # Quadratic term (Interaction between adjacent qubits)
                    interaction = p_op @ c_op
                    penalty = util.get_neighbour_repetition(p_seq, c_seq)
                    h_rep += interaction * penalty

        return h_rep * self.w_rep * (self.n_len ** 2)

    def create_redundant_encoding(self):
        redundant_list = []
        for codon in self.codon_list:
            redundant_list.extend(codon.redundant_indicator_list)

        return self.penalty * sum(redundant_list)

    def create_qubit_op(self):
        """Aggregates all terms into the final Hamiltonian."""
        # Penalty for invalid/redundant bitstrings (One-hot constraints)
        h_total = (
                self.create_usage_term() +
                self.create_gc_term() +
                self.create_repetition_term() +
                self.create_redundant_encoding()
        )
        return h_total.simplify()