import random
from itertools import product


class DNAplore:
    def __init__(self,
                 n: int,
                 m: int,
                 n_base: int,
                 input_sequence: list[str] = None,
                 forbid_bases: list[str] = None,
                 base_ratio=True,
                 base_ratio_min=0.25,
                 base_ratio_max=0.75):
        self.n = n
        self.m = m
        self.n_base = n_base
        self.input_sequence = [] if input_sequence is None else input_sequence
        self.forbid_bases = [] if forbid_bases is None else forbid_bases
        self.available = set()
        self.base_ratio = base_ratio
        self.base_ratio_min = base_ratio_min
        self.base_ratio_max = base_ratio_max

    def generate(self) -> list[str]:
        all_bases = DNAplore.generate_n_permutation(self.n_base)
        unavailable = set()
        for sequence in self.input_sequence:
            bases = DNAplore.get_n_base(sequence, self.n_base)
            unavailable.update(bases)
        self.available = {
            sequence
            for sequence in all_bases
            if sequence not in list(unavailable) + self.forbid_bases
        }

        return self.__generate_sequence([])

    def __generate_sequence(self, seqs: list[str]):
        if len(seqs) >= self.n:
            return seqs
        available_list = list(self.available)
        random.shuffle(available_list)
        for sequence in available_list:
            pair = DNAplore.get_pair_complement(sequence)
            self.available.discard(pair)
            res = self.__find_seq(sequence, seqs)
            if res is not None:
                return res
            self.available.add(pair)

    def __find_seq(self, seq: str, seqs: list[str]):
        if len(seq) >= self.m:
            if self.base_ratio:
                if self.base_ratio_min < DNAplore.calc_base_ratio(seq) < self.base_ratio_max:
                    seqs.append(seq)
                else:
                    return
            else:
                seqs.append(seq)
            return self.__generate_sequence(seqs)
        start_with = seq[-self.n_base + 1:]
        filtered = [sequence for sequence in self.available if sequence.startswith(start_with)]
        random.shuffle(filtered)
        for sequence in filtered:
            new_seq = seq + sequence[self.n_base - 1]
            if any(map(lambda x: new_seq.endswith(x), self.forbid_bases)):
                continue
            pair = DNAplore.get_pair_complement(sequence)
            self.available.discard(pair)
            res = self.__find_seq(seq + sequence[self.n_base - 1], seqs)
            if res is not None:
                return res
            self.available.add(pair)

    # ------------------- utils -------------------
    DNA_BASE = ['A', 'T', 'C', 'G']
    DNA_BASE_PAIR = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
    }

    @staticmethod
    def generate_n_permutation(n: int):
        """
        Generate the sequence from AA...A to GG...G.
        the length of sequence is n.
        order: A -> T -> C -> G
        :param n:
        :return:
        """
        return [''.join(i) for i in product(DNAplore.DNA_BASE, repeat=n)]

    @staticmethod
    def is_sequence_pair(sequence_a: str, sequence_b: str) -> bool:
        """
        Check if two sequences are pair;

        :param sequence_a:
        :param sequence_b:
        :return: is pair
        """
        if len(sequence_a) != len(sequence_b):
            return False
        return all(DNAplore.DNA_BASE_PAIR[a] == b for a, b in zip(sequence_a, sequence_b))

    @staticmethod
    def get_pair_complement(dna_string: str) -> str:
        """
        Calculate the pair complement of a DNA string;
        :param dna_string:
        :return: pair complement
        """
        return ''.join(DNAplore.DNA_BASE_PAIR[base] for base in dna_string)

    @staticmethod
    def get_n_base(seq: str, n: int) -> set[str]:
        """
        Get all n-mers from a sequence;
        :param seq: sequence
        :param n: length of n-mer
        :return: set of n-mers
        """
        return {seq[i:i + n] for i in range(len(seq) - n + 1)}

    @staticmethod
    def calc_base_ratio(sequence: str) -> float:
        """
        calculate (A+T) divide (A+T+C+G)
        :param sequence:
        :return:
        """
        return (sequence.count("A") + sequence.count("T")) / len(sequence)
