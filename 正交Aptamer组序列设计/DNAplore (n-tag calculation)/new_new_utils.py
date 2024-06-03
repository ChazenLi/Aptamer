from itertools import permutations, product
from random import choice
from typing import List, Tuple, Set

from nupack import Model, mfe, Ensemble

DNA_BASE_PAIR = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
}

DNA_BASE = ["A", "T", "C", "G"]

DNA_FORBID = ["AAA", "TTT", "CCC", "GGG"]


def is_sequence_pair(sequence_a: str, sequence_b: str) -> bool:
    """
    Check if two sequences are pair.
    :param sequence_a:
    :param sequence_b:
    :return: is pair
    """
    if len(sequence_a) != len(sequence_b):
        return False
    return all(DNA_BASE_PAIR[a] == b for a, b in zip(sequence_a, sequence_b))


def get_pair_complement(dna_string: str) -> str:
    """
    Calculate the pair complement of a DNA string.
    :param dna_string:
    :return: pair complement
    """
    return ''.join(DNA_BASE_PAIR[base] for base in dna_string)


def generate_random_dna_sequence(n: int) -> str:
    """
    Generate a random DNA sequence.
    :param n: length of DNA sequence
    :return: DNA sequence
    """
    return ''.join(choice(DNA_BASE) for _ in range(n))


def generate_3_permutation():
    """
    Generate the sequence from AAA to GGG.
    order: A -> T -> C -> G
    :return: list of 3-mers
    """
    res = []
    for i in DNA_BASE:
        for j in DNA_BASE:
            res.extend(i + j + k for k in DNA_BASE)
    return res


def generate_n_permutation(n: int):
    """
    Generate the sequence from AA...A to GG...G.
    the length of sequence is n.
    order: A -> T -> C -> G

    :param n:
    :return:
    """
    return [''.join(i) for i in product(DNA_BASE, repeat=n)]


def get_3_base(seq: str) -> set:
    return {seq[i:i + 3] for i in range(len(seq) - 2)}


def get_n_base(seq: str, n: int) -> Set[str]:
    return {seq[i:i + n] for i in range(len(seq) - n + 1)}


def generate_4_permutation():
    res = []
    for i in DNA_BASE:
        for j in DNA_BASE:
            for k in DNA_BASE:
                res.extend(i + j + k + m for m in DNA_BASE)
    return res


def energy_filter(raw_sequence: str, cal_sequence: List[str], rules: List[Tuple[float, int]]) -> bool:
    for s in cal_sequence:
        energy = energy_cal(raw_sequence, s)
        for rule in rules:
            if len(s) <= rule[0]:
                if energy < rule[1]:
                    return False
                else:
                    break
    return True


def calc_base_ratio(sequence: str) -> float:
    """
    calculate (A+T) divide (A+T+C+G)
    :param sequence:
    :return:
    """
    return (sequence.count("A") + sequence.count("T")) / len(sequence)


model = Model(material="dna", ensemble=Ensemble.stacking, celsius=37, sodium=1.0, magnesium=0.0)


# def energy_satisfy(seq, input_seq):
#     mfes = []
#     for sn in input_seq:
#         eng = energy_cal(seq, sn.seq)
#         if eng <= sn.gth:
#             break
#         mfes.append(eng)
#     return None if len(mfes) != len(input_seq) else mfes


def energy_cal(seq_a, seq_b):
    structures = None
    try:
        structures = mfe(strands=[seq_a, seq_b], model=model)
    except Exception as ex:
        print(ex)
    if not structures:
        return -999.0
    energy = structures[0].energy
    return energy
