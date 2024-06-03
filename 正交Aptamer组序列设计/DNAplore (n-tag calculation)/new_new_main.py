import json
import random
from math import inf

from new_new_utils import get_pair_complement, generate_random_dna_sequence, DNA_BASE, DNA_FORBID, \
    generate_3_permutation, \
    get_3_base, calc_base_ratio, generate_n_permutation, get_n_base, energy_filter, energy_cal

n = 5
m = 20
o = 1000
rules = [
    (8.0, -4.5),
    (inf, -8)
]
n_base = 3
forbid_bases = [
    'AAA', 'TTT', 'CCC', 'GGG'
]
base_ratio = True
base_ratio_min = 0.25
base_ratio_max = 0.75
bgyl = True

input_sequence = ["TTACGGATGC", "CGATGGCAAGATATC", "GCCACTA", ]
input_pair_sequence = [get_pair_complement(sequence) for sequence in input_sequence]

available = generate_n_permutation(n_base)

unavailable = set()
for sequence in input_sequence + input_pair_sequence:
    bases = get_n_base(sequence, n_base)
    unavailable.update(bases)
available = {sequence for sequence in available if sequence not in list(unavailable) + forbid_bases}


generate_sequence = []


def find_all_seq(seq: str, seq_list: list):
    if len(seq) >= m:
        seq_list.append(seq)
        return
    start_with = seq[-n_base + 1:]
    filtered = [sequence for sequence in available if sequence.startswith(start_with)]
    random.shuffle(filtered)
    for sequence in filtered:
        new_seq = seq + sequence[n_base - 1]
        if any(map(lambda x: new_seq.endswith(x), forbid_bases)):
            continue
        pair = get_pair_complement(sequence)
        available.discard(pair)
        find_all_seq(seq + sequence[n_base - 1], seq_list)
        available.add(pair)


def find_seq(seq: str, seqs: list):
    if len(seq) >= m:
        ratio = calc_base_ratio(seq)
        if base_ratio and (ratio < base_ratio_min or ratio > base_ratio_max):
            return
        if energy_filter(seq, input_sequence, rules) \
                and energy_filter(seq, seqs, rules) \
                and energy_filter(seq, [seq], rules):
            generate_seq(seqs + [seq])
        return
    start_with = seq[-n_base + 1:]
    filtered = [sequence for sequence in available if sequence.startswith(start_with)]
    random.shuffle(filtered)
    for sequence in filtered:
        new_seq = seq + sequence[n_base - 1]
        if any(map(lambda x: new_seq.endswith(x), forbid_bases)):
            continue
        pair = get_pair_complement(sequence)
        available.discard(pair)
        find_seq(seq + sequence[n_base - 1], seqs)
        available.add(pair)


bgyl_flag = False


def generate_seq(seqs: list):
    global bgyl_flag, bgyl
    if len(generate_sequence) >= o:
        return
    if len(seqs) == 1:
        bgyl_flag = False
    if bgyl_flag and bgyl:
        return
    if len(seqs) >= n:
        # generate_sequence.append(seqs)
        # print(f"Generate: {seqs}")
        res_seq = []
        for seq in seqs:
            res = {
                "seq": seq,
                "self": energy_cal(seq, seq)
            }
            for idx, in_seq in enumerate(input_sequence):
                res[f"in_{idx}"] = energy_cal(seq, in_seq)
            res_seq.append(res)
        generate_sequence.append(res_seq)
        print(json.dumps(res_seq))
        bgyl_flag = True
        return
    available_list = list(available)
    random.shuffle(available_list)
    for sequence in available_list:
        pair = get_pair_complement(sequence)
        available.discard(pair)
        find_seq(sequence, seqs)
        available.add(pair)
    # seq_list = []
    # available_list = list(available)
    # random.shuffle(available_list)
    # for sequence in available_list:
    #     pair = get_pair_complement(sequence)
    #     available.discard(pair)
    #     find_all_seq(sequence, seq_list)
    #     available.add(pair)
    # for seq in seq_list:
    #     if 0.25 < calc_base_ratio(seq) < 0.75:
    #         use_bases = get_n_base(seq, n_base)
    #         unavailable_pair = set(map(get_pair_complement, use_bases))
    #         available.difference_update(unavailable_pair)
    #         generate_seq(seqs + [seq])
    #         available.update(unavailable_pair)


generate_seq([])

for seqs in generate_sequence:
    print(seqs)

# changed = True
# while changed:
#     changed = False
#     for i in range(n):
#         sequence = ''
#         for q in qwq:
#             flag = True
#             for in_squ in input_pair_sequence + generate_sequence[:i]:
#                 if in_squ.find(q) != -1 or in_squ.find(q) != -1:
#                     flag = False
#                     break
#             if flag:
#                 sequence = q
#                 break
#         if sequence == '':
#             print('不存在')
#             exit(0)
#         for j in range(m - 4):
#             while True:
#                 squ = sequence[j + 1:j + 4]
#                 while squ in DNA_FORBID:
#                     squ = squ[:1] + random.choice(DNA_BASE)
#                 squ_pair = get_pair_complement(squ)
#
#                 flag = True
#                 for in_squ in input_pair_sequence + generate_sequence[:i]:
#                     if in_squ.find(squ_pair) != -1 or in_squ.find(squ) != -1:
#                         change_idx = j + 2  # random.randint(0, 2)
#                         old_base = sequence[change_idx]
#                         new_base = random.choice([base for base in DNA_BASE if base != old_base])
#                         sequence = sequence[:change_idx] + new_base + sequence[change_idx + 1:]
#                         # 如果发生变化导致了禁止的序列，则重新生成
#                         l = [sequence[change_idx + i - 2:change_idx + i + 1] for i in range(3)]
#                         while set(l) & set(DNA_FORBID):
#                             new_base = random.choice([base for base in DNA_BASE if base != old_base])
#                             sequence = sequence[:change_idx] + new_base + sequence[change_idx + 1:]
#                             l = [sequence[change_idx + i - 2:change_idx + i + 1] for i in range(3)]
#                         changed = True
#                         flag = False
#                         print(f"[{i}] {change_idx} {old_base} -> {new_base}")
#                 if flag:
#                     break
#         print(f"[{i}]: {sequence}")
