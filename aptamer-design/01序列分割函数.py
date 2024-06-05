def process_dna_sequence(sequence):
    # 将输入的序列存入数组
    sequence_array = list(sequence)

    # 按G分割序列
    split_arrays = []
    current_array = []

    for nucleotide in sequence_array:
        current_array.append(nucleotide)
        if nucleotide == 'G':
            split_arrays.append(current_array)
            current_array = []

    if current_array:
        split_arrays.append(current_array)

    # 计算每个分割数组中的A、T、C、G数量和对应的值
    results = []

    for array in split_arrays:
        a_count = array.count('A')
        t_count = array.count('T')
        c_count = array.count('C')
        g_count = array.count('G')

        value = a_count * 3 + (t_count + c_count) * 2

        results.append({
            'array': array,
            'counts': {'A': a_count, 'T': t_count, 'C': c_count, 'G': g_count},
            'value': value
        })

    num_arrays = len(results)

    return num_arrays, results


# 测试函数
sequence = "ATCGGATCGTGGAACG"
num_arrays, results = process_dna_sequence(sequence)

print(f"分成了 {num_arrays} 个数组")
for i, result in enumerate(results):
    print(
        f"数组 {i + 1}: {result['array']}，A: {result['counts']['A']}，T: {result['counts']['T']}，C: {result['counts']['C']}，G: {result['counts']['G']}，值: {result['value']}")
