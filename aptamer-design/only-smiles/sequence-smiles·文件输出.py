import re

# 定义每个碱基的特征字符模板
base_features = {
    'A': 'OP(=O)(O)OCC1OC(n2cnc3c2ncnc3N)CC1',
    'C': 'OP(=O)(O)OCC1OC(n2ccc(nc2=O)N)CC1',
    'T': 'OP(=O)(O)OCC1OC(N2C=C(C)C(=O)NC2=O)CC1',
    'G': 'OP(=O)(O)OCC1OC(n2cnc3c2nc(N)[nH]c(=O)3)CC1O'
}

def adjust_numbers(base_str, offset):
    """Adjust the numbers in the base feature string by the given offset."""

    def replace_match(match):
        return str(int(match.group()) + offset)

    adjusted_str = re.sub(r'\d+', replace_match, base_str)
    return adjusted_str

def process_sequence(sequence):
    # 将输入的序列存入数组
    sequence_array = list(sequence)

    # 特征描述替换
    result_sequence = []

    for index, nucleotide in enumerate(sequence_array):
        if nucleotide in base_features:
            feature_str = base_features[nucleotide]
            adjusted_feature = adjust_numbers(feature_str, index * 2)
            result_sequence.append(adjusted_feature)

    # 将所有处理后的特征字符串连接起来
    final_sequence = ''.join(result_sequence)

    return final_sequence

# 测试函数
sequence = "ATCG"
processed_sequence = process_sequence(sequence)

print("处理后的特征序列:")
print(processed_sequence)
