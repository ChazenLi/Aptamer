import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor

# 定义每个碱基的特征字符模板
base_features = {
    'A': 'OP(=O)(O)OCC1OC(n2cnc3c2ncnc3N)CC1',
    'C': 'OP(=O)(O)OCC1OC(n2ccc(nc2=O)N)CC1',
    'T': 'OP(=O)(O)OCC1OC(N2C=C(C)C(=O)NC2=O)CC1',
    'G': 'OP(=O)(O)OCC1OC(n2cnc3c2nc(N)[nH]c(=O)3)CC1O'
}

# 环数映射
ring_counts = {
    'A': 3,
    'G': 3,
    'C': 2,
    'T': 2
}


def adjust_numbers(base_str, offset, ring_count):
    """Adjust the numbers in the base feature string by the given offset and ring count."""

    def replace_match(match):
        return str(int(match.group()) + offset * ring_count)

    adjusted_str = re.sub(r'\d+', replace_match, base_str)
    return adjusted_str


def process_sequence(sequence):
    # 特征描述替换
    result_sequence = []

    for index, nucleotide in enumerate(sequence):
        if nucleotide in base_features:
            feature_str = base_features[nucleotide]
            ring_count = ring_counts[nucleotide]
            adjusted_feature = adjust_numbers(feature_str, index, ring_count)
            result_sequence.append(adjusted_feature)

    return ''.join(result_sequence)


def generate_molecule(smiles_str):
    """生成分子，进行加氢、检查等操作，输出图片和PDB文件格式"""
    mol = Chem.MolFromSmiles(smiles_str)
    if mol is None:
        raise ValueError(f"无法从SMILES字符串生成分子: {smiles_str}")

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    # 输出图片
    img = Draw.MolToImage(mol, size=(300, 300))
    img.save("/mnt/data/molecule.png")

    # 输出PDB文件
    Chem.MolToPDBFile(mol, "/mnt/data/molecule.pdb")

    return mol


# 测试函数
sequence = "ATCG"
processed_sequence = process_sequence(sequence)

print("处理后的特征序列:")
print(processed_sequence)

# 生成并检查分子
mol = generate_molecule(processed_sequence)

# 输出分子图片和PDB文件
rdDepictor.Compute2DCoords(mol)
img = Draw.MolToImage(mol, size=(300, 300))
img.show()
Chem.MolToPDBFile(mol, "/mnt/data/molecule.pdb")
