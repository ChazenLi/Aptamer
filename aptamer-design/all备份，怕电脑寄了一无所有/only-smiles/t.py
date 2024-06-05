from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import os
import re

# 定义每个碱基的特征字符模板
base_features = {
    'A': 'OP(=O)(O)O[C@@H]1O[C@@H](n2cnc3c2ncnc3N)[C@H]1',
    'C': 'OP(=O)(O)O[C@@H]1O[C@@H](n2ccc(nc2=O)N)[C@H]1',
    'T': 'OP(=O)(O)O[C@@H]1O[C@@H](N2C=C(C)C(=O)NC2=O)[C@H]1',
    'G': 'OP(=O)(O)O[C@@H]1O[C@@H](n2cnc3c2nc(N)[nH]c(=O)3)[C@H]1'
}

def adjust_numbers(base_str, offset):
    """Adjust the numbers in the base feature string by the given offset."""
    def replace_match(match):
        number = int(match.group())
        if number + offset >= 10:
            return '%' + str(number + offset)
        else:
            return str(number + offset)
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
            # 最后一个手性碳需要[C@]而不是[C@@]
            if index == len(sequence_array) - 1:
                adjusted_feature = adjusted_feature.replace('[C@@H]', '[C@H]', 1)
            result_sequence.append(adjusted_feature)
    # 将所有处理后的特征字符串连接起来
    final_sequence = ''.join(result_sequence)
    return final_sequence

def generate_pdb(sequence, output_file):
    # 处理序列以获得最终的SMILES字符串
    smiles_sequence = process_sequence(sequence)
    print(f"生成的SMILES序列: {smiles_sequence}")
    # 使用RDKit将SMILES字符串转换为分子对象
    mol = Chem.MolFromSmiles(smiles_sequence)
    if mol is None:
        print("无法从SMILES字符串生成分子对象")
        return
    # 如果分子对象成功创建，则进行3D坐标生成并保存为PDB文件
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if result != 0:
        print("分子嵌入失败")
        return
    try:
        AllChem.UFFOptimizeMolecule(mol)
        with open(output_file, 'wb') as f:
            f.write(Chem.MolToPDBBlock(mol).encode('utf-8'))
        print(f"已生成PDB文件: {output_file}")
    except Exception as e:
        print(f"分子优化失败: {e}")

def save_mol_file(sequence, output_file):
    # 处理序列以获得最终的SMILES字符串
    smiles_sequence = process_sequence(sequence)
    print(f"生成的SMILES序列: {smiles_sequence}")
    # 使用RDKit将SMILES字符串转换为分子对象
    mol = Chem.MolFromSmiles(smiles_sequence)
    if mol is None:
        print("无法从SMILES字符串生成分子对象")
        return
    try:
        with open(output_file, 'wb') as f:
            f.write(Chem.MolToMolBlock(mol).encode('utf-8'))
        print(f"已生成MOL文件: {output_file}")
    except Exception as e:
        print(f"保存MOL文件失败: {e}")

def draw_molecule(sequence, image_file, svg_file):
    # 处理序列以获得最终的SMILES字符串
    smiles_sequence = process_sequence(sequence)
    print(f"生成的SMILES序列: {smiles_sequence}")
    # 使用RDKit将SMILES字符串转换为分子对象
    mol = Chem.MolFromSmiles(smiles_sequence)
    if mol is None:
        print("无法从SMILES字符串生成分子对象")
        return
    # 如果分子对象成功创建，则生成2D图像和SVG文件
    AllChem.Compute2DCoords(mol)
    # 生成PNG图像
    img = Draw.MolToImage(mol, size=(300, 300))
    img.save(image_file)
    print(f"已生成PNG图像: {image_file}")
    # 生成SVG图像
    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    with open(svg_file, 'w', encoding='utf-8') as f:
        f.write(svg)
    print(f"已生成SVG图像: {svg_file}")

# 测试函数
sequence = "AGAGACCC"
output_dir = os.getcwd()  # 获取当前工作目录
output_pdb_file = os.path.join(output_dir, "output.pdb")
mol_file = os.path.join(output_dir, "molecule.mol")
image_file = os.path.join(output_dir, "molecule.png")
svg_file = os.path.join(output_dir, "molecule.svg")
generate_pdb(sequence, output_pdb_file)
save_mol_file(sequence, mol_file)
draw_molecule(sequence, image_file, svg_file)
