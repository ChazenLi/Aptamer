from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
import os

# 定义每个碱基的SMILES表示
base_smiles = {
    'A': 'OP(=O)(O)OCC1OC(n2cnc3c2ncnc3N)CC1',
    'C': 'OP(=O)(O)OCC1OC(n2ccc(nc2=O)N)CC1',
    'T': 'OP(=O)(O)OCC1OC(N2C=C(C)C(=O)NC2=O)CC1',
    'G': 'OP(=O)(O)OCC1OC(n2cnc3c2nc(N)[nH]c(=O)3)CC1'
}

# 使用化学反应模板连接单个碱基生成多核苷酸链
def create_dna_sequence(sequence):
    # 重新定义反应模板，确保原子正确映射
    reaction_smarts = '[C:1][O:2][P:3](=[O:4])([O-:5])[O:6].[C:7][O:8][P:9](=[O:10])([O-:11])[O:12]>>[C:1][O:2][P:3](=[O:4])([O-:5])[O:6][C:7][O:8][P:9](=[O:10])([O-:11])[O:12]'
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)

    # 初始碱基
    product = Chem.MolFromSmiles(base_smiles[sequence[0]])

    for base in sequence[1:]:
        reactant = Chem.MolFromSmiles(base_smiles[base])
        if not reactant:
            print(f"无法生成碱基: {base}")
            continue

        products = rxn.RunReactants((product, reactant))
        if len(products) > 0:
            product = products[0][0]
        else:
            print(f"无法连接碱基: {base}")
            return None

    return Chem.MolToSmiles(product)

def generate_3d_structure(smiles, output_pdb_file, output_mol_file, image_file, svg_file):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("无法从SMILES字符串生成分子对象")
        return

    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if result != 0:
        print("分子嵌入失败")
        return

    try:
        AllChem.UFFOptimizeMolecule(mol)

        # 保存PDB文件
        with open(output_pdb_file, 'wb') as f:
            f.write(Chem.MolToPDBBlock(mol).encode('utf-8'))
        print(f"已生成PDB文件: {output_pdb_file}")

        # 保存MOL文件
        with open(output_mol_file, 'wb') as f:
            f.write(Chem.MolToMolBlock(mol).encode('utf-8'))
        print(f"已生成MOL文件: {output_mol_file}")

        # 生成PNG图像
        AllChem.Compute2DCoords(mol)
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

    except Exception as e:
        print(f"分子优化失败: {e}")

# 测试函数
sequence = "AGAGACCCTGACTGCGAACCCATATCGCATTTCCATCCCA"
output_dir = os.getcwd()  # 获取当前工作目录
output_pdb_file = os.path.join(output_dir, "output.pdb")
output_mol_file = os.path.join(output_dir, "molecule.mol")
image_file = os.path.join(output_dir, "molecule.png")
svg_file = os.path.join(output_dir, "molecule.svg")

# 创建DNA序列
smiles_sequence = create_dna_sequence(sequence)
if smiles_sequence:
    print(f"生成的SMILES序列: {smiles_sequence}")

    # 生成3D结构并输出文件
    generate_3d_structure(smiles_sequence, output_pdb_file, output_mol_file, image_file, svg_file)
