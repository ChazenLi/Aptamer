from rdkit import Chem
from rdkit.Chem import AllChem
import os


def convert_and_save(smiles, output_pdb_file, output_mol_file):
    # 将SMILES字符串转换为分子对象
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("无法从SMILES字符串生成分子对象")
        return

    # 生成标准SMILES字符串
    standard_smiles = Chem.MolToSmiles(mol)
    print(f"标准SMILES序列: {standard_smiles}")

    # 使用RDKit将标准SMILES字符串转换为分子对象
    mol = Chem.MolFromSmiles(standard_smiles)
    if mol is None:
        print("无法从标准SMILES字符串生成分子对象")
        return

    # 添加氢原子
    mol = Chem.AddHs(mol)

    # 生成3D坐标
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if result != 0:
        print("分子嵌入失败")
        return

    # 优化分子结构
    try:
        AllChem.UFFOptimizeMolecule(mol)

        # 保存为PDB文件
        with open(output_pdb_file, 'wb') as f:
            f.write(Chem.MolToPDBBlock(mol).encode('utf-8'))
        print(f"已生成PDB文件: {output_pdb_file}")

        # 保存为MOL文件
        with open(output_mol_file, 'wb') as f:
            f.write(Chem.MolToMolBlock(mol).encode('utf-8'))
        print(f"已生成MOL文件: {output_mol_file}")

    except Exception as e:
        print(f"分子优化失败: {e}")


# 测试函数
smiles = "OP(=O)(O)OCC1OC(n2cnc3c2ncnc3N)CC1OP(=O)(O)OCC3OC(n4cnc5c4nc(N)[nH]c(=O)5)CC3OP(=O)(O)OCC5OC(n6cnc7c6ncnc7N)CC5OP(=O)(O)OCC7OC(n8cnc9c8nc(N)[nH]c(=O)9)CC7OP(=O)(O)OCC9OC(n%10cnc%11c%10ncnc%11N)CC9OP(=O)(O)OCC%11OC(n%12ccc(nc%12=O)N)CC%11OP(=O)(O)OCC%13OC(n%14ccc(nc%14=O)N)CC%13OP(=O)(O)OCC%15OC(n%16ccc(nc%16=O)N)CC%15OP(=O)(O)OCC%17OC(N%18C=C(C)C(=O)NC%18=O)CC%17OP(=O)(O)OCC%19OC(n%20cnc%21c%20nc(N)[nH]c(=O)%21)CC%19OP(=O)(O)OCC%21OC(n%22cnc%23c%22ncnc%23N)CC%21OP(=O)(O)OCC%23OC(n%24ccc(nc%24=O)N)CC%23OP(=O)(O)OCC%25OC(N%26C=C(C)C(=O)NC%26=O)CC%25OP(=O)(O)OCC%27OC(n%28cnc%29c%28nc(N)[nH]c(=O)%29)CC%27OP(=O)(O)OCC%29OC(n%30ccc(nc%30=O)N)CC%29OP(=O)(O)OCC%31OC(n%32cnc%33c%32nc(N)[nH]c(=O)%33)CC%31OP(=O)(O)OCC%33OC(n%34cnc%35c%34ncnc%35N)CC%33OP(=O)(O)OCC%35OC(n%36cnc%37c%36ncnc%37N)CC%35OP(=O)(O)OCC%37OC(n%38ccc(nc%38=O)N)CC%37OP(=O)(O)OCC%39OC(n%40ccc(nc%40=O)N)CC%39"
output_dir = os.getcwd()  # 获取当前工作目录
output_pdb_file = os.path.join(output_dir, "toutput.pdb")
output_mol_file = os.path.join(output_dir, "toutput.mol")

convert_and_save(smiles, output_pdb_file, output_mol_file)
