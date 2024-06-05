import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import re

# 核苷酸的标准 SMILES
nucleotides_smiles = {
    'A': 'OP(=O)(O)OC[C@@]1O[C@@](n2cnc3c2ncnc3N)C[C@]1',
    'T': 'OP(=O)(O)OC[C@@]1O[C@@](N2C=C(C)C(=O)NC2=O)C[C@]1',
    'C': 'OP(=O)(O)OC[C@@]1O[C@@](n2ccc(nc2=O)N)C[C@]1',
    'G': 'OP(=O)(O)OC[C@@]1O[C@@](n2cnc3c2nc(N)[nH]c(=O)3)C[C@]1'
}


def update_progress(step):
    print(f"正在进行：{step}")


def smiles_adjust_rings(smiles):
    # 调整环标签
    def replacer(match):
        number = int(match.group(1))
        if number > 99:
            return f"%{number}"
        return match.group(0)

    return re.sub(r"(\d{2,})", replacer, smiles)


def generate_smiles_sequence(sequence):
    update_progress("生成序列的 SMILES 格式")
    smiles_sequence = "".join([nucleotides_smiles[nuc] for nuc in sequence])
    smiles_sequence = smiles_adjust_rings(smiles_sequence)
    return smiles_sequence


def modify_chirality_and_bonding(smiles_sequence):
    # 根据手性构象变化和成键变化修改 SMILES
    # 这里只是一个简单的例子，具体修改需要根据实际情况和需要进行调整
    smiles_sequence = re.sub(r'C@@', 'C@', smiles_sequence)
    smiles_sequence = re.sub(r'N([H])=NC=1C', 'N([H])C=NC1', smiles_sequence)
    return smiles_sequence


def check_valence(mol):
    for atom in mol.GetAtoms():
        valence = atom.GetExplicitValence()
        default_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetSymbol())
        if valence > default_valence:
            print(
                f"原子 {atom.GetIdx()} ({atom.GetSymbol()}) 超过了允许的化合价。当前化合价: {valence}, 默认化合价: {default_valence}")
            return atom.GetIdx(), valence, default_valence
    return None


def smiles_to_3D(smiles, filename_base):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("SMILES 格式不正确，无法生成分子。")

    valence_error = check_valence(mol)
    if valence_error is not None:
        atom_idx, valence, default_valence = valence_error
        raise ValueError(f"原子 {atom_idx} 超过了允许的化合价。当前化合价: {valence}, 默认化合价: {default_valence}")

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    # 保存为 PDB 格式
    pdb_filename = filename_base + ".pdb"
    with open(pdb_filename, 'w') as f:
        f.write(Chem.MolToPDBBlock(mol))
    update_progress(f"3D 结构保存为 PDB 格式: {pdb_filename}")

    # 保存为 MOL 格式
    mol_filename = filename_base + ".mol"
    with open(mol_filename, 'w') as f:
        f.write(Chem.MolToMolBlock(mol))
    update_progress(f"3D 结构保存为 MOL 格式: {mol_filename}")


def save_molecule_image(smiles, filename_base):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol, size=(300, 300))
    img_filename = filename_base + ".png"
    img.save(img_filename)
    update_progress(f"分子结构图片保存为 PNG 格式: {img_filename}")

    svg = Draw.MolToImage(mol, size=(300, 300), kekulize=True, wedgeBonds=True)
    svg_filename = filename_base + ".svg"
    with open(svg_filename, 'w') as f:
        f.write(svg)
    update_progress(f"分子结构图片保存为 SVG 格式: {svg_filename}")


def main():
    sequence = "ATCGATCGATCG"
    update_progress("开始处理")
    smiles_sequence = generate_smiles_sequence(sequence)
    print(f"生成的 SMILES: {smiles_sequence}")

    smiles_sequence = modify_chirality_and_bonding(smiles_sequence)
    print(f"修改后的 SMILES: {smiles_sequence}")

    filename_base = "output_molecule"
    try:
        smiles_to_3D(smiles_sequence, filename_base)
        save_molecule_image(smiles_sequence, filename_base)
    except ValueError as e:
        print(f"错误: {e}")
    update_progress("处理完成")


if __name__ == "__main__":
    main()
