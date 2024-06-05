from rdkit import Chem
from rdkit.Chem import AllChem
import re


def create_reactant_string(sequence, base_features):
    reactant_str = ""
    for index, nucleotide in enumerate(sequence):
        if nucleotide in base_features:
            feature_str = base_features[nucleotide]
            adjusted_feature = adjust_numbers(feature_str, index * 2)
            reactant_str += adjusted_feature
    return reactant_str


def perform_reaction(reactant_str, g_str):
    # 创建RDKit分子对象
    reactant_mol = Chem.MolFromSmiles(reactant_str)
    g_mol = Chem.MolFromSmiles(g_str)

    if not reactant_mol or not g_mol:
        raise ValueError("无效的SMILES字符串用于反应物或G分子。")

    # 创建反应对象 (示例反应)
    rxn = AllChem.ReactionFromSmarts('[Reactant:1].[G:2]>>[Product:1][Product:2]')

    # 执行反应
    products = rxn.RunReactants((reactant_mol, g_mol))

    # 获取产物的SMILES
    product_smiles = [Chem.MolToSmiles(product[0]) for product in products]

    return product_smiles


def adjust_numbers(base_str, offset):
    """调整基序特征字符串中的数字根据给定的偏移量。"""

    def replace_match(match):
        return str(int(match.group()) + offset)

    adjusted_str = re.sub(r'\d+', replace_match, base_str)
    return adjusted_str


# 定义碱基特征
base_features = {
    'A': 'O=P(O)(O)OCC1OC(C(C1O)n2cnc3c2ncnc3N)',
    'C': 'O=P(O)(O)OCC1OC(C(C1O)n2ccc(nc2=O)N)',
    'T': 'O=P(O)(O)OCC1OC(C(C1O)N2C=CC(=O)NC2=O)'
}

# 定义G特征字符串
g_str = 'Nc1nc2c(ncn2C2CC(O)C(COP(=O)(O)O)O2)c(=O)[nH]1'

# 测试函数
sequence = "ATC"
reactant_str = create_reactant_string(sequence, base_features)
print("生成的反应物字符串:", reactant_str)

# 由于反应SMARTS和G分子的复杂性，执行反应部分会比较复杂
# 这里示例的反应SMARTS '[Reactant:1].[G:2]>>[Product:1][Product:2]' 可能需要根据具体化学反应调整

# 打印生成的反应物字符串和反应结果
# product_smiles = perform_reaction(reactant_str, g_str)
# print("反应结果的SMILES表示:")
# for smiles in product_smiles:
#     print(smiles)
