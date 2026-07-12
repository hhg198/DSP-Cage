import os
import numpy as np

# ================= 配置区域 =================
INPUT_PATH = "datas/input/"
FILE_PREFIX = "B"
FILE_RANGE = range(1, 14)

# 阈值设定 (请根据质心算出的新 Diameter 结果再次微调这里)
SMALL_LIMIT = 9.5
LARGE_LIMIT = 10.0

SMALL_THRESH = 0.5
LARGE_THRESH = 4.0

# 原子质量字典 (常见有机元素)
ATOMIC_MASSES = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "P": 30.974,
    "S": 32.06,
    "CL": 35.45,
    "BR": 79.904,
    "I": 126.90,
    "HE": 4.0026  # 有时用氦做dummy atom
}


# ================= 逻辑函数 =================

def get_element_mass(line):
    """
    从 PDB 行解析元素并返回质量
    优先读取 76-78 列的元素符号，如果为空则尝试从原子名称(12-16列)解析
    """
    # 1. 尝试读取标准 Element 字段 (76-78列)
    element = line[76:78].strip().upper()

    # 2. 如果为空，从原子名称解析 (例如 " C1 " -> "C")
    if not element:
        atom_name = line[12:16].strip()
        # 提取字母部分
        element = ''.join([char for char in atom_name if char.isalpha()]).upper()
        # 截取前1-2位 (例如 Ca, Fe, C, N)
        if len(element) > 2:
            element = element[:2]

    return ATOMIC_MASSES.get(element, 12.01)  # 如果找不到，默认按碳原子(12.01)处理，防止报错


def read_pdb_data(filepath):
    """读取 PDB 坐标和对应的质量"""
    positions = []
    masses = []

    if not os.path.exists(filepath):
        return None, None

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    # 读取坐标
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    # 获取质量
                    mass = get_element_mass(line)

                    positions.append([x, y, z])
                    masses.append(mass)
                except ValueError:
                    continue
    return positions, masses


def calculate_approx_inner_diameter_mass(positions, masses):
    """
    计算基于质心的近似内接直径
    """
    if not positions or not masses:
        return 0.0

    pos_array = np.array(positions)
    mass_array = np.array(masses)

    # --- 核心修改：计算质心 (Center of Mass) ---
    # 公式：COM = sum(mi * ri) / sum(mi)

    # mass_array[:, None] 将 (N,) 变为 (N, 1) 以便与 (N, 3) 的坐标数组进行广播乘法
    weighted_positions = pos_array * mass_array[:, None]
    total_mass = np.sum(mass_array)

    center_of_mass = np.sum(weighted_positions, axis=0) / total_mass
    # ----------------------------------------

    # 计算质心到所有原子的距离
    distances = np.linalg.norm(pos_array - center_of_mass, axis=1)

    # 找到最近原子的距离
    min_dist = np.min(distances)

    # 近似内径
    return min_dist * 2.0


def get_suggested_threshold(diameter):
    """根据内径计算建议阈值"""
    if diameter < SMALL_LIMIT:
        return SMALL_THRESH
    elif diameter > LARGE_LIMIT:
        return LARGE_THRESH
    else:
        ratio = (diameter - SMALL_LIMIT) / (LARGE_LIMIT - SMALL_LIMIT)
        return SMALL_THRESH + ratio * (LARGE_THRESH - SMALL_THRESH)


# ================= 主程序 =================

print(f"{'File':<10} | {'COM Inner Dia (Å)':<20} | {'Suggested Threshold'}")
print("-" * 60)

for i in FILE_RANGE:
    filename = f"{FILE_PREFIX}{i}.pdb"
    filepath = os.path.join(INPUT_PATH, filename)

    positions, masses = read_pdb_data(filepath)

    if positions is None:
        print(f"{filename:<10} | {'Not Found':<20} | -")
        continue

    if len(positions) == 0:
        print(f"{filename:<10} | {'Empty':<20} | -")
        continue

    # 使用新的质心计算函数
    inner_dia = calculate_approx_inner_diameter_mass(positions, masses)
    thresh = get_suggested_threshold(inner_dia)

    print(f"{filename:<10} | {inner_dia:<20.2f} | {thresh:.4f}")

print("-" * 60)
print(f"Logic (Center of Mass): Dia < {SMALL_LIMIT} -> {SMALL_THRESH}, Dia > {LARGE_LIMIT} -> {LARGE_THRESH}")