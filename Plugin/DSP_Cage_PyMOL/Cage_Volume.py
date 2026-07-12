import os
import glob
import numpy as np
import scipy.ndimage


# ==========================================
# 1. 核心调优区 (修改这里的参数)
# ==========================================
# [关键] 设为 True 以消除 B7(C60) 约 0.9% 的虚高
DO_EROSION = True

# [关键] 腐蚀深度。
# 1 = 剥离 0.1 埃 (推荐值)
# 0 = 不剥离 (如果结果偏小，请改为 0)
EROSION_ITERATIONS = 1

# [微调] 忽略小于 10.0 A^3 的独立空隙 (论文原话)
VOID_CUTOFF = 10.0

# [基础] 网格精度 (勿动)
GRID_STEP = 0.1

# ==========================================
# 2. 基础配置
# ==========================================
INPUT_PATH = "datas/guests/"
FILE_RANGE = range(1, 14)
PADDING = 2.0
VOID_TARGETS = ["B7", "B8", "B12", "B13"]

# 论文标准半径 (已验证完美，勿动)
PAPER_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "B": 1.85,
    "CO": 1.80, "0": 1.52
}
DEFAULT_RADIUS = 1.70


# ================= 3. 计算引擎 =================
def get_radius(atom_name, element_column):
    elem = element_column.strip().upper()
    if not elem and len(atom_name) <= 2: elem = atom_name.upper()
    return PAPER_RADII.get(elem, PAPER_RADII.get(atom_name.strip().upper(), DEFAULT_RADIUS))


def read_pdb_manual(filepath):
    positions, radii = [], []
    if not os.path.exists(filepath): return np.array([]), np.array([])
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    element_col = line[76:78]
                    atom_name = line[12:16]
                    r = get_radius(atom_name, element_col)
                    positions.append([x, y, z])
                    radii.append(r)
                except ValueError:
                    continue
    return np.array(positions), np.array(radii)


def calculate_optimized_volume(positions, radii, step=0.1, mode="Normal"):
    if len(positions) == 0: return 0.0

    min_coords = positions.min(axis=0) - PADDING
    max_coords = positions.max(axis=0) + PADDING
    dims = np.ceil((max_coords - min_coords) / step).astype(int)
    nx, ny, nz = dims

    grid_shell = np.zeros((nx, ny, nz), dtype=bool)

    x_range = min_coords[0] + np.arange(nx) * step
    y_range = min_coords[1] + np.arange(ny) * step
    z_range = min_coords[2] + np.arange(nz) * step

    for pos, r in zip(positions, radii):
        r_indices = int(np.ceil(r / step))
        idx = np.floor((pos - min_coords) / step).astype(int)

        x_start, x_end = max(0, idx[0] - r_indices), min(nx, idx[0] + r_indices + 1)
        y_start, y_end = max(0, idx[1] - r_indices), min(ny, idx[1] + r_indices + 1)
        z_start, z_end = max(0, idx[2] - r_indices), min(nz, idx[2] + r_indices + 1)

        if x_start >= x_end or y_start >= y_end or z_start >= z_end: continue

        lx = x_range[x_start:x_end]
        ly = y_range[y_start:y_end]
        lz = z_range[z_start:z_end]
        mx, my, mz = np.meshgrid(lx, ly, lz, indexing='ij')

        dist_sq = (mx - pos[0]) ** 2 + (my - pos[1]) ** 2 + (mz - pos[2]) ** 2
        grid_shell[x_start:x_end, y_start:y_end, z_start:z_end] |= (dist_sq <= (r ** 2))

    if mode == "Normal":
        return np.sum(grid_shell) * (step ** 3)

    elif mode == "Smart_Void":
        # A. 全填充
        grid_filled = scipy.ndimage.binary_fill_holes(grid_shell)

        # B. 提取空腔
        grid_void = grid_filled ^ grid_shell

        # C. 过滤微小空隙 (使用全局变量 VOID_CUTOFF)
        labeled_array, num_features = scipy.ndimage.label(grid_void)
        if num_features > 0:
            sizes = scipy.ndimage.sum(grid_void, labeled_array, range(num_features + 1))
            cutoff_voxels = VOID_CUTOFF / (step ** 3)
            mask_keeper = sizes > cutoff_voxels
            grid_void = mask_keeper[labeled_array]

        # D. 腐蚀修整 (使用全局变量 DO_EROSION 和 EROSION_ITERATIONS)
        if DO_EROSION and EROSION_ITERATIONS > 0:
            grid_void = scipy.ndimage.binary_erosion(grid_void, iterations=EROSION_ITERATIONS)

        # E. 合并
        final_grid = grid_shell | grid_void
        return np.sum(final_grid) * (step ** 3)
    return 0.0


# ================= 3. 主程序 (已修改显示 Rebek 分子笼体积) =================
def run():
    print(f"{'ID':<5} | {'Filename':<16} | {'Guest(A^3)':<10} | {'Status':<8} | {'Cage Vol (Rebek / 0.55)':<25}")
    print("-" * 80)

    TARGETS = {
        "B1": 150.32, "B2": 154.51, "B3": 136.63, "B4": 309.11,
        "B5": 49.62, "B6": 52.68, "B7": 519.2, "B8": 511.7,
        "B9": 141.16, "B10": 150.5, "B11": 306.64, "B12": 524.47,
        "B13": 617.55
    }

    for i in FILE_RANGE:
        guest_id = f"B{i}"
        search_pattern = os.path.join(INPUT_PATH, f"{guest_id}-*.pdb")
        found_files = glob.glob(search_pattern)
        if not found_files: continue

        filepath = found_files[0]
        filename = os.path.basename(filepath)
        positions, radii = read_pdb_manual(filepath)

        mode = "Smart_Void" if guest_id in VOID_TARGETS else "Normal"

        # 1. 计算客体体积
        vol = calculate_optimized_volume(positions, radii, step=GRID_STEP, mode=mode)

        # 2. 计算 Rebek 分子笼体积 ( V_host = V_guest / 0.55 )
        cage_volume = vol / 0.55

        # 状态检查
        target = TARGETS.get(guest_id, 0)
        diff_pct = ((vol - target) / target) * 100 if target else 0
        status = "OK"
        if abs(diff_pct) < 1.0: status = "PERFECT"

        # 格式化输出
        # Filename 截断显示以免太长
        fname_display = (filename[:13] + '..') if len(filename) > 15 else filename

        print(f"{guest_id:<5} | {fname_display:<16} | {vol:<10.2f} | {status:<8} | {cage_volume:<10.2f} A^3")


if __name__ == "__main__":
    run()