import numpy as np
import math
import os
import shutil
import glob
import scipy.ndimage
from scipy.spatial import KDTree
from pathlib import Path

from Balloon import Start_Imitation, Init_Vertex
from data import vdw_radii
from Input_methond import read_positions_and_atom_names_from_file


# ================= 1. 客体体积计算引擎 =================
# 当前文件所在目录：既适用于项目源码，也适用于 PyMOL 安装后的插件目录
MODULE_DIR = Path(__file__).resolve().parent

# 开发机上的 DSP-Cage 主项目目录。
# 可通过环境变量 DSP_CAGE_PROJECT_ROOT 覆盖，便于在其他电脑上使用。
PROJECT_ROOT = Path(
    os.environ.get("DSP_CAGE_PROJECT_ROOT", r"D:\daima\DSP_Cage")
).resolve()

# 默认资源目录。实际搜索客体文件时还会根据当前输入文件目录动态推导。
PROJECT_DATA_DIR = PROJECT_ROOT / "datas"
LOCAL_DATA_DIR = MODULE_DIR / "datas"

GUEST_DIR = PROJECT_DATA_DIR / "guests"
PROBE_DIR = PROJECT_DATA_DIR / "PROBE"

# 只确保探针输出目录存在；客体目录不存在时不自动创建空目录，
# 以免掩盖数据文件未复制的问题。
PROBE_DIR.mkdir(parents=True, exist_ok=True)
GUEST_GRID_STEP = 0.1
GUEST_PADDING = 2.0
DO_EROSION = True
EROSION_ITERATIONS = 1
VOID_CUTOFF = 10.0
VOID_TARGETS = ["B7", "B8", "B12", "B13"]

PAPER_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "B": 1.85,
    "CO": 1.80, "0": 1.52
}
DEFAULT_RADIUS = 1.70


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


def calculate_guest_volume(positions, radii, step=0.1, mode="Normal"):
    if len(positions) == 0: return 0.0
    min_coords = positions.min(axis=0) - GUEST_PADDING
    max_coords = positions.max(axis=0) + GUEST_PADDING
    dims = np.ceil((max_coords - min_coords) / step).astype(int)
    nx, ny, nz = dims
    grid_shell = np.zeros((nx, ny, nz), dtype=bool)
    x_range = min_coords[0] + np.arange(nx) * step

    for pos, r in zip(positions, radii):
        r_indices = int(np.ceil(r / step))
        idx = np.floor((pos - min_coords) / step).astype(int)
        x_start, x_end = max(0, idx[0] - r_indices), min(nx, idx[0] + r_indices + 1)
        y_start, y_end = max(0, idx[1] - r_indices), min(ny, idx[1] + r_indices + 1)
        z_start, z_end = max(0, idx[2] - r_indices), min(nz, idx[2] + r_indices + 1)
        if x_start >= x_end or y_start >= y_end or z_start >= z_end: continue
        lx = x_range[x_start:x_end]
        ly = min_coords[1] + np.arange(y_start, y_end) * step
        lz = min_coords[2] + np.arange(z_start, z_end) * step
        mx, my, mz = np.meshgrid(lx, ly, lz, indexing='ij')
        dist_sq = (mx - pos[0]) ** 2 + (my - pos[1]) ** 2 + (mz - pos[2]) ** 2
        grid_shell[x_start:x_end, y_start:y_end, z_start:z_end] |= (dist_sq <= (r ** 2))

    if mode == "Normal":
        return np.sum(grid_shell) * (step ** 3)
    elif mode == "Smart_Void":
        grid_filled = scipy.ndimage.binary_fill_holes(grid_shell)
        grid_void = grid_filled ^ grid_shell
        labeled_array, num_features = scipy.ndimage.label(grid_void)
        if num_features > 0:
            sizes = scipy.ndimage.sum(grid_void, labeled_array, range(num_features + 1))
            cutoff_voxels = VOID_CUTOFF / (step ** 3)
            mask_keeper = sizes > cutoff_voxels
            grid_void = mask_keeper[labeled_array]
        if DO_EROSION and EROSION_ITERATIONS > 0:
            grid_void = scipy.ndimage.binary_erosion(grid_void, iterations=EROSION_ITERATIONS)
        final_grid = grid_shell | grid_void
        return np.sum(final_grid) * (step ** 3)
    return 0.0


# ================= 2. 分子笼计算类 =================
class cavity():
    def __init__(self):
        self.atom_type_list = None
        self.atom_idx_dict = None
        self.dummy_atom_radii = 1
        self.distanceFromCOMFactor = 1
        self.positions = None
        self.atom_names = None
        self.atom_masses = None
        self.atom_vdw = None
        self.n_atoms = 0
        self.filename = None
        self.input_dir = None
        self.INT_MAX = 100005
        self.vdwR_dict = {}
        self.last_frame_data = None

        # [修改点1] 增加变量存储 Rebek 理论体积
        self.cached_rebek_vol = None


    def read_file(self, input_path, filename):
        """读取主体结构，并记录本次输入目录供客体文件检索使用。"""
        self.filename = filename
        self.input_dir = Path(input_path).expanduser().resolve()
        full_path = self.input_dir / filename

        print("[Input] Host file:", full_path)

        self.positions, self.atom_names, self.atom_masses, self.atom_vdw = (
            read_positions_and_atom_names_from_file(str(full_path))
        )
        self.n_atoms = len(self.positions)

    def distance_point_to_ray(self, point, ray_origin, ray_direction):
        point = np.array(point)
        ray_origin = np.array(ray_origin)
        ray_direction = np.array(ray_direction)
        t = np.dot(point - ray_origin, ray_direction) / np.dot(ray_direction, ray_direction)
        if t < 0: return self.INT_MAX
        nearest_point = ray_origin + t * ray_direction
        distance = np.linalg.norm(point - nearest_point)
        return distance

    def calculate_pore_center(self):
        pore_center = np.array(sum(self.positions[i] for i in range(self.n_atoms))) / len(self.atom_masses)
        return pore_center

    def get_guest_based_threshold(self):
        """根据与主体编号对应的客体文件计算自适应边长阈值。"""
        self.cached_rebek_vol = None

        if not self.filename:
            print("[Guest Search] Host filename is empty.")
            return None

        cage_id = Path(self.filename).stem.split("-")[0]

        # 按优先级搜索：
        # 1. 当前输入目录同级的 guests（datas/input -> datas/guests）
        # 2. 主项目 datas/guests
        # 3. 当前插件目录 datas/guests
        candidate_dirs = []

        if self.input_dir is not None:
            candidate_dirs.append(self.input_dir.parent / "guests")

        candidate_dirs.extend([
            GUEST_DIR,
            LOCAL_DATA_DIR / "guests",
        ])

        # 去重，同时保留搜索顺序
        unique_dirs = []
        seen_dirs = set()
        for guest_dir in candidate_dirs:
            resolved_dir = guest_dir.expanduser().resolve()
            key = os.path.normcase(str(resolved_dir))
            if key not in seen_dirs:
                seen_dirs.add(key)
                unique_dirs.append(resolved_dir)

        found_files = []

        print("[Guest Search] Cage ID:", cage_id)
        for guest_dir in unique_dirs:
            search_pattern = str(guest_dir / f"{cage_id}-*.pdb")
            matched = glob.glob(search_pattern)

            print("[Guest Search] Directory:", guest_dir)
            print("[Guest Search] Pattern:", search_pattern)
            print("[Guest Search] Matched:", matched)

            found_files.extend(matched)

        found_files = sorted(set(found_files))
        print("[Guest Search] Found files:", found_files)

        if not found_files:
            print(
                f"[Guest Search] No guest file found for {cage_id}; "
                "falling back to Geometry Mode."
            )
            return None

        guest_path = found_files[0]
        g_pos, g_radii = read_pdb_manual(guest_path)

        if len(g_pos) == 0:
            print(
                f"[Guest Search] Guest file contains no valid ATOM/HETATM "
                f"records: {guest_path}"
            )
            return None

        mode = "Smart_Void" if cage_id in VOID_TARGETS else "Normal"
        guest_vol = calculate_guest_volume(
            g_pos,
            g_radii,
            step=GUEST_GRID_STEP,
            mode=mode,
        )

        rebek_cage_vol = guest_vol / 0.55
        self.cached_rebek_vol = rebek_cage_vol

        if rebek_cage_vol <= 300:
            target_val = 0.4
        elif rebek_cage_vol <= 600:
            target_val = 4
        else:
            target_val = 5

        print(
            f"--- [Guest Mode] ID:{cage_id} "
            f"| File:{guest_path} "
            f"| Guest:{guest_vol:.1f} "
            f"| Rebek:{rebek_cage_vol:.1f} "
            f"| Threshold:{target_val} ---"
        )
        return target_val

    def calculate_volum_by_balloon(self, Path="", centerType=1, times=0, dynamic_threshold=None):
        from Balloon import Init_All_DATA, get_balloon_state, set_max_edge_length

        Init_All_DATA()
        # 💡 【最高优先级拦截】：如果外部传入了独立的阈值，直接使用它，不再运行下面的自动推导逻辑！
        if dynamic_threshold is not None:
            set_max_edge_length(dynamic_threshold)
        else:
            # 否则，才走原来的自动推导逻辑
            guest_threshold = self.get_guest_based_threshold()
            if guest_threshold is not None:
                set_max_edge_length(guest_threshold)
            else:
                # 几何兜底
                if self.positions is not None and len(self.positions) > 0:
                    pos_array = np.array(self.positions)
                    center_point = np.mean(pos_array, axis=0)
                    distances = np.linalg.norm(pos_array - center_point, axis=1)
                    approx_inner_diameter = np.min(distances) * 2.0

                    target_val = 0.8
                    if hasattr(self, 'edge_thresholds') and self.edge_thresholds:
                        if approx_inner_diameter < 5:
                            target_val = self.edge_thresholds["under_5"]
                        elif approx_inner_diameter <= 10:
                            target_val = self.edge_thresholds["5_to_10"]
                        elif approx_inner_diameter <= 15:
                            target_val = self.edge_thresholds["10_to_15"]
                        elif approx_inner_diameter <= 20:
                            target_val = self.edge_thresholds["15_to_20"]
                        else:
                            target_val = self.edge_thresholds["over_20"]
                    else:
                        if approx_inner_diameter <= 5:
                            target_val = 0.05
                        elif approx_inner_diameter <= 10:
                            target_val = 0.25   #####计算O1.pdb时参数采用2
                        elif approx_inner_diameter <= 15:
                            target_val = 0.5   ####默认取值0.5，取值0.45达到最优效果
                        elif approx_inner_diameter <= 20:
                            target_val = 1
                        elif approx_inner_diameter <= 25:
                            target_val = 2.5
                        elif approx_inner_diameter <= 30:
                            target_val = 3.5
                        else:
                            target_val = 15
                    print(
                        f"--- [Geometry Mode] Inner Dia={approx_inner_diameter:.2f}Å -> Threshold={target_val:.4f} ---")
                    set_max_edge_length(target_val)
        '''
        guest_threshold = self.get_guest_based_threshold()
        if guest_threshold is not None:
            set_max_edge_length(guest_threshold)
        else:
            # 几何兜底
            if self.positions is not None and len(self.positions) > 0:
                pos_array = np.array(self.positions)
                center_point = np.mean(pos_array, axis=0)
                distances = np.linalg.norm(pos_array - center_point, axis=1)
                approx_inner_diameter = np.min(distances) * 2.0

                target_val = 0.8
                
                if approx_inner_diameter < 5:
                    target_val = 0.05
                elif approx_inner_diameter > 15.0:
                    target_val = 9
                else:
                    #target_val = 0.3 + ((approx_inner_diameter - 9.5) / 0.5) * (3.0 - 0.3)
                    target_val =0.5
                
                if approx_inner_diameter < 5:
                    target_val = 0.05
                elif approx_inner_diameter <= 10:  # 包含 5 到 15 的区间
                    target_val = 0.25
                elif approx_inner_diameter <= 15:  # 包含 5 到 15 的区间
                    target_val = 0.5
                elif approx_inner_diameter <= 20:  # 包含 15 到 20 的区间
                    target_val = 1
                else:  # 即对应 > 20 的情况
                    target_val = 10
                
                # === 替换开始：使用前端传入的动态阈值，并保留默认 Fallback ===
                if hasattr(self, 'edge_thresholds') and self.edge_thresholds:
                    if approx_inner_diameter < 5:
                        target_val = self.edge_thresholds["under_5"]
                    elif approx_inner_diameter <= 10:
                        target_val = self.edge_thresholds["5_to_10"]
                    elif approx_inner_diameter <= 15:
                        target_val = self.edge_thresholds["10_to_15"]
                    elif approx_inner_diameter <= 20:
                        target_val = self.edge_thresholds["15_to_20"]
                    else:  # 即对应 > 20 的情况
                        target_val = self.edge_thresholds["over_20"]
                else:
                    # Fallback 默认值（当没有从可视化前端传入时使用）
                    if approx_inner_diameter < 5:
                        target_val = 0.05
                    elif approx_inner_diameter <= 10:
                        target_val = 0.3
                    elif approx_inner_diameter <= 15:
                        target_val = 0.5
                    elif approx_inner_diameter <= 20:
                        target_val = 1
                    else:
                        target_val = 10
                # === 替换结束 ===
                print(f"--- [Geometry Mode] Inner Dia={approx_inner_diameter:.2f}Å -> Threshold={target_val:.4f} ---")
                set_max_edge_length(target_val)
        '''
        pore_center_of_mass, pore_radius = self.calculate_center_and_radius()
        self.atom_type_list = []
        self.vdwR_dict = {}
        self.atom_idx_dict = {}

        for atom_idx, cage_name in enumerate(self.atom_names):
            atom_type = cage_name
            pos = self.positions[atom_idx]
            if atom_type not in self.atom_type_list:
                self.atom_type_list.append(atom_type)
                self.vdwR_dict.setdefault(atom_type, []).append(self.atom_vdw[atom_idx])
            self.atom_idx_dict.setdefault(atom_type, []).append(atom_idx)

        pore_center = self.calculate_pore_center()
        balloon_vertex_neighbors = Init_Vertex(times, 1, pore_center_of_mass, pore_center, centerType, self, Path,inherited_state=self.last_frame_data)

        nearest_atom2vertex = []

        for vertex_idx, vertex in enumerate(balloon_vertex_neighbors):
            nearest_atom2vertex.append(self.INT_MAX)
            for atom_idx, cage_name in enumerate(self.atom_names):
                atom_type = cage_name
                pos = self.positions[atom_idx]
                if (self.distance_point_to_ray(pos, pore_center_of_mass,
                                               [vertex.x - pore_center_of_mass[0], vertex.y - pore_center_of_mass[1],
                                                vertex.z - pore_center_of_mass[2]]) < self.vdwR_dict[atom_type][0]):
                    if (nearest_atom2vertex[vertex_idx] == self.INT_MAX or self.distance(vertex, pos) < self.distance(
                            vertex, self.positions[nearest_atom2vertex[vertex_idx]])):
                        nearest_atom2vertex.pop()
                        nearest_atom2vertex.append(atom_idx)
        #vol = Start_Imitation(self.atom_names, self.vdwR_dict, self.positions, nearest_atom2vertex, self.filename)
        #self.last_frame_data = get_balloon_state()
        #return vol
        # 修改调用和接收方式
        result_dict = Start_Imitation(self.atom_names, self.vdwR_dict, self.positions, nearest_atom2vertex,
                                      self.filename)
        # 存储窗口信息供外部查询
        self.last_window_info = result_dict

        self.last_frame_data = get_balloon_state()
        return result_dict["volume"]
###########6.8新增################
    def calculate_volum_by_balloon_dynamic(self, Path="", centerType=1, times=0):
        """
        [新增：动态轨迹专用接口]
        """
        from Balloon import run_trajectory_frame

        # 1. 计算当前帧的各项几何中心
        pore_center_of_mass, pore_radius = self.calculate_center_and_radius()
        pore_center = self.calculate_pore_center()
        if centerType == 1:
            current_center = [pore_center[0], pore_center[1], pore_center[2]]
        elif centerType == 3:
            current_center = [pore_center[0] * 2 - pore_center_of_mass[0], pore_center[1] * 2 - pore_center_of_mass[1],
                              pore_center[2] * 2 - pore_center_of_mass[2]]
        else:
            current_center = [pore_center_of_mass[0], pore_center_of_mass[1], pore_center_of_mass[2]]

        # 2. 构建当前帧的元素范德华字典
        self.vdwR_dict = {}
        for atom_idx, cage_name in enumerate(self.atom_names):
            if cage_name not in self.vdwR_dict:
                self.vdwR_dict[cage_name] = [self.atom_vdw[atom_idx]]

        # 3. 路由分流：直接投入专为轨迹设计的无循环轻量引擎
        result_dict = run_trajectory_frame(
            old_state=self.last_frame_data,
            new_center=current_center,
            atom_names=self.atom_names,
            vdw_dict=self.vdwR_dict,
            atom_positions=self.positions,
            filename=self.filename,
            path=Path
        )

        # 4. 保持链式传递，更新当前帧的状态供下一帧作为基准继承
        self.last_window_info = result_dict
        self.last_frame_data = (result_dict["mesh_neighbors"], current_center)

        return result_dict["volume"]

    def distance(self, vertex_position, atom_position):
        if type(vertex_position) != list:
            x = vertex_position.x - atom_position[0]
            y = vertex_position.y - atom_position[1]
            z = vertex_position.z - atom_position[2]
        else:
            x = vertex_position[0] - atom_position[0]
            y = vertex_position[1] - atom_position[1]
            z = vertex_position[2] - atom_position[2]
        return x * x + y * y + z * z

    def calculate_center_of_mass(self):
        pore_center_of_mass = np.array(sum(self.atom_masses[i] * self.positions[i] for i in range(self.n_atoms))) / sum(
            self.atom_masses)
        return pore_center_of_mass.tolist()
    '''
    def Calculate_Cavity(self, fileName, ball_center_type, divide_times, file_input_path="", file_output_path=""):
        fileName = fileName.split(',')
        volume_list = []
        for file in fileName:
            self.read_file(file_input_path, file)
            vol = self.calculate_volum_by_balloon(file_input_path, ball_center_type, divide_times)
            volume_list.append(vol)
            filepath = file_input_path + "" + file.split('.')[0] + "_cavity.pdb"
            filepath = filepath.replace("\\", "/")
            try:
                shutil.copy(filepath, file_output_path)
                os.remove(filepath)
            except:
                pos = filepath.rfind("/")
                if pos != -1: filepath = filepath[:pos + 1]
                if not file_output_path.endswith("/"): file_output_path += "/"
        return volume_list
    '''
    '''
    # 💡 顺便微调一下原有的 Calculate_Cavity 函数，使其能够根据开关灵活切换调用接口：
    # 💡 【修改】：函数签名增加 dynamic_threshold 参数，默认值为 None
    def Calculate_Cavity(self, fileName, ball_center_type, divide_times, file_input_path="", file_output_path="",
                         use_dynamic_engine=False, dynamic_threshold=None):
        fileName = fileName.split(',')
        volume_list = []
        for file in fileName:
            self.read_file(file_input_path, file)

            # 💡 【核心拦截逻辑】：如果上层传了专属的动态阈值，直接强制设置，跳过引擎内部的自动推导
            if dynamic_threshold is not None:
                from Balloon import set_max_edge_length
                set_max_edge_length(dynamic_threshold)

            # 根据开关选择运行模式
            if use_dynamic_engine and self.last_frame_data is not None:
                vol = self.calculate_volum_by_balloon_dynamic(file_input_path, ball_center_type, divide_times)
            else:
                vol = self.calculate_volum_by_balloon(file_input_path, ball_center_type, divide_times)

            # 💡 【补充安全保护】：如果是普通批量模式（未传动态阈值），内部的 calculate_volum_by_balloon 仍会自动推导，互不干扰。

            volume_list.append(vol)

            filepath = file_input_path + "" + file.split('.')[0] + "_cavity.pdb"
            filepath = filepath.replace("\\", "/")
            try:
                shutil.copy(filepath, file_output_path)
                os.remove(filepath)
            except:
                pass
        return volume_list
    '''

    def Calculate_Cavity(self, fileName, ball_center_type, divide_times, file_input_path="", file_output_path="",
                         use_dynamic_engine=False, dynamic_threshold=None):
        fileName = fileName.split(',')
        volume_list = []
        for file in fileName:
            self.read_file(file_input_path, file)

            # 根据开关选择运行模式
            if use_dynamic_engine and self.last_frame_data is not None:
                # 如果您后续使用了独立的无循环动态引擎，也可以把参数传给它
                vol = self.calculate_volum_by_balloon_dynamic(file_input_path, ball_center_type, divide_times)
            else:
                # 💡 【核心修改】：在这里将 dynamic_threshold 显式透传给 calculate_volum_by_balloon！
                vol = self.calculate_volum_by_balloon(file_input_path, ball_center_type, divide_times,
                                                      dynamic_threshold=dynamic_threshold)

            volume_list.append(vol)

            filepath = file_input_path + "" + file.split('.')[0] + "_cavity.pdb"
            filepath = filepath.replace("\\", "/")
            try:
                shutil.copy(filepath, file_output_path)
                os.remove(filepath)
            except:
                pass
        return volume_list
    def calculate_center_and_radius(self):
        pore_center_of_mass = np.array(sum(self.atom_masses[i] * self.positions[i] for i in range(self.n_atoms))) / sum(
            self.atom_masses)
        kdtxyzAtoms = KDTree(self.positions, leafsize=20)
        distancesFromCOM = kdtxyzAtoms.query(pore_center_of_mass, k=self.n_atoms, p=2)
        pore_radius = distancesFromCOM[0][0] * self.distanceFromCOMFactor
        pore_radius = pore_radius - 1.01 * vdw_radii[
            self.atom_names[distancesFromCOM[1][0]]] - 1.01 * self.dummy_atom_radii
        if pore_radius < 0: pore_radius = 0
        return pore_center_of_mass, pore_radius

