import math
import time
import os
import warnings
import pymeshlab as ml
import numpy as np
from Bio.PDB import *
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import pathlib  # <--- 新增：导入 pathlib
import pywindow as pw
import csv
cavity = []

# 新增：定义触发细分的最大边长（平方，单位 Å^2）
# 1.0Å 的边长是一个合理的开始值，您可以调优这个值
MAX_EDGE_LENGTH_SQUARED = 0.8

# 是否是第一次扩展
first_extension = True
# ... (其他全局变量)
# 是否是第一次扩展
first_extension = True
# 膨胀次数
extension_times = 0
# 当前体积大小
volume = 1
# 球心
balloon_center = []
# 气球顶点邻居位置
balloon_vertex_neighbors = {}
# 气球顶点邻居下标
balloon_index_neighbors = []
# 节点是否仍可扩展
balloon_vertex_extensible = {}
# 距离当前节点最近的节点
balloon_nearest_atom2vertex = []
# 存储已经画出的三角形,防止重复绘制
exist = {}
# 初始状态拷贝
balloon_vertex_copy = {}
# list 状态的 balloon_vertex_neighbors.keys()
list_neighbors = []
# 是否召回过 通过设置
balloon_vertex_recall = []
# 顶点分组
balloon_vertex_group = []

# 是否完成扩展
extension_finished = False
balloon_extension_origin = []
balloon_extension_vector = []

file_name = ""
balloon_influence = []
# 限制条件
balloon_influence_condition = []

# 原子种类
atom_names = []
# 范德华半径
vdwR_dict = []
# 原子位置
atom_positions = []
# 原子球面三角形
Faces = []
# 所有的三角形顶点
npTriangle = []

# 三角形顶点数量
vertex_count = 0
# 球体顶点数量
balloon_vertex_count = 0
# 气球顶点
balloon_triangle = []
# 原子顶点数量
atom_vertex_count = 0
# atom 三角形顶点
atom_triangle = []
# 是否更新球体顶点数量
balloon_vertex_count_update = True
# 是否更新原子顶点数量
atom_vertex_count_update = True
# 对当前节点影响最大的点
balloon_vertex_influence_min = {}
# 最大数量
INT_MAX_COUNT = 100005
# 路径
Path = ""
# ⬇️ 新增：直接锁定绝对路径 ⬇️
BASE_DIR_BALLOON = pathlib.Path(__file__).resolve().parent
PROBE_Path = str(BASE_DIR_BALLOON / "datas" / "PROBE") + os.sep
# [新增] 动态设置函数
def set_max_edge_length(val):
    global MAX_EDGE_LENGTH_SQUARED
    # 注意：refine_mesh 里比较的是距离的平方，但通常我们设置的是边长阈值本身还是它的平方？
    # 原代码变量名是 MAX_EDGE_LENGTH_SQUARED，但在 refine_mesh 中：
    # if edge_len_sq > MAX_EDGE_LENGTH_SQUARED
    # 所以传入的值应该是 "目标边长的平方" (例如边长 1.0 -> 传入 1.0; 边长 0.7 -> 传入 0.49)
    # 或者我们为了方便，约定传入的是平方值。
    MAX_EDGE_LENGTH_SQUARED = val
    print(f">>> Dynamic Threshold Adjusted: MAX_EDGE_LENGTH_SQUARED = {val:.4f}")
def Init_All_DATA():
    global cavity
    global first_extension
    global extension_times
    global volume
    global balloon_center
    global balloon_vertex_neighbors
    global balloon_index_neighbors
    global balloon_vertex_extensible
    global balloon_nearest_atom2vertex
    global exist
    global balloon_vertex_copy
    global list_neighbors
    global balloon_vertex_recall
    global balloon_vertex_group
    global extension_finished
    global balloon_extension_origin
    global balloon_extension_vector
    global file_name
    global balloon_influence
    global balloon_influence_condition
    global atom_names
    global vdwR_dict
    global atom_positions
    global Faces
    global npTriangle
    global vertex_count
    global balloon_vertex_count
    global balloon_triangle
    global atom_vertex_count
    global atom_triangle
    global balloon_vertex_count_update
    global atom_vertex_count_update
    global balloon_vertex_influence_min
    global INT_MAX_COUNT
    global Path
    # ⬇️⬇️ 新增：引入控制开关 global ⬇️⬇️
    global enable_refinement
    cavity = []
    # ⬇️⬇️ 新增：默认允许细分（为了让第一帧正常生长） ⬇️⬇️
    enable_refinement = True
    # 是否是第一次扩展
    first_extension = True
    # 膨胀次数
    extension_times = 0
    # 当前体积大小
    volume = 1
    # 球心
    balloon_center = []
    # 气球顶点邻居位置
    balloon_vertex_neighbors = {}
    # 气球顶点邻居下标
    balloon_index_neighbors = []
    # 节点是否仍可扩展
    balloon_vertex_extensible = {}
    # 距离当前节点最近的节点
    balloon_nearest_atom2vertex = []
    # 存储已经画出的三角形,防止重复绘制
    exist = {}
    # 初始状态拷贝
    balloon_vertex_copy = {}
    # list 状态的 balloon_vertex_neighbors.keys()
    list_neighbors = []
    # 是否召回过 通过设置
    balloon_vertex_recall = []
    # 顶点分组
    balloon_vertex_group = []

    # 是否完成扩展
    extension_finished = False
    balloon_extension_origin = []
    balloon_extension_vector = []

    file_name = ""
    balloon_influence = []
    # 限制条件
    balloon_influence_condition = []

    # 原子种类
    atom_names = []
    # 范德华半径
    vdwR_dict = []
    # 原子位置
    atom_positions = []
    # 原子球面三角形
    Faces = []
    # 所有的三角形顶点
    npTriangle = []

    # 三角形顶点数量
    vertex_count = 0
    # 球体顶点数量
    balloon_vertex_count = 0
    # 气球顶点
    balloon_triangle = []
    # 原子顶点数量
    atom_vertex_count = 0
    # atom 三角形顶点
    atom_triangle = []
    # 是否更新球体顶点数量
    balloon_vertex_count_update = True
    # 是否更新原子顶点数量
    atom_vertex_count_update = True
    # 对当前节点影响最大的点
    balloon_vertex_influence_min = {}
    # 最大数量
    INT_MAX_COUNT = 100005
    # 路径
    Path = ""
    PROBE_Path = str(BASE_DIR_BALLOON / "datas" / "PROBE") + os.sep


class XYZ:
    """
        function: Point结构体。
    """

    def __init__(self, x=0.0, y=0.0, z=0.0):
        if type(x) == list:
            self.x = x[0]
            self.y = x[1]
            self.z = x[2]
        else:
            self.x = x
            self.y = y
            self.z = z

    def __hash__(self):
        return hash((self.x * 10000, self.y * 10000, self.z * 10000))

    def __eq__(self, other):
        return isinstance(other, XYZ) and self.x == other.x and self.y == other.y and self.z == other.z

    def __add__(self, other):
        if type(other) == list:
            try:
                return XYZ(self.x + other[0], self.y + other[1], self.z + other[2])
            except:
                return self
        return XYZ(self.x + other.x, self.y + other.y, self.z + other.z)

    def __distance__(self, other):
        if type(other) == list:
            x = (self.x - other[0])
            y = (self.y - other[1])
            z = (self.z - other[2])
        else:
            x = (self.x - other.x)
            y = (self.y - other.y)
            z = (self.z - other.z)

        return x * x + y * y + z * z

    def __list__(self):
        return list([self.x, self.y, self.z])

# ... (在 XYZ 类定义之后)

def distance_point_to_ray(point, ray_origin, ray_direction):
    """
        计算点到射线的最短距离。
        (从 Cavity_Calculation.py 复制而来)
    """
    # 将输入转换为 NumPy 数组以便进行数学运算
    point = np.array(point)
    ray_origin = np.array(ray_origin)
    ray_direction = np.array(ray_direction)

    # 计算射线上的最近点
    t = np.dot(point - ray_origin, ray_direction) / np.dot(ray_direction, ray_direction)

    if t < 0:
        # 点在射线起点的后面，返回无穷大（或一个非常大的数）
        # 我们只关心前方的碰撞
        global INT_MAX_COUNT
        return INT_MAX_COUNT

    nearest_point = ray_origin + t * ray_direction

    # 计算点到射线上的最近点的距离
    distance = np.linalg.norm(point - nearest_point)

    return distance


# ... (放在 distance_point_to_ray 之后)

def find_nearest_atom_for_vertex_ray(vertex_pos_xyz):
    """
    为新创建的顶点计算其射线方向上最近的原子。
    (逻辑改编自 Cavity_Calculation.py)
    """
    global balloon_center, atom_positions, vdwR_dict, atom_names, INT_MAX_COUNT

    nearest_atom_index = INT_MAX_COUNT
    origin_point = balloon_center

    # 射线方向
    ray_dir = [
        vertex_pos_xyz.x - origin_point[0],
        vertex_pos_xyz.y - origin_point[1],
        vertex_pos_xyz.z - origin_point[2]
    ]
    # 归一化 (防止除零错误)
    norm = np.linalg.norm(ray_dir)
    if norm > 1e-6:
        ray_dir = ray_dir / norm
    else:
        # 如果顶点在中心，随便选个方向（这不应该发生）
        ray_dir = [1, 0, 0]

    min_collision_dist_sq = INT_MAX_COUNT * 1.0

    for atom_idx, cage_name in enumerate(atom_names):
        atom_type = cage_name
        pos = atom_positions[atom_idx]

        # 1. 检查射线是否与原子的 vdW 半径相交
        dist_to_ray = distance_point_to_ray(pos, origin_point, ray_dir)

        if dist_to_ray < vdwR_dict[atom_type][0]:
            # 2. 如果相交，检查它是否是"最近"的相交点
            # 我们关心的是沿射线的*第一个*碰撞

            # 计算原子中心在射线上的投影点
            t = np.dot(np.array(pos) - np.array(origin_point), ray_dir)
            # 碰撞距离约等于 t 减去 vdW 半径（这是一个近似值）
            # 我们用更简单的：只比较到原子中心的距离

            dist_sq_to_center = distance_vertex_atom_center(origin_point, pos)

            if dist_sq_to_center < min_collision_dist_sq:
                min_collision_dist_sq = dist_sq_to_center
                nearest_atom_index = atom_idx

    # 备用方案：如果射线逻辑没有找到（例如，在一个开放的窗口）
    # 我们退回到寻找最近的原子（无论方向）
    if nearest_atom_index == INT_MAX_COUNT:
        min_dist_sq = INT_MAX_COUNT * 1.0
        for atom_idx, pos in enumerate(atom_positions):
            dist_sq = distance_vertex_atom_center(vertex_pos_xyz, pos)
            if dist_sq < min_dist_sq:
                min_dist_sq = dist_sq
                nearest_atom_index = atom_idx

    return nearest_atom_index


# ... (放在 find_nearest_atom_for_vertex_ray 之后)

def rebuild_adjacency_lists(triangles_set):
    """
    根据一个 (i, j, k) 元组的集合，完全重建邻接字典和列表。
    """
    global balloon_vertex_neighbors, balloon_index_neighbors, list_neighbors

    # 1. 创建新的空数据结构
    #    list_neighbors 必须是最新
    list_neighbors = list(balloon_vertex_neighbors.keys())
    num_vertices = len(list_neighbors)

    new_vertex_neighbors = {key: [] for key in list_neighbors}
    new_index_neighbors = [[] for _ in range(num_vertices)]

    # 2. 填充结构
    #    save_Calculation_Result 的逻辑是 (i, neighbors[j], neighbors[j+1])
    #    所以 new_index_neighbors[i] 应该包含 [j, k, l, m, ...]
    #    其中 (i, j, k) 和 (i, l, m) 是三角形

    temp_index_neighbors = [[] for _ in range(num_vertices)]

    for (i, j, k) in triangles_set:
        temp_index_neighbors[i].extend([j, k])
        temp_index_neighbors[j].extend([i, k])
        temp_index_neighbors[k].extend([i, j])

    # 3. 赋值
    balloon_index_neighbors = temp_index_neighbors

    # 4. 填充 vertex_neighbors (字典)
    for i in range(num_vertices):
        v_i = list_neighbors[i]
        for neighbor_idx in balloon_index_neighbors[i]:
            new_vertex_neighbors[v_i].append(list_neighbors[neighbor_idx])

    balloon_vertex_neighbors = new_vertex_neighbors


# ... (放在 rebuild_adjacency_lists 之后)

def refine_mesh():
    """
    检查所有边，细分超过 MAX_EDGE_LENGTH_SQUARED 的边。
    这将就地 (in-place) 修改所有全局顶点列表。
    """
    global balloon_vertex_neighbors, balloon_index_neighbors, list_neighbors
    global balloon_vertex_extensible, balloon_nearest_atom2vertex, balloon_influence
    global balloon_influence_condition, balloon_vertex_recall, balloon_extension_origin
    global balloon_extension_vector, MAX_EDGE_LENGTH_SQUARED, balloon_center

    # 0. 刷新 list_neighbors (重要)
    list_neighbors = list(balloon_vertex_neighbors.keys())

    edges_to_split = set()
    old_triangles = set()
    midpoint_map = {}  # 映射 (idx1, idx2) -> new_idx

    # 1. 识别所有长边和旧三角形
    old_N = len(list_neighbors)
    for i in range(old_N):
        v_i = list_neighbors[i]
        neighbors_indices = balloon_index_neighbors[i]

        # 收集三角形 (i, j, k)
        for j_idx in range(0, len(neighbors_indices), 2):
            try:
                k = neighbors_indices[j_idx]
                l = neighbors_indices[j_idx + 1]
                triangle = tuple(sorted((i, k, l)))
                old_triangles.add(triangle)
            except IndexError:
                continue  # 邻接表可能暂时不完整

        # 识别长边 (i, k)
        for k in neighbors_indices:
            if k < i:  # 只检查一次 (i, k) 而不是 (k, i)
                continue
            v_k = list_neighbors[k]

            # 使用 distance_vertex_atom_center 来计算平方距离
            edge_len_sq = distance_vertex_atom_center(v_i, v_k.__list__())

            if edge_len_sq > MAX_EDGE_LENGTH_SQUARED:
                edges_to_split.add(tuple(sorted((i, k))))

    if not edges_to_split:
        return False  # 没有细分

    # print(f"Refining {len(edges_to_split)} edges...") # 用于调试

    # 2. 创建新顶点并填充其状态
    current_new_idx = old_N
    for (idx1, idx2) in edges_to_split:
        v1 = list_neighbors[idx1]
        v2 = list_neighbors[idx2]

        # a. 位置
        v_new_pos = XYZ((v1.x + v2.x) / 2, (v1.y + v2.y) / 2, (v1.z + v2.z) / 2)

        # b. 状态继承和计算
        extensible = balloon_vertex_extensible.get(idx1, True) and \
                     balloon_vertex_extensible.get(idx2, True)

        nearest_atom = find_nearest_atom_for_vertex_ray(v_new_pos)

        origin = balloon_center
        vector = [v_new_pos.x - origin[0], v_new_pos.y - origin[1], v_new_pos.z - origin[2]]

        # c. 将新顶点及其状态添加到全局列表
        new_idx = current_new_idx
        balloon_vertex_neighbors[v_new_pos] = []  # 邻居稍后重建
        balloon_vertex_extensible[new_idx] = extensible
        balloon_nearest_atom2vertex.append(nearest_atom)
        balloon_influence.append([])  # 重置影响
        balloon_influence_condition.append([])  # 重置影响
        balloon_vertex_recall.append(True)
        balloon_extension_origin.append(origin)
        balloon_extension_vector.append(vector)
        balloon_index_neighbors.append([])  # 邻居稍后重建

        midpoint_map[(idx1, idx2)] = new_idx
        current_new_idx += 1

    # 3. 重建三角形
    new_triangles = set()
    # 必须刷新 list_neighbors，因为它包含了新创建的 midpoints 坐标
    list_neighbors_current = list(balloon_vertex_neighbors.keys())
    for (i, j, k) in old_triangles:
        # 检查这个三角形的三条边
        m_ij = midpoint_map.get(tuple(sorted((i, j))))
        m_jk = midpoint_map.get(tuple(sorted((j, k))))
        m_ki = midpoint_map.get(tuple(sorted((k, i))))
        # 4 种细分情况 (添加前进行三点共线检查)

        if m_ij and m_jk and m_ki:  # 3条边分裂
            tris_to_check = [
                (i, m_ij, m_ki), (j, m_jk, m_ij),
                (k, m_ki, m_jk), (m_ij, m_jk, m_ki)
            ]
        elif m_ij and m_jk:  # 2条边 (i-j, j-k)
            tris_to_check = [
                (i, m_ij, k), (m_ij, j, m_jk),
                (k, m_ij, m_jk)
            ]
        elif m_jk and m_ki:  # 2条边 (j-k, k-i)
            tris_to_check = [
                (j, m_jk, i), (m_jk, k, m_ki),
                (i, m_jk, m_ki)
            ]
        elif m_ki and m_ij:  # 2条边 (k-i, i-j)
            tris_to_check = [
                (k, m_ki, j), (m_ki, i, m_ij),
                (j, m_ki, m_ij)
            ]
        elif m_ij:  # 1条边 (i-j)
            tris_to_check = [(i, m_ij, k), (m_ij, j, k)]
        elif m_jk:  # 1条边 (j-k)
            tris_to_check = [(j, m_jk, i), (m_jk, k, i)]
        elif m_ki:  # 1条边 (k-i)
            tris_to_check = [(k, m_ki, j), (m_ki, i, j)]
        else:  # 0条边
            tris_to_check = [(i, j, k)]
        # 检查并添加所有合法的三角形
        for t_indices in tris_to_check:
            idx_a, idx_b, idx_c = t_indices

            if is_valid_triangle_by_indices(idx_a, idx_b, idx_c, list_neighbors_current):
                new_triangles.add(tuple(sorted(t_indices)))

    # 4. 根据新三角形重建邻接表
    rebuild_adjacency_lists(new_triangles)

    # 5. 刷新 list_neighbors (非常重要)
    list_neighbors = list(balloon_vertex_neighbors.keys())

    return True  # 细分完成


def is_valid_triangle_by_indices(idx1, idx2, idx3, list_neighbors_all, tolerance=1e-6):
    """
    检查由三个顶点索引定义的三角形是否为非退化三角形 (非三点共线)。
    使用叉乘的模长作为面积检查。
    """
    global balloon_center

    # 获取坐标
    p1 = np.array(list_neighbors_all[idx1].__list__())
    p2 = np.array(list_neighbors_all[idx2].__list__())
    p3 = np.array(list_neighbors_all[idx3].__list__())

    # 检查三个点是否重合 (虽然不常见，但可以作为安全检查)
    if np.array_equal(p1, p2) or np.array_equal(p2, p3) or np.array_equal(p1, p3):
        return False

    # 计算向量
    vector12 = p2 - p1
    vector13 = p3 - p1

    # 计算叉乘，其模长是三角形面积的两倍
    cross_product = np.cross(vector12, vector13)

    # 计算叉乘模长的平方 (避免额外的开方运算)
    area_magnitude_sq = np.dot(cross_product, cross_product)

    # 如果面积的平方小于容差的平方，则认为三点共线
    if area_magnitude_sq < tolerance * tolerance:
        return False
    return True
'''
def refine_mesh():
    """
    [修复版] 检查并细分边。
    使用“遗传逻辑”确定新顶点属性：CV + UV = CV。
    这能确保口袋深处的顶点被正确标记为 CV，从而拥有停止条件。
    """
    global balloon_vertex_neighbors, balloon_index_neighbors, list_neighbors
    global balloon_vertex_extensible, balloon_nearest_atom2vertex, balloon_influence
    global balloon_influence_condition, balloon_vertex_recall, balloon_extension_origin
    global balloon_extension_vector, MAX_EDGE_LENGTH_SQUARED, balloon_center, atom_positions
    global INT_MAX_COUNT

    # 0. 刷新 list_neighbors
    list_neighbors = list(balloon_vertex_neighbors.keys())

    edges_to_split = set()
    old_triangles = set()
    midpoint_map = {}

    old_N = len(list_neighbors)

    # 1. 识别长边 (保持不变)
    for i in range(old_N):
        v_i = list_neighbors[i]
        neighbors_indices = balloon_index_neighbors[i]

        # 收集三角形
        for j_idx in range(0, len(neighbors_indices), 2):
            try:
                k = neighbors_indices[j_idx]
                l = neighbors_indices[j_idx + 1]
                triangle = tuple(sorted((i, k, l)))
                old_triangles.add(triangle)
            except IndexError:
                continue

        # 识别长边
        for k in neighbors_indices:
            if k < i: continue
            v_k = list_neighbors[k]
            edge_len_sq = distance_vertex_atom_center(v_i, v_k.__list__())
            if edge_len_sq > MAX_EDGE_LENGTH_SQUARED:
                edges_to_split.add(tuple(sorted((i, k))))

    if not edges_to_split:
        return False

    # 2. 创建新顶点并处理状态 (核心修改部分)
    current_new_idx = old_N
    for (idx1, idx2) in edges_to_split:
        v1 = list_neighbors[idx1]
        v2 = list_neighbors[idx2]

        # a. 位置 (几何中点)
        v_new_pos = XYZ((v1.x + v2.x) / 2, (v1.y + v2.y) / 2, (v1.z + v2.z) / 2)

        # b. 状态继承
        extensible = balloon_vertex_extensible.get(idx1, True) or \
                     balloon_vertex_extensible.get(idx2, True)

        # --- [核心修改：基于遗传的 UV/CV 分类] ---
        # 逻辑：只要有一个父节点是 CV，子节点就是 CV。防止口袋内部出现假性 UV。
        atom1_idx = balloon_nearest_atom2vertex[idx1]
        atom2_idx = balloon_nearest_atom2vertex[idx2]

        nearest_atom_idx = INT_MAX_COUNT

        if atom1_idx != INT_MAX_COUNT and atom2_idx != INT_MAX_COUNT:
            # 两个都是 CV，继承距离更近的那个原子
            d1 = distance_vertex_atom_center(v_new_pos, atom_positions[atom1_idx])
            d2 = distance_vertex_atom_center(v_new_pos, atom_positions[atom2_idx])
            if d1 < d2:
                nearest_atom_idx = atom1_idx
            else:
                nearest_atom_idx = atom2_idx
        elif atom1_idx != INT_MAX_COUNT:
            nearest_atom_idx = atom1_idx  # 继承 v1
        elif atom2_idx != INT_MAX_COUNT:
            nearest_atom_idx = atom2_idx  # 继承 v2
        else:
            nearest_atom_idx = INT_MAX_COUNT  # 只有两个都是 UV 时，才是 UV
        # ---------------------------------------

        origin = balloon_center
        vector = [v_new_pos.x - origin[0], v_new_pos.y - origin[1], v_new_pos.z - origin[2]]

        new_idx = current_new_idx
        balloon_vertex_neighbors[v_new_pos] = []
        balloon_vertex_extensible[new_idx] = extensible
        balloon_nearest_atom2vertex.append(nearest_atom_idx)  # 使用继承的索引
        balloon_influence.append([])
        balloon_influence_condition.append([])
        balloon_vertex_recall.append(True)
        balloon_extension_origin.append(origin)
        balloon_extension_vector.append(vector)
        balloon_index_neighbors.append([])

        midpoint_map[(idx1, idx2)] = new_idx
        current_new_idx += 1

    # 3. 重建三角形 (保持不变)
    new_triangles = set()
    for (i, j, k) in old_triangles:
        m_ij = midpoint_map.get(tuple(sorted((i, j))))
        m_jk = midpoint_map.get(tuple(sorted((j, k))))
        m_ki = midpoint_map.get(tuple(sorted((k, i))))

        if m_ij and m_jk and m_ki:
            new_triangles.add(tuple(sorted((i, m_ij, m_ki))))
            new_triangles.add(tuple(sorted((j, m_jk, m_ij))))
            new_triangles.add(tuple(sorted((k, m_ki, m_jk))))
            new_triangles.add(tuple(sorted((m_ij, m_jk, m_ki))))
        elif m_ij and m_jk:
            new_triangles.add(tuple(sorted((i, m_ij, k))))
            new_triangles.add(tuple(sorted((m_ij, j, m_jk))))
            new_triangles.add(tuple(sorted((k, m_ij, m_jk))))
        elif m_jk and m_ki:
            new_triangles.add(tuple(sorted((j, m_jk, i))))
            new_triangles.add(tuple(sorted((m_jk, k, m_ki))))
            new_triangles.add(tuple(sorted((i, m_jk, m_ki))))
        elif m_ki and m_ij:
            new_triangles.add(tuple(sorted((k, m_ki, j))))
            new_triangles.add(tuple(sorted((m_ki, i, m_ij))))
            new_triangles.add(tuple(sorted((j, m_ki, m_ij))))
        elif m_ij:
            new_triangles.add(tuple(sorted((i, m_ij, k))))
            new_triangles.add(tuple(sorted((m_ij, j, k))))
        elif m_jk:
            new_triangles.add(tuple(sorted((j, m_jk, i))))
            new_triangles.add(tuple(sorted((m_jk, k, i))))
        elif m_ki:
            new_triangles.add(tuple(sorted((k, m_ki, j))))
            new_triangles.add(tuple(sorted((m_ki, i, j))))
        else:
            new_triangles.add(tuple(sorted((i, j, k))))

    rebuild_adjacency_lists(new_triangles)
    list_neighbors = list(balloon_vertex_neighbors.keys())
    return True
'''

'''
def refine_mesh():
    """
    [实验版] 检查并细分边。
    逻辑修改：UV + CV = UV (UV 优先/窗口优先)。
    只要父节点中有一个是 UV，新节点就被视为 UV。
    """
    global balloon_vertex_neighbors, balloon_index_neighbors, list_neighbors
    global balloon_vertex_extensible, balloon_nearest_atom2vertex, balloon_influence
    global balloon_influence_condition, balloon_vertex_recall, balloon_extension_origin
    global balloon_extension_vector, MAX_EDGE_LENGTH_SQUARED, balloon_center, atom_positions
    global INT_MAX_COUNT

    # 0. 刷新 list_neighbors
    list_neighbors = list(balloon_vertex_neighbors.keys())

    edges_to_split = set()
    old_triangles = set()
    midpoint_map = {}

    old_N = len(list_neighbors)

    # 1. 识别长边 (保持不变)
    for i in range(old_N):
        v_i = list_neighbors[i]
        neighbors_indices = balloon_index_neighbors[i]

        for j_idx in range(0, len(neighbors_indices), 2):
            try:
                k = neighbors_indices[j_idx]
                l = neighbors_indices[j_idx + 1]
                triangle = tuple(sorted((i, k, l)))
                old_triangles.add(triangle)
            except IndexError:
                continue

        for k in neighbors_indices:
            if k < i: continue
            v_k = list_neighbors[k]
            edge_len_sq = distance_vertex_atom_center(v_i, v_k.__list__())
            if edge_len_sq > MAX_EDGE_LENGTH_SQUARED:
                edges_to_split.add(tuple(sorted((i, k))))

    if not edges_to_split:
        return False

    # 2. 创建新顶点并处理状态
    current_new_idx = old_N
    for (idx1, idx2) in edges_to_split:
        v1 = list_neighbors[idx1]
        v2 = list_neighbors[idx2]

        # a. 位置
        v_new_pos = XYZ((v1.x + v2.x) / 2, (v1.y + v2.y) / 2, (v1.z + v2.z) / 2)

        # b. 状态继承
        extensible = balloon_vertex_extensible.get(idx1, True) or \
                     balloon_vertex_extensible.get(idx2, True)

        # --- [核心修改：UV 优先逻辑] ---
        atom1_idx = balloon_nearest_atom2vertex[idx1]
        atom2_idx = balloon_nearest_atom2vertex[idx2]

        nearest_atom_idx = INT_MAX_COUNT

        # 只有当两个父节点 *都* 是 CV 时，新节点才是 CV
        if atom1_idx != INT_MAX_COUNT and atom2_idx != INT_MAX_COUNT:
            d1 = distance_vertex_atom_center(v_new_pos, atom_positions[atom1_idx])
            d2 = distance_vertex_atom_center(v_new_pos, atom_positions[atom2_idx])
            if d1 < d2:
                nearest_atom_idx = atom1_idx
            else:
                nearest_atom_idx = atom2_idx
        else:
            # 只要有一个是 UV (INT_MAX_COUNT)，或者两个都是 UV
            # 新节点就“感染”为 UV
            nearest_atom_idx = INT_MAX_COUNT
        # ---------------------------------------

        origin = balloon_center
        vector = [v_new_pos.x - origin[0], v_new_pos.y - origin[1], v_new_pos.z - origin[2]]

        new_idx = current_new_idx
        balloon_vertex_neighbors[v_new_pos] = []
        balloon_vertex_extensible[new_idx] = extensible
        balloon_nearest_atom2vertex.append(nearest_atom_idx)
        balloon_influence.append([])
        balloon_influence_condition.append([])
        balloon_vertex_recall.append(True)
        balloon_extension_origin.append(origin)
        balloon_extension_vector.append(vector)
        balloon_index_neighbors.append([])

        midpoint_map[(idx1, idx2)] = new_idx
        current_new_idx += 1

    # 3. 重建三角形 (保持不变)
    new_triangles = set()
    for (i, j, k) in old_triangles:
        m_ij = midpoint_map.get(tuple(sorted((i, j))))
        m_jk = midpoint_map.get(tuple(sorted((j, k))))
        m_ki = midpoint_map.get(tuple(sorted((k, i))))

        if m_ij and m_jk and m_ki:
            new_triangles.add(tuple(sorted((i, m_ij, m_ki))))
            new_triangles.add(tuple(sorted((j, m_jk, m_ij))))
            new_triangles.add(tuple(sorted((k, m_ki, m_jk))))
            new_triangles.add(tuple(sorted((m_ij, m_jk, m_ki))))
        elif m_ij and m_jk:
            new_triangles.add(tuple(sorted((i, m_ij, k))))
            new_triangles.add(tuple(sorted((m_ij, j, m_jk))))
            new_triangles.add(tuple(sorted((k, m_ij, m_jk))))
        elif m_jk and m_ki:
            new_triangles.add(tuple(sorted((j, m_jk, i))))
            new_triangles.add(tuple(sorted((m_jk, k, m_ki))))
            new_triangles.add(tuple(sorted((i, m_jk, m_ki))))
        elif m_ki and m_ij:
            new_triangles.add(tuple(sorted((k, m_ki, j))))
            new_triangles.add(tuple(sorted((m_ki, i, m_ij))))
            new_triangles.add(tuple(sorted((j, m_ki, m_ij))))
        elif m_ij:
            new_triangles.add(tuple(sorted((i, m_ij, k))))
            new_triangles.add(tuple(sorted((m_ij, j, k))))
        elif m_jk:
            new_triangles.add(tuple(sorted((j, m_jk, i))))
            new_triangles.add(tuple(sorted((m_jk, k, i))))
        elif m_ki:
            new_triangles.add(tuple(sorted((k, m_ki, j))))
            new_triangles.add(tuple(sorted((m_ki, i, j))))
        else:
            new_triangles.add(tuple(sorted((i, j, k))))

    rebuild_adjacency_lists(new_triangles)
    list_neighbors = list(balloon_vertex_neighbors.keys())
    return True
'''




def init_atom(_atom_names, _vdwR_dict, _positions):
    """
    function: 初始化元素信息
    :param _atom_names: 元素名称列表，例如 ['H', 'O', 'C']
    :param _vdwR_dict: 范德华半径字典，例如 {'H': 1.2, 'O': 1.52, 'C': 1.7}
    :param _positions: 原子三维坐标列表，例如 [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    :return: None
    """
    global atom_names
    global vdwR_dict
    global atom_positions

    # 初始化全局变量
    atom_names = _atom_names
    vdwR_dict = _vdwR_dict
    atom_positions = _positions


def update_vertex_neighbors(triangle):
    """
    function: 设计顶点邻接表。
    :param triangle: 三角面片信息，包含三个顶点
    :return: None
    """
    global balloon_vertex_neighbors, balloon_index_neighbors

    # 遍历三角形的每个顶点
    for v in triangle:
        # if v not in balloon_vertex_neighbors :
        if not balloon_vertex_neighbors.__contains__(v):
            balloon_vertex_neighbors[v] = []  # 初始化邻接表

    # 更新邻接关系
    for i in range(3):
        # 将相邻顶点添加到集合中
        balloon_vertex_neighbors[triangle[i]].append(triangle[(i + 1) % 3])
        balloon_vertex_neighbors[triangle[i]].append(triangle[(i + 2) % 3])


def find_neighborIndex(vertex_neighbor):
    """
    function: 根据位置获取顶点索引
    :param vertex_neighbor: 顶点位置，类型为 XYZ 或者其他可以比较的类型
    :return: 顶点索引，类型为 int，如果未找到则返回 None
    """
    # 遍历所有的邻接顶点
    for i, v in enumerate(balloon_vertex_neighbors):
        # 检查当前顶点是否与给定的邻接顶点相同
        if v == vertex_neighbor:
            return i  # 返回找到的索引

    return None  # 如果未找到，返回 None


def update_index_neighbors():
    """
        function: 更新邻居下标信息。
        :param: None
        :return: None
    """
    global balloon_index_neighbors, balloon_vertex_neighbors
    index = 0
    vertex_dic = {}
    for i, v in enumerate(balloon_vertex_neighbors.keys()):
        vertex_dic[v] = i
    #“对象邻接关系” 到 “索引邻接关系” 的转换工具
    for v, neighbors in balloon_vertex_neighbors.items():
        balloon_index_neighbors.append([])
        for vertex_neighbor in neighbors:
            balloon_index_neighbors[index].append(vertex_dic[vertex_neighbor])
        index += 1


def distSquare(a, b):
    """
    function: 欧氏距离计算
    :param a: 顶点 a，类型为 XYZ
    :param b: 顶点 b，类型为 XYZ
    :return: 距离，类型为 float
    """
    # 计算每个坐标轴上的差值
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z

    # 计算并返回欧氏距离
    return math.sqrt(dx * dx + dy * dy + dz * dz)


# 球体中位点
def midArcPoint(a, b):
    """
    function: 计算球面中点。
    :param a: 顶点 a，类型为 XYZ
    :param b: 顶点 b，类型为 XYZ
    :return: 球面中点，类型为 XYZ
    """
    c = XYZ(a.x + b.x, a.y + b.y, a.z + b.z)
    mod = math.sqrt(c.x * c.x + c.y * c.y + c.z * c.z)
    c.x /= mod
    c.y /= mod
    c.z /= mod
    return c


def normalization(point):
    """
    function: 对三维点进行归一化处理
    :param point: 三维坐标点，格式为列表或数组
    :return: 归一化后的三维坐标点
    """
    length = np.linalg.norm(point)
    return [point[0] / length, point[1] / length, point[2] / length]


def distance_vertex_atom_center(vertex, atom_center):
    """
        function: 计算距离
        :param vertex: 顶点坐标，可以是 XYZ 对象或坐标列表
        :param atom_center: 原子中心坐标，可以是 XYZ 对象或坐标列表
        :return: 顶点到原子中心的平方距离
    """
    if type(vertex) == XYZ and type(atom_center) != XYZ:
        x = vertex.x - atom_center[0]
        y = vertex.y - atom_center[1]
        z = vertex.z - atom_center[2]
    if type(vertex) == XYZ and type(atom_center) == XYZ:
        x = vertex.x - atom_center.x
        y = vertex.y - atom_center.y
        z = vertex.z - atom_center.z
    if type(vertex) != XYZ and type(atom_center) == XYZ:
        x = vertex[0] - atom_center.x
        y = vertex[1] - atom_center.y
        z = vertex[2] - atom_center.z
    if type(vertex) != XYZ and type(atom_center) != XYZ:
        x = vertex[0] - atom_center[0]
        y = vertex[1] - atom_center[1]
        z = vertex[2] - atom_center[2]

    return x * x + y * y + z * z


def collision_influence(influence_point, atom_center):
    """
        function: 碰撞限制扩散。
        :param influence_point: 碰撞顶点下标
        :param atom_center: 最近元素位置
        :return: None
    """
    global balloon_index_neighbors, balloon_nearest_atom2vertex
    global balloon_influence, balloon_influence_condition
    global list_neighbors
    list_neighbors = list(balloon_vertex_neighbors.keys())

    # 传递序列，无线距离的碰撞影响会进行传递的
    influence_arry = []

    # 对当前的相邻节点进行标记
    for index in balloon_index_neighbors[influence_point]:
        result = np.array(balloon_influence_condition[index].__contains__(list_neighbors[influence_point].__list__()))
        if not result.any():
            # 使用顶点进行召回
            balloon_influence[index].append(2)
            balloon_influence_condition[index].append(list_neighbors[influence_point].__list__())

        # 将碰撞的原子中心加入影响因素之中
        exis = False
        for j, tmp_list in enumerate(balloon_influence_condition[index]):
            # 影响是1 并且 已经保存过当前的节点
            if balloon_influence[index][j] == 1 and list(atom_center) == list(tmp_list):
                exis = True
                break
        if exis == False:
            # 由于邻接节点碰撞到导致的限制 使用原子球进行限制
            balloon_influence[index].append(1)
            balloon_influence_condition[index].append(list(atom_center))

            # 如果当前被影响的节点会传递到无限远处，那么这个限制应该被传递
            if balloon_nearest_atom2vertex[index] == INT_MAX_COUNT:
                if influence_arry.__contains__(index) == False:
                    influence_arry.append(index)

    # 向所有的无限远扩散点进行限制
    for index in influence_arry:
        for i in balloon_index_neighbors[index]:
            result = np.array(balloon_influence_condition[i].__contains__(list_neighbors[influence_point].__list__()))
            if not result.any():
                # 使用顶点进行召回
                balloon_influence[i].append(2)
                balloon_influence_condition[i].append(list_neighbors[influence_point].__list__())
            # 如果当前节点影响过
            exis = False
            for tmp_list in balloon_influence_condition[i]:
                if list(atom_center) == list(tmp_list):
                    exis = True
                    break
            if exis == False:
                # 添加影响因素以及因素点
                balloon_influence[i].append(1)
                balloon_influence_condition[i].append(list(atom_center))

                # 将当前影响因素添加进这个内容之中，如果用__contain__方法会更快一些 not in 会有卡顿
                if balloon_nearest_atom2vertex[i] == INT_MAX_COUNT:
                    if influence_arry.__contains__(i) == False:
                        influence_arry.append(i)


'''
def collision_influence(influence_point, atom_center):
    """
        function: 碰撞限制扩散。
        [方案二修改]: 将原子中心沿径向向外推 COLLISION_SHIFT_OFFSET。
        :param influence_point: 碰撞顶点下标
        :param atom_center: 最近元素位置
        :return: None
    """
    global balloon_index_neighbors, balloon_nearest_atom2vertex
    global balloon_influence, balloon_influence_condition
    global list_neighbors, balloon_center

    # 定义偏移量：原子中心虚拟向外推的距离（例如：每次膨胀步长的一半）
    COLLISION_SHIFT_OFFSET = 0.05

    # 1. 计算虚拟的原子中心位置 C'_A
    np_atom_center = np.array(atom_center)
    np_center = np.array(balloon_center)

    # 向量：从球心 C_ball 指向原子中心 C_A
    vector_to_atom = np_atom_center - np_center
    norm = np.linalg.norm(vector_to_atom)

    if norm < 1e-6:
        # 如果原子中心与球心重合，不进行偏移
        shifted_atom_center = np_atom_center
    else:
        # 单位向量
        unit_vector = vector_to_atom / norm
        # 沿单位向量方向，向外推 COLLISION_SHIFT_OFFSET
        shifted_atom_center = np_atom_center + unit_vector * COLLISION_SHIFT_OFFSET

    # 转换回列表格式用于存储
    shifted_atom_center_list = list(shifted_atom_center)

    list_neighbors = list(balloon_vertex_neighbors.keys())

    # 传递序列，无线距离的碰撞影响会进行传递的
    influence_arry = []

    # 对当前的相邻节点进行标记
    for index in balloon_index_neighbors[influence_point]:
        result = np.array(balloon_influence_condition[index].__contains__(list_neighbors[influence_point].__list__()))
        if not result.any():
            # 使用顶点进行召回
            balloon_influence[index].append(2)
            balloon_influence_condition[index].append(list_neighbors[influence_point].__list__())

        # 将碰撞的原子中心加入影响因素之中
        exis = False
        for j, tmp_list in enumerate(balloon_influence_condition[index]):
            # 影响是1 并且 已经保存过当前的节点
            # ⚠️ 注意：这里必须使用原始的 atom_center 来检查是否已经存在该原子墙的影响
            if balloon_influence[index][j] == 1 and list(atom_center) == list(tmp_list):
                exis = True
                break

        # ⬇️ 核心修改 1：传递 C'_A (shifted_atom_center_list) ⬇️
        # 如果当前原子墙的影响不存在，则添加新的限制
        if exis == False:
            # 由于邻接节点碰撞到导致的限制 使用原子球进行限制
            balloon_influence[index].append(1)
            balloon_influence_condition[index].append(shifted_atom_center_list)  # <--- 使用虚拟中心

            # 如果当前被影响的节点会传递到无限远处，那么这个限制应该被传递
            if balloon_nearest_atom2vertex[index] == INT_MAX_COUNT:
                if influence_arry.__contains__(index) == False:
                    influence_arry.append(index)

    # 向所有的无限远扩散点进行限制
    for index in influence_arry:
        for i in balloon_index_neighbors[index]:
            result = np.array(balloon_influence_condition[i].__contains__(list_neighbors[influence_point].__list__()))
            if not result.any():
                # 使用顶点进行召回
                balloon_influence[i].append(2)
                balloon_influence_condition[i].append(list_neighbors[influence_point].__list__())

            # 如果当前节点影响过
            exis = False
            for tmp_list in balloon_influence_condition[i]:
                # ⚠️ 注意：这里必须使用原始的 atom_center 来检查是否已经存在该原子墙的影响
                if list(atom_center) == list(tmp_list):
                    exis = True
                    break

            # ⬇️ 核心修改 2：传递 C'_A (shifted_atom_center_list) ⬇️
            if exis == False:
                # 添加影响因素以及因素点
                balloon_influence[i].append(1)
                balloon_influence_condition[i].append(shifted_atom_center_list)  # <--- 使用虚拟中心

                # 将当前影响因素添加进这个内容之中，如果用__contain__方法会更快一些 not in 会有卡顿
                if balloon_nearest_atom2vertex[i] == INT_MAX_COUNT:
                    if influence_arry.__contains__(i) == False:
                        influence_arry.append(i)
'''
def extension_influence_judgement(vertex_index, new_vertex, old_vertex):
    """
        function: 判断顶点是否可以进行扩展
        :param vertex_index: 顶点索引
        :param new_vertex: 扩展后顶点的位置
        :param old_vertex: 扩展前顶点的位置
        :return: 是否允许扩展
    """
    global balloon_influence, balloon_influence_condition
    global balloon_index_neighbors, balloon_vertex_extensible, balloon_vertex_neighbors
    global balloon_center
    global balloon_nearest_atom2vertex, INT_MAX_COUNT
    judgement_result = True

    min_distance = INT_MAX_COUNT
    neighbors_distance_sum = 0
    step_len = distance_vertex_atom_center(new_vertex, old_vertex)

    for index, condition in enumerate(balloon_influence[vertex_index]):
        # ⬇️ ----------------- 在这里插入您的新代码 ----------------- ⬇️

        # 检查我是否是一个 UV (窗口顶点)
        if balloon_nearest_atom2vertex[vertex_index] == INT_MAX_COUNT:
            # 我是一个 UV。我的工作是"逃逸"出去，
            # 以便 "Elastic Regression" (recall) 稍后能找到我。
            # 我将忽略所有从遥远原子传递过来的"条件1"约束。
            if condition == 1:
                continue  # 跳过这个约束，继续膨胀
            # (我仍然会遵守来自其他 UV 的"条件2"约束)
        # ⬆️ ----------------- 插入结束 ----------------- ⬆️
        # 根据每一个限制条件进行
        # 传递过来的限制条件
        if condition == 1:
            atom_center = balloon_influence_condition[vertex_index][index]

            # 这是对斜边的整个的影响
            distance_newvertex_atom = distance_vertex_atom_center(new_vertex, atom_center)
            distance_oldvertex_atom = distance_vertex_atom_center(old_vertex, atom_center)
            distance = distance_newvertex_atom - distance_oldvertex_atom

            neighbors_distance_sum += distance
            if distance < min_distance:
                balloon_vertex_influence_min[vertex_index] = atom_center

            judgement_result = judgement_result and distance < 0
            #“扩展能力 vs 原子总限制力” 的安全校验 ？？？？？？？？？？
            if step_len < neighbors_distance_sum:
                # judgement_result = False
                return False

    sum_radis = 0
    num = 0
    avg = INT_MAX_COUNT
    list_neighbors = list(balloon_vertex_neighbors.keys())
    for index in balloon_index_neighbors[vertex_index]:
        if balloon_vertex_extensible.keys().__contains__(index):
            if balloon_vertex_extensible[index] == False:
                sum_radis += distance_vertex_atom_center(list_neighbors[index], balloon_center)
                num += 1
        else:
            num = -INT_MAX_COUNT
            break
    # 如果大于当前的一半的点
    if num > len(balloon_index_neighbors[vertex_index]) / 2:
       return False
    return True
    # 邻接约束：如果一半以上的邻居停了，我也得停。
'''
    if num > len(balloon_index_neighbors[vertex_index]) / 2:
        # [交叉逻辑]
        if balloon_nearest_atom2vertex[vertex_index] == INT_MAX_COUNT:
            # 如果是 UV (窗口)，必须遵守邻接约束。这防止 UV 无限膨胀 (Infinite Loop Fix)
            return False
        else:
            # 如果是 CV (口袋/墙壁)，则忽略此约束。
            # 这样 CV 就能钻入 B12 的狭窄通道，直到撞到墙为止。
            pass
'''



'''
def extension_influence_judgement(vertex_index, new_vertex, old_vertex):
    """
    判断顶点是否可以进行扩展。
    [修正配置]：
    1. 必须让所有顶点（包括UV）遵守传递约束 (condition==1)，以防止无限循环。
       (B7窗口的UV会停在边缘，这没关系，recall函数会修正它们)
    2. 必须禁用邻接约束 (num > len/2)，以防止B12口袋过早停止。
       (B1的泄漏问题通过 refine_mesh 的高密度网格来缓解)
    """
    global balloon_influence, balloon_influence_condition
    global balloon_index_neighbors, balloon_vertex_extensible, balloon_vertex_neighbors
    global balloon_center
    global balloon_nearest_atom2vertex, INT_MAX_COUNT

    judgement_result = True
    min_distance = INT_MAX_COUNT
    neighbors_distance_sum = 0
    step_len = distance_vertex_atom_center(new_vertex, old_vertex)

    for index, condition in enumerate(balloon_influence[vertex_index]):
        # 传递约束：所有顶点都必须遵守，这是程序能停止的保证。
        if condition == 1:
            atom_center = balloon_influence_condition[vertex_index][index]
            distance_newvertex_atom = distance_vertex_atom_center(new_vertex, atom_center)
            distance_oldvertex_atom = distance_vertex_atom_center(old_vertex, atom_center)
            distance = distance_newvertex_atom - distance_oldvertex_atom

            neighbors_distance_sum += distance
            if distance < min_distance:
                balloon_vertex_influence_min[vertex_index] = atom_center

            judgement_result = judgement_result and distance < 0

            # 如果这一步会撞到邻居传来的原子墙，停止。
            if step_len < neighbors_distance_sum:
                return False

    # --- 邻接约束 ---
    # 之前这个约束导致 B12 卡住。我们现在在高密度动态网格下，可以安全地禁用它。
    # 只要传递约束(上面的代码)在工作，我们就不会无限循环。

    # sum_radis = 0
    # num = 0
    # list_neighbors = list(balloon_vertex_neighbors.keys())
    # for index in balloon_index_neighbors[vertex_index]:
    #     if balloon_vertex_extensible.keys().__contains__(index):
    #         if balloon_vertex_extensible[index] == False:
    #             sum_radis += distance_vertex_atom_center(list_neighbors[index], balloon_center)
    #             num += 1
    #     else:
    #         num = -INT_MAX_COUNT
    #         break

    # if num > len(balloon_index_neighbors[vertex_index]) / 2:
    #    return False
    return True
'''



def calculate_centroid(vertices):
    """
        function: 计算顶点数组中心点
        :param: 顶点数组
        :return: 中心点坐标
    """
    sum_x, sum_y, sum_z = 0.0, 0.0, 0.0
    N = len(vertices)

    # 累加顶点坐标
    for vertex in vertices:
        sum_x += vertex[0]
        sum_y += vertex[1]
        sum_z += vertex[2]

    # 计算中心顶点的坐标
    centroid_x = sum_x / N
    centroid_y = sum_y / N
    centroid_z = sum_z / N

    return (centroid_x, centroid_y, centroid_z)


def find_plane_from_points(P1, P2):
    """
        function: 计算平面方程。
        :param: 顶点，平面
        :return: （x,y,z）
    """
    x1, y1, z1 = P1
    x2, y2, z2 = P2

    # 计算法向量
    A = x2 - x1
    B = y2 - y1
    C = z2 - z1

    # 计算 D
    D = - (A * x1 + B * y1 + C * z1)

    return (A, B, C, D)


def get_plane(vertices):
    """
        function: 获取回归平面方程。
        :param: 周围顶点位置信息数组
        :return: 平面方程
    """
    global balloon_center
    center = calculate_centroid(vertices)
    plane = find_plane_from_points(center, balloon_center)
    return plane


def get_group(index):
    """
        function: 获取顶点分组。
        :param index: 顶点的下标
        :return: 顶点的分组索引
    """
    global balloon_vertex_group

    if balloon_vertex_group[index] == -1:
        return -1

    if index != balloon_vertex_group[index]:
        return get_group(balloon_vertex_group[index])
    return index


def project_point_to_plane(C, plane):
    """
        function: 求出顶点到平面的映射。
        :param C: 顶点坐标 (x, y, z)
        :param plane: 平面的参数 (A, B, C, D)，表示平面方程 Ax + By + Cz + D = 0
        :return: 投影点的坐标 (x_p, y_p, z_p)
    """
    x_c, y_c, z_c = C
    A, B, C, D = plane

    # 平面法向量
    normal = (A, B, C)

    # 计算平面法向量的模
    normal_magnitude = math.sqrt(A ** 2 + B ** 2 + C ** 2)

    # 计算点到平面的距离
    d = (A * x_c + B * y_c + C * z_c + D) / normal_magnitude

    # 计算投影点坐标
    x_p = x_c - d * (A / normal_magnitude)
    y_p = y_c - d * (B / normal_magnitude)
    z_p = z_c - d * (C / normal_magnitude)

    return (x_p, y_p, z_p)
'''
# 对于每一组数据咱们采用多顶点拟合平面的方式
def recall():
    """
        function: 球面顶点召回
        :param: None
        :return: None
    """
    global balloon_vertex_neighbors
    global list_neighbors
    global balloon_center
    global balloon_influence
    global balloon_nearest_atom2vertex
    global balloon_vertex_group
    global atom_positions

    balloon_vertex_group = []
    list_neighbors = list(balloon_vertex_neighbors.keys())

    # 判断当前是否进行过分组，如果 分组过则为 True 没分组为 False
    grouped_list = []

    # 初始化默认值
    for i, index in enumerate(balloon_nearest_atom2vertex):
        if index == INT_MAX_COUNT:
            balloon_vertex_group.append(i)
            grouped_list.append(False)
        else:
            balloon_vertex_group.append(-1)
            grouped_list.append(True)

    # 存储所有inf节点的下标
    inf_vertex_list = []

    for i, index in enumerate(balloon_nearest_atom2vertex):
        if balloon_vertex_group[i] != -1:
            inf_vertex_list.append(i)

    # 进行分组
    for i, index in enumerate(inf_vertex_list):
        # 如果是inf的
        if not grouped_list[index]:
            # 把当前节点加入list之中，然后开始扩散
            list_extension = [index]
            grouped_list[index] = True

            # 从下标为 j 这个点进行扩散
            for j in list_extension:
                neighbor_index_list = balloon_index_neighbors[j]

                # 扩散到 k
                for k in neighbor_index_list:
                    if not grouped_list[k]:
                        list_extension.append(k)
                        grouped_list[k] = True
                        balloon_vertex_group[k] = get_group(j)
    # 分组没问题
    # 分组数组
    groups = set()
    # 分组中心
    group_dic = {}
    # 分组数量
    groups_count = {}
    for i in balloon_vertex_group:
        if i != -1:
            groups.add(i)

    for i in groups:
        group_dic[i] = [0, 0, 0]
        groups_count[i] = 0

    for index, i in enumerate(balloon_vertex_group):
        if i != -1:
            group_dic[i] = [group_dic[i][0] + list_neighbors[index].x, group_dic[i][1] + list_neighbors[index].y,
                            group_dic[i][2] + list_neighbors[index].z] #####出错出错出错！！！
            groups_count[i] += 1

    for i in groups:
        group_dic[i] = [group_dic[i][0] / groups_count[i], group_dic[i][1] / groups_count[i],
                        group_dic[i][2] / groups_count[i]]

    # 每个点到中心的距离
    groups_distance = {}
    # 每个点的下标编号
    groups_index = {}
    # 临近的节点坐标的List
    groups_near = {}
    # 临近的原子
    group_near_atom = {}
    # 根据分组找他们最近的三个碰撞的节点

    for i in groups:
        groups_distance[i] = []
        groups_index[i] = []
        groups_near[i] = []
        group_near_atom[i] = set()

    # 召回后的key数值
    recall_list = []

    for j, vertex in enumerate(list_neighbors):
        recall_list.append(vertex.__list__())

    for i in groups:
        for j, vertex in enumerate(list_neighbors):
            if balloon_vertex_group[j] == -1:
                is_near = False
                for k in balloon_index_neighbors[j]:
                    if get_group(k) == i:
                        is_near = True
                        break
                if is_near:
                    groups_near[i].append(vertex.__list__())
                    group_near_atom[i].add(XYZ(list(atom_positions[balloon_nearest_atom2vertex[j]])))
                    groups_index[i].append(j)

        plane = get_plane(groups_near[i])
        cou = 0
        for j, vertex in enumerate(list_neighbors):
            if get_group(j) == i:
                try:
                    recall_list[j] = list(project_point_to_plane(recall_list[j], plane))
                except:
                    cou += 1

    new_neighbors = {}
    for index, vertex in enumerate(recall_list):
        if distance_vertex_atom_center(list_neighbors[index], balloon_center) < distance_vertex_atom_center(vertex,
                                                                                                            balloon_center):
            XYZ_vertex = list_neighbors[index]
        else:
            XYZ_vertex = XYZ(vertex[0], vertex[1], vertex[2])
        while new_neighbors.keys().__contains__(XYZ_vertex):
            XYZ_vertex = XYZ_vertex.__add__(XYZ(0.00005, 0.00005, 0.00005))
        new_neighbors[XYZ_vertex] = []

    key_list = list(new_neighbors.keys())
    for i, index_list in enumerate(balloon_index_neighbors):
        for j in index_list:
            new_neighbors[key_list[i]].append(key_list[j])

    balloon_vertex_neighbors = new_neighbors
    list_neighbors = list(balloon_vertex_neighbors.keys())
'''

'''
# 对于每一组数据咱们采用多顶点拟合平面的方式
def recall():
    """
        function: 球面顶点召回
        (重写版本：修复了召回位置未被应用的bug)
        :param: None
        :return: None
    """
    global balloon_vertex_neighbors
    global list_neighbors
    global balloon_center
    global balloon_influence
    global balloon_nearest_atom2vertex
    global balloon_vertex_group
    global atom_positions
    global balloon_index_neighbors  # <-- 确保 balloon_index_neighbors 是全局的

    balloon_vertex_group = []
    list_neighbors = list(balloon_vertex_neighbors.keys())

    # 判断当前是否进行过分组，如果 分组过则为 True 没分组为 False
    grouped_list = []

    # 初始化默认值
    for i, index in enumerate(balloon_nearest_atom2vertex):
        if index == INT_MAX_COUNT:
            balloon_vertex_group.append(i)
            grouped_list.append(False)
        else:
            balloon_vertex_group.append(-1)
            grouped_list.append(True)

    # 存储所有inf节点的下标
    inf_vertex_list = []

    for i, index in enumerate(balloon_nearest_atom2vertex):
        if balloon_vertex_group[i] != -1:
            inf_vertex_list.append(i)

    # 进行分组
    for i, index in enumerate(inf_vertex_list):
        # 如果是inf的
        if not grouped_list[index]:
            # 把当前节点加入list之中，然后开始扩散
            list_extension = [index]
            grouped_list[index] = True

            # 从下标为 j 这个点进行扩散
            for j in list_extension:
                neighbor_index_list = balloon_index_neighbors[j]

                # 扩散到 k
                for k in neighbor_index_list:
                    if not grouped_list[k]:
                        list_extension.append(k)
                        grouped_list[k] = True
                        balloon_vertex_group[k] = get_group(j)
    # 分组没问题
    # 分组数组
    groups = set()
    # 分组中心
    group_dic = {}
    # 分组数量
    groups_count = {}
    for i in balloon_vertex_group:
        if i != -1:
            groups.add(i)

    for i in groups:
        group_dic[i] = [0, 0, 0]
        groups_count[i] = 0

    for index, i in enumerate(balloon_vertex_group):
        if i != -1:
            group_dic[i] = [group_dic[i][0] + list_neighbors[index].x, group_dic[i][1] + list_neighbors[index].y,
                            group_dic[i][2] + list_neighbors[index].z]  #####出错出错出错！！！
            groups_count[i] += 1

    for i in groups:
        group_dic[i] = [group_dic[i][0] / groups_count[i], group_dic[i][1] / groups_count[i],
                        group_dic[i][2] / groups_count[i]]

    # 每个点到中心的距离
    groups_distance = {}
    # 每个点的下标编号
    groups_index = {}
    # 临近的节点坐标的List
    groups_near = {}
    # 临近的原子
    group_near_atom = {}
    # 根据分组找他们最近的三个碰撞的节点

    for i in groups:
        groups_distance[i] = []
        groups_index[i] = []
        groups_near[i] = []
        group_near_atom[i] = set()

    # 召回后的key数值
    recall_list = []

    for j, vertex in enumerate(list_neighbors):
        recall_list.append(vertex.__list__())

    for i in groups:
        for j, vertex in enumerate(list_neighbors):
            if balloon_vertex_group[j] == -1:
                is_near = False
                for k in balloon_index_neighbors[j]:
                    if get_group(k) == i:
                        is_near = True
                        break
                if is_near:
                    groups_near[i].append(vertex.__list__())
                    group_near_atom[i].add(XYZ(list(atom_positions[balloon_nearest_atom2vertex[j]])))
                    groups_index[i].append(j)

        plane = get_plane(groups_near[i])
        cou = 0
        for j, vertex in enumerate(list_neighbors):
            if get_group(j) == i:
                try:
                    recall_list[j] = list(project_point_to_plane(recall_list[j], plane))
                except:
                    cou += 1

    # --- [BUG修复] ---
    # 旧的代码到此为止计算正确，但未能应用 recall_list。
    # 下面的代码是新的，它将 recall_list 中的坐标
    # 正确地应用回 balloon_vertex_neighbors。

    new_balloon_vertex_neighbors = {}
    new_key_list = []  # 存储新的 XYZ 对象，按原始索引顺序

    # 1. 从 recall_list (它是 list 列表) 创建新的 XYZ 对象 (字典键)
    for i in range(len(list_neighbors)):
        # 获取计算出的新坐标
        recalled_coords = recall_list[i]

        # 检查原始顶点 (CV) 是否应该保留其原始位置
        # (project_point_to_plane 只修改了 UV 顶点)
        # 我们比较新旧坐标的距离，以防万一
        original_coords = list_neighbors[i].__list__()

        # 弹性回归可能会轻微移动CV，我们只取更接近中心的那个
        if distance_vertex_atom_center(recalled_coords, balloon_center) > \
                distance_vertex_atom_center(original_coords, balloon_center):
            # 召回的点 "逃逸" 得更远了，使用原始点
            new_xyz_key = list_neighbors[i]
        else:
            # 召回的点更近或未改变，使用召回的点
            new_xyz_key = XYZ(recalled_coords[0], recalled_coords[1], recalled_coords[2])

        # (关键) 防止哈希碰撞 (如果两个顶点被召回到完全相同的位置)
        while new_xyz_key in new_balloon_vertex_neighbors:
            new_xyz_key = new_xyz_key.__add__(XYZ(0.00005, 0.00005, 0.00005))

        new_balloon_vertex_neighbors[new_xyz_key] = []
        new_key_list.append(new_xyz_key)

    # 2. 使用原始的邻接拓扑 (balloon_index_neighbors)
    #    来重建新字典的邻接列表
    for i in range(len(balloon_index_neighbors)):

        # 获取与索引 i 对应的 *新* XYZ 键
        current_new_key = new_key_list[i]

        # 获取索引 i 的 *原始* 邻居索引列表
        original_neighbor_indices = balloon_index_neighbors[i]

        for neighbor_idx in original_neighbor_indices:
            # 获取那个邻居索引对应的 *新* XYZ 键
            new_neighbor_key = new_key_list[neighbor_idx]

            # 将这个新的键 (XYZ对象) 添加到邻接列表中
            new_balloon_vertex_neighbors[current_new_key].append(new_neighbor_key)

    # 3. (关键) 用包含正确召回坐标的新字典替换旧字典
    balloon_vertex_neighbors = new_balloon_vertex_neighbors
    list_neighbors = list(balloon_vertex_neighbors.keys())

    # --- [BUG修复结束] ---
'''


def recall():
    """
        function: 球面顶点召回
        [最终优化版]：
        1. 引入 DELTA_R_OFFSET 恢复体积损失。
        2. 投影到球面盖子上 (Spherical Bulge)，消除平面封口。
        3. 修复了旧版本未应用召回位置的bug。
        :param: None
        :return: None
    """
    global balloon_vertex_neighbors
    global list_neighbors
    global balloon_center
    global balloon_influence
    global balloon_nearest_atom2vertex
    global balloon_vertex_group
    global atom_positions
    global balloon_index_neighbors  # 确保 balloon_index_neighbors 是全局的

    # [新增] 几何偏移因子 (Bulge Factor): 用来恢复体积损失
    DELTA_R_OFFSET = 0.8

    balloon_vertex_group = []
    list_neighbors = list(balloon_vertex_neighbors.keys())

    # 判断当前是否进行过分组，如果 分组过则为 True 没分组为 False
    grouped_list = []

    # 初始化默认值
    for i, index in enumerate(balloon_nearest_atom2vertex):
        if index == INT_MAX_COUNT:
            balloon_vertex_group.append(i)
            grouped_list.append(False)
        else:
            balloon_vertex_group.append(-1)
            grouped_list.append(True)

    # 存储所有inf节点的下标 (即所有非CV节点的下标)
    inf_vertex_list = []

    for i, index in enumerate(balloon_nearest_atom2vertex):
        if balloon_vertex_group[i] != -1:
            inf_vertex_list.append(i)

    # 进行分组
    for i, index in enumerate(inf_vertex_list):
        # 如果是 inf 的
        if not grouped_list[index]:
            # 把当前节点加入list之中，然后开始扩散
            list_extension = [index]
            grouped_list[index] = True

            # 从下标为 j 这个点进行扩散
            for j in list_extension:
                neighbor_index_list = balloon_index_neighbors[j]

                # 扩散到 k
                for k in neighbor_index_list:
                    if not grouped_list[k]:
                        list_extension.append(k)
                        grouped_list[k] = True
                        balloon_vertex_group[k] = get_group(j)
    # 分组没问题
    groups = set()
    for i in balloon_vertex_group:
        if i != -1:
            groups.add(i)

    # (省略了原来计算 group_dic 的代码，因为后面用不到)

    groups_near = {}

    for i in groups:
        groups_near[i] = []

    # 召回后的key数值
    recall_list = []

    for j, vertex in enumerate(list_neighbors):
        recall_list.append(vertex.__list__())

    for i in groups:
        # 1. 收集每组的边界邻居 (CV)
        for j, vertex in enumerate(list_neighbors):
            if balloon_vertex_group[j] == -1:  # 找到 CV 节点
                is_near = False
                for k in balloon_index_neighbors[j]:
                    if get_group(k) == i:
                        is_near = True
                        break
                if is_near:
                    groups_near[i].append(vertex.__list__())

        # --- [核心修改 1: 计算边界的平均半径并施加偏移] ---
        avg_radius = 0.0
        count_near = len(groups_near[i])

        if count_near > 0:
            for cv_coord in groups_near[i]:
                # 计算 CV 坐标到球心的距离（非平方）
                d_sq = (cv_coord[0] - balloon_center[0]) ** 2 + \
                       (cv_coord[1] - balloon_center[1]) ** 2 + \
                       (cv_coord[2] - balloon_center[2]) ** 2
                avg_radius += math.sqrt(d_sq)
            avg_radius /= count_near

            # 施加几何偏移因子，恢复体积损失
            avg_radius += DELTA_R_OFFSET
        else:
            # 只有 UV 没有 CV 邻居，理论上不该发生，给个默认值
            avg_radius = 10.0
            # ------------------------------------

        plane = get_plane(groups_near[i])

        for j, vertex in enumerate(list_neighbors):
            if get_group(j) == i:  # 遍历当前分组中的所有 UV 节点
                try:
                    # 1. 先投影到平面 (获取正确的方向分布)
                    flat_point = list(project_point_to_plane(recall_list[j], plane))

                    # --- [核心修改 2: 径向推挤 (鼓起)] ---
                    # 计算平面点相对于球心的向量
                    vec_x = flat_point[0] - balloon_center[0]
                    vec_y = flat_point[1] - balloon_center[1]
                    vec_z = flat_point[2] - balloon_center[2]

                    dist_flat = math.sqrt(vec_x ** 2 + vec_y ** 2 + vec_z ** 2)

                    if dist_flat > 1e-6:
                        # 将点沿径向向外推，直到达到 avg_radius
                        scale = avg_radius / dist_flat
                        bulge_point = [
                            balloon_center[0] + vec_x * scale,
                            balloon_center[1] + vec_y * scale,
                            balloon_center[2] + vec_z * scale
                        ]
                        recall_list[j] = bulge_point
                    else:
                        recall_list[j] = flat_point
                    # ------------------------------------

                except:
                    pass

    # --- [BUG修复] ---
    # 应用结果回 balloon_vertex_neighbors (保持之前的修复逻辑)

    new_balloon_vertex_neighbors = {}
    new_key_list = []  # 存储新的 XYZ 对象，按原始索引顺序

    for i in range(len(list_neighbors)):
        recalled_coords = recall_list[i]

        # 默认使用召回后的新位置
        new_xyz_key = XYZ(recalled_coords[0], recalled_coords[1], recalled_coords[2])

        # 如果是 CV 节点，我们保留它的原始位置 (因为 CV 已经正确停止在原子壁了)
        if balloon_vertex_group[i] == -1:
            new_xyz_key = list_neighbors[i]

        # (关键) 防止哈希碰撞 (如果两个顶点被召回到完全相同的位置)
        while new_xyz_key in new_balloon_vertex_neighbors:
            new_xyz_key = new_xyz_key.__add__(XYZ(0.00005, 0.00005, 0.00005))

        new_balloon_vertex_neighbors[new_xyz_key] = []
        new_key_list.append(new_xyz_key)

    # 2. 重建邻接列表
    for i in range(len(balloon_index_neighbors)):
        current_new_key = new_key_list[i]
        original_neighbor_indices = balloon_index_neighbors[i]

        for neighbor_idx in original_neighbor_indices:
            new_neighbor_key = new_key_list[neighbor_idx]
            new_balloon_vertex_neighbors[current_new_key].append(new_neighbor_key)

    # 3. 替换全局字典
    balloon_vertex_neighbors = new_balloon_vertex_neighbors
    list_neighbors = list(balloon_vertex_neighbors.keys())

def extension_sphere():
    """
        function: 球面顶点扩散
        :param: None
        :return: None
    """
    # volume 指的是膨胀的倍数
    global balloon_vertex_neighbors
    global balloon_nearest_atom2vertex
    global balloon_vertex_copy
    global balloon_center
    global volume
    global balloon_vertex_count_update, balloon_index_neighbors
    global list_neighbors
    global extension_finished
    global balloon_extension_origin

    if extension_finished:
        return
    # 需要更新球体节点
    # print("balloon_center = "+str(balloon_center[0])+" "+str(balloon_center[1])+" "+str(balloon_center[2])+" ")
    balloon_vertex_count_update = True
    new_neighbors = {}
    stop_extension_count = 0
    # 判断当前点下次是否可以扩展
    for i, vertex in enumerate(balloon_vertex_neighbors.keys()):
        x = balloon_extension_origin[i][0]
        y = balloon_extension_origin[i][1]
        z = balloon_extension_origin[i][2]

        nx = x + balloon_extension_vector[i][0] * volume
        ny = y + balloon_extension_vector[i][1] * volume
        nz = z + balloon_extension_vector[i][2] * volume

        new_vertex = XYZ(nx, ny, nz)
        atom_center = []

        # 获取当前距离最近的节点
        if balloon_nearest_atom2vertex[i] != INT_MAX_COUNT:
            atom_center = atom_positions[balloon_nearest_atom2vertex[i]]

        index = balloon_nearest_atom2vertex[i]

        # 如果当前是 无限点 或者 下一步会距离节点更近
        if balloon_nearest_atom2vertex[i] == INT_MAX_COUNT or len(balloon_vertex_extensible) == i \
                or distance_vertex_atom_center(new_vertex, atom_center) < distance_vertex_atom_center(vertex,
                                                                                                      atom_center) \
                and vdwR_dict[atom_names[index]][0] * vdwR_dict[atom_names[index]][0] < distance_vertex_atom_center(
            new_vertex, atom_positions[index]):
            # 当前节点为 无限点 and 当前节点无法再次扩展
            if balloon_vertex_extensible.keys().__contains__(i) and balloon_vertex_extensible[i] == False:
                stop_extension_count += 1
                continue
            # 根据 周围节点的限制 以及 距离限制
            balloon_vertex_extensible[i] = extension_influence_judgement(i, new_vertex, vertex)

        #     如果上一次可以扩展，这次不可以扩展了，那么我们进行一次扩散
        elif balloon_vertex_extensible[i] == True:

            # 当前是因为碰撞不可以在扩展了
            balloon_vertex_extensible[i] = False
            # 当前点不可以被召回
            balloon_vertex_recall[i] = False
            # 散播碰撞信息
            collision_influence(i, atom_center)

        if balloon_vertex_extensible[i] == False:
            stop_extension_count += 1

    # 如果所有节点全部都受到了
    if stop_extension_count == len(balloon_index_neighbors):
        extension_finished = True
        recall()
        return False

    vertex_index = 0
    # 扩大球体，更新邻居信息
    for vertex, neighbors in balloon_vertex_neighbors.items():
        new_vertex = []
        if balloon_vertex_extensible[vertex_index]:

            x = balloon_extension_origin[vertex_index][0]
            y = balloon_extension_origin[vertex_index][1]
            z = balloon_extension_origin[vertex_index][2]

            nx = x + balloon_extension_vector[vertex_index][0] * volume
            ny = y + balloon_extension_vector[vertex_index][1] * volume
            nz = z + balloon_extension_vector[vertex_index][2] * volume
            # 新的节点位置
            new_vertex = XYZ(nx, ny, nz)

            # 在这里先不管相邻节点,让相邻节点自己去更新,等整个节点更新完了,在通过索引更新相邻节点
        else:
            new_vertex = vertex

        # 添加相邻节点
        '''
        while new_vertex in new_neighbors:
        #添加微小偏移（避免坐标完全重合）
            new_vertex = XYZ(new_vertex.x + 1e-8,new_vertex.y + 1e-8,new_vertex.z + 1e-8)
        '''
        new_neighbors[new_vertex] = []
        tmp_List = list(neighbors)
        for i, neighbor in enumerate(tmp_List):
            neighbor_index = balloon_index_neighbors[vertex_index][i]

            if balloon_vertex_extensible[neighbor_index]:
                x = balloon_extension_origin[i][0]
                y = balloon_extension_origin[i][1]
                z = balloon_extension_origin[i][2]

                nx = x + balloon_extension_vector[neighbor_index][0] * volume
                ny = y + balloon_extension_vector[neighbor_index][1] * volume
                nz = z + balloon_extension_vector[neighbor_index][2] * volume
                new_neighbor = XYZ(nx, ny, nz)

            else:
                new_neighbor = neighbor
            new_neighbors[new_vertex].append(new_neighbor)
        vertex_index += 1

    balloon_vertex_neighbors = new_neighbors
    list_neighbors = list(balloon_vertex_neighbors.keys())
    return True

def Init_balloon_influence():
    """
    function: 初始化气球的影响力扩散信息
    :return: None
    """
    global balloon_influence, balloon_vertex_neighbors
    global balloon_extension_origin
    global balloon_extension_vector
    global balloon_center
    list_neighbors = balloon_vertex_neighbors.keys()

    for i, item in enumerate(list_neighbors):
        balloon_influence.append([])
        balloon_influence_condition.append([])

        # 添加扩散起点
        balloon_extension_origin.append([balloon_center[0], balloon_center[1], balloon_center[2]])
        x = item.x - balloon_center[0]
        y = item.y - balloon_center[1]
        z = item.z - balloon_center[2]
        # 添加扩散方向
        balloon_extension_vector.append([x, y, z])





# 初始化球的顶点
def Init_Vertex(time, r, mess_center, center, balloon_center_type, cage_cavity, path,inherited_state=None):
    """
    function: 初始化球的顶点
    :param time: 初始三角形递归细分的次数
    :param r: 球的半径
    :param mess_center: 质心坐标
    :param center: 球心坐标
    :param balloon_center_type: 球心类型
    :param cage_cavity: 空腔
    :return: None
    """

    # Dictionary to store neighboring vertices
    global balloon_vertex_neighbors
    global balloon_vertex_copy
    global list_neighbors
    global cavity
    global balloon_center
    global Path
    # [新增] 引入需要的全局变量
    global balloon_index_neighbors
    # ⬇️⬇️ 新增：引入开关 global ⬇️⬇️
    global enable_refinement

    Path = path
    if not path.endswith('/') and not path.endswith('\\'):
        path += '/'
    if not os.path.exists(PROBE_Path):
        os.makedirs(PROBE_Path)

    cavity = cage_cavity
    # 三角形数列
    triangles = []
    '''
    # 球心
    if balloon_center_type == 1:
        balloon_center = [center[0], center[1], center[2]]
    # 对称点
    elif balloon_center_type == 3:

        balloon_center = [center[0] * 2 - mess_center[0], center[1] * 2 - mess_center[1],
                          center[2] * 2 - mess_center[2]]
    # 质心
    else:
        balloon_center = [mess_center[0], mess_center[1], mess_center[2]]
    '''
    # 确定当前帧的球心
    if balloon_center_type == 1:
        current_center = [center[0], center[1], center[2]]
    elif balloon_center_type == 3:
        current_center = [center[0] * 2 - mess_center[0], center[1] * 2 - mess_center[1],
                          center[2] * 2 - mess_center[2]]
    else:
        current_center = [mess_center[0], mess_center[1], mess_center[2]]

    balloon_center = current_center

    # === [插入修改] 继承逻辑开始 ===
    if inherited_state is not None:
        print("Inheriting vertex data from previous frame...")
        old_neighbors, old_center = inherited_state

        # 计算位移向量：新中心 - 旧中心
        shift_vector = [
            current_center[0] - old_center[0],
            current_center[1] - old_center[1],
            current_center[2] - old_center[2]
        ]

        new_neighbors_dict = {}

        # 遍历旧的顶点，将其平移到新位置
        # 注意：这里需要处理 XYZ 对象的运算
        for old_v, old_n_list in old_neighbors.items():
            # 1. 平移顶点
            new_x = old_v.x + shift_vector[0]
            new_y = old_v.y + shift_vector[1]
            new_z = old_v.z + shift_vector[2]

            # 2. [可选] 冲突检测：如果平移后的点在原子内部，需要处理
            # 简单策略：如果发生冲突，暂时不处理，依赖后续 run_process 中的逻辑
            # 或者将其强制拉回中心附近（略复杂，暂略）

            new_v_obj = XYZ(new_x, new_y, new_z)

            # 3. 平移邻居列表
            new_n_list_objs = []
            for n in old_n_list:
                nx = n.x + shift_vector[0]
                ny = n.y + shift_vector[1]
                nz = n.z + shift_vector[2]
                new_n_list_objs.append(XYZ(nx, ny, nz))

            new_neighbors_dict[new_v_obj] = new_n_list_objs

        # 更新全局变量
        balloon_vertex_neighbors = new_neighbors_dict

        # 必须重新生成 list_neighbors 和 index_neighbors，因为对象变了
        list_neighbors = list(balloon_vertex_neighbors.keys())
        balloon_index_neighbors = []  # 需清空重算
        update_index_neighbors()  # 利用现有的 update_index_neighbors 重建索引

        # 初始化辅助变量 (copy, recall 等)
        for i, vertex in enumerate(list_neighbors):
            balloon_vertex_recall.append(True)
            # 更新 copy 偏移量
            balloon_vertex_copy[i] = [vertex.x - balloon_center[0], vertex.y - balloon_center[1],
                                      vertex.z - balloon_center[2]]

        # 初始化影响因子
        Init_balloon_influence()
        # ================= ⬇️ 新增这两行 ⬇️ =================
        global first_extension, volume
        first_extension = False  # 告诉程序：这不是第一次，不要重置我的体积！
        volume = 0.3  # 设置初始倍率为 0.9 (下一帧+0.1变回1.0，即保持当前大小开始检查)
        # 3. ⬇️⬇️ 核心修改：禁止后续帧加点！ ⬇️⬇️
        enable_refinement = False
        # ================= ⬆️ 新增结束 ⬆️ =================
        return balloon_vertex_neighbors  # 继承完成，直接返回
    # === [插入修改] 继承逻辑结束 ===
    #  初始化数值
    triangles.append([XYZ(r, 0, 0), XYZ(0, r, 0), XYZ(0, 0, r), int(time)])

    if len(balloon_vertex_neighbors) != 0:
        return
    while len(triangles) > 0:
        t = triangles[0]
        triangles.pop(0)
        if t[0].__eq__(t[1]) or t[1].__eq__(t[2]) or t[2].__eq__(t[0]):
            continue

        if t[3] != 0:
            d = midArcPoint(t[0], t[1])
            e = midArcPoint(t[1], t[2])
            f = midArcPoint(t[2], t[0])
            triangles.append([t[0], f, d, t[3] - 1])
            triangles.append([t[1], d, e, t[3] - 1])
            triangles.append([t[2], e, f, t[3] - 1])
            triangles.append([d, e, f, t[3] - 1])
        else:
            vertex_point_1 = XYZ(t[0].x, t[0].y, t[0].z).__add__(balloon_center)
            vertex_point_2 = XYZ(t[1].x, t[1].y, t[1].z).__add__(balloon_center)
            vertex_point_3 = XYZ(t[2].x, t[2].y, t[2].z).__add__(balloon_center)
            update_vertex_neighbors([vertex_point_1, vertex_point_2, vertex_point_3])

            # 第二象限
            vertex_point_1 = XYZ(-t[0].x, t[0].y, t[0].z).__add__(balloon_center)
            vertex_point_2 = XYZ(-t[1].x, t[1].y, t[1].z).__add__(balloon_center)
            vertex_point_3 = XYZ(-t[2].x, t[2].y, t[2].z).__add__(balloon_center)
            update_vertex_neighbors([vertex_point_1, vertex_point_2, vertex_point_3])

            # 第三象限
            vertex_point_1 = XYZ(-t[0].x, -t[0].y, t[0].z).__add__(balloon_center)
            vertex_point_2 = XYZ(-t[1].x, -t[1].y, t[1].z).__add__(balloon_center)
            vertex_point_3 = XYZ(-t[2].x, -t[2].y, t[2].z).__add__(balloon_center)
            update_vertex_neighbors([vertex_point_1, vertex_point_2, vertex_point_3])

            #
            # # 第四象限
            vertex_point_1 = XYZ(t[0].x, -t[0].y, t[0].z).__add__(balloon_center)
            vertex_point_2 = XYZ(t[1].x, -t[1].y, t[1].z).__add__(balloon_center)
            vertex_point_3 = XYZ(t[2].x, -t[2].y, t[2].z).__add__(balloon_center)
            update_vertex_neighbors([vertex_point_1, vertex_point_2, vertex_point_3])

            # # 第五象限
            vertex_point_1 = XYZ(t[0].x, -t[0].y, -t[0].z).__add__(balloon_center)
            vertex_point_2 = XYZ(t[1].x, -t[1].y, -t[1].z).__add__(balloon_center)
            vertex_point_3 = XYZ(t[2].x, -t[2].y, -t[2].z).__add__(balloon_center)
            update_vertex_neighbors([vertex_point_1, vertex_point_2, vertex_point_3])

            # # 第六象限
            vertex_point_1 = XYZ(-t[0].x, -t[0].y, -t[0].z).__add__(balloon_center)
            vertex_point_2 = XYZ(-t[1].x, -t[1].y, -t[1].z).__add__(balloon_center)
            vertex_point_3 = XYZ(-t[2].x, -t[2].y, -t[2].z).__add__(balloon_center)
            update_vertex_neighbors([vertex_point_1, vertex_point_2, vertex_point_3])

            # # 第七象限
            vertex_point_1 = XYZ(-t[0].x, t[0].y, -t[0].z).__add__(balloon_center)
            vertex_point_2 = XYZ(-t[1].x, t[1].y, -t[1].z).__add__(balloon_center)
            vertex_point_3 = XYZ(-t[2].x, t[2].y, -t[2].z).__add__(balloon_center)
            update_vertex_neighbors([vertex_point_1, vertex_point_2, vertex_point_3])

            # # 第八象限
            vertex_point_1 = XYZ(t[0].x, t[0].y, -t[0].z).__add__(balloon_center)
            vertex_point_2 = XYZ(t[1].x, t[1].y, -t[1].z).__add__(balloon_center)
            vertex_point_3 = XYZ(t[2].x, t[2].y, -t[2].z).__add__(balloon_center)
            update_vertex_neighbors([vertex_point_1, vertex_point_2, vertex_point_3])

    for i, vertex in enumerate(list(balloon_vertex_neighbors.keys())):
        # 当前点可以被召回
        balloon_vertex_recall.append(True)
        if (type(vertex) == list):
            vertex = XYZ(vertex[0], vertex[1], vertex[2])
        #记录顶点相对中心的偏移（为膨胀方向提供基准）
        balloon_vertex_copy[i] = [vertex.x - balloon_center[0], vertex.y - balloon_center[1],
                                  vertex.z - balloon_center[2]]

    list_neighbors = list(balloon_vertex_neighbors.keys())
    # 更新邻居信息
    update_index_neighbors()
    # 初始化影响
    Init_balloon_influence()
    return balloon_vertex_neighbors


def read_obj_file(file_path):
    """
    function: 从 OBJ 文件中读取顶点和面片数据。
    :param file_path: 文件路径
    :return: 顶点数组和面片数组
    """
    vertices = []
    faces = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('v '):
                # 读取顶点
                parts = line.split()
                vertex = [float(parts[1]), float(parts[2]), float(parts[3])]
                vertices.append(vertex)
            elif line.startswith('f '):
                # 读取面片
                parts = line.split()
                # OBJ 文件的面片索引从1开始，因此需要减去1
                face = [int(part.split('/')[0]) - 1 for part in parts[1:]]
                faces.append(face)

    return np.array(vertices), np.array(faces)

'''
没使用
def Init_atom_vertex(r, position):
    """
    function: 初始化顶点定点信息
    :param time: 迭代次数，用于细分三角形
    :param r: 气球半径
    :param mess_center: 质心位置 (x, y, z)
    :param center: 球心位置 (x, y, z)
    :param balloon_center_type: 气球中心类型 (1: 球心, 2: 质心, 3: 对称点)
    :param cage_cavity: 笼空腔信息
    :return: None
    """
    floor_num = 4
    pointCount = 6
    Vertexs = []
    step = 2 * r / floor_num
    nowfloor = 0
    angle_step = 360 / pointCount

    if len(Vertexs) == 0:
        while nowfloor < floor_num:
            nowFloor_dep = (nowfloor + 1) * step - r
            r_floor = math.sqrt(r * r - nowFloor_dep * nowFloor_dep)
            angle = 0
            Vertexs.append([])
            # 制造顶点
            while angle < 360:
                if (360 - angle < 0.05):
                    break
                x = r_floor * math.sin(angle * math.pi / 180)
                z = r_floor * math.cos(angle * math.pi / 180)
                Vertexs[nowfloor].append((x + position[0], nowFloor_dep + position[1], z + position[2]))
                angle += angle_step
            nowfloor += 1

    for i, vertex in enumerate(Vertexs):
        length = len(vertex)
        for j, vertex_point in enumerate(vertex):
            # 与最下面的顶点链接
            if i == 0:
                Faces.append([vertex_point, vertex[(j + 1) % length], [position[0], -r + position[1], position[2]]])
            else:
                Faces.append([vertex_point, Vertexs[i - 1][j], Vertexs[i - 1][(j + 1) % length]])
                Faces.append([vertex_point, vertex[(j + 1) % length], vertex[(j + 2) % length]])
            # 顶层的顶点链接
            if i == len(Vertexs) - 1:
                Faces.append([vertex_point, vertex[(j + 1) % length], [position[0], r + position[1], position[2]]])
'''

def compute_normal_vector(point1, point2, point3):
    """
        function: 计算由三个点定义的平面的法向量
        :param point1: 第一个点坐标
        :param point2: 第二个点坐标
        :param point3: 第三个点坐标
        :return: 单位法向量
    """
    vector1 = np.array(point2) - np.array(point1)
    vector2 = np.array(point3) - np.array(point1)
    normal_vector = np.cross(vector1, vector2)
    sum = abs(normal_vector[0]) + abs(normal_vector[1]) + abs(normal_vector[2])
    return [normal_vector[0] / sum, normal_vector[1] / sum, normal_vector[2] / sum]


def calculate_angle(a, b):
    """
        function: 计算两个向量之间的角度
        :param a: 第一个向量
        :param b: 第二个向量
        :return: 两个向量之间的角度（度数）
    """
    a_norm = np.sqrt(np.sum(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]))
    b_norm = np.sqrt(np.sum(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]))
    cos_value = np.dot(a, b) / (a_norm * b_norm)
    arc_value = np.arccos(cos_value)
    angle = arc_value * 180 / np.pi

    if angle > 270:
        angle -= 360

    return angle


def get_extent_vector(point_1, point_2, point_3, origin_point, now_point, result_type=list):
    """
        function: 计算法向量是否与扩散方向相同
        :param: 三角形平面a,b,c,射线起点，射线重点
        :return: True/False
    """
    point_1 = point_1.__list__()
    point_2 = point_2.__list__()
    point_3 = point_3.__list__()
    normal_vector = compute_normal_vector(point_1, point_2, point_3)
    ray_vector = [now_point[0] - origin_point[0], now_point[1] - origin_point[1], now_point[2] - origin_point[2]]
    angle = calculate_angle(normal_vector, ray_vector)
    if angle <= 90 and angle >= -90:
        if result_type == bool:
            return True
        return normal_vector
    if result_type == bool:
        return False
    return [normal_vector[0] * -1, normal_vector[1] * -1, normal_vector[2] * -1]


def save_Calculation_Result():
    """
        function: 计算结果
        :param: None
        :return: None
    """
    global file_name
    global list_neighbors
    global balloon_index_neighbors
    global balloon_center
    global balloon_vertex_neighbors
    global PROBE_Path


    list_neighbors = list(balloon_vertex_neighbors.keys())
    # --- ⬇️ 插入新的函数调用 ⬇️ ---
    pdb_probe_path = PROBE_Path + file_name[:-4] + "_PROBE.pdb"
    write_final_probe_pdb(list_neighbors, pdb_probe_path, balloon_center)
    print(f"Final probe vertex count saved to: {pdb_probe_path}")
    # --- ⬆️ 插入结束 ⬆️ ---
    # --- 在这里添加代码 ---
    print(f"Final probe vertex count: {len(list_neighbors)}")
    # --- 添加结束 ---
    faces = []

    list_neighbor = []
    for i in list_neighbors:
        point = i.__list__()
        list_neighbor.append([point[0] - balloon_center[0], point[1] - balloon_center[1], point[2] - balloon_center[2]])
    vertices = np.array(list_neighbor)

    index = 0

    all_vertices = []
    all_faces = []
    # 遍历每个原子
    new_vertices = []

    for i, vertex in enumerate(vertices):
        vx, vy, vz = vertex
        new_vertices.append([vx, vy, vz])

    all_vertices.extend(new_vertices)

    face_set = set()
    for index, neighbors_list in enumerate(balloon_index_neighbors):
        for j, neighbors_index in enumerate(neighbors_list):
            if j % 2 == 0:
                tmp = [index, neighbors_list[j], neighbors_list[j + 1]]
                tmp.sort()
                tmp_XYZ = XYZ(tmp)
                if not face_set.__contains__(tmp_XYZ):  # tmp_XYZ not in face_set
                    if get_extent_vector(list_neighbors[index], list_neighbors[neighbors_list[j]],
                                         list_neighbors[neighbors_list[j + 1]], balloon_center,
                                         list_neighbors[index].__list__(), bool):
                        faces.append([index, neighbors_list[j], neighbors_list[j + 1]])
                    else:
                        faces.append([neighbors_list[j + 1], neighbors_list[j], index])
                    face_set.add(tmp_XYZ)

    # 创建新的网格对象并设置顶点和面
    tmp_name = file_name[:-4]
    pdb_out_path = Path
    obj_out_path = Path + "/OBJ/"
    obj_out_path = ""

    vertices = np.asarray(new_vertices, dtype=np.float32)
    faces = np.asarray(faces, dtype=np.int32)

    ms = ml.MeshSet()
    mesh = ml.Mesh(vertices, faces)
    ms.add_mesh(mesh, 'my_mesh')
    ms.save_current_mesh(obj_out_path + tmp_name + '.obj')

    name = tmp_name
    obj_file_path = obj_out_path + name + '.obj'
    pdb_file_path = pdb_out_path + name + '_Cavity.pdb'

    vertices = parse_obj(obj_file_path)
    write_pdb(vertices, pdb_file_path)

    vertices, faces = read_obj_file(obj_out_path + name + '.obj')
    #print(f"volume: {volume_of_mesh(vertices, faces)}")
    final_volume = volume_of_mesh(vertices, faces)  # 1. 这里定义变量
    print(f"volume: {final_volume}")  # 2. 这里打印变量

    os.remove(obj_out_path + name + '.obj')
    # [新增] 返回体积值
    return final_volume
    # 创建pml文件
    # name = name+"_Cavity"
    # pml_name = Path+tmp_name + '.pml'
    # with open(pml_name, 'w') as file:
    #     file.write(f"load " + tmp_name + ".pdb\n")
    #     file.write(f"load " + name +'.pdb' + "\n")
    #
    #     file.write(f"select " + name + "\n")
    #     file.write(f"hide everything, " + name + "\n")
    #     file.write(f"show surface, " + name + "\n")
    #     file.write(f"cmd.color_deep(\"white\", '" + name + "', 0)\n")
    #     file.write(f"util.cba(33,\"" + tmp_name + "\",_self=cmd)\n")


def parse_obj(obj_file_path):
    """
        function: 解析obj文件。
        :param: 问及那路径
        :return: None
    """
    global balloon_center
    vertices = []

    with open(obj_file_path, 'r') as file:
        for line in file:
            if line.startswith('v '):
                parts = line.strip().split()
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                nor = normalization([x, y, z])
                vertices.append((x + balloon_center[0] - 1.7 * nor[0], y + balloon_center[1] - 1.7 * nor[1],
                                 z + balloon_center[2] - 1.7 * nor[2]))
    return vertices


def write_pdb(vertices, pdb_file_path):
    global balloon_center
    with open(pdb_file_path, 'w') as file:
        for i, (x, y, z) in enumerate(vertices, start=1):
            file.write(f"ATOM  {i:5d}  C   UNK     1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")


def write_final_probe_pdb(vertices_list, pdb_file_path, center_offset):
    """
    function: 将最终召回后的探针顶点位置坐标输出为 PDB 文件。
    :param vertices_list: 最终的 XYZ 对象列表 (即 list_neighbors)。
    :param pdb_file_path: PDB 文件保存路径。
    :param center_offset: 球心坐标 (用于参考，但坐标本身是绝对坐标)。
    :return: None
    """
    with open(pdb_file_path, 'w') as file:
        file.write("REMARK CMCC FINAL PROBE VERTICES\n")
        file.write("REMARK This file contains the final probe vertex positions used for volume calculation.\n")

        # 定义 PDB 记录的固定字段
        record_name = "HETATM"
        res_name = "SPB"  # 探针残基名 (Surface Probe)
        element = "P"  # 元素符号
        chain_id = "A"  # 链标识符

        for i, xyz_vertex in enumerate(vertices_list, start=1):
            # 提取 XYZ 坐标
            x = xyz_vertex.x
            y = xyz_vertex.y
            z = xyz_vertex.z

            # 格式化 PDB 字段
            atom_num = f"{i:5d}"
            atom_name = "SP"
            res_seq = "1"

            # 坐标字段 (31-38, 39-46, 47-54) 必须是 8.3f 格式，以确保 PyMOL 对齐
            x_coord = f"{x:8.3f}"
            y_coord = f"{y:8.3f}"
            z_coord = f"{z:8.3f}"

            occupancy = "1.00"
            temp_factor = "0.00"

            # 写入 PDB 行 (确保各字段对齐，尤其是坐标)
            file.write(
                f"{record_name:<6}{atom_num:>5}  {atom_name:<3} {res_name:<3} {chain_id:<1}{res_seq:>4}    "
                f"{x_coord}{y_coord}{z_coord}{occupancy:>6}{temp_factor:>6}          {element:>2}\n"
            )
        file.write("END\n")


def calulate_volme():
    """
        function: 球面单次扩展。
        :param: None
        :return: None
    """
    global volume
    global balloon_nearest_atom2vertex
    global list_neighbors
    global balloon_center

    distance = INT_MAX_COUNT
    min_atom_index = 0

    for i, vertex in enumerate(list_neighbors):
        # 获取当前距离最近的节点
        if balloon_nearest_atom2vertex[i] != INT_MAX_COUNT:
            atom_center = atom_positions[balloon_nearest_atom2vertex[i]]
            tmp_dis = distance_vertex_atom_center(balloon_center, atom_center)
            if math.sqrt(tmp_dis) < math.sqrt(distance):
                distance = tmp_dis
                min_atom_index = balloon_nearest_atom2vertex[i]

    atom_type = atom_names[min_atom_index]
    volume = math.sqrt(distance) - vdwR_dict[atom_type][0] - 0.1
    volume = "{:.1f}".format(volume)
    volume = float(volume)
    # try:
    extension_sphere()
    # except:
    #     volume = 1.1
    #     extension_sphere()


def signed_volume_of_triangle(v1, v2, v3):
    """
    function: 计算三角面片体积。
    :param: 顶点1，2，3
    :return: 三角面体积
    """
    return np.dot(np.cross(v1, v2), v3) / 6.0


def volume_of_mesh(vertices, faces):
    """
        function: 计算由顶点和面片定义的封闭三角网格的体积。
        :param: 顶点坐标，面片索引
        :return: 体积
    """
    volume = 0.0
    for face in faces:
        if len(face) == 4:
            # 将四边形分解为两个三角形
            v1 = vertices[face[0]]
            v2 = vertices[face[1]]
            v3 = vertices[face[2]]
            v4 = vertices[face[3]]
            volume += abs(signed_volume_of_triangle(v1, v2, v3))
            volume += abs(signed_volume_of_triangle(v1, v3, v4))
        elif len(face) == 3:
            # 处理三角形
            v1 = vertices[face[0]]
            v2 = vertices[face[1]]
            v3 = vertices[face[2]]
            volume += abs(signed_volume_of_triangle(v1, v2, v3))
        else:
            raise ValueError("Unsupported face with more than 4 vertices")
    return volume


def run_process():
    """
    function: 计算由顶点和面片定义的封闭三角网格的体积。
    :param: None
    :return: None
    """
    global first_extension
    global extension_times
    global extension_finished
    global volume
    # ⬇️ [新增] 设置最大迭代次数，防止无限循环 ⬇️
    MAX_ITERATIONS = 100
    # ⬆️ ------------------------------------- ⬆️
    # ⬇️⬇️ 新增：引入开关 global ⬇️⬇️
    global enable_refinement

    start_time = time.time()
    while extension_finished == False:
        # ⬇️ [新增] 紧急刹车检查 ⬇️
        if extension_times > MAX_ITERATIONS:
            print(f"WARNING: Force stopping at {MAX_ITERATIONS} iterations!")
            extension_finished = True
            recall()  # 强制召回
            break
        # ⬆️ ------------------------- ⬆️
        if first_extension:
            calulate_volme()
            first_extension = False
            # break
        else:
            volume = volume + 0.1
            extension_sphere()
            if extension_finished == False:
                extension_times += 1
        # --- 新增的细分步骤 ---
        # 每 5 次膨胀检查一次，并且在非首次膨胀后
        if enable_refinement and not extension_finished and not first_extension and extension_times % 1 == 0:
        # print(f"Checking for refinement at step {extension_times}...") # 调试
            refine_mesh()
        # refine_mesh() 会就地更新所有列表，
        # 下一个 extension_sphere() 循环将自动使用新顶点。
        # --- 细分步骤结束 ---

    #save_Calculation_Result()
    # [修改] 接收并返回结果
    final_volume = save_Calculation_Result()

    end_time = time.time()
    execution_time = end_time - start_time
    print("iteration times: " + str(extension_times))
    print(f"execution time: {execution_time} Second")
    print()
    return final_volume

''''
#没用
def read_pdb_file(file_name):
    """
        function: 纠正PDB位置信息。
        :param: 文件名称
        :return: None
        """
    global Path

    file_name = Path + file_name
    atom_info = []
    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_name = line[12:16].strip()
                element = line[76:78].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                atom_info.append({
                    'name': atom_name,
                    'element': element,
                    'coord': (x, y, z)
                })

    return atom_info
'''
'''
def Start_Imitation(atom_names, vdwR_dict, positions, nearest_atom2vertex, fileName):
    global balloon_nearest_atom2vertex
    global file_name
    file_name = fileName
    balloon_nearest_atom2vertex = nearest_atom2vertex

    if file_name != "":
        init_atom(atom_names, vdwR_dict, positions)
        # for atom_idx, cage_name in enumerate(atom_names):
        #     atom_type = cage_name
        #     # Init_atom_vertex(vdwR_dict[atom_type][0], positions[atom_idx])
    #run_process()
    # [修改] 接收体积值
    vol = run_process()
    #Init_All_DATA()
    # [新增] 返回体积值
    return vol
'''
def get_balloon_state():
    """
    新增函数：获取当前气球的状态（顶点邻接表和球心）
    """
    global balloon_vertex_neighbors, balloon_center
    # 返回深拷贝或当前状态的引用（由于下一帧会重建，这里返回引用即可，但为了安全建议复制逻辑在接收端处理）
    # 注意：这里我们返回 balloon_vertex_neighbors 字典和 balloon_center 列表
    return balloon_vertex_neighbors, balloon_center

warnings.simplefilter('ignore', PDBConstructionWarning)


def analyze_windows():
    """
    分析气球逃逸产生的窗口数量和每个窗口的近似直径。
    逻辑：对所有 UV 顶点（无最近原子的顶点）进行聚类。
    """
    global balloon_nearest_atom2vertex, balloon_index_neighbors, list_neighbors, INT_MAX_COUNT

    # 1. 筛选出所有 UV (窗口/逃逸) 顶点的索引
    uv_indices = [i for i, atom_idx in enumerate(balloon_nearest_atom2vertex) if atom_idx == INT_MAX_COUNT]

    if not uv_indices:
        return 0, []

    # 2. 聚类：使用 BFS 寻找连通分量 (每个簇代表一个窗口)
    visited = set()
    window_clusters = []

    for idx in uv_indices:
        if idx not in visited:
            cluster = []
            queue = [idx]
            visited.add(idx)
            while queue:
                curr = queue.pop(0)
                cluster.append(curr)
                # 在邻居中寻找同样是 UV 且未访问的点
                for neighbor_idx in balloon_index_neighbors[curr]:
                    if neighbor_idx in uv_indices and neighbor_idx not in visited:
                        visited.add(neighbor_idx)
                        queue.append(neighbor_idx)
            # 过滤噪声：点数太少的簇可能是网格破碎产生的伪窗口
            if len(cluster) > 3:
                window_clusters.append(cluster)

    window_count = len(window_clusters)
    window_details = []

    # 3. 计算每个窗口的几何特征
    for i, cluster in enumerate(window_clusters):
        # 获取该窗口所有顶点的 3D 坐标
        coords = np.array([list_neighbors[idx].__list__() for idx in cluster])

        # 计算质心
        centroid = np.mean(coords, axis=0)

        # 计算直径：寻找该簇中距离质心最远的两个点之间的距离（或边缘最远两点）
        # 这里采用一种稳健的估计：计算所有点到质心平均距离的 2 倍
        dist_to_centroid = np.linalg.norm(coords - centroid, axis=1)
        avg_radius = np.mean(dist_to_centroid)
        max_dist = 0

        # 精确直径：计算簇内任意两点间的最大距离
        if len(coords) < 200:  # 采样计算以保证性能
            for p1 in coords:
                dists = np.linalg.norm(coords - p1, axis=1)
                max_dist = max(max_dist, np.max(dists))
        else:
            max_dist = avg_radius * 2  # 大规模簇使用平均直径估算

        window_details.append({
            "id": i + 1,
            "vertex_count": len(cluster),
            "diameter": max_dist,
            "centroid": centroid.tolist()
        })

    return window_count, window_details


def Start_Imitation(atom_names, vdwR_dict, positions, nearest_atom2vertex, fileName):
    """
    [修改版] 返回体积和窗口信息
    """
    global balloon_nearest_atom2vertex, file_name
    file_name = fileName
    balloon_nearest_atom2vertex = nearest_atom2vertex

    if file_name != "":
        init_atom(atom_names, vdwR_dict, positions)

    vol = run_process()

    # --- 新增：分析窗口 ---
    win_count, win_details = analyze_windows()

    # 返回包含体积和窗口数据的字典
    return {
        "volume": vol,
        "window_count": win_count,
        "windows": win_details
    }

'''
def run_trajectory_frame(old_state, new_center, atom_names, vdw_dict, atom_positions, filename, path):
    """
    [新增：动态轨迹专用引擎]
    1. 点数和拓扑结构 100% 不变。
    2. 仅根据新一帧的原子坐标微调探针位置，防止原子穿透。
    3. 彻底跳过 while 膨胀循环，从根本上杜绝体积叠加 Bug！
    """
    import math
    import numpy as np

    old_neighbors, old_center = old_state

    # 计算当前帧球心相对于上一帧球心的平移位移
    shift = [new_center[0] - old_center[0], new_center[1] - old_center[1], new_center[2] - old_center[2]]

    old_keys = list(old_neighbors.keys())
    new_keys = []
    key_mapping = {}

    # 1. 1对1映射顶点位置，保证点数绝不发生任何变化
    for v in old_keys:
        # 先随整体结构做刚体平移
        new_v = XYZ(v.x + shift[0], v.y + shift[1], v.z + shift[2])

        # 【表面贴合微调】：寻找当前点在新一帧里距离最近的原子中心
        min_dist_sq = 100005.0
        nearest_atom_idx = 0
        for idx, pos in enumerate(atom_positions):
            d_sq = (new_v.x - pos[0]) ** 2 + (new_v.y - pos[1]) ** 2 + (new_v.z - pos[2]) ** 2
            if d_sq < min_dist_sq:
                min_dist_sq = d_sq
                nearest_atom_idx = idx

        vdw_r = vdw_dict[atom_names[nearest_atom_idx]][0]
        curr_dist = math.sqrt(min_dist_sq)

        # 如果由于分子热运动，该顶点深深地陷入了新一帧的原子范德华内部
        if curr_dist < vdw_r:
            # 计算从原子中心指向该顶点的排斥向量
            vx = new_v.x - atom_positions[nearest_atom_idx][0]
            vy = new_v.y - atom_positions[nearest_atom_idx][1]
            vz = new_v.z - atom_positions[nearest_atom_idx][2]
            v_len = math.sqrt(vx * vx + vy * vy + vz * vz)
            if v_len > 1e-6:
                # 将其优雅地推回到该原子范德华表面外侧 0.02 Å 处，实现位置稍微改变
                new_v.x = atom_positions[nearest_atom_idx][0] + (vx / v_len) * (vdw_r + 0.02)
                new_v.y = atom_positions[nearest_atom_idx][1] + (vy / v_len) * (vdw_r + 0.02)
                new_v.z = atom_positions[nearest_atom_idx][2] + (vz / v_len) * (vdw_r + 0.02)

        new_keys.append(new_v)
        key_mapping[v] = new_v

    # 2. 完美复制拓扑邻接字典
    new_neighbors_dict = {}
    for v in old_keys:
        new_neighbors_dict[key_mapping[v]] = [key_mapping[n] for n in old_neighbors[v]]

    # 3. 重新绑定全局变量，直接调用已有的构面与体积计算引擎
    global balloon_vertex_neighbors, balloon_center, file_name, Path, balloon_index_neighbors, list_neighbors,PROBE_Path
    balloon_vertex_neighbors = new_neighbors_dict
    balloon_center = new_center
    file_name = filename
    Path = path
    if not os.path.exists(PROBE_Path):
        os.makedirs(PROBE_Path)
    list_neighbors = list(balloon_vertex_neighbors.keys())

    # 刷新索引映射表
    balloon_index_neighbors = []
    vertex_dic = {v: i for i, v in enumerate(list_neighbors)}
    for v, neighbors in balloon_vertex_neighbors.items():
        balloon_index_neighbors.append([vertex_dic[n] for n in neighbors])

    # 🔥 【核心避坑点】：直接保存并计算当前网格体积，由于不走 while 循环，体积绝不叠加翻倍！
    final_volume = save_Calculation_Result()

    return {
        "volume": final_volume,
        "window_count": 0,
        "windows": [],
        "mesh_neighbors": new_neighbors_dict
    }
'''
def run_trajectory_frame(old_state, new_center, atom_names, vdw_dict, atom_positions, filename, path):
    """
    [新增：动态轨迹专用引擎]
    1. 点数和拓扑结构 100% 不变。
    2. 仅根据新一帧的原子坐标微调探针位置，防止原子穿透。
    3. 彻底跳过 while 膨胀循环，从根本上杜绝体积叠加 Bug！
    """
    import math
    import numpy as np

    old_neighbors, old_center = old_state

    # 计算当前帧球心相对于上一帧球心的平移位移
    shift = [new_center[0] - old_center[0], new_center[1] - old_center[1], new_center[2] - old_center[2]]

    old_keys = list(old_neighbors.keys())
    new_keys = []
    key_mapping = {}

    # 1. 1对1映射顶点位置，保证点数绝不发生任何变化
    for v in old_keys:
        # 先随整体结构做刚体平移
        new_v = XYZ(v.x + shift[0], v.y + shift[1], v.z + shift[2])

        # 【表面贴合微调】：寻找当前点在新一帧里距离最近的原子中心
        min_dist_sq = 100005.0
        nearest_atom_idx = 0
        for idx, pos in enumerate(atom_positions):
            d_sq = (new_v.x - pos[0]) ** 2 + (new_v.y - pos[1]) ** 2 + (new_v.z - pos[2]) ** 2
            if d_sq < min_dist_sq:
                min_dist_sq = d_sq
                nearest_atom_idx = idx

        vdw_r = vdw_dict[atom_names[nearest_atom_idx]][0]
        curr_dist = math.sqrt(min_dist_sq)

        # 如果由于分子热运动，该顶点深深地陷入了新一帧的原子范德华内部
        if curr_dist < vdw_r:
            # 计算从原子中心指向该顶点的排斥向量
            vx = new_v.x - atom_positions[nearest_atom_idx][0]
            vy = new_v.y - atom_positions[nearest_atom_idx][1]
            vz = new_v.z - atom_positions[nearest_atom_idx][2]
            v_len = math.sqrt(vx * vx + vy * vy + vz * vz)
            if v_len > 1e-6:
                # 将其优雅地推回到该原子范德华表面外侧 0.02 Å 处，实现位置稍微改变
                new_v.x = atom_positions[nearest_atom_idx][0] + (vx / v_len) * (vdw_r + 0.02)
                new_v.y = atom_positions[nearest_atom_idx][1] + (vy / v_len) * (vdw_r + 0.02)
                new_v.z = atom_positions[nearest_atom_idx][2] + (vz / v_len) * (vdw_r + 0.02)

        new_keys.append(new_v)
        key_mapping[v] = new_v

    # 2. 完美复制拓扑邻接字典
    new_neighbors_dict = {}
    for v in old_keys:
        new_neighbors_dict[key_mapping[v]] = [key_mapping[n] for n in old_neighbors[v]]

    # 3. 重新绑定全局变量，直接调用已有的构面与体积计算引擎
    global balloon_vertex_neighbors, balloon_center, file_name, Path, balloon_index_neighbors, list_neighbors,PROBE_Path
    balloon_vertex_neighbors = new_neighbors_dict
    balloon_center = new_center
    file_name = filename
    Path = path
    if not os.path.exists(PROBE_Path):
        os.makedirs(PROBE_Path)
    list_neighbors = list(balloon_vertex_neighbors.keys())

    # 刷新索引映射表
    balloon_index_neighbors = []
    vertex_dic = {v: i for i, v in enumerate(list_neighbors)}
    for v, neighbors in balloon_vertex_neighbors.items():
        balloon_index_neighbors.append([vertex_dic[n] for n in neighbors])

    # 🔥 保存并计算当前网格体积
    final_volume = save_Calculation_Result()

    # 🔴 【核心修复】：在这里挂载窗口检测函数！
    # 因为拓扑顶点索引没有变化，它会完美识别出原本是窗口的点，并计算它们在新一帧的新直径。
    win_count, win_details = analyze_windows()

    return {
        "volume": final_volume,
        "window_count": win_count,    # 替换掉写死的 0
        "windows": win_details,       # 替换掉写死的 []
        "mesh_neighbors": new_neighbors_dict
    }

