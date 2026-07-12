import pathlib
import os
import sys
import logging
import matplotlib.pyplot as plt
import numpy as np  # 用于计算 MAPE
import time  # 确保在文件开头引入 time 模块
import json
# 在 Main_Process.py 的头部导入区域添加这一行
from generate_pymol_maps import write_mapped_pdb
BASE_DIR = pathlib.Path(__file__).resolve().parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))
os.chdir(BASE_DIR)
os.environ.setdefault("BABEL_DATADIR", str(BASE_DIR))

from Cavity_Calculation import cavity
import traj_tool
from calculate_ESP_MHP import (
    calculate_esp_for_files,
    calculate_mhp_for_files,
    detect_main_metal_charge,
)

logging.getLogger("hydrophobicity").setLevel(logging.ERROR)

# ================= ⚙️ 配置区域 ⚙️ =================
IS_TRAJECTORY_MODE = False
ENABLE_INHERITANCE = True  # 💡 新增控制开关：True 代表开启时间继承加速 / False 代表关闭继承（完全独立计算）
Path = str(BASE_DIR / "datas" / "input") + os.sep
Output_Path = str(BASE_DIR / "datas" / "output") + os.sep
PROBE_Path = str(BASE_DIR / "datas" / "PROBE") + os.sep
os.makedirs(Output_Path, exist_ok=True)
os.makedirs(PROBE_Path, exist_ok=True)  # <--- 新增：确保运行时文件夹被物理创建
ball_center_type = 2
subdivision_time = 4
charge_method = "eem"
hydroph_method = "Ghose"
distance_function = "Fauchere"
# 💡 【新增】：动态轨迹专用边长阈值变量（单独调节旋钮）
# 填入你期望的独立阈值（例如 0.49 代表目标边长约 0.7Å）
DYNAMIC_EDGE_THRESHOLD = 1
cav = cavity()
fileName = ""
volume_data = []
plot_title = ""
x_label = ""


def calculate_fields_for_cage(filename):
    molecule_path = os.path.join(Path, filename)
    stem, _ = os.path.splitext(filename)
    cavity_path = os.path.join(Output_Path, f"{stem}_cavity.pdb")

    if not os.path.exists(cavity_path):
        print(f"[Field Warning] Probe file not found for {filename}: {cavity_path}")
        return {"status": "missing_probe", "esp": None, "mhp": None}

    metal, metal_charge = detect_main_metal_charge(molecule_path)
    if metal and metal_charge is not None:
        print(f"[Field] {filename}: metal charge override {metal} = {metal_charge}")
    else:
        print(f"[Field] {filename}: metal charge override none")

    try:
        esp_result = calculate_esp_for_files(
            molecule_path,
            cavity_path,
            method=charge_method,
            metal_name=metal,
            metal_charge=metal_charge,
        )
    except Exception as exc:
        print(f"[Field Error] ESP calculation failed for {filename}: {type(exc).__name__}: {exc}")
        return {"status": "esp_failed", "esp": None, "mhp": None}

    esp_values = np.asarray(esp_result["esp"], dtype=float)
    print(f"Mean ESP: {np.mean(esp_values):.6f}")
    print(f"Min/Max ESP: {np.min(esp_values):.6f} / {np.max(esp_values):.6f}")

    # ================= ✨ 更改位置 A：即使 MHP 失败，也先为 ESP 生成可视化 PDB =================
    esp_vis_path = os.path.join(Output_Path, f"{stem}_ESP_vis.pdb")
    try:
        write_mapped_pdb(cavity_path, esp_vis_path, esp_values, field_name="ESP")
    except Exception as vis_exc:
        print(f"[Vis Warning] Failed to write ESP PDB: {vis_exc}")
    # =======================================================================================

    try:
        mhp_result = calculate_mhp_for_files(
            molecule_path,
            cavity_path,
            method=hydroph_method,
            distance_function=distance_function,
        )
    except Exception as exc:
        print(f"[Field Warning] MHP calculation failed for {filename}; ESP result printed only.")
        print(f"[Field Warning] MHP error: {type(exc).__name__}: {exc}")
        return {
            "status": "esp_only",
            "esp": float(np.mean(esp_values)) if len(esp_values) else None,
            "mhp": None,
            "hi": None,
        }

    mhp_values = np.asarray(mhp_result["mhp"], dtype=float)
    mhp_pos = mhp_values[mhp_values > 0]
    mhp_neg = mhp_values[mhp_values < 0]
    hi_denominator = np.sum(mhp_pos) - np.sum(mhp_neg)
    hydrophobic_index = np.sum(mhp_pos) / hi_denominator if hi_denominator != 0 else np.nan

    print(f"Mean MHP: {np.mean(mhp_values):.6f}")
    print(f"Min/Max MHP: {np.min(mhp_values):.6f} / {np.max(mhp_values):.6f}")
    print(f"Hydrophobic index (HI): {hydrophobic_index:.6f}")

    # ================= ✨ 更改位置 B：MHP 成功，生成 MHP 的可视化 PDB =================
    mhp_vis_path = os.path.join(Output_Path, f"{stem}_MHP_vis.pdb")
    try:
        write_mapped_pdb(cavity_path, mhp_vis_path, mhp_values, field_name="MHP")
    except Exception as vis_exc:
        print(f"[Vis Warning] Failed to write MHP PDB: {vis_exc}")
    # =======================================================================================

    return {
        "status": "esp_mhp",
        "esp": float(np.mean(esp_values)) if len(esp_values) else None,
        "mhp": float(np.mean(mhp_values)) if len(mhp_values) else None,
        "hi": float(hydrophobic_index),
    }

# ================= 🚀 逻辑处理 🚀 =================

if IS_TRAJECTORY_MODE:
    # 💡 优化 1：让提示信息根据你的开关动态切换
    status_str = "🚀已启用时间继承/加速版" if ENABLE_INHERITANCE else "🐢无继承/独立计算计时版"
    print(f">>>以此模式运行: [动态轨迹分析 - {status_str}]")

    my_history_file = os.path.join(Path, "HISTORY_singlemol_short")
    pdb_list = traj_tool.split_trajectory_to_pdbs(
        input_filepath=my_history_file,
        output_directory=Path,
        swap_atoms={"he": "H"},
        forcefield="opls"
    )
    print(f"成功切分并获取了 {len(pdb_list)} 帧文件。")
    file_name_list = [pathlib.Path(p).name for p in pdb_list]

    volume_data = []
    time_data = []  # ⌚ 用于存储每一帧的耗时

    print("\n" + "=" * 30 + " 开始计算 " + "=" * 30)
    total_start_time = time.time()  # 记录总开始时间

    for f in file_name_list:
        frame_start_time = time.time()  # 记录当前帧开始时间

        # 💡 修改 2：只有当【未开启继承】时，才在循环内部强制清空上一帧的数据，打断继承
        if not ENABLE_INHERITANCE and hasattr(cav, 'last_frame_data'):
            cav.last_frame_data = None

        # 逐帧进行计算（如果 ENABLE_INHERITANCE 为 True，cav.last_frame_data 会保留并传递给下一帧）
        v_list = cav.Calculate_Cavity(f, ball_center_type, subdivision_time, Path, Output_Path, use_dynamic_engine=ENABLE_INHERITANCE,dynamic_threshold=DYNAMIC_EDGE_THRESHOLD)

        if v_list:
            volume_data.append(v_list[0])
        else:
            volume_data.append(0.0)

        frame_end_time = time.time()  # 记录当前帧结束时间
        frame_cost = frame_end_time - frame_start_time
        time_data.append(frame_cost)  # 存入列表供表格打印

        print(f"⏱️ 进度: [{f}] 处理完成 (耗时: {frame_cost:.2f} s)")

    total_end_time = time.time()
    print(f"\n✅ 轨迹分析全部完成！总耗时: {total_end_time - total_start_time:.2f} 秒\n")
    plot_title = 'Cavity Volume Change (Dynamic Trajectory)'
    x_label = 'Frame Number'

    # ================= 📊 结果展示 (表格输出) 📊 =================
    table_width = 50
    print("\n\n" + "=" * table_width)
    print(f"{'Framename':<15} | {'Volume (Å³)':<15} | {'Time (s)':<12}")
    print("-" * table_width)

    if volume_data:
        # 使用 zip 同时遍历 文件名、体积、耗时 三个列表
        for fname, vol, cost in zip(file_name_list, volume_data, time_data):
            print(f"{fname.split('.')[0]:<15} | {vol:<15.2f} | {cost:<12.2f}")
    else:
        print("未获取到体积数据，请检查前置解析步骤。")

    print("=" * table_width)

    # 3. 绘图
    plt.figure(figsize=(10, 6))
    plt.plot(volume_data, marker='o', linestyle='-', color='b', label='Cavity Volume')
    plt.title(plot_title)
    plt.xlabel(x_label)
    plt.ylabel('Volume (Å³)')
    plt.grid(True)
    plt.legend()
    save_name = "volume_dynamic.png"
    plt.savefig(os.path.join(Output_Path, save_name))

else:
    print(">>> 以此模式运行: [普通 PDB 批量计算]")

    # 你的文件列表
    #target_files = ["B1.pdb", "B2.pdb", "B3.pdb", "B4.pdb", "B5.pdb", "B6.pdb", "B7.pdb", "B8.pdb", "B9.pdb", "B10.pdb",
                   #"B11.pdb", "B12.pdb", "B13.pdb"]
    #target_files = [ "B4.pdb", "B11.pdb"]
    #target_files = [ "B7.pdb", "B8.pdb","B12.pdb"]
    #target_files = ["B1.pdb", "B2.pdb", "B3.pdb",  "B5.pdb", "B6.pdb","B9.pdb", "B10.pdb"]
    #target_files = ["B1.pdb", "B2.pdb", "B3.pdb"]
    #target_files =["cage_1_JACS_2006_128_14120.mol2"]
    #target_files = ["frame_0.pdb","frame_1.pdb"]
    #target_files = [f"B{i}.pdb" for i in range(1, 14)]
    target_files = ["A1.pdb","C1.pdb","F1.pdb","F2.pdb","H1.pdb","N1.pdb","O1.pdb","W1.pdb","O2.pdb","C13.pdb", "C14.pdb", "C16.pdb"]


    target_files =["B1.pdb"]
    #target_files = ["C1.pdb"]
    #target_files = ["O1.pdb","H1.pdb"]
    #target_files = ["B2.pdb"]
    #target_files = ["C7_new.pdb"]
    #target_files = ["C16.pdb"]
    #target_files = ["C13.pdb"]
    #target_files = ["C14.pdb"]
    #target_files = ["C13.pdb", "C14.pdb", "C16.pdb"]
    print(f"待处理文件列表: {len(target_files)} 个文件")

    full_volume_list = []

    # [关键修改] 用于存储最终表格数据的列表
    summary_table_data = []
    error_values = []
    window_data_dict = {}  # ✨ 26.6.8新增：初始化用于存放窗口数据的空字典

    print("\n" + "=" * 30 + " Calculation Started " + "=" * 30)

    for f in target_files:
        # ⚠️ 强制清空状态
        if hasattr(cav, 'last_frame_data'): cav.last_frame_data = None

        # 计算
        v_list = cav.Calculate_Cavity(f, ball_center_type, subdivision_time, Path, Output_Path)
        field_result = calculate_fields_for_cage(f)

        if v_list:
            vol = v_list[0]
            full_volume_list.append(vol)
            # 获取窗口信息
            win_info = cav.last_window_info

            print(f"\n>>> 结果分析 [{f}]:")
            print(f"    空腔体积: {vol:.2f} Å³")
            print(f"    窗口数量: {win_info['window_count']}")

            for win in win_info['windows']:
                print(f"    - 窗口 #{win['id']}: 直径 = {win['diameter']:.2f} Å, 顶点数 = {win['vertex_count']}")
            # ✨ 26.6.8新增：提取该笼子所有窗口的直径并存入字典
            cage_name = f.split('.')[0]
            diameters = [win['diameter'] for win in win_info['windows']]
            window_data_dict[cage_name] = diameters


            # [关键] 从 cav 对象中读取刚才算出的 Rebek 理论体积
            rebek_vol = cav.cached_rebek_vol

            # 计算误差
            diff_str = "N/A"
            rebek_str = "N/A"

            if rebek_vol is not None:
                diff_pct = ((vol - rebek_vol) / rebek_vol) * 100
                error_values.append(abs(diff_pct))  # 存绝对值算 MAPE

                diff_str = f"{diff_pct:+.2f}%"
                rebek_str = f"{rebek_vol:.2f}"

            # 收集数据到列表，暂不打印
            summary_table_data.append({
                'name': f.split('.')[0],
                'calc': vol,
                'rebek': rebek_str,
                'diff': diff_str,
                'esp': field_result.get('esp'),
                'mhp': field_result.get('mhp'),
                'hi': field_result.get('hi'),
                'field_status': field_result.get('status'),
                'win_count': win_info['window_count']  # ✨ 新增：存储窗口数量

            })
    # ================= ✨ 在这里添加 JSON 导出代码 ✨ =================
    import json

    data_save_path = os.path.join(Output_Path, "window_data.json")
    with open(data_save_path, 'w', encoding='utf-8') as f:
        json.dump(window_data_dict, f, indent=4)

    print(f"✅ 窗口数据已导出至: {data_save_path}，可用于独立绘图。")
    # =================================================================

    volume_data = full_volume_list
    plot_title = 'Cavity Volume Comparison'
    x_label = 'File Index'


    # ================= 📊 结果展示 📊 =================
    '''
    # 1. 统一打印最终表格
    print("\n\n" + "=" * 75)
    print(f"{'Filename':<15} | {'Calc Vol':<12} | {'Rebek Vol':<12} | {'Diff %':<10}")
    print("-" * 75)
    
    
    
    for row in summary_table_data:
        print(f"{row['name']:<15} | {row['calc']:<12.2f} | {row['rebek']:<12} | {row['diff']:<10}")
    
    
    print("-" * 75)
    '''
    # 1. 统一打印最终表格 (调整了宽度以容纳 Windows 列)
    column_width = 85  # 增加总宽度
    header = f"{'Filename':<12} | {'Calc Vol':<12} | {'Rebek Vol':<12} | {'Diff %':<12} | {'Windows':<10}"
    column_width = 145
    header = (
        f"{'Filename':<12} | {'Calc Vol':<12} | {'Rebek Vol':<12} | "
        f"{'Diff %':<12} | {'Windows':<10} | {'ESP':<14} | {'MHP':<14} | {'HI':<14}"
    )
    print("\n\n" + "=" * column_width)
    print(header)
    print("-" * column_width)

    for row in summary_table_data:
        # 格式化输出每一行
        # 注意：row['rebek'] 如果是字符串 "N/A" 需要特殊处理格式
        rebek_val = f"{row['rebek']:<12}" if isinstance(row['rebek'], str) else f"{float(row['rebek']):<12.2f}"

        esp_val = "ESP_FAILED" if row.get('esp') is None else f"{row['esp']:<14.6f}"
        if row.get('mhp') is None:
            mhp_val = "MHP_FAILED" if row.get('field_status') == "esp_only" else "N/A"
        else:
            mhp_val = f"{row['mhp']:<14.6f}"
        if row.get('hi') is None:
            hi_val = "HI_FAILED" if row.get('field_status') == "esp_only" else "N/A"
        else:
            hi_val = f"{row['hi']:<14.6f}"

        print(
            f"{row['name']:<12} | {row['calc']:<12.2f} | {rebek_val} | "
            f"{row['diff']:<12} | {row['win_count']:<10} | {esp_val:<14} | {mhp_val:<14} | {hi_val:<14}"
        )

    print("-" * column_width)
    # 2. 打印平均误差 (MAPE)
    if len(error_values) > 0:
        mape = np.mean(error_values)
        print(f"Final Summary: Processed {len(summary_table_data)} cages.")
        print(f"Mean Absolute Percentage Error (MAPE): {mape:.2f}%")
    else:
        print("No guest data available for error analysis.")
    print("=" * column_width + "\n")

    # 3. 绘图 (可选)
    plt.figure(figsize=(10, 6))
    plt.plot(volume_data, marker='o', linestyle='-', color='b', label='Cavity Volume')
    plt.title(plot_title)
    plt.xlabel(x_label)
    plt.ylabel('Volume (Å³)')
    plt.grid(True)
    plt.legend()
    save_name = "volume_independent.png"
    plt.savefig(os.path.join(Output_Path, save_name))
