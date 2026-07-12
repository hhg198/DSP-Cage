import os
import sys
import numpy as np


def write_mapped_pdb(input_probe_pdb, output_pdb, values, field_name="FIELD"):
    """
    读取原始探针 PDB 文件，将计算出的数值写入到 B-factor 字段，生成用于可视化的 PDB。
    """
    if not os.path.exists(input_probe_pdb):
        raise FileNotFoundError(f"找不到原始探针 PDB 文件: {input_probe_pdb}")

    lines_written = 0
    val_idx = 0

    with open(input_probe_pdb, "r") as infile, open(output_pdb, "w") as outfile:
        for line in infile:
            if line.startswith(("ATOM", "HETATM")):
                if val_idx >= len(values):
                    # 如果数值不够，后面补0
                    val = 0.0
                else:
                    val = values[val_idx]

                # PDB 格式严格切片与拼接：
                # B-factor 位于第 61-66 列（Python 索引 60:66），要求宽度为 6，保留 2 位小数
                # Occupancy 位于第 55-60 列（Python 索引 54:60），这里保持原样或格式化
                before_b = line[:60]
                after_b = line[66:]

                # 格式化数值，防止溢出 6 位宽度（例如：-99.99 到 999.99）
                # 如果数值极大或极小，进行安全截断
                if val > 999.99: val = 999.99
                if val < -99.99: val = -99.99

                b_factor_str = f"{val:6.2f}"
                new_line = f"{before_b}{b_factor_str}{after_b}"

                outfile.write(new_line)
                val_idx += 1
                lines_written += 1
            else:
                outfile.write(line)

    print(f"[成功] 已将 {field_name} 写入到 B-factor。")
    print(f"       输入探针: {input_probe_pdb}")
    print(f"       输出可视化文件: {output_pdb}")
    print(f"       成功映射了 {lines_written} 个顶点的数值。")


# 示例：如何嵌入到你现有的 main() 流程中
def add_visualization_to_pipeline(molecule_path, probe_path, output_dir, molecule_filename, esp_values,
                                  mhp_values=None):
    stem, _ = os.path.splitext(molecule_filename)

    # 1. 映射 ESP 到 PDB
    esp_pdb_path = os.path.join(output_dir, f"{stem}_ESP_vis.pdb")
    write_mapped_pdb(probe_path, esp_pdb_path, esp_values, field_name="ESP")

    # 2. 如果 MHP 存在，映射 MHP 到 PDB
    if mhp_values is not None:
        mhp_pdb_path = os.path.join(output_dir, f"{stem}_MHP_vis.pdb")
        write_mapped_pdb(probe_path, mhp_pdb_path, mhp_values, field_name="MHP")