"""Backend bridge between PyMOL and the DSP-Cage calculation engine."""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path
from typing import Any

from pymol import cmd


PLUGIN_DIR = Path(__file__).resolve().parent

# 兼容你当前算法文件中的：
# from Balloon import ...
# from data import ...
# 这类非相对导入写法
if str(PLUGIN_DIR) not in sys.path:
    sys.path.insert(0, str(PLUGIN_DIR))

from Cavity_Calculation import cavity  # noqa: E402


def make_safe_name(name: str) -> str:
    """Convert a PyMOL selection name into a safe filename."""
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    safe_name = safe_name.strip("._")

    return safe_name or "pymol_selection"


def run_calculation(
    selection: str,
    output_dir: str = "",
    center_type: int = 2,
    subdivision: int = 4,
) -> dict[str, Any]:
    """Calculate the cavity for a PyMOL object or selection."""

    selection = selection.strip()

    if not selection:
        raise ValueError("PyMOL 对象或选择不能为空。")

    atom_count = cmd.count_atoms(selection)

    if atom_count <= 0:
        raise ValueError(
            f"选择 '{selection}' 中没有检测到原子，请检查对象名称。"
        )

    if output_dir:
        work_dir = Path(output_dir).expanduser().resolve()
    else:
        work_dir = Path.home() / "DSP_Cage_results"

    input_dir = work_dir / "input"
    result_dir = work_dir / "output"
    probe_dir = work_dir / "PROBE"

    input_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)
    probe_dir.mkdir(parents=True, exist_ok=True)

    safe_name = make_safe_name(selection)
    input_file = input_dir / f"{safe_name}.pdb"

    # 将 PyMOL 当前对象或 selection 导出为 PDB
    cmd.save(str(input_file), selection)

    if not input_file.exists():
        raise RuntimeError(f"PyMOL 未能导出结构文件：{input_file}")

    engine = cavity()

    # 防止前一次动态轨迹状态影响本次计算
    if hasattr(engine, "last_frame_data"):
        engine.last_frame_data = None

    input_path_string = str(input_dir) + os.sep
    output_path_string = str(result_dir) + os.sep

    volume_list = engine.Calculate_Cavity(
        input_file.name,
        int(center_type),
        int(subdivision),
        input_path_string,
        output_path_string,
    )

    if not volume_list:
        raise RuntimeError(
            "DSP-Cage 没有返回体积结果，请检查输入结构和计算参数。"
        )

    volume = float(volume_list[0])

    cavity_file = result_dir / f"{input_file.stem}_cavity.pdb"

    result = {
        "selection": selection,
        "atom_count": atom_count,
        "volume": volume,
        "input_file": str(input_file),
        "output_dir": str(result_dir),
        "cavity_file": str(cavity_file) if cavity_file.exists() else None,
        "window_info": getattr(engine, "last_window_info", None),
        "rebek_volume": getattr(engine, "cached_rebek_vol", None),
    }

    if cavity_file.exists():
        load_cavity_into_pymol(
            cavity_file=cavity_file,
            source_name=safe_name,
        )

    return result


def load_cavity_into_pymol(
    cavity_file: Path,
    source_name: str,
) -> None:
    """Load the generated cavity probe file into PyMOL."""

    object_name = f"DSP_Cage_{source_name}_cavity"

    # 避免重复运行后产生重名对象
    cmd.delete(object_name)
    cmd.load(str(cavity_file), object_name)

    cmd.hide("everything", object_name)
    cmd.show("spheres", object_name)

    cmd.set("sphere_scale", 0.16, object_name)
    cmd.set("sphere_transparency", 0.35, object_name)
    cmd.color("yellow", object_name)

    cmd.group("DSP_Cage_Results", object_name)
    cmd.zoom(f"({source_name}) or ({object_name})")