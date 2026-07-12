import pathlib
import logging
import pywindow as pw
from typing import List  # 引入类型提示
import sys

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    stream=sys.stdout  # <--- 强制让日志也走标准输出通道
)
logger = logging.getLogger(__name__)

def split_trajectory_to_pdbs(
        input_filepath: str,
        output_directory: str,
        swap_atoms: dict = None,
        forcefield: str = "opls",
        filename_prefix: str = "frame"
) -> List[str]:  # <--- 修改点1：返回值类型改为字符串列表
    """
    加载轨迹文件，将每一帧单独保存为 PDB 文件，并返回包含所有生成的 PDB 文件路径的列表。
    """

    # 1. 路径处理
    input_path = pathlib.Path(input_filepath).resolve()
    out_dir = pathlib.Path(output_directory).resolve()

    if not out_dir.exists():
        out_dir.mkdir(parents=True)
        logger.info(f"Created output directory: {out_dir}")

    # 2. 加载轨迹
    logger.info(f"Loading trajectory from: {input_path}")
    try:
        traj = pw.DLPOLY(input_path)
    except Exception as e:
        logger.error(f"Failed to load trajectory: {e}")
        return []  # 如果失败，返回空列表

    total_frames = traj.no_of_frames
    logger.info(f"Total frames found: {total_frames}")

    # 3. 逐帧处理并重命名
    if swap_atoms is None:
        swap_atoms = {"he": "H"}

    temp_base_name = "temp_worker.pdb"
    temp_filepath = out_dir / temp_base_name

    # <--- 修改点2：初始化一个空列表用来存路径
    generated_files = []

    for i in range(total_frames):
        try:
            # A. 保存临时文件
            traj.save_frames(
                frames=[i],
                filepath=temp_filepath,
                forcefield=forcefield,
                swap_atoms=swap_atoms
            )

            # B. 修正文件名
            ugly_filename = f"{temp_base_name}_{i}.pdb"
            ugly_path = out_dir / ugly_filename

            pretty_filename = f"{filename_prefix}_{i}.pdb"
            pretty_path = out_dir / pretty_filename

            # C. 重命名并记录路径
            if ugly_path.exists():
                if pretty_path.exists():
                    pretty_path.unlink()

                ugly_path.rename(pretty_path)

                # <--- 修改点3：将生成的文件的绝对路径添加到列表中
                # 使用 str(pretty_path) 转为字符串，方便后续通用处理
                generated_files.append(str(pretty_path))

                if (i + 1) % 10 == 0 or (i + 1) == total_frames:
                    logger.info(f"Saved: {pretty_filename}")
            else:
                logger.warning(f"Expected file not found: {ugly_filename}")

        except Exception as e:
            logger.error(f"Error processing frame {i}: {e}")

    logger.info("All frames have been successfully exported.")

    # <--- 修改点4：返回列表
    return generated_files