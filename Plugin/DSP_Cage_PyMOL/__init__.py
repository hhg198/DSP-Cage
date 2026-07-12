"""DSP-Cage PyMOL plugin entry."""

_dialog = None


def __init_plugin__(app=None):
    """PyMOL loads this function when the plugin starts."""
    from pymol import cmd
    from pymol.plugins import addmenuitemqt

    # 添加 PyMOL 插件菜单
    addmenuitemqt("DSP-Cage", show_dialog)

    # 添加 PyMOL 命令行命令
    cmd.extend("dsp_cage", dsp_cage_command)


def show_dialog():
    """Open or reactivate the DSP-Cage dialog."""
    global _dialog

    from .dialog import DSPCageDialog

    if _dialog is None:
        _dialog = DSPCageDialog()

    _dialog.refresh_objects()
    _dialog.show()
    _dialog.raise_()
    _dialog.activateWindow()


def dsp_cage_command(
    selection="all",
    output_dir="",
    center_type=2,
    subdivision=4,
):
    """Run DSP-Cage from the PyMOL command line.

    Example:
        dsp_cage cage1, center_type=2, subdivision=4
    """
    from .backend import run_calculation

    result = run_calculation(
        selection=selection,
        output_dir=output_dir,
        center_type=int(center_type),
        subdivision=int(subdivision),
    )

    print("=" * 50)
    print("DSP-Cage calculation completed")
    print(f"Cavity volume: {result['volume']:.3f} Å³")
    print(f"Cavity file: {result.get('cavity_file', 'Not generated')}")

    window_info = result.get("window_info") or {}
    print(f"Window count: {window_info.get('window_count', 0)}")
    print("=" * 50)

    return result
 