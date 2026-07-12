import argparse
import csv
import os
import sys
from generate_pymol_maps import write_mapped_pdb


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
if BASE_DIR not in sys.path:
    sys.path.insert(0, BASE_DIR)
os.environ.setdefault("BABEL_DATADIR", BASE_DIR)

import numpy as np

from electrostatics import calculate_esp_for_files
from hydrophobicity import calculate_mhp_for_files


COMMON_METAL_CHARGES = {
    "Pd": 2,
    "Pt": 2,
    "Zn": 2,
    "Cu": 2,
    "Co": 2,
    "Fe": 2,
    "Ru": 2,
    "Ag": 1,
    "Ni": 2,
    "Ti": 4,
    "Ga": 3,
}


def detect_main_metal_charge(pdb_path):
    transition_metals = {metal.upper() for metal in COMMON_METAL_CHARGES}
    counts = {}

    with open(pdb_path, "r") as file:
        for line in file:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            element = line[76:78].strip()
            if not element:
                atom_name = line[12:16].strip()
                element = "".join(char for char in atom_name if char.isalpha())

            element = element.capitalize()
            if element.upper() in transition_metals:
                counts[element] = counts.get(element, 0) + 1

    if not counts:
        return None, None

    metal = max(counts, key=counts.get)
    return metal, COMMON_METAL_CHARGES.get(metal)


def build_default_probe_path(input_dir, filename):
    stem, _ = os.path.splitext(filename)
    # 获取 input_dir 的上一级目录，并拼接出 datas/output 路径
    base_datas_dir = os.path.dirname(input_dir)
    output_dir = os.path.join(base_datas_dir, "output")

    return os.path.join(output_dir, f"{stem}_cavity.pdb")



def build_default_output_path(output_dir, filename):
    stem, _ = os.path.splitext(filename)
    return os.path.join(output_dir, f"{stem}_ESP_MHP.csv")


def build_default_esp_output_path(output_dir, filename):
    stem, _ = os.path.splitext(filename)
    return os.path.join(output_dir, f"{stem}_ESP.csv")


def write_combined_csv(output_path, points, esp_values, mhp_values):
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(output_path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["x", "y", "z", "ESP", "MHP"])
        for point, esp, mhp in zip(points, esp_values, mhp_values):
            writer.writerow([point[0], point[1], point[2], esp, mhp])


def write_esp_csv(output_path, points, esp_values):
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(output_path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["x", "y", "z", "ESP"])
        for point, esp in zip(points, esp_values):
            writer.writerow([point[0], point[1], point[2], esp])


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate ESP and MHP values on CMCC *_PROBE.pdb vertices."
    )
    parser.add_argument(
        "filename",
        nargs="?",
        default="B7.pdb",
        help="Input cage file name under --input-dir, default: B7.pdb",
    )
    parser.add_argument("--input-dir", default="datas/input", help="Input directory.")
    parser.add_argument("--output-dir", default="datas/output", help="Output directory.")
    parser.add_argument("--probe", default=None, help="Explicit *_PROBE.pdb path.")
    parser.add_argument("--output", default=None, help="Explicit output CSV path.")
    parser.add_argument("--charge-method", default="eem", help="OpenBabel charge model.")
    parser.add_argument("--metal", default=None, help="Metal element name, e.g. Ga.")
    parser.add_argument("--metal-charge", type=int, default=None, help="Metal formal charge.")
    parser.add_argument("--hydroph-method", default="Ghose", choices=["Ghose", "Crippen"])
    parser.add_argument(
        "--distance-function",
        default="Fauchere",
        choices=["Audry", "Fauchere", "Fauchere2", "OnlyValues"],
    )
    return parser.parse_args()


def main():
    os.chdir(BASE_DIR)
    args = parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    if not os.path.isabs(input_dir):
        input_dir = os.path.join(BASE_DIR, input_dir)
    if not os.path.isabs(output_dir):
        output_dir = os.path.join(BASE_DIR, output_dir)

    if os.path.isabs(args.filename):
        molecule_path = args.filename
    elif os.path.dirname(args.filename):
        molecule_path = os.path.join(BASE_DIR, args.filename)
    else:
        molecule_path = os.path.join(input_dir, args.filename)

    molecule_filename = os.path.basename(molecule_path)
    molecule_dir = os.path.dirname(molecule_path)
    probe_path = args.probe or build_default_probe_path(molecule_dir, molecule_filename)
    output_path = args.output or build_default_output_path(output_dir, molecule_filename)
    if args.probe and not os.path.isabs(probe_path):
        probe_path = os.path.join(BASE_DIR, probe_path)
    if args.output and not os.path.isabs(output_path):
        output_path = os.path.join(BASE_DIR, output_path)

    if not os.path.exists(molecule_path):
        raise FileNotFoundError(f"Input molecule file not found: {molecule_path}")
    if not os.path.exists(probe_path):
        raise FileNotFoundError(
            f"Probe file not found: {probe_path}. Run cavity calculation first."
        )

    metal = args.metal
    metal_charge = args.metal_charge
    if metal is None and metal_charge is None:
        metal, metal_charge = detect_main_metal_charge(molecule_path)

    print(f"Molecule: {molecule_path}")
    print(f"Probe:    {probe_path}")
    if metal and metal_charge is not None:
        print(f"Metal charge override: {metal} = {metal_charge}")
    else:
        print("Metal charge override: none")

    esp_result = calculate_esp_for_files(
        molecule_path,
        probe_path,
        method=args.charge_method,
        metal_name=metal,
        metal_charge=metal_charge,
    )
    esp_values = np.asarray(esp_result["esp"], dtype=float)

    try:
        mhp_result = calculate_mhp_for_files(
            molecule_path,
            probe_path,
            method=args.hydroph_method,
            distance_function=args.distance_function,
        )
    except Exception as exc:
        esp_output_path = output_path
        if args.output is None:
            esp_output_path = build_default_esp_output_path(output_dir, molecule_filename)
        write_esp_csv(esp_output_path, esp_result["points"], esp_result["esp"])
        print("WARNING: MHP calculation failed; saved ESP values only.")
        print(f"MHP error: {type(exc).__name__}: {exc}")
        print(f"Saved: {esp_output_path}")
        print(f"Point count: {len(esp_values)}")
        print(f"Mean ESP: {np.mean(esp_values):.6f}")
        print(f"Min/Max ESP: {np.min(esp_values):.6f} / {np.max(esp_values):.6f}")
        return

    write_combined_csv(
        output_path,
        esp_result["points"],
        esp_result["esp"],
        mhp_result["mhp"],
    )

    mhp_values = np.asarray(mhp_result["mhp"], dtype=float)
    mhp_pos = mhp_values[mhp_values > 0]
    mhp_neg = mhp_values[mhp_values < 0]
    hi_denominator = np.sum(mhp_pos) - np.sum(mhp_neg)
    hydrophobic_index = np.sum(mhp_pos) / hi_denominator if hi_denominator != 0 else np.nan

    print(f"Saved: {output_path}")
    print(f"Point count: {len(esp_values)}")
    print(f"MHP assignment: {mhp_result.get('hydrophobic_assignment', args.hydroph_method)}")
    print(f"Mean ESP: {np.mean(esp_values):.6f}")
    print(f"Min/Max ESP: {np.min(esp_values):.6f} / {np.max(esp_values):.6f}")
    print(f"Mean MHP: {np.mean(mhp_values):.6f}")
    print(f"Min/Max MHP: {np.min(mhp_values):.6f} / {np.max(mhp_values):.6f}")
    print(f"Hydrophobic index (HI): {hydrophobic_index:.6f}")


if __name__ == "__main__":
    main()
