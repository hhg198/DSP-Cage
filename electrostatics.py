"""
DSP-Cage

This file contains modified components adapted from:

CageCavityCalc
Vicente Marti-Centelles

Copyright (c) 2021 VicenteMartiCentelles

Licensed under the MIT License.

See LICENSE_C3 for details.

The present implementation additionally includes
new Dynamic Splitting of Probes algorithms.
"""
import csv
import os
import shutil
from tempfile import mkdtemp

import numpy as np


ESP_FACTOR = (8.987551792e9) * (1.602176634e-19) * (1e10)
MIN_DISTANCE = 0.1


def print_to_pdb_file(filename, positions, atom_names):
    """Write a minimal PDB file for OpenBabel charge assignment."""
    with open(filename, "w") as file:
        for idx, (pos, atom_name) in enumerate(zip(positions, atom_names), start=1):
            element = str(atom_name).strip().upper()
            if len(element) > 1:
                element = element[0] + element[1:].lower()
            file.write(
                f"HETATM{idx:5d} {str(atom_name).upper():>4s} MOL     1"
                f"    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                f"  1.00  0.00          {element:>2s}\n"
            )
        file.write("END\n")


def read_pdb_points(pdb_path):
    """Read coordinates from ATOM/HETATM records."""
    points = []
    with open(pdb_path, "r") as file:
        for line in file:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    points.append(
                        [
                            float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54]),
                        ]
                    )
                except ValueError:
                    continue
    return np.asarray(points, dtype=float)


def calculate_partial_charges_using_ob(positions, atom_names, method='eem'):
    '''
    From Openbabel description:

    eem    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Bultinck B3LYP/6-31G*/MPA
    eem2015ba    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf B3LYP/6-311G/AIM
    eem2015bm    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf B3LYP/6-311G/MPA
    eem2015bn    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf B3LYP/6-311G/NPA
    eem2015ha    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf HF/6-311G/AIM
    eem2015hm    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf HF/6-311G/MPA
    eem2015hn    Assign Electronegativity Equilization Method (EEM) atomic partial charges. Cheminf HF/6-311G/NPA
    eqeq    Assign EQEq (charge equilibration) partial charges.
    fromfile    Assign charges from file containing {'atom-name', charge} pairs
    gasteiger    Assign Gasteiger-Marsili sigma partial charges
    mmff94       Assign MMFF94 partial charges
    none    Clear all partial charges
    qeq    Assign QEq (charge equilibration) partial charges (Rappe and Goddard, 1991)
    qtpie    Assign QTPIE (charge transfer, polarization and equilibration) partial charges (Chen and Martinez, 2007)
    '''
    # openbabel
    try:
        from openbabel import openbabel
    except ImportError as exc:
        raise ImportError("OpenBabel is required to calculate partial charges.") from exc

    tmpdir_path = mkdtemp()
    try:
        print_to_pdb_file(os.path.join(tmpdir_path, "temp.pdb"), positions, atom_names)

        ob_conversion = openbabel.OBConversion()
        ob_conversion.SetInAndOutFormats("pdb", "mol2")

        # We firstly change pdb -> mol2, because charges are not assigned correctly when used with pdb
        mol = openbabel.OBMol()
        ob_conversion.ReadFile(mol, os.path.join(tmpdir_path, "temp.pdb"))
        ob_conversion.WriteFile(mol, os.path.join(tmpdir_path, "temp.mol2"))
        ob_conversion.SetInAndOutFormats("mol2", "mol2")

        mol = openbabel.OBMol()
        ob_conversion.ReadFile(mol, os.path.join(tmpdir_path, "temp.mol2"))

        ob_charge_model = openbabel.OBChargeModel.FindType(method)
        if ob_charge_model is None:
            raise ValueError(f"OpenBabel charge model not found: {method}")

        is_calculated = ob_charge_model.ComputeCharges(mol)

        if is_calculated:
            return list(ob_charge_model.GetPartialCharges())
        return None
    finally:
        shutil.rmtree(tmpdir_path, ignore_errors=True)


def calculate_partial_charges(positions, atom_names, method, metal_name=None, metal_charge=None):
    positions = np.asarray(positions, dtype=float)
    atom_names = np.asarray(atom_names)

    if metal_name is not None:
        list_of_metals = [a for a, atom_name in enumerate(atom_names) if atom_name.title() == metal_name.title()]
        positions = np.delete(positions, list_of_metals, axis=0)
        atom_names = np.delete(atom_names, list_of_metals)
        partial_charges = list(calculate_partial_charges_using_ob(positions, atom_names, method))
        for index in list_of_metals:
            partial_charges.insert(index, metal_charge)
    else:
        partial_charges = list(calculate_partial_charges_using_ob(positions, atom_names, method))
    #print("Sum", np.sum(partial_charges))
    return partial_charges


def calculate_esp_on_points(
    atom_positions,
    atom_names,
    grid_points,
    method="gasteiger",
    metal_name=None,
    metal_charge=None,
    partial_charges=None,
    chunk_size=5000,
    min_distance=MIN_DISTANCE,
):
    """
    Calculate electrostatic potential on CMCC probe/cavity points.

    Parameters
    ----------
    atom_positions, atom_names
        Cage atom coordinates and element names.
    grid_points
        Probe vertices, preferably read from the CMCC *_PROBE.pdb file.
    partial_charges
        Optional precomputed charge list. If omitted, OpenBabel is used.
    """
    atom_positions = np.asarray(atom_positions, dtype=float)
    grid_points = np.asarray(grid_points, dtype=float)

    if len(atom_positions) == 0:
        raise ValueError("No cage atoms were provided for ESP calculation.")
    if len(grid_points) == 0:
        return np.asarray([], dtype=float), []

    if partial_charges is None:
        partial_charges = calculate_partial_charges(
            atom_positions,
            atom_names,
            method=method,
            metal_name=metal_name,
            metal_charge=metal_charge,
        )
    if partial_charges is None:
        raise RuntimeError(f"Failed to calculate partial charges with method: {method}")

    partial_charges = np.asarray(partial_charges, dtype=float)
    if len(partial_charges) != len(atom_positions):
        raise ValueError(
            f"Charge count ({len(partial_charges)}) does not match atom count ({len(atom_positions)})."
        )

    esp_values = np.empty(len(grid_points), dtype=float)
    for start in range(0, len(grid_points), chunk_size):
        stop = min(start + chunk_size, len(grid_points))
        diff = grid_points[start:stop, None, :] - atom_positions[None, :, :]
        distances = np.linalg.norm(diff, axis=2)
        distances = np.where(distances < min_distance, min_distance, distances)
        esp_values[start:stop] = (ESP_FACTOR / distances).dot(partial_charges)

    return esp_values, partial_charges.tolist()


def calculate_esp_for_files(
    molecule_path,
    probe_path,
    method="gasteiger",
    metal_name=None,
    metal_charge=None,
    partial_charges=None,
    chunk_size=5000,
):
    """Convenience wrapper for molecule PDB/MOL2 + CMCC *_PROBE.pdb files."""
    try:
        from Input_methond import read_positions_and_atom_names_from_file
    except ImportError:
        from .Input_methond import read_positions_and_atom_names_from_file

    atom_positions, atom_names, _, _ = read_positions_and_atom_names_from_file(molecule_path)
    grid_points = read_pdb_points(probe_path)
    esp_values, charges = calculate_esp_on_points(
        atom_positions,
        atom_names,
        grid_points,
        method=method,
        metal_name=metal_name,
        metal_charge=metal_charge,
        partial_charges=partial_charges,
        chunk_size=chunk_size,
    )
    return {
        "points": grid_points,
        "esp": esp_values,
        "partial_charges": charges,
    }


def write_field_csv(csv_path, points, esp_values=None):
    """Write ESP values on probe points for downstream plotting or analysis."""
    with open(csv_path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["x", "y", "z", "ESP"])
        esp_values = np.zeros(len(points)) if esp_values is None else esp_values
        for point, esp in zip(points, esp_values):
            writer.writerow([point[0], point[1], point[2], esp])
