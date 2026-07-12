import csv
import logging
import os
import shutil
from tempfile import mkdtemp

import numpy as np

try:
    from rdkit import Chem
except ImportError:
    Chem = None

try:
    from data import hydrophValuesGhose1998, hydrophValuesCrippen1999
except ImportError:
    from .data import hydrophValuesGhose1998, hydrophValuesCrippen1999


logger = logging.getLogger(__name__)


HYDROPHOBICITY_METHODS = {
    "Ghose": hydrophValuesGhose1998,
    "Crippen": hydrophValuesCrippen1999,
}



def _require_rdkit():
    if Chem is None:
        raise ImportError("RDKit is required to calculate molecular hydrophobic potential.")


def convert_pdb_to_mol2_with_openbabel(pdb_path):
    """Convert PDB to a temporary MOL2 file for RDKit SMARTS atom typing."""
    try:
        from openbabel import openbabel
    except ImportError as exc:
        raise ImportError(
            "OpenBabel is required to convert PDB to MOL2 for MHP calculation."
        ) from exc

    tmpdir_path = mkdtemp()
    mol2_path = os.path.join(tmpdir_path, "mhp_topology.mol2")

    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("pdb", "mol2")

    mol = openbabel.OBMol()
    if not ob_conversion.ReadFile(mol, pdb_path):
        shutil.rmtree(tmpdir_path, ignore_errors=True)
        raise ValueError(f"OpenBabel failed to read PDB file: {pdb_path}")
    if not ob_conversion.WriteFile(mol, mol2_path):
        shutil.rmtree(tmpdir_path, ignore_errors=True)
        raise ValueError(f"OpenBabel failed to write temporary MOL2 file: {mol2_path}")

    return mol2_path, tmpdir_path


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


def read_rdkit_mol(molecule_path, sanitize=True):
    """Read a cage molecule with RDKit, converting PDB to MOL2 first when needed."""
    _require_rdkit()
    suffix = os.path.splitext(molecule_path)[1].lower()
    cleanup_dir = None
    rdkit_path = molecule_path

    if suffix == ".pdb":
        rdkit_path, cleanup_dir = convert_pdb_to_mol2_with_openbabel(molecule_path)

    try:
        if os.path.splitext(rdkit_path)[1].lower() == ".mol2":
            mol = Chem.MolFromMol2File(rdkit_path, sanitize=sanitize, removeHs=False)
        else:
            mol = Chem.MolFromPDBFile(rdkit_path, removeHs=False, sanitize=False)
            if mol is not None and sanitize:
                Chem.SanitizeMol(mol)

        if mol is None:
            raise ValueError(f"RDKit failed to load molecule: {rdkit_path}")
        return mol
    finally:
        if cleanup_dir is not None:
            shutil.rmtree(cleanup_dir, ignore_errors=True)

# Assign the hydrophobic values to cage atoms from the library
def assignHydrophobicValuesToCageAtoms(rdkit_cage, hydrophValues, printLevel=2):
    _require_rdkit()
    listGlobalOut = list()
    listGlobalOut2 = list()
    for i in hydrophValues:
        if isinstance(hydrophValues[i][2], str) or hydrophValues[i][2] == "undefined":
            listGlobalOut2.append([])
            continue
        match = rdkit_cage.GetSubstructMatches(Chem.MolFromSmarts(hydrophValues[i][1]), False, False)
        listOut = list()
        for x in match:
            listOut.append(x[0] + 1)
        listGlobalOut2.append(list(set(listOut)))
        if listOut:
            listGlobalOut.extend(list(set(listOut)))

    numberOfAtoms = rdkit_cage.GetNumAtoms()

    atomTypesList = []
    for i in range(0, numberOfAtoms):
        atomTypesList.append([])

    atomTypesListHydrophValues = []
    for i in range(0, numberOfAtoms):
        atomTypesListHydrophValues.append([])

    atomTypesValuesListHydrophValues = []
    atomTypesMeanListHydrophValues = []
    atomTypesInfoAtomSymbol = []
    atomTypesInfoAtomGlobalIndex = []
    atomTypesAssignemet = []

    for i in range(0, len(listGlobalOut2)):
        if listGlobalOut2[i]:
            for j in listGlobalOut2[i]:
                atomTypesList[j - 1].append(i + 1)

    for i in range(0, numberOfAtoms):
        atomSymbol = rdkit_cage.GetAtomWithIdx(i).GetSymbol()
        valuesList = []
        if len(atomTypesList[i]) > 0:
            valuesList = []
            for j in atomTypesList[i]:
                valuesList.append(hydrophValues[j][2])
        else:
            valuesList.append(0)
            logger.warning(f"WARNING: atom {atomSymbol:}, {i + 1:}, not found. Assinged 0 as hydropobicity factor.")

        atomTypesListHydrophValues.append(valuesList)
        meanValuestList = np.mean(valuesList)
        atomTypesMeanListHydrophValues.append(meanValuestList)
        atomTypesValuesListHydrophValues.append(valuesList)
        atomTypesInfoAtomSymbol.append(atomSymbol.upper()) # we need upper case
        atomTypesInfoAtomGlobalIndex.append(i + 1)
        atomTypesAssignemet.append(atomTypesList[i])
        logger.info(f"{atomSymbol:}, {i + 1:}, {atomTypesList[i]:}, {valuesList:}, {meanValuestList:}")

    return atomTypesMeanListHydrophValues, atomTypesValuesListHydrophValues, atomTypesInfoAtomSymbol, atomTypesInfoAtomGlobalIndex, atomTypesAssignemet


def assign_hydrophobic_values_from_mol(rdkit_cage, method="Ghose"):
    """Return one mean hydrophobic contribution value for each cage atom."""
    if method not in HYDROPHOBICITY_METHODS:
        raise ValueError(f"Unknown hydrophobicity method: {method}")
    values, _, _, _, _ = assignHydrophobicValuesToCageAtoms(
        rdkit_cage,
        HYDROPHOBICITY_METHODS[method],
        printLevel=0,
    )
    return np.asarray(values, dtype=float)


def assign_hydrophobic_values_from_pdb(pdb_path, method="Ghose", expected_count=None):
    """Read a PDB file and assign C3/Ghose-style atom hydrophobic values."""
    mol = read_rdkit_mol(pdb_path)
    if expected_count is not None and mol.GetNumAtoms() != expected_count:
        logger.warning(
            "RDKit atom count (%s) differs from coordinate atom count (%s) for %s.",
            mol.GetNumAtoms(),
            expected_count,
            pdb_path,
        )
    values = assign_hydrophobic_values_from_mol(mol, method=method)
    if expected_count is not None and len(values) != expected_count:
        fixed = np.zeros(expected_count, dtype=float)
        fixed[: min(expected_count, len(values))] = values[: min(expected_count, len(values))]
        return fixed
    return values


def calc_single_hydrophobicity(distance,Hydroph_Value, distance_function):
    if distance_function == "Audry":
        return Hydroph_Value / (1 + distance)
    elif distance_function == "Fauchere":
        return Hydroph_Value * np.exp(-1 * distance)
    elif distance_function == "Fauchere2":
        return Hydroph_Value * np.exp(-1 / 2 * distance)
    elif distance_function == "OnlyValues":
        return Hydroph_Value
    return 0


def calculate_mhp_on_points(
    atom_positions,
    hydrophobic_values,
    grid_points,
    distance_function="Fauchere",
    cutoff=20.0,
    chunk_size=5000,
):
    """
    Calculate molecular hydrophobic potential on CMCC probe/cavity points.

    The distance attenuation follows the original C3 hydrophobicity functions.
    """
    atom_positions = np.asarray(atom_positions, dtype=float)
    hydrophobic_values = np.asarray(hydrophobic_values, dtype=float)
    grid_points = np.asarray(grid_points, dtype=float)

    if len(atom_positions) == 0:
        raise ValueError("No cage atoms were provided for MHP calculation.")
    if len(atom_positions) != len(hydrophobic_values):
        raise ValueError(
            f"Hydrophobic value count ({len(hydrophobic_values)}) does not match atom count ({len(atom_positions)})."
        )
    if len(grid_points) == 0:
        return np.asarray([], dtype=float)

    mhp_values = np.empty(len(grid_points), dtype=float)
    for start in range(0, len(grid_points), chunk_size):
        stop = min(start + chunk_size, len(grid_points))
        diff = grid_points[start:stop, None, :] - atom_positions[None, :, :]
        distances = np.linalg.norm(diff, axis=2)
        within_cutoff = distances < cutoff

        if distance_function == "Audry":
            contributions = hydrophobic_values[None, :] / (1.0 + distances)
        elif distance_function == "Fauchere":
            contributions = hydrophobic_values[None, :] * np.exp(-distances)
        elif distance_function == "Fauchere2":
            contributions = hydrophobic_values[None, :] * np.exp(-0.5 * distances)
        elif distance_function == "OnlyValues":
            contributions = np.broadcast_to(hydrophobic_values, distances.shape)
        else:
            raise ValueError(f"Unknown hydrophobicity distance function: {distance_function}")

        mhp_values[start:stop] = np.where(within_cutoff, contributions, 0.0).sum(axis=1)

    return mhp_values


def calculate_mhp_for_files(
    molecule_path,
    probe_path,
    method="Ghose",
    distance_function="Fauchere",
    cutoff=20.0,
    chunk_size=5000,
):
    """Convenience wrapper for molecule PDB/MOL2 + CMCC *_PROBE.pdb files."""
    try:
        from Input_methond import read_positions_and_atom_names_from_file
    except ImportError:
        from .Input_methond import read_positions_and_atom_names_from_file

    atom_positions, atom_names, _, _ = read_positions_and_atom_names_from_file(molecule_path)
    hydrophobic_values = assign_hydrophobic_values_from_pdb(
        molecule_path,
        method=method,
        expected_count=len(atom_positions),
    )
    grid_points = read_pdb_points(probe_path)
    mhp_values = calculate_mhp_on_points(
        atom_positions,
        hydrophobic_values,
        grid_points,
        distance_function=distance_function,
        cutoff=cutoff,
        chunk_size=chunk_size,
    )
    return {
        "points": grid_points,
        "mhp": mhp_values,
        "atom_hydrophobic_values": hydrophobic_values,
        "atom_names": atom_names,
        "hydrophobic_assignment": method,
    }


def write_field_csv(csv_path, points, mhp_values=None):
    """Write MHP values on probe points for downstream plotting or analysis."""
    with open(csv_path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["x", "y", "z", "MHP"])
        mhp_values = np.zeros(len(points)) if mhp_values is None else mhp_values
        for point, mhp in zip(points, mhp_values):
            writer.writerow([point[0], point[1], point[2], mhp])
