#!/usr/bin/env python
"""Protein-ligand PDB -> qcelemental dimer -> AP-Net2 interaction energy.

Two calls do all the work:

    dimer = pdb_to_dimer("pdbs/cdk2_17.pdb")
    energies = predict_ap2([dimer])

`pdb_to_dimer` assumes the PDB is already prepared: every hydrogen present
(protein and ligand), one protonation state per residue, and the ligand in the
same file as the protein.  It does no protonation and no geometry repair.  It
cuts a pocket of whole residues around the ligand, caps every severed bond with
a hydrogen, reads each residue's formal charge from the hydrogens actually
present, and returns a two-fragment Molecule (fragment 0 = pocket,
fragment 1 = ligand), the dimer layout AP2 expects.

The pocket rules follow ~/data/TAQO22/apnet/build_pockets.py; the hydrogen-cap
bond lengths follow ~/gits/protein-prep (src/cdsg_prep/geometry.py).  The
bundled pdbs/cdk2_*.pdb files are protein-prep's 16 minimized CDK2 complexes;
cdk2_ligands.csv supplies their ligand charges and experimental binding data.

Run with the qcml environment:

    conda activate qcml
    export QCMLFORGE_AUTO_DOWNLOAD_PRETRAINED=1   # first run only
    python demo.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import qcelemental as qcel

HERE = Path(__file__).resolve().parent

WATER = {"HOH", "WAT", "TIP3", "SOL"}
IONS = {"ZN": 2, "MG": 2, "CA": 2, "MN": 2, "FE": 3, "NA": 1, "K": 1, "CL": -1}
NEUTRAL = {
    "ALA", "ASN", "GLN", "GLY", "ILE", "LEU", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "VAL", "ACE", "NME", "NHE", "ASH", "GLH", "LYN", "CYX",
    "HID", "HIE",
}
CAPS = {"ACE", "NME", "NHE", "NH2"}
CARBOXY_TERMINAL = {"OXT", "OT1", "OT2", "OC1", "OC2"}
# Elements AP2 is parameterized for (apnet_pt.constants.z_to_ind).
SUPPORTED_ELEMENTS = {"H", "B", "C", "N", "O", "F", "Na", "P", "S", "Cl", "Br"}
# Hydrogen-cap X-H lengths (angstrom) for each bond the pocket cut can sever.
CAP_LENGTH = {("N", "C"): 1.01, ("C", "N"): 1.09, ("SG", "SG"): 1.34}
CLASH = 1.2  # any intermolecular contact shorter than this is a bad pose


# --------------------------------------------------------------------------
# PDB parsing and residue charges
# --------------------------------------------------------------------------
def _read_pdb(path: Path) -> dict:
    """Group ATOM/HETATM records into residues, dropping water and alt sites."""
    residues: dict[tuple, dict] = {}
    for line in path.read_text().splitlines():
        if line[:6] not in ("ATOM  ", "HETATM") or line[16] not in " A":
            continue
        resname = line[17:20].strip()
        if resname in WATER:
            continue
        element = line[76:78].strip().capitalize()
        if not element:
            raise ValueError(f"{path}: atom without element column: {line!r}")
        charge_field = line[78:80].strip()  # e.g. "1+" or "1-"
        formal = int(charge_field[::-1]) if charge_field else 0
        key = (line[21], int(line[22:26]), line[26])
        res = residues.setdefault(key, {
            "resname": resname, "hetatm": line.startswith("HETATM"),
            "names": [], "elements": [], "xyz": [], "formal": 0})
        res["names"].append(line[12:16].strip())
        res["elements"].append(element)
        res["xyz"].append([float(line[30 + 8 * i:38 + 8 * i]) for i in range(3)])
        res["formal"] += formal
    for res in residues.values():
        res["xyz"] = np.asarray(res["xyz"])
        res["heavy"] = res["xyz"][[e != "H" for e in res["elements"]]]
    return residues


def _residue_charge(resname: str, names: set[str], disulfide: bool) -> int:
    """Formal charge implied by the hydrogens actually present in a residue."""
    if resname in IONS:
        return IONS[resname]
    rules = {
        "ASP": 0 if "HD2" in names else -1,
        "GLU": 0 if "HE2" in names else -1,
        "LYS": 1 if {"HZ1", "HZ2", "HZ3"} <= names else 0,
        "ARG": 1, "HIP": 1, "CYM": -1,
        "HIS": 1 if {"HD1", "HE2"} <= names else 0,
        "CYS": 0 if ("HG" in names or disulfide) else -1,
        "TYR": 0 if "HH" in names else -1,
    }
    if resname in rules:
        charge = rules[resname]
    elif resname in NEUTRAL:
        charge = 0
    else:
        raise ValueError(f"no charge rule for residue {resname}")
    if resname not in CAPS:
        # Free ammonium N-terminus: OpenMM writes H/H2/H3, Amber H1/H2/H3.
        if {"H2", "H3"} <= names:
            charge += 1
        if names & CARBOXY_TERMINAL:
            charge -= 1
    return charge


def _atom(res: dict, name: str) -> np.ndarray:
    return res["xyz"][res["names"].index(name)]


def _protein_bonds(protein: dict) -> list[tuple]:
    """Inter-residue bonds a pocket cut can sever: peptide C-N and S-S."""
    bonds = []
    keys = sorted(protein)
    for left, right in zip(keys, keys[1:]):
        a, b = protein[left], protein[right]
        if left[0] == right[0] and "C" in a["names"] and "N" in b["names"]:
            if np.linalg.norm(_atom(a, "C") - _atom(b, "N")) < 1.8:
                bonds.append((left, "C", right, "N"))
    sulfurs = [k for k in keys if protein[k]["resname"] in ("CYS", "CYX")
               and "SG" in protein[k]["names"]]
    for i, a in enumerate(sulfurs):
        for b in sulfurs[i + 1:]:
            gap = _atom(protein[a], "SG") - _atom(protein[b], "SG")
            if np.linalg.norm(gap) < 2.5:
                bonds.append((a, "SG", b, "SG"))
    return bonds


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------
def pdb_to_dimer(
    pdb_path: str | Path,
    ligand_resname: str | None = None,
    ligand_charge: int | None = None,
    pocket_radius: float = 12.0,
) -> qcel.models.Molecule:
    """Build a pocket+ligand qcelemental dimer from a prepared complex PDB.

    Parameters
    ----------
    pdb_path
        Protonated protein-ligand complex with element columns filled in.
        Waters are dropped; ions stay in the pocket.
    ligand_resname
        Residue name of the ligand.  If None, the single HETATM residue that
        is neither water nor a known ion is used.
    ligand_charge
        Net ligand charge.  If None, the sum of the PDB formal-charge columns
        (79-80) is used, which is 0 when those columns are blank.
    pocket_radius
        Whole residues with any heavy atom within this distance (angstrom) of
        any ligand heavy atom form the pocket.

    Returns
    -------
    qcelemental.models.Molecule
        Fragment 0 is the hydrogen-capped pocket, fragment 1 the ligand.
    """
    pdb_path = Path(pdb_path)
    residues = _read_pdb(pdb_path)
    if ligand_resname is None:
        candidates = {r["resname"] for r in residues.values()
                      if r["hetatm"] and r["resname"] not in IONS}
        if len(candidates) != 1:
            raise ValueError(f"pass ligand_resname; HETATM residues: {candidates}")
        ligand_resname = candidates.pop()
    lig_keys = [k for k, r in residues.items() if r["resname"] == ligand_resname]
    if not lig_keys:
        raise ValueError(f"no residue named {ligand_resname} in {pdb_path}")
    ligand = {
        "elements": [e for k in lig_keys for e in residues[k]["elements"]],
        "xyz": np.vstack([residues[k]["xyz"] for k in lig_keys]),
    }
    if ligand_charge is None:
        ligand_charge = sum(residues[k]["formal"] for k in lig_keys)
    protein = {k: r for k, r in residues.items() if k not in lig_keys}

    bonds = _protein_bonds(protein)
    bridged = {k for a, an, b, bn in bonds if an == "SG" for k in (a, b)}
    lig_heavy = ligand["xyz"][[e != "H" for e in ligand["elements"]]]
    pocket = []
    for key, res in sorted(protein.items()):
        if len(res["heavy"]) == 0:
            continue
        gap = np.linalg.norm(res["heavy"][:, None] - lig_heavy[None], axis=-1)
        if gap.min() <= pocket_radius:
            pocket.append(key)
    chosen = set(pocket)

    elements, xyz = [], []
    pocket_charge = 0
    for key in pocket:
        res = protein[key]
        elements += res["elements"]
        xyz.append(res["xyz"])
        pocket_charge += _residue_charge(
            res["resname"], set(res["names"]), key in bridged)
    caps = []
    for a, a_name, b, b_name in bonds:
        if (a in chosen) == (b in chosen):
            continue
        (inside, i_name), (outside, o_name) = (
            ((a, a_name), (b, b_name)) if a in chosen else ((b, b_name), (a, a_name)))
        origin = _atom(protein[inside], i_name)
        direction = _atom(protein[outside], o_name) - origin
        length = CAP_LENGTH[(i_name, o_name)]
        caps.append(origin + length * direction / np.linalg.norm(direction))
    elements += ["H"] * len(caps)
    pocket_xyz = np.vstack(xyz + ([np.asarray(caps)] if caps else []))

    _check(elements, pocket_xyz, pocket_charge, ligand, ligand_charge)
    n_pocket = len(elements)
    n_total = n_pocket + len(ligand["elements"])
    print(f"{pdb_path.name}: ligand {ligand_resname} ({len(ligand['elements'])} "
          f"atoms, charge {ligand_charge:+d}); pocket {len(pocket)} residues, "
          f"{n_pocket} atoms incl. {len(caps)} H caps, charge {pocket_charge:+d}")
    return qcel.models.Molecule(
        name=f"{pdb_path.stem}",
        symbols=elements + ligand["elements"],
        geometry=np.vstack([pocket_xyz, ligand["xyz"]])
        / qcel.constants.bohr2angstroms,
        fragments=[list(range(n_pocket)), list(range(n_pocket, n_total))],
        fragment_charges=[float(pocket_charge), float(ligand_charge)],
        fragment_multiplicities=[1, 1],
        molecular_charge=float(pocket_charge + ligand_charge),
        molecular_multiplicity=1,
        fix_com=True, fix_orientation=True,
    )


def _check(pocket_elements, pocket_xyz, pocket_charge, ligand, ligand_charge):
    """Reject inputs AP2 cannot score meaningfully."""
    unsupported = set(pocket_elements + ligand["elements"]) - SUPPORTED_ELEMENTS
    if unsupported:
        raise ValueError(f"elements outside AP2's training set: {unsupported}")
    for name, symbols, charge in (("pocket", pocket_elements, pocket_charge),
                                  ("ligand", ligand["elements"], ligand_charge)):
        electrons = sum(qcel.periodictable.to_Z(s) for s in symbols) - charge
        if electrons % 2:
            raise ValueError(
                f"{name} has an odd electron count with charge {charge:+d}; "
                "check protonation or pass ligand_charge")
    closest = np.linalg.norm(
        pocket_xyz[:, None] - ligand["xyz"][None], axis=-1).min()
    if closest < CLASH:
        raise ValueError(f"pocket-ligand contact {closest:.2f} A: clashing pose")


def predict_ap2(
    dimers: list[qcel.models.Molecule],
    weights: str = "ap2_tf_paper",
    batch_size: int = 1,
) -> np.ndarray:
    """AP-Net2 SAPT interaction energies (kcal/mol) for a list of dimers.

    Returns an array of shape (N, 5) with columns
    [E_total, E_elst, E_exch, E_ind, E_disp], averaged over the ensemble.
    """
    import torch

    from apnet_pt.pretrained_models import apnet2_model_predict

    with torch.inference_mode():
        return apnet2_model_predict(
            dimers, compile=False, batch_size=batch_size, weights=weights)


COMPONENTS = ["E_total", "E_elst", "E_exch", "E_ind", "E_disp"]


def rank_series(
    table: str | Path,
    batch_size: int = 32,
    weights: str = "ap2_tf_paper",
):
    """Score every complex in a ligand table and rank them against experiment.

    `table` is a CSV with columns ligand, pdb (relative to the CSV),
    ligand_charge and dg_exp_kcal_mol.  Returns (components, ranking)
    DataFrames: AP2 components sorted by E_total, and AP2 vs experiment ranks
    (rank 1 = most negative E_total / tightest binder).
    `weights` selects the AP2 ensemble; the default is the paper's weights.
    """
    import pandas as pd

    table = Path(table)
    ligands = pd.read_csv(table)
    dimers = [pdb_to_dimer(table.parent / row.pdb,
                           ligand_charge=int(row.ligand_charge))
              for row in ligands.itertuples()]
    energies = predict_ap2(dimers, batch_size=batch_size, weights=weights)

    components = pd.DataFrame(energies, columns=COMPONENTS)
    components.insert(0, "ligand", ligands["ligand"])
    components = components.sort_values("E_total", ignore_index=True)

    ranking = pd.DataFrame({
        "ligand": ligands["ligand"],
        "E_total_ap2": energies[:, 0],
        "dG_exp": ligands["dg_exp_kcal_mol"],
    })
    ranking["rank_ap2"] = ranking["E_total_ap2"].rank(method="min").astype(int)
    ranking["rank_exp"] = ranking["dG_exp"].rank(method="min").astype(int)
    ranking = ranking.sort_values("rank_exp", ignore_index=True)
    return components, ranking


if __name__ == "__main__":
    import pandas as pd
    from scipy.stats import kendalltau, spearmanr

    pd.set_option("display.float_format", "{:.2f}".format)
    components, ranking = rank_series(HERE / "cdk2_ligands.csv")

    print("\nAP2 (ap2_tf_paper) SAPT components, kcal/mol, sorted by E_total:")
    print(components.to_string(index=False))

    print("\nCDK2 ligand ranking: AP2 E_total vs experimental dG (kcal/mol):")
    print(ranking.to_string(index=False))
    tau = kendalltau(ranking["E_total_ap2"], ranking["dG_exp"]).statistic
    rho = spearmanr(ranking["E_total_ap2"], ranking["dG_exp"]).statistic
    print(f"\nKendall tau = {tau:.2f}   Spearman rho = {rho:.2f}   "
          f"(n = {len(ranking)})")
