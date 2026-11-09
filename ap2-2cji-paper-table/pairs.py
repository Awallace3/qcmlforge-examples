"""Validate ap2_tf_paper AP-Net2 predictions against the paper's matched-pair table.

For each protein-ligand matched pair in PL-matched-pairs/, compute the
protein...ligand interaction energy with AP-Net2 (one or more weight sets), read the
reference SAPT0 total from the Psi4 outputs, and form
    ddE = IE(A) - IE(B)
with (A, B) = (Cl, Me), (Br, CN) for 2O7N, and (2Y5G, 2Y5H). This sign
convention reproduces the SAPT0 column of the paper table exactly.
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import qcelemental as qcel

import apnet_pt

ROOT = Path(__file__).resolve().parent / "PL-matched-pairs"
GEOM_DIR = ROOT / "geoms"
SAPT_DIR = ROOT / "SAPT0_outputs"

# Paper table: PDB -> (SAPT0 ddEint, AP-Net ddEint), kcal/mol
PAPER = {
    "2CJI": (-2.07, -1.67),
    "2O7N": (+4.16, +4.26),
    "2PR3": (-1.03, -1.42),
    "2UZT": (+0.17, +0.91),
    "2W26": (-2.00, -1.00),
    "2Y5(G/H)": (-3.47, -3.65),
    "3ENS": (-4.95, -2.94),
    "4YFF": (+0.86, +0.38),
    "4YHT": (-2.23, -2.97),
}

# PDB -> (system A, system B); ddE = IE(A) - IE(B). Names are xyz stems
# without the -prot/-ligand suffix.
PAIRS = {
    "2CJI": ("2CJI-cl-s5-v2", "2CJI-me-s5-v2"),
    "2O7N": ("2O7N-br-s5", "2O7N-cn-s5"),
    "2PR3": ("2PR3-cl-s5-v2", "2PR3-me-s5-v2"),
    "2UZT": ("2UZT-cl-s5", "2UZT-me-s5"),
    "2W26": ("2W26-cl-s5-v2", "2W26-me-s5-v2"),
    "2Y5(G/H)": ("2Y5G-s5-v2", "2Y5H-s5-v2"),
    "3ENS": ("3ens-cl-s5-v2", "3ens-me-s5-v2"),
    "4YFF": ("4YFF-cl-s5", "4YFF-me-s5"),
    "4YHT": ("4YHT-cl-s5", "4YHT-me-s5"),
}


def load_charges():
    props = pd.read_csv(ROOT / "Properties.csv")
    return {
        t.lower(): int(q) for t, q in zip(props["Title"], props["Molecular charge"])
    }


def xyz_block(path):
    lines = path.read_text().splitlines()
    natom = int(lines[0])
    return "\n".join(lines[2 : 2 + natom])


def build_dimer(system, charges):
    prot = GEOM_DIR / f"{system}-prot.xyz"
    lig = GEOM_DIR / f"{system}-ligand.xyz"
    q_prot = charges[f"{system}-prot".lower()]
    q_lig = charges[f"{system}-ligand".lower()]
    return qcel.models.Molecule.from_data(
        f"{q_prot} 1\n{xyz_block(prot)}\n--\n{q_lig} 1\n{xyz_block(lig)}\n"
        "units angstrom\nno_com\nno_reorient\n"
    )


def sapt0_total(system):
    pdb, rest = system.split("-", 1)  # outputs use upper-case PDB IDs (3ENS vs 3ens)
    out = SAPT_DIR / f"{pdb.upper()}-{rest}-prot-ligand.out"
    m = re.search(r"Total SAPT0\s+\S+ \[mEh\]\s+(\S+) \[kcal/mol\]", out.read_text())
    return float(m.group(1))


def evaluate(weights, systems, mols):
    """Per-system IEs and per-pair ddE for one AP-Net2 weight set."""
    preds = apnet_pt.pretrained_models.apnet2_model_predict(
        mols,
        compile=False,
        batch_size=2,
        weights=weights,
    )

    cols = ["IE", "Elst", "Exch", "Ind", "Disp"]
    ie = pd.DataFrame(preds, index=systems, columns=[f"AP2 {c}" for c in cols])
    ie.insert(0, "SAPT0 IE", [sapt0_total(s) for s in systems])
    ie["AP2 - SAPT0"] = ie["AP2 IE"] - ie["SAPT0 IE"]
    ie.index.name = "system"

    rows = []
    for pdb, (a, b) in PAIRS.items():
        sapt_ref, ap2_ref = PAPER[pdb]
        sapt_dd = ie.at[a, "SAPT0 IE"] - ie.at[b, "SAPT0 IE"]
        ap2_dd = ie.at[a, "AP2 IE"] - ie.at[b, "AP2 IE"]
        rows.append(
            {
                "PDB": pdb,
                "pair (A - B)": f"{a.split('-s5')[0]} - {b.split('-s5')[0]}",
                "SAPT0 ddE (paper)": sapt_ref,
                "SAPT0 ddE (outputs)": sapt_dd,
                "AP-Net ddE (paper)": ap2_ref,
                "AP-Net ddE (ours)": ap2_dd,
                "AP-Net diff (ours - paper)": ap2_dd - ap2_ref,
                "SAPT0 diff (ours - paper)": sapt_dd - sapt_ref,
            }
        )
    dd = pd.DataFrame(rows).set_index("PDB")
    # Agreement at the table's reported precision (0.01 kcal/mol, rounded).
    tol = 0.005 + 1e-9
    dd["AP-Net matches"] = dd["AP-Net diff (ours - paper)"].abs() <= tol
    dd["SAPT0 matches"] = dd["SAPT0 diff (ours - paper)"].abs() <= tol
    dd["sign agrees w/ SAPT0"] = np.sign(dd["AP-Net ddE (ours)"]) == np.sign(
        dd["SAPT0 ddE (outputs)"]
    )

    print(f"\n===== weights = {weights} =====")
    print("\nInteraction energies (kcal/mol)")
    print(ie)
    print("\nMatched-pair ddEint (kcal/mol), ddE = IE(A) - IE(B)")
    print(dd)

    d = dd["AP-Net diff (ours - paper)"]
    print(
        f"\nAP-Net vs paper: MAE {d.abs().mean():.2f}, max |diff| {d.abs().max():.2f}, "
        f"{dd['AP-Net matches'].sum()}/{len(dd)} match to 0.01"
    )
    e = dd["AP-Net ddE (ours)"] - dd["SAPT0 ddE (outputs)"]
    print(f"AP-Net (ours) vs SAPT0 ddE: MAE {e.abs().mean():.2f}")
    return ie, dd


def summary_table(dd_by_weights):
    """Side-by-side ddE: SAPT0, paper AP-Net, and each weight set, plus MAE rows."""
    first = next(iter(dd_by_weights.values()))
    summary = pd.DataFrame(
        {
            "SAPT0": first["SAPT0 ddE (outputs)"],
            "Paper AP-Net": first["AP-Net ddE (paper)"],
        }
    )
    for w, dd in dd_by_weights.items():
        summary[w] = dd["AP-Net ddE (ours)"]
    model_cols = ["Paper AP-Net", *dd_by_weights]
    mae_paper = {c: (summary[c] - summary["Paper AP-Net"]).abs().mean() for c in model_cols}
    max_paper = {c: (summary[c] - summary["Paper AP-Net"]).abs().max() for c in model_cols}
    mae_sapt = {c: (summary[c] - summary["SAPT0"]).abs().mean() for c in model_cols}
    sign_ok = {
        c: f"{(np.sign(summary[c]) == np.sign(summary['SAPT0'])).sum()}/{len(summary)}"
        for c in model_cols
    }
    stats = pd.DataFrame(
        [mae_paper, max_paper, mae_sapt],
        index=["MAE vs paper AP-Net", "max |diff| vs paper AP-Net", "MAE vs SAPT0"],
    )
    table = pd.concat([summary, stats])
    table.index.name = "PDB"
    signs = pd.Series(sign_ok, name="sign agrees w/ SAPT0")
    return table, signs


def main(weights_list=("ap2_tf_paper", "qcmlforge")):
    pd.set_option("display.float_format", "{:.2f}".format)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", None)

    charges = load_charges()
    systems = [s for pair in PAIRS.values() for s in pair]
    mols = [build_dimer(s, charges) for s in systems]

    results = {w: evaluate(w, systems, mols) for w in weights_list}
    table, signs = summary_table({w: dd for w, (_, dd) in results.items()})

    print("\n===== Summary: matched-pair ddEint (kcal/mol) =====")
    print(table.to_string(na_rep=""))
    print("\nsign of ddE agrees with SAPT0:")
    print(signs.to_string())
    return results, table


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--weights",
        nargs="+",
        default=["ap2_tf_paper", "qcmlforge"],
        help="AP-Net2 weight set(s) to compare (qcmlforge = package default)",
    )
    results, summary = main(parser.parse_args().weights)
