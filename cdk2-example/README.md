# AP-Net2 from protein-ligand PDBs: ranking CDK2 ligands

Turn prepared (protonated) protein-ligand complex PDBs into qcelemental
pocket+ligand dimers, predict their SAPT interaction energies (IEs) with the
`ap2_tf_paper` AP-Net2 ensemble, and rank the 16 CDK2 ligands against
experimental binding free energies.

```bash
conda activate qcml
export QCMLFORGE_AUTO_DOWNLOAD_PRETRAINED=1   # first run downloads the weights
CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=8 python demo.py
```

Single complex:

```python
from demo import pdb_to_dimer, predict_ap2

dimer = pdb_to_dimer("pdbs/cdk2_17.pdb")   # ligand found automatically
energies = predict_ap2([dimer])            # (N, 5): total, elst, exch, ind, disp
```

Whole series (what `python demo.py` runs):

```python
from demo import rank_series

components, ranking = rank_series("cdk2_ligands.csv")
```

## Files

| Path | Contents |
|---|---|
| `demo.py` | `pdb_to_dimer`, `predict_ap2`, `rank_series` |
| `cdk2_ligands.csv` | ligand, PDB path, ligand charge, experimental dG (kcal/mol) |
| `pdbs/cdk2_*.pdb` | 16 protonated, minimized CDK2 complexes |

## Output

The first table lists AP2 components sorted by E_total, the second the ranking
against experiment (rank 1 = most negative E_total or dG):

```text
AP2 (ap2_tf_paper) SAPT components, kcal/mol, sorted by E_total:
ligand  E_total  E_elst  E_exch  E_ind  E_disp
    32  -132.99  -97.42  122.10 -47.84 -109.84
  1h1s  -132.68  -96.62  116.51 -45.98 -106.60
    28  -131.01  -93.35  119.34 -48.66 -108.34
    29  -123.41  -76.98   96.52 -41.49 -101.46
    30  -113.40  -71.08  103.02 -40.76 -104.58
    20  -107.67  -75.10  105.97 -40.72  -97.82
  1oiu  -105.12  -65.07   92.56 -38.79  -93.82
  1oiy  -104.34  -62.59   85.03 -33.96  -92.81
    22  -104.15  -62.29   95.22 -37.49  -99.59
    26  -103.96  -62.78   85.22 -35.01  -91.40
  1oi9  -100.71  -61.93   82.97 -33.96  -87.79
    21  -100.48  -59.48   85.79 -35.28  -91.51
    31   -94.29  -51.74   84.99 -34.51  -93.03
    17   -94.25  -52.70   80.00 -32.59  -88.96
  1h1r   -93.00  -52.44   79.76 -32.05  -88.27
  1h1q   -91.69  -53.32   79.50 -31.51  -86.36

CDK2 ligand ranking: AP2 E_total vs experimental dG (kcal/mol):
ligand  E_total_ap2  dG_exp  rank_ap2  rank_exp
  1h1s      -132.68  -11.25         2         1
    28      -131.01  -11.11         3         2
    29      -123.41   -9.88         4         3
    30      -113.40   -9.81         5         4
  1oiy      -104.34   -9.79         8         5
    32      -132.99   -9.75         1         6
  1oi9      -100.71   -9.74        11         7
    31       -94.29   -9.54        13         8
  1oiu      -105.12   -9.08         7         9
    20      -107.67   -8.72         6        10
    26      -103.96   -8.43        10        11
  1h1q       -91.69   -8.18        16        12
    22      -104.15   -7.86         9        13
    21      -100.48   -7.83        12        14
  1h1r       -93.00   -7.67        15        15
    17       -94.25   -7.04        14        16

Kendall tau = 0.62   Spearman rho = 0.79   (n = 16)
```

An IE is not a binding free energy: it is a single point at the PDB geometry,
with no solvation, entropy, ligand strain or conformational sampling. Its
magnitude is roughly 10x dG; only the ordering within a congeneric series is
meaningful.

## Validation against TAQO22

The full CPU run was checked against
`~/data/TAQO22/campaigns/cdsg-cdk2-v1/`. All 16 bundled PDBs are
byte-identical to the campaign's minimized complexes, and all experimental
energies and ligand charges match its `crosswalk.csv`.
Using the campaign's selected residues, the demo's charge rules reproduce
all 16 stored pocket charges; severed-bond/cap counts also match, and cap
coordinates agree within 0.002 Å (the input PDBs round coordinates).
The upstream pinned-energy test
`tests/test_ap2_tf_paper_route.py::test_paper_ensemble_reproduces_reference_interaction_energy`
also passes in the same `qcml` environment, checking all five ensemble
output columns against its reference with `atol=1e-5` kcal/mol.

| Comparison against experiment (16 ligands) | Kendall tau | Spearman rho |
|---|---:|---:|
| This demo, `ap2_tf_paper` | 0.62 | 0.79 |
| Stored TAQO22 AP2 five-member ensemble | 0.62 | 0.79 |

The metrics agree to two decimals, not exactly: the demo's rho is
0.794118 versus 0.791176 in TAQO22. The two predicted energy rankings
themselves have rho = 0.93 and tau = 0.83.

This is a consistency check, **not numerical reproduction**: the campaign
used the older QCMLForge `am_ensemble`/`ap2_ensemble` checkpoints, not
`ap2_tf_paper`, and fixed the residue selection before minimization around
the reference ligand's mass-weighted centroid, rather than selecting within
12 Å of any minimized ligand heavy atom. For ligand 17, the stored dimer has
933 atoms versus 1398 here. Mean absolute total-energy difference is
2.75 kcal/mol (maximum 7.41), calculated from this demo's printed energies
and the mean of the five stored `predictions/ap2-member-*/cdk2.csv` files.
The models and pockets differ, so individual ranks need not match.

### Current QCMLForge weights on identical demo pockets

An additional CPU run used the corrected five-member `qcmlforge` ensemble
on the same 16 PDBs, 12 Å ligand-contact pockets and batch size 4 as the
paper-weight run:

| Ensemble | Kendall tau | Spearman rho |
|---|---:|---:|
| `ap2_tf_paper` | 0.62 | 0.79 |
| `qcmlforge` | 0.58 | 0.75 |

The paper ensemble performs slightly better on this particular series;
16 ligands alone do not establish a general accuracy advantage.
These current weights are distinct from the older stored TAQO22 AP2 weights.
To select them without changing the demo's paper-weight default:

```python
components, ranking = rank_series(
    "cdk2_ligands.csv", batch_size=4, weights="qcmlforge"
)
```

## What `pdb_to_dimer` assumes and does

- **Input:** all hydrogens present, one protonation state per residue, element
  columns filled in, ligand in the same file as the protein. Nothing is
  protonated or repaired. Waters are dropped; ions stay in the pocket.
- **Pocket:** whole residues with any heavy atom within `pocket_radius`
  (default 12 Å) of any ligand heavy atom. Severed peptide (C–N) and disulfide
  (S–S) bonds are capped with hydrogens.
- **Charges:** each residue's charge comes from the hydrogens it carries (ASP
  without HD2 is −1, HIS with HD1 and HE2 is +1, termini, and so on). The ligand
  charge is `ligand_charge` if you pass it, otherwise the sum of the PDB
  formal-charge columns (0 when those are blank, as they are in OpenMM output).
  **Pass `ligand_charge` for a charged ligand**; `rank_series` takes it from
  the CSV.
- **Checks:** it raises `ValueError` for elements AP2 was not trained on, an
  odd electron count in either fragment (usually a wrong ligand charge), or a
  pocket-ligand contact under 1.2 Å.

The pocket rules come from `~/data/TAQO22/apnet/build_pockets.py`, and the cap
bond lengths from `~/gits/protein-prep`.

## Where the complexes come from

The 16 complexes are the CDK2 set of the Schrödinger JACS FEP benchmark,
prepared and minimized by protein-prep for the TAQO22 campaign
(`~/data/TAQO22/campaigns/cdsg-cdk2-v1/prepared/cdk2/*/complex_minimized.pdb`);
experimental dG and ligand charges come from that campaign's `crosswalk.csv`
(the source SDF's `r_exp_dg`). protein-prep protonated the receptor with
OpenMM's pH 7.4 rules and mapped the phosphothreonine (TPO) to THR. With an
effectively infinite radius, the demo's charge rules give the whole receptor a
net charge of 0, the same value protein-prep's force-field charges give.
