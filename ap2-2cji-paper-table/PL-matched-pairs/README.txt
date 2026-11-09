README.txt
============

This directory contains protein-ligand "matched pairs", in which we have
a fragment of a protein and a fragment of a ligand (truncated to keep
the model systems of practical size), and a substituent on the ligand
is changed from X to Y (usually Methyl to Chlorine, but sometimes Br
to CN). In one case (the 2Y5G, 2Y5H pair), the change is the flip of
a heterocycle rather than a chemical substitution.

Some of these systems are for factor Xa and have been previously
studied in "The Surprising Importance of Peptide Bond Contacts
in Drug-Protein Interactions," R. M. Parrish, D. F. Sitkoff,
D. L. Cheney, and C. D. Sherrill, Chem. Eur. J. 23, 7887 (2017) (doi:
10.1002/chem.201701031) at the SAPT0/jun-cc-pVDZ level of theory.
These are: 2W26, 3ENS, 2CJI, and 2PR3.  There may be small differences
in the coordinates in these files vs the 2017 study.

Model systems differ from those examined in "CLIFF: A Component-Based,
Machine-Learned, Intermolecular Force Field," J. B. Schriber,
D. R. Nascimento, A. Koutsoukas, S. A. Spronk, D. L. Cheney,
and C. D. Sherrill, J. Chem. Phys. 154, 184110 (2021) (doi:
10.1063/5.0042989). CLIFF has been limited to neutral systems only,
so the model systems in that study used neutralized sidechains.

The subdirectory 'geoms/' contains XYZ files of the truncated
proteins and fragment ligands.  The XYZ files are named acoording to:
[PDB_ID]-[substituent]-s5(-v2)-[prot or ligand].xyz. The subdirectory
'SAPT0_outputs/' contains Psi4 output files for each complex, and these
files are named similarly to the XYZs.

Each file in 'geoms/' is named by the PDB for the original crystal
structure, the substituent on the ligand (usually Cl or Me, but sometimes
Br or CN), either "s5" or "s5-v2" for a truncated protein model, and
"prot" or "ligand" to specify that the file has coordinates for the
protein or the ligand fragment, respectively.

Properties.csv contains information about the model systems. Rows are
named consistently with geometry XYZ filenames given above, except
that it additionally contains entries XXX-"cx", where "cx" denotes a
protein-ligand complex.  The file contains the number of molecules per
model system (e.g., number of chains + water molecules), number of atoms,
PDB filenames and resolutions, and SMILES strings for ligands.

In the present model systems, all ligands are neutral.  The truncated
protein may be charged, and charges are given in the "Molecular charge"
column of Properties.csv.

2Y5G and 2Y5H comprise a factor Xa matched pair with ligand heterocycles
of different dipolar natures. Because there is no chemical substitution,
these are named without [substituent].

Proteins were subject to protein preparation by Dan Cheney as outlined below. 

Whole protein preparation was conducted using the Protein Prep (Ref
1) utility in Schrodingers Maesro 2023-3 program (Ref 2), as in the
following steps:

1. Assignment of bond orders, charges and hydrogens 
2. Retention of the highest quality structure in the unit cell, where
multiple copies were present
3. Filling in missing atoms of side chains and backbones
4. Modification of the ligand in each of the 10 complexes to afford
a matched molecular pair (MMP), or, in the instance of 3ENS, a
matched triplet. In most cases, this entailed a single modification
of the crystallographic ligand, for example replacing a methyl with
chlorine in 2PR3.  In some cases, multiple, small modifications of
the crystallographic ligand were required to obtain a matched pair for
which binding data was available.
5. Retention of waters with a minimum of 1 hydrogen bond to non-water
molecules
6. Ionization of side chains and terminal amino acids (pH = 7.4 +/- 1)
7. Hydrogen bond networks were optimized, including rotation of water
molecules and flips of Asn and Gln side chains amides.
8. An optimization constrained to 0.3 RMSD was conducted with OPLS4.
9. Cosolvent molecules and ions incompatible with AP-Net were removed.

(1) Truncated protein + ligand fragment:

a. Ligands are truncated similar to what was previously described
(Ref 4, Ref 5, Ref 6).

b. In "s5" systems, proteins are truncated in a manner similar to what
was previously descibed (Ref 4, Ref 5, Ref 6). Terminal methyls created
as a result of protein cuts were observed to adopt arbitrary rotations,
misaligning their hydrogens with corresponding hydrogens in the protein
model from which they were derived.  Hydrogen alignments were recaptured
by bond rotations.

c. Note that in factor Xa whole protein pairs, one will not be able to
rank order entire set by energies as the number and positions (near
the ligand) of waters vary from pair to pair.

----------References----------
Ref 1. Maestro, Schrödinger, LLC, New York, NY, 2023
Ref 2. Protein Preparation Wizard; Epik,  Schrödinger, LLC, New York, NY, 2023; Impact, Schrödinger, LLC, New York, NY; Prime, Schrödinger, LLC, New York, NY, 2023
Ref 3. Spronk, Steven A., et al. "A quantum chemical interaction energy dataset for accurately modeling protein-ligand interactions." Scientific Data 10.1 (2023): 619.
Ref 4. The Surprising Importance of Peptide Bond Contacts in Drug–Protein Interactions
Dr. Robert M. Parrish,  Dr. Doree F. Sitkoff,  Dr. Daniel L. Cheney,  Prof. C. David Sherrill
Ref 5. Jeffrey B. Schriber, Daniel R. Nascimento, Alexios Koutsoukas, Steven A. Spronk, Daniel L. Cheney, C. David Sherrill; CLIFF: A component-based, machine-learned, intermolecular force field. J. Chem. Phys. 14 May 2021; 154 (18): 184110. https://doi.org/10.1063/5.0042989

----------Additional Details-----
Quality control checks were made on the basis of MMP's
- having similar numbers of atoms, residues, molecules and molecular charges
- visual inspection to insure OH and SH containing side chains had similar torsions

- 2O7N-s0-whole exhibits a slight H-H clash with Val223 (2.11 A) which could not alleviated with repeated constrained minimizations. This persisted for the -CN analogue (nitrile) in sets S1, S2, S3
- 2cji-s0-whole exhibits a slight H-H clash with Trp215 cbetaH which could not be alleviated with repeated constrained minimization. This persists for both Cl and Me analogues in sets S1,S2,S3
- 2cji water near Asp189 was determined to be a Ca2+ cation (based on raw xtal data from other fXa structures).
- In the 2pr3 MMP of Sets 1,2,3 there exists for both analogues a slight steric clash between the terminal ligand phenyl and the Trp215 indole ring, which could not be resolved with repeated constrained minimizations.

2o7n LFA-1.
Resolution = 1.7.
PDB: Bioorganic & Medicinal Chemistry Letters 17 (2007) 1908–1911
Biological data: Bioorganic & Medicinal Chemistry Letters 17 (2007) 1908–1911
Compound 15 (2o7n-br)   106 nM IC50 racemate, HeLa/ HSB binding
Compound 18 (2o7n-cn)   52 nM " "
Unit cell: 1 copy

3ens factor Xa
PDB: J. Med. Chem. 2008, 51, 7541–7551 (compound 18)
Biological data: J. Med. Chem. 2008, 51, 7541–7551
Compound 19 (3ens-me)   118 nM IC50 FXa utilizing the cleavage of a synthetic substrate S-2222
Compound 20 (3ens-cl)   2.4 nM " "
Compound 21 (3ens-br)   9 nM " "
Unit cell: 2 copies; Chains A/B selected as being more complete around binding site

2cji factor Xa
PDB: J. Med. Chem. 2007, 50, 1546-1557 (compound 1a)
Biological data: J. Med. Chem. 2007, 50, 1546-1557
Compound  1a (2cji-cl)  6 nM Ki Fluorogenic Factor Xa Assay
Compound  1g (2cji-me)  109 nM " "
Unit cell: 1 copy

2w26 factor Xa
PDB: J. Med. Chem. 2005, 48, 5900-5908 (compound 5)
Biological data: J. Med. Chem. 2005, 48, 5900-5908
Compound: 5 (2w26-cl) 0.7 nM IC50
Compound:  18 (2w26-me) 4.2 nM " "
Unit cell: 1 copy

2prj factor Xa
PDB: (compound 18)
Biological data:
Compound: 23 (2prj-cl) 36 nM IC50
Compound: 24 (2prj-me) 389 nM IC50
Unit cell: 1 copy

Note for following two structures, 2y5g and 2y5h:  Cryst. structures were solved for both analogues, which differed in water content.  Thus for sets s1-s5, water content was omitted except for the two conserved waters in the vicinity of Asp189.
2y5g factor Xa
PDB: ChemMedChem 2013, 8, 397 – 404 (ChemMedChem 2013, 8, 397 – 404)
Biological data: ChemMedChem 2013, 8, 397 – 404
Compound:  2y5g  146 nM IC50?
Unit cell:

2y5h factor Xa
PDB: ChemMedChem 2013, 8, 397 – 404 (ChemMedChem 2013, 8, 397 – 404)
Biological data: ChemMedChem 2013, 8, 397 – 404
Compound:  2y5h 1620 nM IC50?

2uzt Protein Kinase A
Bovine protein used in crystallography and models -- is nearly identical to human; binding site residues are identical.
PDB: J. Med. Chem. 2007, 50, 2990-3003  (compound 7)
Biological data:  Bioorganic & Medicinal Chemistry 15 (2007) 2441–2452
Compound: 9j (4yht-cl) 59  nM IC50
Compound: 9i (4yht-me) 25  nM IC50
Unit cell: 1 copy

4yht B-RAF (Kinase)
PDB:  BIOORG.MED.CHEM.LETT.         V.  26  3355 2016
Resolution: 3.1
Biological data:  Bioorganic & Medicinal Chemistry Letters 26 (2016) 3355–3358 (also: J. Med. Chem. 2015, 58, 7431−7448)
Compound: 4 (4yht-me) 63 nM
Compound: 6 (4yht-cl) 12 nM
Unit cell:  2 copies; Chain B selected

4yff TNNI3K
PDB: Bioorganic & Medicinal Chemistry Letters 26 (2016) 3355–3358 (compound 10)
Resolution: 3.1
Biological data: Bioorganic & Medicinal Chemistry Letters 26 (2016) 3355–3358 (also: J. Med. Chem. 2015, 58, 7431−7448)
Compound: 4 (4yff-me) 40 nM IC50
Compound:  6 (4yff-cl) 32 nM IC50
Unit cell: 4 copies; Chain A selected

