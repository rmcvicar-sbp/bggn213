Structure based drug design
================

## 

Me: First, when we import the protein structure from PDB (it will be
bound to a compound), we want to separate the protein structure from
that molecule.

Barry: Here we download and clean up the HIV-Pr structure(PDB: 1HSG)
from PDB. Make a separate set of “protein only” and “ligand only” files.

``` r
library(bio3d)
file.name<- get.pdb("1HSG")
```

    ## Warning in get.pdb("1HSG"): ./1HSG.pdb exists. Skipping download

``` r
file.name
```

    ## [1] "./1HSG.pdb"

We’ll use `read.pdb`, `atom.select()`, and write `write.pdb()` functions
to make our separate “protein only” and “ligand only” file

``` r
read.pdb("1HSG")
```

    ##   Note: Accessing on-line PDB file

    ## 
    ##  Call:  read.pdb(file = "1HSG")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
hiv<- read.pdb(file.name)
?atom.select
```

``` r
hiv.protein.only<- atom.select(hiv, "protein", value=TRUE)
hiv.ligand.only<- atom.select(hiv, "ligand", value=TRUE)
```

``` r
?write.pdb
write.pdb(hiv.protein.only, file="1HSG_protein.pdb")
write.pdb(hiv.ligand.only, file="1HSG_ligand.pdb")
```

The next step is to use AutoDockTools to setup protein docking input
Need hydrogens to calculate energies when docking, done through
“AutoDockTools” charge and atom-type is the difference between .pdbqt
vs .pdb files (charge & type)

``` r
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```

``` r
library(bio3d)
res <- read.pdb("all.pdbqt", multi= TRUE)
write.pdb(res, "results.pdb")
```
