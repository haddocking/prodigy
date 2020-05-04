# PRODIGY / Binding Affinity Prediction

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1193244.svg)](https://doi.org/10.5281/zenodo.1193244)

Collection of scripts to predict binding affinity values
for protein-protein complexes from atomic structures.

The online version of PRODIGY predictor can be found here:
* [PRODIGY](http://milou.science.uu.nl/services/PRODIGY/)

Details of the binding affinity predictor implemented in PRODIGY can be found here:
* [Contacts-based model](http://www.ncbi.nlm.nih.gov/pubmed/26193119)

# Requirements

Python 3

# Installation

```
git clone http://github.com/haddocking/prodigy
cd prodigy
git checkout python3_package
pip install .

# Have fun!
```

# Usage

```bash
prodigy <pdb file> [--selection <chain1><chain2>]
```

To get a list of all the possible options.
```bash
prodigy --help 
```

# Information about dependencies
The scripts rely on [Biopython](www.biopython.org) to validate the PDB structures and calculate
interatomic distances. [freesasa](https://github.com/mittinatten/freesasa), with the parameter
set used in NACCESS ([Chothia, 1976](http://www.ncbi.nlm.nih.gov/pubmed/994183)), is also
required for calculating the buried surface area.

**DISCLAIMER**: given the different software to calculate solvent accessiblity, predicted
values might differ (very slightly) from those published in the reference implementations.
The correlation of the actual atomic accessibilities is over 0.99, so we expect these
differences to be very minor.

To install and use the scripts, just clone the git repository or download the tarball zip
archive. Make sure `freesasa` and Biopython are accessible to the Python scripts
through the appropriate environment variables ($PYTHONPATH).

# License
These utilities are open-source and licensed under the Apache License 2.0. For more information
read the LICENSE file.

# Citing us
If our predictive model or any scripts are useful to you, consider citing them in your
publications:

**Xue L, Rodrigues J, Kastritis P, Bonvin A.M.J.J, Vangone A.**: PRODIGY: a web server for predicting the binding affinity of protein-protein complexes. *Bioinformatics* (2016) ([link](http://bioinformatics.oxfordjournals.org/content/early/2016/08/27/bioinformatics.btw514))

**Anna Vangone and Alexandre M.J.J. Bonvin**: Contacts-based prediction of binding affinity in protein-protein complexes. *eLife*, e07454 (2015) ([link](http://www.ncbi.nlm.nih.gov/pubmed/26193119))

**Panagiotis L. Kastritis , João P.G.L.M. Rodrigues, Gert E. Folkers, Rolf Boelens, Alexandre M.J.J. Bonvin**: Proteins Feel More Than They See: Fine-Tuning of Binding Affinity by Properties of the Non-Interacting Surface. *Journal of Molecular Biology*, 14, 2632–2652 (2014). ([link](http://www.ncbi.nlm.nih.gov/pubmed/24768922))

# Contact
For questions about PRODIGY usage, please contact the team at: prodigy.bonvinlab@gmail.com
