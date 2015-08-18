#Binding Affinity Prediction Tools
Collection of scripts to calculate predictive indexes of binding affinity values
for protein-protein complexes from atomic structures.

The scripts implement several binding affinity predictors:
* [Non-Interacting Surface (NIS) model](http://www.ncbi.nlm.nih.gov/pubmed/24768922)
* [Contacts-based model](http://www.ncbi.nlm.nih.gov/pubmed/26193119)

#Quick & Dirty Installation
```bash
git clone http://github.com/biopython/biopython.git
cd biopython
sudo python setup.py install # Alternatively, install locally but fix $PYTHONPATH

git clone http://github.com/JoaoRodrigues/freesasa
cd freesasa
git checkout with-configure
./configure && make
export FREESASA_SRC="$PWD"

git clone http://github.com/JoaoRodrigues/binding_affinity
git checkout refactor
# Have fun!
```

#Usage
* __Non-Interacting Surface (NIS) model__  
```bash
python predict_NIS.py <pdb file>
```

* __Contacts-based model__  
```bash
python predict_IC.py <pdb file>
```

#Installation & Dependencies
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
through the appropriate environment variables ($FREESASA_SRC and $PYTHONPATH). $FREESASA_SRC is
required by the prediction scripts to locate the freesasa executable and atomic radii
configuration file. It should point to the directory where freesasa is cloned into, containing
both the src/ and config/ sub-directories.

#License
These utilities are open-source and licensed under the Apache License 2.0. For more information
read the LICENSE file.

#Citing us
If any of the predictive models or scripts are useful to you, consider citing them in your
publications:

**Anna Vangone and Alexandre M.J.J. Bonvin**: Contacts-based prediction of binding affinity in protein-protein complexes. Revision in eLife (2015) ([link](http://www.ncbi.nlm.nih.gov/pubmed/26193119))

**Panagiotis L. Kastritis , João P.G.L.M. Rodrigues, Gert E. Folkers, Rolf Boelens, Alexandre M.J.J. Bonvin**: Proteins Feel More Than They See: Fine-Tuning of Binding Affinity by Properties of the Non-Interacting Surface. *Journal of Molecular Biology*, 14, 2632–2652 (2014). ([link](http://www.ncbi.nlm.nih.gov/pubmed/24768922))
