# PRODIGY / Binding Affinity Prediction 

![PyPI - License](https://img.shields.io/pypi/l/prodigy-prot)
![PyPI - Status](https://img.shields.io/pypi/status/prodigy-prot)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/prodigy-prot)
[![ci](https://github.com/haddocking/prodigy/actions/workflows/ci.yml/badge.svg)](https://github.com/haddocking/prodigy/actions/workflows/ci.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/98180cbac27d4a5aaf46a3dd72c3174d)](https://www.codacy.com/gh/haddocking/prodigy/dashboard?utm_source=github.com&utm_medium=referral&utm_content=haddocking/prodigy&utm_campaign=Badge_Grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/98180cbac27d4a5aaf46a3dd72c3174d)](https://www.codacy.com/gh/haddocking/prodigy/dashboard?utm_source=github.com&utm_medium=referral&utm_content=haddocking/prodigy&utm_campaign=Badge_Coverage)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)

[![SQAaaS badge shields.io](https://img.shields.io/badge/sqaaas%20software-gold-yellow)](https://api.eu.badgr.io/public/assertions/w8HcpcH4Svi3-UZ93LHHMA "SQAaaS gold badge achieved")




[![SQAaaS badge](https://github.com/EOSC-synergy/SQAaaS/raw/master/badges/badges_150x116/badge_software_gold.png)](https://api.eu.badgr.io/public/assertions/w8HcpcH4Svi3-UZ93LHHMA "SQAaaS gold badge achieved")

* * *

PRODIGY is also available as a web service @ [wenmr.science.uu.nl/prodigy](https://wenmr.science.uu.nl/prodigy/)

## Installation

```text
pip install prodigy-prot
```

If you want to develop PRODIGY, check [DEVELOPMENT](DEVELOPMENT.md) for more details.

## Usage

```bash
prodigy <pdb file> [--selection <chain1><chain2>]
```

To get a list of all the possible options.

```bash
$ prodigy --help
usage: prodigy [-h] [--distance-cutoff DISTANCE_CUTOFF] [--acc-threshold ACC_THRESHOLD] [--temperature TEMPERATURE]
               [--contact_list] [--pymol_selection] [-q] [-V] [--selection A B [A,B C ...]]
               structf

Binding affinity predictor based on Intermolecular Contacts (ICs).

Anna Vangone and Alexandre M.J.J. Bonvin,
Contacts-based prediction of binding affinity in protein-protein complexes.
eLife (2015)

positional arguments:
  structf               Structure to analyse in PDB or mmCIF format

options:
  -h, --help            show this help message and exit
  --distance-cutoff DISTANCE_CUTOFF
                        Distance cutoff to calculate ICs
  --acc-threshold ACC_THRESHOLD
                        Accessibility threshold for BSA analysis
  --temperature TEMPERATURE
                        Temperature (C) for Kd prediction
  --contact_list        Output a list of contacts
  --pymol_selection     Output a script to highlight the interface (pymol)
  -q, --quiet           Outputs only the predicted affinity value
  -V, --version         Print the version and exit.

Selection Options:

      By default, all intermolecular contacts are taken into consideration,
      a molecule being defined as an isolated group of amino acids sharing
      a common chain identifier. In specific cases, for example
      antibody-antigen complexes, some chains should be considered as a
      single molecule.

      Use the --selection option to provide collections of chains that should
      be considered for the calculation. Separate by a space the chains that
      are to be considered _different_ molecules. Use commas to include multiple
      chains as part of a single group:

      --selection A B => Contacts calculated (only) between chains A and B.
      --selection A,B C => Contacts calculated (only) between chains A and C; and B and C.
      --selection A B C => Contacts calculated (only) between chains A and B; B and C; and A and C.


  --selection A B [A,B C ...]
```

## Example

Download the PDB [3BZD](https://www.rcsb.org/structure/3bzd) and run PRODIGY on it.

```bash
$ curl -o 3bzd.pdb https://files.rcsb.org/download/3BZD.pdb
$ prodigy 3bzd.pdb
[+] Reading structure file: /Users/rvhonorato/dbg/3bzd.pdb
[+] Parsed structure file 3bzd (2 chains, 343 residues)
[+] No. of intermolecular contacts: 51
[+] No. of charged-charged contacts: 4
[+] No. of charged-polar contacts: 7
[+] No. of charged-apolar contacts: 6
[+] No. of polar-polar contacts: 7
[+] No. of apolar-polar contacts: 15
[+] No. of apolar-apolar contacts: 12
[+] Percentage of apolar NIS residues: 29.48
[+] Percentage of charged NIS residues: 29.48
[++] Predicted binding affinity (kcal.mol-1):     -9.4
[++] Predicted dissociation constant (M) at 25.0˚C:  1.3e-07
```

Details of the binding affinity predictor implemented in PRODIGY can be found at [10.7554/elife.07454](https://doi.org/10.7554/elife.07454)

## Citing us

If our tool is useful to you, please cite PRODIGY in your publications:

- **Xue L, Rodrigues J, Kastritis P, Bonvin A.M.J.J, Vangone A.**: PRODIGY: a web server for predicting the binding affinity of protein-protein complexes. _Bioinformatics_ (2016) ([10.1093/bioinformatics/btw514](https://doi.org/10.1093/bioinformatics/btw514))

- **Anna Vangone and Alexandre M.J.J. Bonvin**: Contacts-based prediction of binding affinity in protein-protein complexes. _eLife_, e07454 (2015) ([10.7554/eLife.07454](https://doi.org/10.7554/elife.07454))

- **Panagiotis L. Kastritis , João P.G.L.M. Rodrigues, Gert E. Folkers, Rolf Boelens, Alexandre M.J.J. Bonvin**: Proteins Feel More Than They See: Fine-Tuning of Binding Affinity by Properties of the Non-Interacting Surface. _Journal of Molecular Biology_, 14, 2632–2652 (2014). ([10.1016/j.jmb.2014.04.017](https://doi.org/10.1016/j.jmb.2014.04.017))

## Contact

For questions about PRODIGY usage, please reach out the team at [ask.bioexcel.eu](https://ask.bioexcel.eu/)

## Information about dependencies

The scripts rely on [Biopython](www.biopython.org) to validate the PDB structures and calculate
interatomic distances. [freesasa](https://github.com/mittinatten/freesasa), with the parameter
set used in NACCESS ([Chothia, 1976](http://www.ncbi.nlm.nih.gov/pubmed/994183)), is also
required for calculating the buried surface area.

**DISCLAIMER**: given the different software to calculate solvent accessiblity, predicted
values might differ (very slightly) from those published in the reference implementations.
The correlation of the actual atomic accessibilities is over 0.99, so we expect these
differences to be very minor.

---
