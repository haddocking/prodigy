# binding_affinity
Repository containing various scripts to predict the binding affinity of protein-protein complexes from structure

# About
These scripts allow to calculate properties such as inter-residue distances (ICs) and non-interacting surface properties (NIC) in protein-protein complexes. Such values can be combined with specific weights to predict the binding affinity of the protein complex.

#Usage:
ICs: Copy ICs-parameters.sh, icsdistance.py and filter_property.py in the same directory.
     Launch the program as:

         ./ICs-parameters.sh <pdb-file>
          Example: ./ICs-parameter.sh 2OOB.pdb

NIS: Launch the program as:

     ./parameter_0.5.sh <pdb-file>

Please note that the complex has to contain only TWO interacting chains, labelled as A and B.

# Requirement
ICs_parameters.sh requires python 2.7 version.

NISparameters_0.5.sh requires naccess to be installed. See:[www.bioinf.manchester.ac.uk/naccess](www.bioinf.manchester.ac.uk/naccess).

=======

# Licensing

- If ICs calcualtion scripts package was useful to your research, please cite us

**Anna Vangone and Alexandre M.J.J. Bonvin**: Contacts-based prediction of binding affinity in protein-protein complexes. Revision in eLife (2015).



- If NIS calculation script was useful to your research, please cite us

**Panagiotis L. Kastritis , João P.G.L.M. Rodrigues, Gert E. Folkers, Rolf Boelens, Alexandre M.J.J. Bonvin**: Proteins Feel More Than They See: Fine-Tuning of Binding Affinity by Properties of the Non-Interacting Surface. *Journal of Molecular Biology*, 14, 2632–2652 (2014). 

Apache Licence 2.0

