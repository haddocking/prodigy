##########################################################################
##    Calculate inter-residue contacts in protein-protein complexes     ##
##                        Author: Anna Vangone                          ##
##########################################################################
## Calculate the inter-protein distances divided by                     ##
## polar/apolar/charged characteristics of the amino acids              ## 
## in a protein-protein complex.                                        ## 
## Proteins as to be labelleed as chain A and chain B.                  ##
## Deault distance threshold: 5.5 Angstrom                              ##
## usage: ICs-parameters.py <pdb file>                                  ##
## Contacts: a.vangone@uu.nl, a.m.j.j.bonvin@uu.nl                      ##
##########################################################################

egrep ^ATOM $1  > $1.temp
python ./icsdistance.py $1.temp
python ./filter_propery.py
rm $1.temp
rm distance.temp
