##########################################################################
##      calculate atom-atom distances in protein complexes              ##
##                        Author: Anna Vangone                          ##
##########################################################################

"""
Calculate the inter-protein distances in a 
protein-protein complex, which protein should be 
labelled as chain A and chain B.
Deault distance threshold: 5.5 Angstrom
usage: icsdistance.py <pdb file>
This script is part of the ICs-parameters.sh program
"""

import string
import sys
from math import *
from scipy import *

# define class for protein1
class PDB1:
      """ Atoms characteristics stored here"""
      def __init__ (self,atom1=None,res1=None,chain1=None,Nres1=None,x1=None,y1=None,z1=None):
          self.res=res1
          self.atom=atom1
          self.chain=chain1
          self.Nres=Nres1
          self.x=x1
          self.y=y1
          self.z=z1
# define class for protein2
class PDB2:
      """ Atoms characteristics stored here"""
      def __init__ (self,atom2=None,res2=None,chain2=None,Nres2=None,x2=None,y2=None,z2=None):
          self.res=res2
          self.atom=atom2
          self.chain=chain2
          self.Nres=Nres2
          self.x=x2
          self.y=y2
          self.z=z2

# define a function to calculate distance between atoms
def dist(i,j):
    dx = Protein2[j].x - Protein1[i].x
    dy = Protein2[j].y - Protein1[i].y
    dz = Protein2[j].z - Protein1[i].z
    d = dx*dx + dy*dy + dz*dz
   #d = sqrt(d)
    return d

#Setting input and output files
input = file(sys.argv[1], "r")  
name = sys.argv[1].split(".")
output = file("distance.temp", "w")

#Reading the PDB file
Protein1=[]
Protein2=[]
lenght_prot1=[]
lenght_prot2=[]

for charact in input.readlines():
    F = charact.split()
    if len(charact) > 53:
       record = charact[0:6]
       if record == "ATOM  ":
          atom_num = eval(charact[6:11])
          atom = charact[13:16]
      #   alternate_character = charact[16]
          res_name = charact[17:20]
          chain = charact[21]
        #  chain = charact[72]
          res_num = eval(charact[22:26])
      #   code_for_insertio_residues = charact[26]
          x_coor = charact[30:38]
          y_coor = charact[38:46]
          z_coor = charact[46:54]
          occupancy = charact[54:59]
          if record == "ATOM  ":
            if chain == "A":
               Protein1.append(PDB1(atom,res_name,chain,res_num,eval(x_coor),eval(y_coor),eval(z_coor)))
               if atom == "CA ":
                  lenght_prot1.append(res_num)
            if chain == "B":
                Protein2.append(PDB2(atom,res_name,chain,res_num,eval(x_coor),eval(y_coor),eval(z_coor)))
                if atom == "CA ":
                   lenght_prot2.append(res_num)


#Calculation of distances and list od ICs
list_of_distances=[]

for i in range(0,len(Protein1)):
 for j in range(0,len(Protein2)):
      if dist(i,j) <= 30.25:
         list_of_distances.append('%-5s %5s %3s %5s %5s %3s    %5.2f' % (Protein1[i].Nres, Protein1[i].res, "A", Protein2[j].Nres, Protein2[j].res,"B", sqrt(dist(i,j))))


LEN= len(list_of_distances)

# Writing temporary output
for i in range(0,LEN):
   print >> output, list_of_distances[i]

