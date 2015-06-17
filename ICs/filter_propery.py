##########################################################################
##      count residue-residue shortest distances based on aa property   ##
##                        Author: Anna Vangone                          ##
##########################################################################

"""
Count iand classfied the number of 
inter-residue contacts between 
two proteins labelled as chain A and chain B.
This script is part of the ICs-parameters.sh program
"""

import string
import sys

#defining properties
property  = {  'HIS' : 'charged',
               'ARG' : 'charged',
               'LYS' : 'charged',
               'ASP' : 'charged',
               'GLU' : 'charged',
               'SER' : 'polar',
               'THR' : 'polar',
               'ASN' : 'polar',
               'GLN' : 'polar',
               'CYS' : 'apolar',
               'GLY' : 'apolar',
               'PRO' : 'apolar',
               'ALA' : 'apolar',
               'VAL' : 'apolar',
               'ILE' : 'apolar',
               'MET' : 'apolar',
               'LEU' : 'apolar',
               'PHE' : 'apolar',
               'TYR' : 'apolar',
               'TRP' : 'apolar'  }

#reading input
input = file("distance.temp", "r")

#setting dictionaries
dictio1={}
charged_charged={}
dictio2={}
charged_polar={}
dictio3={}
charged_apolar={}
dictio4={}
polar_polar={}
dictio5={}
polar_apolar={}
dictio6={}
apolar_apolar={}

#reading input file and filterg based on properties
for line in input.readlines():
   ar=line.split()

   if property[ar[1]] == "charged" and property[ar[4]] == "charged":
      key=ar[0]+ar[2]+ar[3]+ar[5]
      if dictio1.has_key(key):
         if eval(dictio1.get(key))>eval(ar[6]):
            dictio1[key]=ar[6]
            charged_charged[key]=line
      else:
           dictio1.setdefault(key,ar[6])
           charged_charged[key]=line

   if property[ar[1]] == "charged" and property[ar[4]] == "polar" or property[ar[1]] == "polar" and property[ar[4]] == "charged": 
      key=ar[0]+ar[2]+ar[3]+ar[5]
      if dictio2.has_key(key):
         if eval(dictio2.get(key))>eval(ar[6]):
            dictio2[key]=ar[6]
            charged_polar[key]=line
      else:
           dictio2.setdefault(key,ar[6])
           charged_polar[key]=line

   if property[ar[1]] == "charged" and property[ar[4]] == "apolar" or property[ar[1]] == "apolar" and property[ar[4]] == "charged":
      key=ar[0]+ar[2]+ar[3]+ar[5]
      if dictio3.has_key(key):
         if eval(dictio3.get(key))>eval(ar[6]):
            dictio3[key]=ar[6]
            charged_apolar[key]=line
      else:
           dictio3.setdefault(key,ar[6])
           charged_apolar[key]=line

   if property[ar[1]] == "polar" and property[ar[4]] == "polar":
      key=ar[0]+ar[2]+ar[3]+ar[5]
      if dictio4.has_key(key):
         if eval(dictio4.get(key))>eval(ar[6]):
            dictio4[key]=ar[6]
            polar_polar[key]=line
      else:
           dictio4.setdefault(key,ar[6])
           polar_polar[key]=line

   if property[ar[1]] == "polar" and property[ar[4]] == "apolar" or property[ar[1]] == "apolar" and property[ar[4]] == "polar":
      key=ar[0]+ar[2]+ar[3]+ar[5]
      if dictio5.has_key(key):
         if eval(dictio5.get(key))>eval(ar[6]):
            dictio5[key]=ar[6]
            polar_apolar[key]=line
      else:
           dictio5.setdefault(key,ar[6])
           polar_apolar[key]=line

   if property[ar[1]] == "apolar" and property[ar[4]] == "apolar":
      key=ar[0]+ar[2]+ar[3]+ar[5]
      if dictio6.has_key(key):
         if eval(dictio6.get(key))>eval(ar[6]):
            dictio6[key]=ar[6]
            apolar_apolar[key]=line
      else:
           dictio4.setdefault(key,ar[6])
           apolar_apolar[key]=line

#print output

print "ICs_charged-charged ", len(charged_charged)
#print len(charged_polar)
print "ICs_charged-apolar  ", len(charged_apolar)
print "ICs_polar-polar     ", len(polar_polar)
print "ICs_polar-apolar    ", len(polar_apolar)
#print len(apolar_apolar)

#for key in charged_charged:
#    print charged_charged[key],
