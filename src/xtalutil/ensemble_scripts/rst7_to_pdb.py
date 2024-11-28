#!/usr/bin/python

import sys

rst7_file = sys.argv[1]
template_pdb = sys.argv[2]
outpdb = sys.argv[3]

'''
Combine information in {rst7, formatted} with that in the template_pdb
   file (usually 4phenix_xxxx.pdb) to create outpdb 

Note:  the atom order in outpdb will be the "Amber" atom order; use phenix
routines after this to convert to phenix atom order.

Further note: the template_pdb has to have exactly the same order of
atoms as does the rst7 file.  This should/will be the case when
les_builder is creating the template_pdb file, but is potentially
fragile: this routine just inserts the coordinates from the rst7
file into consecutive atoms in the template pdb file.
'''

rst7h = open(rst7_file,"r")
line = rst7h.readline()   # title
line = rst7h.readline()   # natom
nat = int(line[0:6])
x = []
y = []
z = []
for iat in range(0,nat,2):
   line = rst7h.readline()
   x.append(line[ 0: 8])
   y.append(line[12:20])
   z.append(line[24:32])
   x.append(line[36:44])
   y.append(line[48:56])
   z.append(line[60:68])
if nat%2 == 1:    # grab the last coordinate:
   line = rst7h.readline()
   x.append(line[ 0: 8])
   y.append(line[12:20])
   z.append(line[24:32])

rst7h.close()

# transfer minimized coordinates to the 4phenix_xxxx.min.pdb file:

iat = 0
th = open(template_pdb,"r")
outh = open(outpdb, "w" )
for line in th:
   if line[0:4] == "ATOM" or line[0:6] == "HETATM":
      outh.write(line[0:30] + x[iat] + y[iat] + z[iat] + line[54:])
      iat = iat + 1
   else:
      outh.write(line)

th.close()
outh.close()

