#  Script like UnitCell: expands an asu to the full unit cell.
#    Unlike UnitCell, does not require the REAMRK 290 matrices to be present
#
#   Usage:  phenix.python expand_to_p1.py  <asu.pdb>  <uc.pdb>
#
import sys
import iotbx.pdb
import cctbx

if len(sys.argv) != 2:
   print "Usage: phenix.python expand_to_p1.py  <asu.pdb>  <uc.pdb>"
   sys.exit(1)

pdb_inp = iotbx.pdb.input(sys.argv[1])
cs = pdb_inp.crystal_symmetry()
ph = pdb_inp.construct_hierarchy()
ph_p1 = ph.expand_to_p1(crystal_symmetry=cs)
abc = cs.unit_cell().parameters()[:6]
cs_p1 = cctbx.crystal.symmetry(abc, "P 1")
ph_p1.write_pdb_file(sys.argv[2], crystal_symmetry=cs_p1)
