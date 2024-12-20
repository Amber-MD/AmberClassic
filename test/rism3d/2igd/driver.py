#!/usr/bin/env python3

import sander
from parmed.amber.readparm import AmberParm, Rst7

# Initialize the topology object with coordinates
parm = AmberParm("2igd.parm7")
rst = Rst7.open("2igd.rst7")

# Set up input options 
inp = sander.pme_input()
inp.irism = 1
inp.grdspc1 = 0.9
# inp.rism_verbose = 1
# inp.mdiis_del = 0.8

sander.setup(parm, rst.coordinates, rst.box, inp)

# Compute the energies and forces
ene, frc = sander.energy_forces()

# Do whatever you want with the energies and forces
print( "Energy    = %11.4f" % ene.tot )
print( "  vdw     = %11.4f" % ene.vdw )
print( "  elec    = %11.4f" % ene.elec )
print( "  vdw_14  = %11.4f" % ene.vdw_14 )
print( "  elec_14 = %11.4f" % ene.elec_14 )
print( "  bond    = %11.4f" % ene.bond )
print( "  angle   = %11.4f" % ene.angle )
print( "  dihed   = %11.4f" % ene.dihedral )
print( "  rism    = %11.4f" % ene.rism )

# Free up our memory
sander.cleanup()
