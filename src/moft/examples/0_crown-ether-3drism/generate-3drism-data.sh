
# run rism3d

$AMBERHOME/bin/rism3d.snglpnt --pdb bc5-k.pdb --prmtop bc5-k.parm7 --xvv KCl-aq-0.2M-pse3.xvv --closure pse1 pse2 pse3 --tolerance 1e-03 1e-06 --solvcut 999999  --ng 192,192,192 --solvbox 96,96,96 --solvcut 999999 --buffer -1  --mdiis_del 0.5 --mdiis_nvec 10 --mdiis_restart 20 --verbose 2 --npropagate 0 --maxstep 1000 --guv g  > rism.out

# compress volumetric files

bzip2 -f *dx
