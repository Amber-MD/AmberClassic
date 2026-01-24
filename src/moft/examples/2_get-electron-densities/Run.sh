metatwist=$AMBERHOME/bin/metatwist


$metatwist --dx ../0_crown-ether-3drism/g.K+.1.dx.bz2 --species K+ --odx rho.K+.1.dx --map rhoel --bulkdens 0.2
$metatwist --dx ../0_crown-ether-3drism/g.Cl-.1.dx.bz2 --species Cl- --odx rho.Cl-.1.dx --map rhoel --bulkdens 0.2

$metatwist --dx ../0_crown-ether-3drism/g.O.1.dx.bz2 --species O2- --odx rho.O.1.dx --map rhoel --bulkdens 55.55

# assembly of all densities; for now metatwist will average over multiple input densities.
$metatwist --dx rho.Cl-.1.dx   rho.K+.1.dx    rho.O.1.dx --odx rho.dx --species none

bzip2 --force *dx
