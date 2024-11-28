metatwist=$AMBERHOME/bin/metatwist 
# (1) get the laplacian map
# 
$metatwist --dx ../0_crown-ether-3drism/g.K+.1.dx.bz2 --odx lp-K+.dx --sigma 1.0 --convolve 4 \
	   --species K+ K+ --bulk 0.2 > get-laplacian.out


# (2) place ions
#  
$metatwist --dx ../0_crown-ether-3drism/g.K+.1.dx.bz2 --ldx convolution-lp-K+.dx  --map blobsper --thresh 0.1 \
	   --bulkdens 0.2 --species K+ K+ > get-blobs.out



