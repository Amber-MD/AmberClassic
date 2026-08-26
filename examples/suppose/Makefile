include ../../src/config.h

suppose$(SFX): suppose.nab
	nab -o suppose$(SFX) suppose.nab

PCA_suppose$(SFX): PCA_suppose.nab housev.c qrevec.c eigen.c trnm.c
	nab -o PCA_suppose$(SFX) PCA_suppose.nab housev.c qrevec.c eigen.c trnm.c

clean:
	-/bin/rm suppose$(SFX) PCA_suppose$(SFX) suppose.c PCA_suppose.c *.o
