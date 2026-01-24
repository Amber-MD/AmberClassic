
int	belly_parm( int natom, int* frozen )
{
	int	i, j, k, l, ib, nfrozen, nb, ibig, ismall, ka, la;
	int *iptmp, *dihtmp;

	nfrozen = 0;
	for( i=0; i<natom; ++i ){
        if( frozen[i] ) nfrozen++;
    }
    if( nfrozen == 0 ) return( nfrozen );

/*  remove all of the internal coordinates for frozen-frozen pairs:  */

    /*  remove frozen-frozen bonds:  */
    for( nb=0, ib=0; ib<prm->Nbonh; ib++ ){
        i = prm->BondHAt1[ib]/3;
        j = prm->BondHAt2[ib]/3;
        if( !frozen[i] || !frozen[j] ){
            prm->BondHAt1[nb] = prm->BondHAt1[ib];
            prm->BondHAt2[nb] = prm->BondHAt2[ib];
            prm->BondHNum[nb] = prm->BondHNum[ib];
            nb++;
        }
    }
    prm->Nbonh = nb;
    for( nb=0, ib=0; ib<prm->Nbona; ib++ ){
        i = prm->BondAt1[ib]/3;
        j = prm->BondAt2[ib]/3;
        if( !frozen[i] || !frozen[j] ){
            prm->BondAt1[nb] = prm->BondAt1[ib];
            prm->BondAt2[nb] = prm->BondAt2[ib];
            prm->BondNum[nb] = prm->BondNum[ib];
            nb++;
        }
    }
    prm->Nbona = nb;
    prm->Mbona = nb;

    /*  remove frozen-frozen angles:  */
    for( nb=0, ib=0; ib<prm->Ntheth; ib++ ){
        i = prm->AngleHAt1[ib]/3;
        j = prm->AngleHAt2[ib]/3;
        k = prm->AngleHAt3[ib]/3;
        if( !frozen[i] || !frozen[j] || !frozen[k] ){
            prm->AngleHAt1[nb] = prm->AngleHAt1[ib];
            prm->AngleHAt2[nb] = prm->AngleHAt2[ib];
            prm->AngleHAt3[nb] = prm->AngleHAt3[ib];
            prm->AngleHNum[nb] = prm->AngleHNum[ib];
            nb++;
        }
    }
    prm->Ntheth = nb;
    for( nb=0, ib=0; ib<prm->Ntheta; ib++ ){
        i = prm->AngleAt1[ib]/3;
        j = prm->AngleAt2[ib]/3;
        k = prm->AngleAt3[ib]/3;
        if( !frozen[i] || !frozen[j] || !frozen[k] ){
            prm->AngleAt1[nb] = prm->AngleAt1[ib];
            prm->AngleAt2[nb] = prm->AngleAt2[ib];
            prm->AngleAt3[nb] = prm->AngleAt3[ib];
            prm->AngleNum[nb] = prm->AngleNum[ib];
            nb++;
        }
    }
    prm->Ntheta = nb;

    /*  remove frozen dihedrals and 1-4's:  */
#define ABS(x) x > 0 ? x : -x 
    for( i=0; i<prm->Natom; i++ ) prm->N14pairs[i] = 0;
    iptmp = ivector( 0, 12*prm->Natom );
    dihtmp = ivector( 0, 12*prm->Natom );
    for( nb=0, ib=0; ib<prm->Nphih; ib++ ){
        i = prm->DihHAt1[ib]/3;
        j = prm->DihHAt2[ib]/3;
        k = prm->DihHAt3[ib]/3;
        l = prm->DihHAt4[ib]/3;
        ka = ABS( k );
        la = ABS( l );
        if( !frozen[i] || !frozen[j] || !frozen[ka] || !frozen[la] ){
            prm->DihHAt1[nb] = prm->DihHAt1[ib];
            prm->DihHAt2[nb] = prm->DihHAt2[ib];
            prm->DihHAt3[nb] = prm->DihHAt3[ib];
            prm->DihHAt4[nb] = prm->DihHAt4[ib];
            prm->DihHNum[nb] = prm->DihHNum[ib];
            nb++;
            if( k >= 0 && l >= 0 ){
                ismall = i < l ? i : l;
                ibig   = i > l ? i : l;
                iptmp[12*ismall + prm->N14pairs[ismall]] = ibig;
                dihtmp[12*ismall + prm->N14pairs[ismall]++] = prm->DihHNum[ib] - 1;
            }
        }
    }
    prm->Nphih = nb;
    for( nb=0, ib=0; ib<prm->Nphia; ib++ ){
        i = prm->DihAt1[ib]/3;
        j = prm->DihAt2[ib]/3;
        k = prm->DihAt3[ib]/3;
        l = prm->DihAt4[ib]/3;
        ka = ABS( k );
        la = ABS( l );
        if( !frozen[i] || !frozen[j] || !frozen[ka] || !frozen[la] ){
            prm->DihAt1[nb] = prm->DihAt1[ib];
            prm->DihAt2[nb] = prm->DihAt2[ib];
            prm->DihAt3[nb] = prm->DihAt3[ib];
            prm->DihAt4[nb] = prm->DihAt4[ib];
            prm->DihNum[nb] = prm->DihNum[ib];
            nb++;
            if( k >= 0 && l >= 0 ){
                ismall = i < l ? i : l;
                ibig   = i > l ? i : l;
                iptmp[12*ismall + prm->N14pairs[ismall]] = ibig;
                dihtmp[12*ismall + prm->N14pairs[ismall]++] = prm->DihNum[ib] - 1;
            }
        }
    }
    prm->Nphia = nb;
    prm->Mphia = nb;
    j = 0;
    for( i=0; i<prm->Natom - 1; i++ ){
        for( k=0; k<prm->N14pairs[i]; k++ ) {
             prm->N14pairlist[j] = iptmp[12*i + k];
             prm->N14sceelist[j] = prm->Scee[dihtmp[12 * i + k]];
             prm->N14scnblist[j++] = prm->Scnb[dihtmp[12 * i + k]];
        }
    }
    free_ivector( iptmp, 0, 12*prm->Natom );
    free_ivector( dihtmp, 0, 12*prm->Natom );

	return( nfrozen );
}
