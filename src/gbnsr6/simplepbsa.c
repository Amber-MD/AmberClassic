#include "copyright_c.h"
// Description:
//	A program to test C interface for libpbsa
// Author:
// 	Mengjuei Hsieh, University of California Irvine
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prm.h"
#include "interface.h"
#include "gopt.h"

void private_getx_(int*,int*,const char*,REAL_T*);

int main(int argc, const char **argv){
    int    ipb;
    REAL_T e_pb,evdw,eelt,esurf,edisp;
    const char *prm_fpath, *crd_fpath, *tmpstring;
    int    pathlen;
    REAL_T *x,*f;
    PARMSTRUCT_T *prm;
    PBOPTSTRUCT_T *opt;

    ipb = 1;
    e_pb = 0.0;

    void *options= gopt_sort( & argc, argv, gopt_start(
	gopt_option( 900, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "ipb" )),
	gopt_option( 901, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "inp" )),
	gopt_option( 902, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "nfocus" )),
	gopt_option( 903, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "smoothopt" )),
	gopt_option( 904, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "fillratio" )),
	gopt_option( 905, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "istrng" )),
	gopt_option( 906, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "dprob" )),
	gopt_option( 907, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "dbfopt" )),
	gopt_option( 908, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "cutnb" )),
	gopt_option( 909, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "maxitn" )),
	gopt_option( 910, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "sprob" )),
	gopt_option( 911, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "frcopt" )),
	gopt_option( 912, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "cutsa" )),
	gopt_option( 914, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "offx" )),
	gopt_option( 915, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "offy" )),
	gopt_option( 916, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "offz" )),
	gopt_option( 917, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "accept" )),
	gopt_option( 918, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "radiopt" )),
	gopt_option( 919, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "bcopt" )),
	gopt_option( 920, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "eneopt" )),
	gopt_option( 921, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "solvopt" )),
	gopt_option( 922, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "fscale" )),
	gopt_option( 923, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "npbopt" )),
	gopt_option( 924, GOPT_ARG, gopt_shorts( 0 ), gopt_longs( "space" )),
        gopt_option( 'p', GOPT_ARG, gopt_shorts( 'p' ), gopt_longs( "prmtop" )),
        gopt_option( 'c', GOPT_ARG, gopt_shorts( 'c' ), gopt_longs( "inpcrd" )),
        gopt_option( 'h', 0, gopt_shorts( 'h', '?' ), gopt_longs( "help", "HELP" ))
	));
    if( gopt_arg( options, 'p', & prm_fpath ) < 1 || !strcmp( prm_fpath, "-" ) ){
       printf( "-p filename is required\n" );exit(1);
    }
    if( gopt_arg( options, 'c', & crd_fpath ) < 1 || !strcmp( crd_fpath, "-" ) ){
       printf( "-c filename is required\n" );exit(1);
    }
    if( gopt( options, 'h' ) ){
      /*
       * if any of the help options was specified
      */
      fprintf( stdout, "help text\n" );
      exit( EXIT_SUCCESS );
    }

    if ((nabout = fopen("/dev/null","w"))==NULL){
        printf("having problem opening nabout\n");
        exit (1);
    }

    prm = (PARMSTRUCT_T *) malloc(sizeof(PARMSTRUCT_T));
    prm = rdparm(prm_fpath); fclose(nabout);

    //to import x from inpcrd
    x = (REAL_T *) malloc(sizeof(REAL_T)*3*prm->Natom);
    pathlen=strlen(crd_fpath);
    private_getx_(&prm->Natom,&pathlen,crd_fpath,x);

    pboptinit(opt);
    opt->ipb=ipb;
    if( gopt_arg( options, 900, &tmpstring ) > 0 )
       opt->ipb=atoi(tmpstring);
    if( gopt_arg( options, 901, &tmpstring ) > 0 )
       opt->inp=atoi(tmpstring);
    if( gopt_arg( options, 902, &tmpstring ) > 0 )
       opt->nfocus=atoi(tmpstring);
    if( gopt_arg( options, 903, &tmpstring ) > 0 )
       opt->smoothopt=atoi(tmpstring);
    if( gopt_arg( options, 904, &tmpstring ) > 0 )
       opt->fillratio=atof(tmpstring);
    if( gopt_arg( options, 905, &tmpstring ) > 0 )
       opt->istrng=atof(tmpstring);
    if( gopt_arg( options, 906, &tmpstring ) > 0 )
       opt->dprob=atof(tmpstring);
    if( gopt_arg( options, 907, &tmpstring ) > 0 )
       opt->dbfopt=atoi(tmpstring);
    if( gopt_arg( options, 908, &tmpstring ) > 0 )
       opt->cutnb=atof(tmpstring);
    if( gopt_arg( options, 909, &tmpstring ) > 0 )
       opt->maxitn=atoi(tmpstring);
    if( gopt_arg( options, 910, &tmpstring ) > 0 )
       opt->sprob=atof(tmpstring);
    if( gopt_arg( options, 911, &tmpstring ) > 0 )
       opt->frcopt=atoi(tmpstring);
    if( gopt_arg( options, 912, &tmpstring ) > 0 )
       opt->cutsa=atof(tmpstring);
    if( gopt_arg( options, 914, &tmpstring ) > 0 )
       opt->offx=atof(tmpstring);
    if( gopt_arg( options, 915, &tmpstring ) > 0 )
       opt->offy=atof(tmpstring);
    if( gopt_arg( options, 916, &tmpstring ) > 0 )
       opt->offz=atof(tmpstring);
    if( gopt_arg( options, 917, &tmpstring ) > 0 )
       opt->accept=atof(tmpstring);
    if( gopt_arg( options, 918, &tmpstring ) > 0 )
       opt->radiopt=atoi(tmpstring);
    if( gopt_arg( options, 919, &tmpstring ) > 0 )
       opt->bcopt=atoi(tmpstring);
    if( gopt_arg( options, 920, &tmpstring ) > 0 )
       opt->eneopt=atoi(tmpstring);
    if( gopt_arg( options, 921, &tmpstring ) > 0 )
       opt->solvopt=atoi(tmpstring);
    if( gopt_arg( options, 922, &tmpstring ) > 0 )
       opt->fscale=atoi(tmpstring);
    if( gopt_arg( options, 923, &tmpstring ) > 0 )
       opt->npbopt=atoi(tmpstring);
    if( gopt_arg( options, 924, &tmpstring ) > 0 )
       opt->space=atof(tmpstring);
    
    f = (REAL_T *) malloc(sizeof(REAL_T)*3*prm->Natom);

    e_pb = epbsa(opt,prm,x,
		 f,&evdw,&eelt,&esurf,
		 &edisp, 0, 0);
    printf("VDWAALS = %12.4f\n",evdw);
    printf("EPB     = %12.4f\n",e_pb);
    printf("EELEC   = %12.4f\n",eelt);
    printf("ECAVITY = %12.4f\n",esurf);
    printf("EDISPER = %12.4f\n",edisp);

    free(f);
    free(opt);
    free(x);
    free(prm);
    return (0);
}
