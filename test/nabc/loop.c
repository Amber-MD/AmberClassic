#include    "nabc.h"

FILE *nabout;    /* can this go into nabc.h?  */

int main( int arg, char *argv[] )
{

    MOLECULE_T *m;
    STRAND_T *sp;    
    RESIDUE_T *rp;
    ATOM_T *ap;

    m = getpdb("d01.pdb", NULL);
    upd_molnumbers( m );

#if 0
    //  strands are stored in a simple linked-list:
    for (sp = m->m_strands; sp; sp = sp->s_next) {

        printf( "strand %d:  %s\n", sp->s_strandnum, sp->s_strandname );

        // but residues are stored in a residue array inside each strand:
        for (int rn = 0; rn < sp->s_nresidues; rn++) {
            rp = sp->s_residues[rn];

            printf( "    residue %d: %s\n", rp->r_resnum, rp->r_resname );

            // atoms are stored in an atom array inside each residue:
            for (int an = 0; an < rp->r_natoms; an++) {
                ap = &rp->r_atoms[an];

                NAB_arc( ap, "fullname" );  // needed to update a_fullname
                // printf( "        atom %d: %s\n", an, ap->a_atomname );

            }                   /* end loop over atoms in this residue  */

        }                       /* end loop over residues in this strand  */

    }                           /* end loop over strands  */

    // second version of the above, like what the nab compiler would produce:

#define ForRinM(r,m) for( r = NULL; r = NAB_rnext( m, r ); )
#define ForAinR(a,r) for( a = NULL; a = NAB_anext( r, a ); )

    // for( rp = NULL; rp = NAB_rnext( m, rp ); ){
    ForRinM( rp, m ){

        // chimera synatx for residue names:
        printf( "residue:  %s.%s\n", rp->r_resname, rp->r_strand->s_strandname );

        //for( ap = NULL; ap = NAB_anext( rp, ap ); ){
        ForAinR( ap, rp ){
            if( NAB_aematch( ap, "::*P*" ) )
                printf( "    atom:  %s\n",  ap->a_atomname );
        }
    }
#endif

    // third version, just look over all atoms:

#define ForAinM(a,m) for( a = NULL; a = NAB_mnext( m, a ); )
#define A_fullname(a)  NAB_arc( a, "fullname" )

    ForAinM( ap, m ){
        if( NAB_aematch( ap, "::*P*" ) ){
            // NAB_arc( ap, "fullname" );  // needed to update a_fullname
            // printf( "    atom:  %s\n",  NAB_arc( ap, "fullname" ) );
            printf( "    atom:  %s\n",  A_fullname( ap ) );
        }
    }

}

