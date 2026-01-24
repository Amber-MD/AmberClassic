
//  simple driver script to test:

#include "nabc.h"

FILE* nabout;

int main(){

    MOLECULE_T  *m;
    char  *seq =  "ggcaattcgc";
    char  *cseq = "ccgttaagcg";

    m = wc_helix( seq, "dna", cseq, "dna", 2.25, -4.96, 36.0, 3.38, "s5a5s3a3");
    putpdb( "wch.pdb", m, NULL );
}
