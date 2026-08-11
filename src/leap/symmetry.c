#include "basics.h"
#include "symmetry.h"

/* ------------------------------------------------------------------ */
/*  parse_symop of the form Y+1/2,X-1/2,Z  Case insensitive           */
/*  matrices are integer form with translation in 1/12 units          */
/* ------------------------------------------------------------------ */
int parse_symop(const char *s, int rot[3][3], int trans[3])
{
    memset(rot,   0, 9 * sizeof(int));
    memset(trans, 0, 3 * sizeof(int));
    int row = 0, sign = +1;

    while (*s && row < 3) {
        while (isspace((unsigned char)*s)) s++;
        if (!*s) break;

        while (*s && *s != ',') {
            while (isspace((unsigned char)*s)) s++;

            if      (*s == '+') { sign = +1; s++; }
            else if (*s == '-') { sign = -1; s++; }

            while (isspace((unsigned char)*s)) s++;
            if (!*s || *s == ',') break;

            if (*s == 'x' || *s == 'X')      { rot[row][0] = sign; s++; sign=+1; }
            else if (*s == 'y' || *s == 'Y') { rot[row][1] = sign; s++; sign=+1; }
            else if (*s == 'z' || *s == 'Z') { rot[row][2] = sign; s++; sign=+1; }
            else if (isdigit((unsigned char)*s)) {
                int p = 0;
                while (isdigit((unsigned char)*s)) p = p*10 + (*s++ - '0');
                int q = 1;
                if (*s == '/') {
                    s++; q = 0;
                    if (!isdigit((unsigned char)*s)) return -1;
                    while (isdigit((unsigned char)*s)) q = q*10 + (*s++ - '0');
                    if (q == 0 || 12 % q != 0) return -1;
                }
                trans[row] += sign * p * (12 / q);
                sign = +1;
            } else return -1;
        }
        if (*s == ',') s++;
        row++;
    }
    return (row == 3) ? 0 : -1;
}

/* ------------------------------------------------------------------ */
/*  Helpers                                                            */
/* ------------------------------------------------------------------ */

/* Strip leading/trailing whitespace in-place, return pointer */
static char *trim(char *s)
{
    while (isspace((unsigned char)*s)) s++;
    char *e = s + strlen(s);
    while (e > s && isspace((unsigned char)e[-1])) *--e = '\0';
    return s;
}

/* Extract the first single-quoted token from src into dst (max dstlen).
 * Returns 1 on success, 0 if no quoted token found. */
static int extract_quoted(const char *src, char *dst, int dstlen)
{
    const char *p = strchr(src, '\'');
    if (!p) return 0;
    p++;
    const char *q = strchr(p, '\'');
    if (!q) return 0;
    int len = (int)(q - p);
    if (len >= dstlen) len = dstlen - 1;
    strncpy(dst, p, len);
    dst[len] = '\0';
    /* trim internal leading/trailing spaces */
    char *t = trim(dst);
    if (t != dst) memmove(dst, t, strlen(t)+1);
    return 1;
}

/* ------------------------------------------------------------------ */
/*  Main parser                                                        */
/* ------------------------------------------------------------------ */

/*
 * Parse CCP4 syminfo.lib for a space group matching number (if > 0)
 * or name (matched against xHM and old symbol fields).
 * Fills *sg on success.  Returns 0 on success, -1 on error/not found.
 */
int parse_spacegroup_file(int number, const char *name, SPACEGROUPt *sg)
{
    STRING sFilename;
    char *CLIBD = getenv("CLIBD");
    char *PHENIX = getenv("PHENIX");
    if (CLIBD) {
        strcpy(sFilename,CLIBD);
        strcat(sFilename,"/syminfo.lib");
    } else if (PHENIX) {
        strcpy(sFilename,PHENIX);
        strcat(sFilename,"modules/ccp4io/libccp4/data/syminfo.lib");
    } else
        strcpy(sFilename,"syminfo.lib");

    FILE *f = FOPENCOMPLAIN(sFilename,"r");
    if (!f) { return -1; }

    memset(sg, 0, sizeof(*sg));

    char line[256];
    int  in_block = 0;
    SYMOPt cenops[MAX_SYMOPS];
    int   n_cenops = 0;

    while (fgets(line, sizeof(line), f)) {
        char *t = trim(line);
        if (!*t || *t == '#') continue;

        /* --- block boundaries --- */
        if (strncmp(t, "begin_spacegroup", 16) == 0) {
            in_block = 1;
            n_cenops = 0;
            memset(sg, 0, sizeof(*sg));
            continue;
        }
        if (strncmp(t, "end_spacegroup", 14) == 0) {
            if ( (number > 0 && sg->number == number && sg->old[0]) ||
                 (name && *name &&
                     (strcmp(sg->xHM, name) == 0 ||
                      strcmp(sg->old, name) == 0)) ) {
                /* expand symops by centering translations before returning */
                if (n_cenops > 1) {
                    int base = sg->n_symops;
                    for (int c = 1; c < n_cenops; c++) {
                        for (int s = 0; s < base; s++) {
                            if (sg->n_symops >= MAX_SYMOPS) {
                                VPFATAL("Warning: MAX_SYMOPS reached, truncating\n");
                                fclose(f);
                                return -1;
                            }
                            SYMOPt *dst = &sg->symops[sg->n_symops++];
                            *dst = sg->symops[s];
                            for (int i = 0; i < 3; i++)
                                dst->trans[i] = (sg->symops[s].trans[i]
                                               + cenops[c].trans[i]) % 12;
                        }
                    }
                }
                fclose(f);
                return 0;
            }
            in_block = 0;
        }
        if (!in_block) continue;

        /* --- keywords --- */
        if (strncmp(t, "number", 6) == 0) {
            sg->number = atoi(t + 6);
        } else if (strncmp(t, "symbol xHM", 10) == 0) {
            extract_quoted(t + 11, sg->xHM, sizeof(sg->xHM));
        } else if (strncmp(t, "symbol old", 10) == 0) {
            extract_quoted(t + 11, sg->old, sizeof(sg->old));
        } else if (strncmp(t, "symbol pgrp", 10) == 0) {
            extract_quoted(t + 12, sg->pgrp, sizeof(sg->pgrp));
            extract_quoted(t + 14 + strlen(sg->pgrp), sg->pgrp, sizeof(sg->pgrp));
        } else if (strncmp(t, "symop", 5) == 0) {
            char *ops = trim(t + 5);
            if (sg->n_symops < MAX_SYMOPS) {
                SYMOPt *op = &sg->symops[sg->n_symops];
                if (parse_symop(ops, op->rot, op->trans) == 0)
                    sg->n_symops++;
                else
                    VPWARN("failed to parse symop '%s'\n", ops);
            }
        } else if (strncmp(t, "cenop", 5) == 0) {
            char *ops = trim(t + 5);
            if (n_cenops < MAX_SYMOPS) {
                SYMOPt *op = &cenops[n_cenops];
                if (parse_symop(ops, op->rot, op->trans) == 0)
                    n_cenops++;
                else
                    VPWARN("failed to parse cenop '%s'\n", ops);
            }
        }
        /* all other keywords silently ignored */
    }

    /* Fell off end of file without matching end_spacegroup */
    fclose(f);
    return -1;
}

/* -----------------------------------------------------------------------
 * BuildSymopMatrices
 * Converts an array of SYMOPt (fractional integer rot + trans/12)
 * to Cartesian 4x4 matrices using:  Mcart = M * Mfrac * Mi
 * The translation (trans/12 fractional) is embedded in column 3 of
 * Mfrac so it is carried through automatically by the multiply.
 * --------------------------------------------------------------------- */
void
BuildSymopMatrices( UNIT uUnit,
                     SYMOPt *symmops, int nSymops,
                     MATRIX **maPSymops )
{
MATRIX  M, Mi, mFrac, mTmp;
MATRIX  *maSymops;
int     i, j;

    BuildFractionalTransforms( uUnit, M, Mi );

    maSymops = (MATRIX *)malloc( nSymops * sizeof(MATRIX) );
    *maPSymops = maSymops;

    for ( i = 0; i < nSymops; i++ ) {
        /* --- build fractional symop as 4x4 ---
         * rotation in upper-left 3x3, translation in column 3 */
        for ( j = 0; j < 4; j++ )
            mFrac[j][0] = mFrac[j][1] = mFrac[j][2] = mFrac[j][3] = 0.0;
        mFrac[3][3] = 1.0;

        int r, c;
        for ( r = 0; r < 3; r++ )
            for ( c = 0; c < 3; c++ )
                mFrac[c][r] = (double)symmops[i].rot[r][c];  /* col-major */

        mFrac[3][0] = (double)symmops[i].trans[0] / 12.0;
        mFrac[3][1] = (double)symmops[i].trans[1] / 12.0;
        mFrac[3][2] = (double)symmops[i].trans[2] / 12.0;

        /* Mcart = M * Mfrac * Mi */
        MatrixMultiply( mTmp,           M,    mFrac );
        MatrixMultiply( maSymops[i],    mTmp, Mi    );
    }
}

/* -----------------------------------------------------------------------
 * BuildFractionalTransforms
 * Builds frac->cart (M) and cart->frac (Mi) as 4x4 matrices.
 * bOrtho is auto-detected from cell angles.
 * --------------------------------------------------------------------- */
void
BuildFractionalTransforms( UNIT uUnit, MATRIX M, MATRIX Mi )
{
double  a, b, c, alpha, beta, gamma;
int     i, j, bOrtho;

    a     = uUnit->dXWidth;
    b     = uUnit->dYWidth;
    c     = uUnit->dZWidth;
    alpha = uUnit->dAlpha;
    beta  = uUnit->dBeta;
    gamma = uUnit->dGamma;

    bOrtho = ( fabs(alpha - 90.0*DEGTORAD) < 1e-6 &&
               fabs(beta  - 90.0*DEGTORAD) < 1e-6 &&
               fabs(gamma - 90.0*DEGTORAD) < 1e-6 );

    /* zero everything first */
    for ( i = 0; i < 4; i++ )
        for ( j = 0; j < 4; j++ )
            M[i][j] = Mi[i][j] = 0.0;
    M[3][3] = Mi[3][3] = 1.0;

    if ( bOrtho ) {
        M[0][0]=a;    M[1][1]=b;    M[2][2]=c;
        Mi[0][0]=1.0/a; Mi[1][1]=1.0/b; Mi[2][2]=1.0/c;
    } else {
        double ar = alpha;
        double br = beta;
        double gr = gamma;
        double ca = cos(ar), cb = cos(br), cg = cos(gr), sg = sin(gr);

        double cx = cb;
        double cy = (ca - cb*cg) / sg;
        double cz = sqrt(1.0 - cx*cx - cy*cy);

        /* frac->cart, column-major: M[col][row] */
        M[0][0]=a;
        M[1][0]=b*cg;  M[1][1]=b*sg;
        M[2][0]=c*cx;  M[2][1]=c*cy;  M[2][2]=c*cz;

        /* analytic inverse of upper-triangular 3x3, stored column-major */
        double m00=M[0][0];
        double m10=M[1][0], m11=M[1][1];
        double m20=M[2][0], m21=M[2][1], m22=M[2][2];

        Mi[0][0] = 1.0/m00;
        Mi[1][0] = -m10/(m00*m11);
        Mi[1][1] = 1.0/m11;
        Mi[2][0] = (m10*m21 - m20*m11)/(m00*m11*m22);
        Mi[2][1] = -m21/(m11*m22);
        Mi[2][2] = 1.0/m22;
    }
}

