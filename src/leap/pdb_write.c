/*
 *      subroutine for writing PDB format files
 */

/*      Modifications induced by the implementation of the savemol2 command
 *      Christine Cezard (2007) 
 *      Universite de Picardie - Jules Verne, Amiens
 *      http://q4md-forcefieldtools.org
 *
 *      declaration of a case MOL2_ATOM line 372
 */    
  
 
/* LINTLIBRARY */

# include       <stdio.h>
# include       <ctype.h>
# include       <assert.h>
# include       "pdb_int.h"
# include       "unit.h"
# include       "defaults.h"

/*
 *      For each pdb record type, there is a format for reading the
 *      record values and for printing them.
 *
 *      The actual format of a line written, is the print format
 *      followed by blank padding to 72 characters, followed by
 *      8 characters of file and line information.
 */

void
pdb_write_record(FILE *f, pdb_record *r)
{
        char            buffer[PDB_BUFSIZ];
        char            *fmt;
        register struct pdb_sheet       *sh;
        register pdb_residue    *shr0, *shr1, *sha0, *sha1;

        /* convert C structure to pdb record */
        assert(r->record_type <= PDB_NUM_R + 1 );
        fmt = pdb_record_format[r->record_type].print_format;
        switch (r->record_type) {

        
		
		case PDB_UNKNOWN:
                pdb_sprintf(buffer, fmt, r->pdb.unknown.junk);
                break;

        case PDB_ANISOU:
        case PDB_SIGUIJ:
                pdb_sprintf(buffer, fmt, r->pdb.anisou.serial_num,
                        r->pdb.anisou.name, r->pdb.anisou.alt_loc,
                        r->pdb.anisou.residue.name,
                        r->pdb.anisou.residue.chain_id,
                        r->pdb.anisou.residue.seq_num,
                        r->pdb.anisou.residue.insert_code,
                        r->pdb.anisou.u[0], r->pdb.anisou.u[1],
                        r->pdb.anisou.u[2], r->pdb.anisou.u[3],
                        r->pdb.anisou.u[4], r->pdb.anisou.u[5]);
                break;

        case PDB_ATOM:
        case PDB_HETATM:
        case PDB_SIGATM:
/*
"ATOM  %5d %-4s%C%3s %C%4d%C   %8.3f%8.3f%8.3f%6.2f%6.2f %3D"
ATOM      1  N   GLY   204      16.233  16.706  11.826  1.00  0.00
ATOM      6  HA2 GLY   204      15.123  14.988  11.463  1.00  0.00

ATOM  150015 Na+  Na+  1443     -17.714  -6.411   7.833  1.00  0.00

"ATOM  %6d%-4s%C%3s %C%4d%C   %8.3f%8.3f%8.3f%6.2f%6.2f %3D";
ATOM  150009Na+  Na+  1439     -23.630   8.766 -18.151  1.00  0.00

ATOM  41150  H1  WAT  8081      21.760   0.196  43.791  1.00  0.00
ATOM  51150  H1  WAT  10581      21.703 -29.969  47.175  1.00  0.00

hybrid36 and chainId2 fields:
ATOM  A3AE9  CB  ALAAk 300     300.704 253.529 174.255  1.00 91.98           C
ATOM  A3AEA  N   TYRBP -45     351.130 172.789 312.980  1.00 43.62           N

*/
        {
                if (!GDefaults.bPdbHybrid36) {
                        int atom_overflow = (r->pdb.atom.serial_num > 99999);
                        int residue_overflow = (r->pdb.atom.residue.seq_num > 9999);
                        if (atom_overflow  &&  !residue_overflow) {
                                fmt = "ATOM  %6d%-4s%C%3s%.2s%4d%C   %8.3f%8.3f%8.3f%6.2f%6.2f         %2s%2s";
                        } else if (!atom_overflow  &&  residue_overflow) {
                                fmt = "ATOM  %5d %-4s%C%3s%.2s%5d%C  %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s";
                        } else if (atom_overflow  &&  residue_overflow) {
                                fmt = "ATOM  %6d%-4s%C%3s%.2s%5d%C  %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s";
                        } else {
                                fmt = "ATOM  %5d %-4s%C%3s%.2s%4d%C   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s";
                        }
                }
                pdb_sprintf(buffer, fmt, r->pdb.atom.serial_num,
                        r->pdb.atom.name, r->pdb.atom.alt_loc,
                        r->pdb.atom.residue.name,
                        r->pdb.atom.residue.chain_id,
                        r->pdb.atom.residue.seq_num,
                        r->pdb.atom.residue.insert_code,
                        r->pdb.atom.x, r->pdb.atom.y, r->pdb.atom.z,
                        r->pdb.atom.occupancy, r->pdb.atom.temp_factor,
                        r->pdb.atom.element,r->pdb.atom.fcharge);
                break;
        }
        case PDB_AUTHOR:
        case PDB_COMPND:
        case PDB_JRNL:
        case PDB_SOURCE:
                pdb_sprintf(buffer, fmt, r->pdb.author.continuation,
                        r->pdb.author.data);
                break;

        case PDB_CONECT:
                pdb_sprintf(buffer, fmt, r->pdb.conect.serial_num,
                        r->pdb.conect.covalent[0], r->pdb.conect.covalent[1], 
						r->pdb.conect.covalent[2]);
                break;

        case PDB_CRYST1:
                pdb_sprintf(buffer, fmt, r->pdb.cryst1.a, r->pdb.cryst1.b,
                        r->pdb.cryst1.c, r->pdb.cryst1.alpha,
                        r->pdb.cryst1.beta, r->pdb.cryst1.gamma,
                        &(r->pdb.cryst1.space_grp), r->pdb.cryst1.z);
                break;

        case PDB_END:
                pdb_sprintf(buffer, fmt);
                break;

        case PDB_FORMUL:
                pdb_sprintf(buffer, fmt, r->pdb.formul.component,
                        r->pdb.formul.het_id, r->pdb.formul.continuation,
                        r->pdb.formul.exclude, r->pdb.formul.formula);
                break;

        case PDB_FTNOTE:
                pdb_sprintf(buffer, fmt, r->pdb.ftnote.num, r->pdb.ftnote.text);
                break;

        case PDB_HEADER:
                pdb_sprintf(buffer, fmt, r->pdb.header.class,
                        r->pdb.header.date, r->pdb.header.id);
                break;

        case PDB_HELIX:
                pdb_sprintf(buffer, fmt, r->pdb.helix.serial_num,
                        r->pdb.helix.id,
                        r->pdb.helix.residues[0].name,
                        r->pdb.helix.residues[0].chain_id,
                        r->pdb.helix.residues[0].seq_num,
                        r->pdb.helix.residues[0].insert_code,
                        r->pdb.helix.residues[1].name,
                        r->pdb.helix.residues[1].chain_id,
                        r->pdb.helix.residues[1].seq_num,
                        r->pdb.helix.residues[1].insert_code,
                        r->pdb.helix.class, r->pdb.helix.comment);
                break;

        case PDB_HET:
                pdb_sprintf(buffer, fmt, r->pdb.het.het_grp.name,
                        r->pdb.het.het_grp.chain_id, r->pdb.het.het_grp.seq_num,
                        r->pdb.het.het_grp.insert_code, r->pdb.het.num_atoms,
                        r->pdb.het.text);
                break;

        case PDB_HETNAM:
                pdb_sprintf(buffer, fmt, r->pdb.hetnam.serial_num, r->pdb.hetnam.het_id,
                        r->pdb.hetnam.desc, r->pdb.hetnam.extra);
                break;

        case PDB_HETSYN:
                pdb_sprintf(buffer, fmt, r->pdb.hetsyn.serial_num, r->pdb.hetsyn.het_id,
                        r->pdb.hetsyn.desc,"");
                break;

        case PDB_MASTER:
                pdb_sprintf(buffer, fmt, r->pdb.master.num_remark,
                        r->pdb.master.num_ftnote, r->pdb.master.num_het,
                        r->pdb.master.num_helix, r->pdb.master.num_sheet,
                        r->pdb.master.num_turn, r->pdb.master.num_site,
                        r->pdb.master.num_transform,
                        r->pdb.master.num_coordinate, r->pdb.master.num_ter,
                        r->pdb.master.num_conect, r->pdb.master.num_seqres);
                break;

        case PDB_MTRIX:
                pdb_sprintf(buffer, fmt, r->pdb.mtrix.row_num,
                        r->pdb.mtrix.serial_num, r->pdb.mtrix.m1,
                        r->pdb.mtrix.m2, r->pdb.mtrix.m3, r->pdb.mtrix.v,
                        r->pdb.mtrix.given);
                break;

        case PDB_OBSLTE:
                pdb_sprintf(buffer, fmt, r->pdb.obslte.continuation,
                        r->pdb.obslte.date, r->pdb.obslte.old_id,
                        r->pdb.obslte.id_map[0], r->pdb.obslte.id_map[1],
                        r->pdb.obslte.id_map[2], r->pdb.obslte.id_map[3],
                        r->pdb.obslte.id_map[4], r->pdb.obslte.id_map[2],
                        r->pdb.obslte.id_map[6], r->pdb.obslte.id_map[7]);
                break;

        case PDB_ORIGX:
                pdb_sprintf(buffer, fmt, r->pdb.origx.row_num, r->pdb.origx.o1,
                        r->pdb.origx.o2, r->pdb.origx.o3, r->pdb.origx.t);
                break;

        case PDB_REMARK:
                pdb_sprintf(buffer, fmt, r->pdb.remark.num, r->pdb.remark.text);
                break;

        case PDB_REVDAT:
                pdb_sprintf(buffer, fmt, r->pdb.revdat.modification,
                        r->pdb.revdat.continuation, r->pdb.revdat.date,
                        r->pdb.revdat.id, r->pdb.revdat.mod_type,
                        r->pdb.revdat.corrections);
                break;

        case PDB_SCALE:
                pdb_sprintf(buffer, fmt, r->pdb.scale.row_num, r->pdb.scale.s1,
                        r->pdb.scale.s2, r->pdb.scale.s3, r->pdb.scale.u);
                break;

        case PDB_SEQRES:
                pdb_sprintf(buffer, fmt, r->pdb.seqres.serial_num,
                        r->pdb.seqres.chain_id, r->pdb.seqres.count,
                        r->pdb.seqres.names[0], r->pdb.seqres.names[1],
                        r->pdb.seqres.names[2], r->pdb.seqres.names[3],
                        r->pdb.seqres.names[4], r->pdb.seqres.names[5],
                        r->pdb.seqres.names[6], r->pdb.seqres.names[7],
                        r->pdb.seqres.names[8], r->pdb.seqres.names[9],
                        r->pdb.seqres.names[10], r->pdb.seqres.names[11],
                        r->pdb.seqres.names[12]);
                break;

        case PDB_SHEET:
                sh = &r->pdb.sheet;
                shr0 = &sh->residues[0];
                shr1 = &sh->residues[1];
                sha0 = &sh->atoms[0].residue;
                sha1 = &sh->atoms[1].residue;
                pdb_sprintf(buffer, fmt, sh->strand_num,
                        sh->id, sh->count,
                        shr0->name, shr0->chain_id, shr0->seq_num,
                        shr0->insert_code,
                        shr1->name, shr1->chain_id, shr1->seq_num,
                        shr1->insert_code,
                        sh->sense,
                        sh->atoms[0].name,
                        sha0->name, sha0->chain_id, sha0->seq_num,
                        sha0->insert_code,
                        sh->atoms[1].name,
                        sha1->name, sha1->chain_id, sha1->seq_num,
                        sha1->insert_code);
                break;

        case PDB_SITE:
                shr0 = &r->pdb.site.residues[0];
                shr1 = &r->pdb.site.residues[1];
                sha0 = &r->pdb.site.residues[2];
                sha1 = &r->pdb.site.residues[3];
                pdb_sprintf(buffer, fmt, r->pdb.site.seq_num,
                        r->pdb.site.id, r->pdb.site.count,
                        shr0->name, shr0->chain_id, shr0->seq_num,
                        shr0->insert_code,
                        shr1->name, shr1->chain_id, shr1->seq_num,
                        shr1->insert_code,
                        sha0->name, sha0->chain_id, sha0->seq_num,
                        sha0->insert_code,
                        sha1->name, sha1->chain_id, sha1->seq_num,
                        sha1->insert_code);
                break;

        case PDB_SPRSDE:
                pdb_sprintf(buffer, fmt, r->pdb.sprsde.continuation,
                        r->pdb.sprsde.date, r->pdb.sprsde.id,
                        r->pdb.sprsde.supersede[0], r->pdb.sprsde.supersede[1],
                        r->pdb.sprsde.supersede[2], r->pdb.sprsde.supersede[3],
                        r->pdb.sprsde.supersede[4], r->pdb.sprsde.supersede[5],
                        r->pdb.sprsde.supersede[6], r->pdb.sprsde.supersede[7]);
                break;

        case PDB_SSBOND:
                pdb_sprintf(buffer, fmt, r->pdb.ssbond.seq_num,
                        r->pdb.ssbond.residues[0].name,
                        r->pdb.ssbond.residues[0].chain_id,
                        r->pdb.ssbond.residues[0].seq_num,
                        r->pdb.ssbond.residues[0].insert_code,
                        r->pdb.ssbond.residues[1].name,
                        r->pdb.ssbond.residues[1].chain_id,
                        r->pdb.ssbond.residues[1].seq_num,
                        r->pdb.ssbond.residues[1].insert_code,
                        r->pdb.ssbond.symop[0],
                        r->pdb.ssbond.symop[1],
                        r->pdb.ssbond.distance);
                break;

        case PDB_LINK:
                pdb_sprintf(buffer, fmt, 
                        r->pdb.link.name[0],
                        r->pdb.link.alt_loc[0],
                        r->pdb.link.residues[0].name,
                        r->pdb.link.residues[0].chain_id,
                        r->pdb.link.residues[0].seq_num,
                        r->pdb.link.residues[0].insert_code,
                        r->pdb.link.name[1],
                        r->pdb.link.alt_loc[1],
                        r->pdb.link.residues[1].name,
                        r->pdb.link.residues[1].chain_id,
                        r->pdb.link.residues[1].seq_num,
                        r->pdb.link.residues[1].insert_code,
                        r->pdb.link.symop[0],
                        r->pdb.link.symop[1],
                        r->pdb.link.distance);
                break;

        case PDB_TER: {
/*
"TER   %5d      %-3s %C%4d%C"
TER    3269      GLN   203
TER   150016      WAT  35296
*/

#ifdef ELABORATE_PDB_TER
                if (!GDefaults.bPdbHybrid36) {
                        int atom_overflow = (r->pdb.atom.serial_num > 99999);
                        int residue_overflow = (r->pdb.atom.residue.seq_num > 9999);
                        if (atom_overflow  &&  !residue_overflow) {
                                fmt = "TER   %6d     %3s%%.2s%4d%C";
                        } else if (!atom_overflow  &&  residue_overflow) {
                                fmt = "TER   %5d      %3s%.2s%5d%C";
                        } else if (atom_overflow  &&  residue_overflow) {
                                fmt = "TER   %6d     %3s%.2s%4d%C";
                        }
                }
                pdb_sprintf(buffer, fmt, r->pdb.ter.serial_num,
                        r->pdb.ter.residue.name,
                        r->pdb.ter.residue.chain_id,
                        r->pdb.ter.residue.seq_num,
                        r->pdb.ter.residue.insert_code);
#else
                fmt = "TER   ";
                pdb_sprintf(buffer, fmt );
#endif
                break;
        }

        case PDB_TURN:
                pdb_sprintf(buffer, fmt, r->pdb.turn.seq_num,
                        r->pdb.turn.id,
                        r->pdb.turn.residues[0].name,
                        r->pdb.turn.residues[0].chain_id,
                        r->pdb.turn.residues[0].seq_num,
                        r->pdb.turn.residues[0].insert_code,
                        r->pdb.turn.residues[1].name,
                        r->pdb.turn.residues[1].chain_id,
                        r->pdb.turn.residues[1].seq_num,
                        r->pdb.turn.residues[1].insert_code,
                        r->pdb.turn.comment);
                break;

        case PDB_TVECT:
                pdb_sprintf(buffer, fmt, r->pdb.tvect.serial_num,
                        r->pdb.tvect.t1, r->pdb.tvect.t2, r->pdb.tvect.t3,
                        r->pdb.tvect.comment);
                break;

        case PDB_USER:
                pdb_sprintf(buffer, fmt, r->pdb.user.subtype, r->pdb.user.text);
                break;

        case PDB_ENDMDL:
                pdb_sprintf(buffer, fmt);
                break;

        case PDB_MODEL:
                pdb_sprintf(buffer, fmt, r->pdb.model.num);
                break;
				
        case MOL2_ATOM:
                pdb_sprintf(buffer, fmt, r->pdb.atom.serial_num,
                        r->pdb.atom.name, 
                        r->pdb.atom.x, r->pdb.atom.y, r->pdb.atom.z,
                        r->pdb.atom.type_at, r->pdb.atom.residue.seq_num,
						r->pdb.atom.residue.name, r->pdb.atom.temp_factor);
                break;
       			
				
				
        default:
                (void) sprintf(buffer, "unknown pdb record #%d",
                                                                r->record_type);
                break;
        }

        register        char    *s, *t;
        /* Do not shorten TER/END cards */
        if (r->record_type == PDB_TER || r->record_type == PDB_END)
                (void) fprintf(f, "%s\n", buffer);
        else {
                /* find last non-blank in buffer, and shorten it */
                t = NULL;
                for (s = buffer; *s != '\0'; s++)
                        if (!isspace(*s))
                                t = s + 1;
                if (t == NULL)  /* this should never happen, but ... */
                        t = buffer;
                *t = '\0';
                (void) fprintf(f, "%s\n", buffer);
        }
}


#include <string.h>
#include <stdlib.h>

/*
 * pdb_write_multiline()
 *
 * Splits content on semicolons (and within segments if needed) into
 * chunks of at most field_len chars, writing each as a separate record.
 * Calls pdb_write_record(f, r) to emit each line; that routine knows
 * the full struct layout. We only touch field[] and *serial_num.
 *
 * Split priority for oversized segments:
 *   1. Last space  (split before it)
 *   2. Last dash   (split after it)
 *   3. Hard cut at field_len
 */

/* Forward declaration: implemented elsewhere, knows the full struct layout. */
void pdb_write_record(FILE *f, pdb_record *r);

static int best_split(const char *s, int max)
{
    int i, last_space = -1, last_dash = -1;
    for (i = 0; i < max; i++) {
        if (s[i] == ' ') last_space = i;
        if (s[i] == '-') last_dash  = i;
    }
    if (last_space > 0) return last_space;
    if (last_dash  > 0) return last_dash + 1;
    return max;
}

void pdb_write_multiline(FILE *f, pdb_record *r, int *serial_num,
                         char *field, size_t field_len,
                         const char *content)
{
    int flen = 0;
    field[0] = '\0';
    const char *seg = content;
    while (seg && *seg) {
        /* Trim leading spaces at the start of a line */
        if (flen == 0) while (*seg == ' ') seg++;
        if (!*seg) break;

        /* Find end of current semicolon segment */
        const char *semi = strchr(seg, ';');
        int slen = semi ? (int)(semi - seg + 1)  /* include ';' */
                        : (int)strlen(seg);

        if (flen + slen <= field_len) {
            /* Whole segment fits — append and advance */
            memcpy(field + flen, seg, slen);
            flen += slen;
            field[flen] = '\0';
            seg = semi ? semi + 1 : NULL;
        } else if (flen > 0) {
            /* Doesn't fit and line has content — flush and retry same seg */
            pdb_write_record(f, r);
            (*serial_num)++;
            flen = 0;
            field[0] = '\0';
            /* (seg unchanged — retry from same position) */
        } else {
            /* Line is empty but segment is still too long — must split */
            int cut = best_split(seg, field_len);
            memcpy(field, seg, cut);
            flen = cut;
            field[flen] = '\0';
            pdb_write_record(f, r);
            (*serial_num)++;
            flen = 0;
            field[0] = '\0';
            seg += cut;
            /* Don't advance past semi — remainder of segment loops back */
        }
    }

    /* Flush any remaining content */
    if (flen > 0) pdb_write_record(f, r);
}
