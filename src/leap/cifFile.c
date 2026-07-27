
#include        "basics.h"
#include        "pdbFile.h"
#include        "defaults.h"
#include        "matrix.h"
#include        "cifparse.h"

// can't use flexible size struct and static initialization
// (GCC does allow it in some contexts, but it's not standard C)
#define CIF_MAXCOLUMNS 20
typedef struct _CIFCATEGORY {
    char *sName;
    NdbCifCategoryFormat *pCategory; // lookup
    struct {
        char *name, *alt_name;
        bool bOptional;
        int iColumn;                 // lookup
    } field[CIF_MAXCOLUMNS];
} CIFCATEGORYt;

static bool
zbCifLookup(CIFCATEGORYt *cifCat, int iBlock, bool bAltAtomName) {
    int iCategory = get_category_index(iBlock, cifCat->sName);
    if (iCategory < 0) {
        VPFATAL("CIF data block %s does not contain data category %s\n",
                cifFiles.datablocks[iBlock].datablockName, cifCat->sName);
        return FALSE;
    }
    cifCat->pCategory = & cifFiles.datablocks[iBlock].categories[iCategory];
    for (int i=0; i<CIF_MAXCOLUMNS && cifCat->field[i].name; i++) {
        char *name = bAltAtomName ? cifCat->field[i].alt_name : cifCat->field[i].name;
        char *alt_name = bAltAtomName ? cifCat->field[i].name : cifCat->field[i].alt_name;
        cifCat->field[i].iColumn = get_column_index(iBlock, iCategory, name);
        if (cifCat->field[i].iColumn >= 0) continue;
        cifCat->field[i].iColumn = get_column_index(iBlock, iCategory, alt_name);
        if (cifCat->field[i].iColumn >= 0) continue;
        if (!cifCat->field[i].bOptional) {
            VPFATAL("CIF parse error: Column %s.%s not found\n",cifCat->sName,cifCat->field[i].name);
                return FALSE;
        }
    }
    return TRUE;
}

static const char *
zcPCifGetItem(CIFCATEGORYt *cifCat, int iRow, int iColumn) {
    static const char *null = "";
    if (cifCat->field[iColumn].iColumn < 0) return null;
    return cifCat->pCategory->rows[iRow].columns[cifCat->field[iColumn].iColumn];
}

void
CifReadFile( PDBREADt *prPRead )
{
int             iPdbSequence=0; // previous resSeq
bool            bLastReadPdbRecordWasTer = FALSE;
bool            bNewChain = TRUE, bNewRes;
RESIDUENAMEt    rnName;
ATOMNAMEt       anAtom;
int             iTerm, iLast, iSerialNumMax=0;
char            sCurrChain[3]="", cInsertionCode=' ';
int             iMultipleResName=0;
int             iCurrentModel=0;

CIFCATEGORYt cifAtoms = {
    "atom_site", NULL,
    {
        { "id" },  /* atomSerial */
        { "label_atom_id", "auth_atom_id" },  /* name */
        { "label_alt_id", NULL, TRUE },   /* altLoc */
        { "label_comp_id", "auth_comp_id" },  /* resName */
        { "label_asym_id", "auth_asym_id" },  /* chainID */
        { "label_seq_id", "auth_seq_id" },   /* resSeq */
        { "pdbx_PDB_ins_code", NULL, TRUE },/* iCode */
        { "Cartn_x" }, /* x-coord */
        { "Cartn_y" }, /* y-coord */
        { "Cartn_z" }, /* z-coord */
        //{ "occupancy" },
        //{ "B_iso_or_equiv" },
        { "type_symbol" }, // Element name, *not* force field type
                                   // Note: FF type is .type_energy
        //{ "pdbx_formal_charge" }, /* formal charge */
        { "pdbx_PDB_model_num", NULL, TRUE },
        //{ "group_PDB", NULL, TRUE }, -> "ATOM", "HETATM"
        {NULL}
    }
}, cifConn = {
    "struct_conn", NULL,
    {
        { "conn_type_id" },            // enum: covale, disulf, metalc, hydrog
        { "ptnr1_label_asym_id", "ptnr1_auth_asym_id" },     // chainID, alt=ptnr1_auth_asym_id
        { "ptnr1_label_comp_id", "ptnr1_auth_comp_id" },     // resName, alt=ptnr1_auth_comp_id
        { "ptnr1_label_seq_id", "ptnr1_auth_seq_id" },      // resSeq, alt=ptnr1_auth_seq_id
        { "ptnr1_label_atom_id" },     // name
        { "pdbx_ptnr1_label_alt_id", NULL, TRUE }, // altLoc (optional)
        { "pdbx_ptnr1_PDB_ins_code", NULL, TRUE }, // iCode  (optional)
        { "ptnr1_symmetry", NULL, TRUE },          //<op>_<dx+5><dy+5><dz+5> (optional)
        { "ptnr2_label_asym_id", "ptnr2_auth_asym_id" },     // alt=ptnr2_auth_asym_id
        { "ptnr2_label_comp_id", "ptnr2_auth_comp_id" },     // alt=ptnr2_auth_comp_id
        { "ptnr2_label_seq_id", "ptnr2_auth_seq_id" },      // alt=ptnr2_auth_seq_id
        { "ptnr2_label_atom_id" },
        { "pdbx_ptnr2_label_alt_id", NULL, TRUE },
        { "pdbx_ptnr2_PDB_ins_code", NULL, TRUE },
        { "ptnr2_symmetry", NULL, TRUE },
        { "pdbx_dist_value", NULL, TRUE },
        {NULL}
    }
}, cifMtrix = {
    // converted to "pdbx_struct_oper_list" if GDefaults.iPdbReadBioMT
    "struct_ncs_oper", NULL,
    {
        { "id" },
        { "code" }, // converted to "type" if GDefaults.iPdbReadBioMT
        { "matrix[1][1]" },
        { "matrix[2][1]" },
        { "matrix[3][1]" },
        { "vector[1]" },
        { "matrix[1][2]" },
        { "matrix[2][2]" },
        { "matrix[3][2]" },
        { "vector[2]" },
        { "matrix[1][3]" },
        { "matrix[2][3]" },
        { "matrix[3][3]" },
        { "vector[3]" },
        {NULL}
    }
};
    VPTRACEENTER("zCifReadFile" );

    // NOTE: 1990s simple NDB CIF reader, modified to give caller direct access via globals
    // Hacked for antechamber use, still seems to work OK for mmCIF coordinate files
    ndb_cif_init();
    int nBlocks = ndb_cif_read_file(prPRead->fPdbFile);
    if (!nBlocks) {
        VPFATAL("No datablocks read from CIF file\n");
        return;
    }
    if (nBlocks>1) {
        VPWARN("Multiple data blocks in CIF file\n");
    }
    int iCifBlock = 0;
    VP0("Reading first CIF datablock: %s\n",cifFiles.datablocks[iCifBlock].datablockName);

    bNewChain = TRUE;

    if (GDefaults.iPdbReadBioMT) {
        cifMtrix.sName = "pdbx_struct_oper_list";
        cifMtrix.field[1].name = "type";
/* FIXME: we need to parse this data to process correctly
_pdbx_struct_assembly_gen.assembly_id       1
_pdbx_struct_assembly_gen.oper_expression   1
_pdbx_struct_assembly_gen.asym_id_list      A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,AA,BA,CA,DA,EA,FA,GA,HA,IA,JA
               Operation expressions may have the forms:

                (1)        the single operation 1
                (1,2,5)    the operations 1, 2, 5
                (1-4)      the operations 1,2,3 and 4
                (1,2)(3,4) the combinations of operations
                           3 and 4 followed by 1 and 2 (i.e.
                           the cartesian product of parenthetical
                           groups applied from right to left)
example 6O3H.cif data:
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '(1-60)'           A,B,C,D,E,F,G,H,I,J,K,L,M,N
2 1                  A,B,C,D,E,F,G,H,I,J,K,L,M,N
3 '(1-5)'            A,B,C,D,E,F,G,H,I,J,K,L,M,N
4 '(1,2,6,10,23,24)' A,B,C,D,E,F,G,H,I,J,K,L,M,N
5 P                  A,B,C,D,E,F,G,H,I,J,K,L,M,N

*/
    }
    if (!( zbCifLookup(&cifAtoms,iCifBlock,GDefaults.bCIFReadAuth)
             && zbCifLookup(&cifConn,iCifBlock,GDefaults.bCIFReadAuth)
             && zbCifLookup(&cifMtrix,iCifBlock,GDefaults.bCIFReadAuth) ) ) {
        return;
    }
    // Collect atom-site data corresponding to ATOM, HETATM
    for (int i = 0; i < cifAtoms.pCategory->numRow; i++) {
        int iCol=0;
        anAtom.iAtomSerial = atoi(zcPCifGetItem(&cifAtoms, i, iCol++));
        memcpy(anAtom.sName, zcPCifGetItem(&cifAtoms, i, iCol++),sizeof(anAtom.sName));
        anAtom.sName[sizeof(anAtom.sName)-1]=0;
        char altLoc = zcPCifGetItem(&cifAtoms, i, iCol++)[0];
        if (altLoc != 0 && altLoc != GDefaults.cPdbAltLocSelect) continue;
        const char *resName = zcPCifGetItem(&cifAtoms, i, iCol++);
        char sChainId[3];
        memcpy(sChainId,zcPCifGetItem(&cifAtoms, i, iCol++),sizeof(sChainId));
        sChainId[2]=0;
        if (!strcmp(sChainId," ")) strcpy(sChainId, "");
        int resSeq = atoi(zcPCifGetItem(&cifAtoms, i, iCol++));
        char iCode = zcPCifGetItem(&cifAtoms, i, iCol++)[0];
        if (!iCode) iCode = ' ';
        anAtom.x = atof(zcPCifGetItem(&cifAtoms, i, iCol++));
        anAtom.y = atof(zcPCifGetItem(&cifAtoms, i, iCol++));
        anAtom.z = atof(zcPCifGetItem(&cifAtoms, i, iCol++));
        const char *element = zcPCifGetItem(&cifAtoms, i, iCol++);
        if (*element) anAtom.iElement = iPdbElementNumber(element);
        else anAtom.iElement = iElementNumberFromAmber(anAtom.sName);
        iCurrentModel = atoi(zcPCifGetItem(&cifAtoms, i, iCol++));

        if (anAtom.iAtomSerial > iSerialNumMax) iSerialNumMax = anAtom.iAtomSerial;
        if (iCurrentModel) {
            if (!GDefaults.iPdbReadModel) GDefaults.iPdbReadModel = iCurrentModel;
            if (GDefaults.iPdbReadModel > 0 && GDefaults.iPdbReadModel != iCurrentModel) continue;
            if (GDefaults.iPdbReadModel < 0 && strlen(sChainId)<2) {
                int j = iCurrentModel-1;
                if (j < CHAINID_LIST_LEN) {
                    char c = (sChainId[0] == 0) ? '_':sChainId[0];
                    sChainId[0] = GsChainIdList[j];
                    sChainId[1] = c;
                    sChainId[2] = 0;
                } else VPFATAL("ChainID overflow for iPDB_Read_Model < 0:"
                               " Unable to relabel MODEL=%d ChainID='%s'\n",iCurrentModel, sChainId);
            }
        }

        iTerm = NOEND;
        bNewChain = bNewChain || strcmp(sCurrChain,sChainId);
        if ( bNewChain ) {
            VP2(" (starting new molecule for chain \"%.2s\")\n", sChainId );
            iTerm = FIRSTEND;
            iLast = iVarArrayElementCount( prPRead->vaResidues );
            if (iLast >= 0) PVAI( (prPRead->vaResidues), RESIDUENAMEt,iLast-1)->iTerminator = LASTEND;
            strcpy(sCurrChain,sChainId);
        }
        //TODO bLastReadPdbRecordWasTer == change in entity_id
        bNewRes = (bNewChain || resSeq != iPdbSequence ||
                    iCode != cInsertionCode ||
                    bLastReadPdbRecordWasTer );
        if (!bNewRes && strcmp(rnName.sName, resName) ) {
            /* First detection of a residue with the same sequence
             * number and insertion code but a different name.
             */
            VPWARN("Name change in pdb file residue %.2s %d%c;\n"
                "this residue is split into %s and %s.\n",
                sCurrChain, iPdbSequence, cInsertionCode, rnName.sName, resName);
            bNewRes = TRUE;
            iMultipleResName++;
        }
        if (bNewRes) {
            VPTRACE("Detected a new residue.\n" );
            rnName.iTerminator = iTerm;
            rnName.iPdbSequence = resSeq;
            strncpy(rnName.sChainId, sChainId, sizeof(rnName.sChainId));
            rnName.sChainId[sizeof(rnName.sChainId)-1]=0;
            memcpy(rnName.sName, resName, sizeof(rnName.sName));
            rnName.sName[sizeof(rnName.sName)-1]=0;
            rnName.iCode = iCode;
            rnName.iFirstAtom = iVarArrayElementCount(prPRead->vaAtomRecs); // zero based array

            anAtom.iResNameIndex = iVarArrayElementCount( prPRead->vaResidues );
            VarArrayAdd( prPRead->vaResidues, (GENP)&rnName );
            bLastReadPdbRecordWasTer = FALSE;

            MESSAGE("Reading residue: <%s>\n", rnName.sName );
            iPdbSequence = resSeq;

        }
        bNewChain = FALSE;

        VarArrayAdd( prPRead->vaAtomRecs, (GENP)&anAtom );
    }
    if (cifMtrix.pCategory) {
        // Convert struct_ncs_oper to MTRIXn form
        for (int i = 0; i < cifMtrix.pCategory->numRow; i++) {
            int iCol=0;
            const char *pItem = zcPCifGetItem(&cifMtrix, i, iCol++);
            // Crude filtering: assume all numerically named transforms applied to all chains
            // And assumes first entry is identity matrix !!FIXME!!
            if (!isdigit((unsigned char)*pItem)) continue;
            int iSerial = atoi(pItem);
            //pItem = zcPCifGetItem(&cifMtrix, i, iCol++); // ncs "code": given or identity = self, always 1?
                                                         // or BioMT "type"
            // FIXME: how to use code/type info
            // e.g. 6O3H P='transform to point frame', 1='identity operation', 2='point symmetry operation'
            iCol++;
            PDBMATRIXt matrix = { (iSerial != 1), {
                 {
                     atof(zcPCifGetItem(&cifMtrix, i, iCol++)),
                     atof(zcPCifGetItem(&cifMtrix, i, iCol++)),
                     atof(zcPCifGetItem(&cifMtrix, i, iCol++)),
                     atof(zcPCifGetItem(&cifMtrix, i, iCol++)),
                 },{
                     atof(zcPCifGetItem(&cifMtrix, i, iCol++)),
                     atof(zcPCifGetItem(&cifMtrix, i, iCol++)),
                     atof(zcPCifGetItem(&cifMtrix, i, iCol++)),
                     atof(zcPCifGetItem(&cifMtrix, i, iCol++)),
                 },{
                     atof(zcPCifGetItem(&cifMtrix, i, iCol++)),
                     atof(zcPCifGetItem(&cifMtrix, i, iCol++)),
                     atof(zcPCifGetItem(&cifMtrix, i, iCol++)),
                     atof(zcPCifGetItem(&cifMtrix, i, iCol++)),
                 },{
                     0.0,0.0,0.0,1.0
                 }
            }};
            VarArrayAdd( prPRead->vaMatrices, (GENP)&matrix );
        }
    }

    if (cifConn.pCategory && GDefaults.bPdbUseLinkRecords) {
        // Convert struct_conn to LINK record format
        for (int i = 0; i < cifConn.pCategory->numRow; i++) {
            int iCol=0;
            const char *type = zcPCifGetItem(&cifConn, i, iCol++);
            if (type[0] != 'c' && type[0] != 'd') continue;
            struct pdb_link pdbLink = {0};
            for (int j=0; j<2; j++) {
                strncpy(pdbLink.residues[j].chain_id,zcPCifGetItem(&cifConn, i, iCol++), sizeof(pdbLink.residues[j].chain_id));
                pdbLink.residues[j].chain_id[sizeof(pdbLink.residues[j].chain_id)-1]=0;
                strcpy(pdbLink.residues[j].name, zcPCifGetItem(&cifConn, i, iCol++) );
                pdbLink.residues[j].seq_num = atoi( zcPCifGetItem(&cifConn, i, iCol++) );
                strcpy(pdbLink.name[j], zcPCifGetItem(&cifConn, i, iCol++) );
                pdbLink.alt_loc[j] = zcPCifGetItem(&cifConn, i, iCol++)[0];
                if (pdbLink.alt_loc[j] == 0 ) pdbLink.alt_loc[j]=' ';
                pdbLink.residues[j].insert_code = zcPCifGetItem(&cifConn, i, iCol++)[0];
                if (pdbLink.residues[j].insert_code == 0) pdbLink.residues[j].insert_code = ' ';
                strcpy(pdbLink.symop[j], zcPCifGetItem(&cifConn, i, iCol++) );
            }
            pdbLink.distance = atof(zcPCifGetItem(&cifConn, i, iCol++));
            VarArrayAdd( prPRead->vaLinkRecs, (GENP)&pdbLink );
        }
    }

                /* Make the last residue a LASTEND */

    if ( (iLast = iVarArrayElementCount( prPRead->vaResidues )) ) {
        iLast--;        /* could be 0th element */
        PVAI( (prPRead->vaResidues), RESIDUENAMEt,iLast)->iTerminator = LASTEND;
        prPRead->iMaxSerialNum = iSerialNumMax+1;
    }

    if ( iMultipleResName ) {
        // another place conflation of PDB ResId and residue Sequence #
        VPNOTE("%d residues had naming warnings.\n"
            "Thus, there are split residues;\n"
            "residue sequence numbers will not correspond to those in the pdb.\n",
            iMultipleResName );
    }

}


