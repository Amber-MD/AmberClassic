#ifndef ParseHeadings
#define ParseHeadings

#include "MatrixDS.h"
#include "TopologyDS.h"
#include "CrdManipDS.h"
#include "ParamFitDS.h"
#include "ChargeFitDS.h"
#include "ConfigSampDS.h"
#include "SinglePointEvalDS.h"
#include "ParseDS.h"
#include "PeptideDS.h"

int CountWords(char* line);

char ToUpper(char c);

char ToLower(char c);

void StringToUpper(char* s);

void StringToLower(char* s);

void Type4Char(char* a, char* t);

int CaselessStrcmp(char* s1, char* s2, int nc);

double RealXpYf(char* word, int X, int Y);

int WordIsInteger(char* word);

int WordIsNumber(char* word);

int WordIsAtomType(char* word);

cmat ParseWords(char* line);

int LocateWordInLine(char* w, char* s, int minpos);

void RemoveWhiteSpace(char* a, int tail, char ws);

void EqualSpace(char* line);

void RemoveComments(char* line);

void NixCommaCarriage(char* line);

void AddToGlossary(cmat *glossary, int nargs, ...);

int WordInGlossary(cmat *glossary, char* word);

int AdvanceToSegment(FILE *inp, char* segname, int scan0);

int DetectNamelistEnd(char* line, char* errmsg);

int ReadNamelistLine(char* line, cmat *lwords, char* callfunc, FILE *inp);

void SeekString(cmat L, char* val, char* sname, char* salias);

void SeekMultiString(cmat L, cmat *val, char* sname, char* salias,
                     cmat *glossary, int append);

void SeekStringInc(cmat L, cmat *val, char* sname, char* salias, int *counter);

void SeekSSR(cmat L, char* val1, char* val2, double *val3, char* sname,
             char* salias, int *counter);

void SeekS3R(cmat L, char* val1, char* val2, char* val3, double *val4,
	     char* sname, char* salias, int *counter);

void SeekStringPlusVal(cmat L, char* val1, double *val2, char* sname,
		       char* salias, int *counter);

void SeekRecord(cmat L, cmat *C, char* sname, char* salias, int *counter);

cmat SeekNString(cmat L, cmat* val, int* fspec, char* sname, char* salias);

void SeekFixQ(cmat L, char* val1, double *val2, char* val3, double *val4,
              char* sname, char* salias, int *counter);

void SeekAmberLibrary(cmat L, fset *myfit, cmat *glossary);

void SeekManipulator(cmat L, configs *cfsinp, int *maxops, char* sname,
                     char* salias, cmat *glossary);

void SeekPeptideSystem(cmat L, pepcon *ppctrl, int *maxsys, char* sname,
		       char* salias, cmat *glossary);

void SeekOpsCombo(cmat L, configs *cfsinp, int *maxops, char* sname,
                  char* salias, cmat *glossary);

void SeekSPItem(cmat L, spdata *spelist, int *maxitems, char* sname,
                char* salias, cmat *glossary);

void SeekNumberedName(cmat L, char* base, char* suffix, char* sname,
                      char* salias);

void SeekReal(cmat L, double *val, char* sname, char* salias);

void SeekNReal(cmat L, double* val, char* sname, char* salias, int maxidx);

void SeekInt(cmat L, int *val, char* sname, char* salias);

void SeekLLInt(cmat L, long long int *val, char* sname, char* salias);

int* ParseAmbMask(char* maskstr, prmtop *tp, coord *crd);

FILE* FOpenSafe(char* fname, int ovrwrt);

void SeekBondTermID(cmat L, prmset *mp, char* sname, char* salias,
                    int *maxhadj, int order);

void SeekRecast(cmat L, prmset *mp, char* sname, char* salias, int *maxhold,
                int specinst);

void TouchupAtomName(char* atn);

void SeekSinglePoint(cmat L, prmset *mp, char* sname, char* salias,
                     int *maxconf, cmat *glossary);

void SeekSpecParmRest(cmat L, prmset *mp, char* sname, char* salias, int order,
                      int *maxsr, cmat *glossary);

void SeekGeomRest(cmat L, prmset *mp, char* sname, char* salias, int *maxgeom);

void SeekSpecReq(cmat L, prmset *mp, char* sname, char* salias, int *maxsr,
		 cmat *glossary);

long long int ReadNumericalShorthand(char* numstr);

int DetectBinaryFile(char* fname);

int FindFileType(const char* path);

int Alphabetical(const void *a, const void *b);

int* AlphabetizeCmat(cmat *C);

void JoinPath(char* newpath, char* path1, char slash, char* path2);

cmat RegExpFileSearch(char* query);

cmat DirectoryFileSearch(char* mydir);

int DoubleFromLine(char* line, int widx, double *val);

int ValidRealData(char* fname);

dmat ReadRealData(char* fname);

nmlgroup ReadNamelistsFromFile(char* fname, int nargs, ...);

#endif
