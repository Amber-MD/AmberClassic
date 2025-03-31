/*
 *      File:   database.h
 *
 *      Description:
 *              Routines to manage an object oriented database.
 *              The database is an ASCII file with
 *              very little internal structure, making
 *              it human readible, editable, and portable.
 */

#ifndef DATABASE_H
#define DATABASE_H


#define MAXDATALINELEN  1000
#define MAXPREFIXSTACK  10

typedef	char	Bool;
typedef	char	String[256];

typedef struct  {
                long            lFileOffset;
                String          sName;
                int             iType;
                int             iRows;
                } ENTRYt;

typedef ENTRYt* ENTRY;

#define	DB_READ		1
#define	DB_WRITE	2

#define	DB_RANDOM_ACCESS	1
#define	DB_SEQUENTIAL_ACCESS	2

typedef struct  {
		int		iAccessMode;
                FILE*           fDataBase;
                String          sFileName;
		int		iOpenMode;
                int             iPrefix;
                String          saPrefixStack[MAXPREFIXSTACK];
                Bool            bCompactFileAtClose;
                void		*dEntries;
                int             iCurrentLine;
                char            sLookAhead[MAXDATALINELEN];
		int		iLastSequentialOperation;

			/* Loop over entries with a certian prefix */
		String		sLoopPrefix;
		int		dlEntryLoop;

                } DATABASEt;

typedef DATABASEt*      DATABASE;

#define	LENGTH_NOT_KNOWN	-1

#define ENTRYTYPE       0x0000000F
#define ENTRYINTEGER    0x00000001
#define ENTRYDOUBLE     0x00000002
#define ENTRYSTRING     0x00000003

#define ENTRYMODIFIER   0x000000F0
#define ENTRYSINGLE     0x00000010
#define ENTRYARRAY      0x00000020
#define ENTRYTABLE      0x00000040


/*
 *	DATABASE open modes.
 *
 *	File can be opened for read-only or read-write
 */

#define	OPENREADONLY	1
#define	OPENREADWRITE	2


/*
 *	DATABASE errors.
 */

extern	int	GiDBLastError;
#define	DB_ERROR_NONE			0
#define	DB_ERROR_INVALID_FILE		1
#define	DB_ERROR_INVALID_DATABASE	2


/*
 *--------------------------------------------------------------------
 *
 *	Routines for random access files
 */

extern  DATABASE dbDBRndOpen( char*, int );
extern  Bool    bDBRndDeleteEntry( DATABASE, char* );

extern	void	DBRndLoopEntryWithPrefix( DATABASE, char* );
extern	Bool	bDBRndNextEntryWithPrefix( DATABASE, char* );

/*
 *--------------------------------------------------------------------
 *
 *	Routines for accessing DB_SEQUENTIAL_ACCESS files
 */


extern	DATABASE dbDBSeqOpen();		/* ( char*, int ) */
extern	void	DBSeqRewind();		/* ( DATABASE ) */
extern	void	DBSeqSkipData();	/* ( DATABASE ) */
extern	Bool	bDBSeqEof();		/* ( DATABASE ) */

extern	long	lDBSeqCurPos();		/* ( DATABASE ) */
extern	void	DBSeqGoto();		/* ( DATABASE, long ) */


/*
 *--------------------------------------------------------------------
 *
 *	Routines used on both DB_RANDOM_ACCESS and
 *	DB_SEQUENTIAL_ACCESS files.
 */


#define	sDBName(d)	(d->sName)

extern  Bool    bDBGetType( DATABASE, char*, int*, int* );
extern  Bool    bDBGetValue();          /* ( DATABASE, char*, int*, GENP, int) */
extern  void    DBPutValue();           /* ( DATABASE, char*, int, int,
                                                GENP, int ) */
extern  Bool    bDBGetTableType();      /* ( DATABASE, char*, int*, int*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char*,
                                                int*, char* ) */

extern  Bool    bDBGetTable();          /* ( DATABASE, char*, int*,
                                                int, GENP, int,
                                                int, GENP, int,
                                                int, GENP, int,
                                                int, GENP, int,
                                                int, GENP, int,
                                                int, GENP, int,
Each line:                                      int, GENP, int,
 column, data, datainc                          int, GENP, int,
                                                int, GENP, int,
                                                int, GENP, int,
                                                int, GENP, int,
                                                int, GENP, int,
                                                int, GENP, int,
                                                int, GENP, int,
                                                int, GENP, int,
                                                int, GENP, int,
                                                int, GENP, int ) */


extern  void    DBPutTable();           /* ( DATABASE, char*, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int,
Each line:                                      int, char*, GENP, int,
  column, name, data, datainc                   int, char*, GENP, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int,
                                                int, char*, GENP, int ) */

extern  void    DBClose( DATABASE* );

extern  void    DBPushPrefix( DATABASE, char* );
extern  void    DBPopPrefix( DATABASE );
extern  void    DBZeroPrefix( DATABASE );
extern  void    DBPushZeroPrefix( DATABASE, char* );



#define	iDBLastError()	(GiDBLastError)

#endif  /* ifdef DATABASE_H */
