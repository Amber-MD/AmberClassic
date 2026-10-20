/*
 *      File:   library.h
 *
 ************************************************************************
 *                                                                      *
 *     Designed by:    Christian Schafmeister                           *
 *     Author:         Christian Schafmeister                           *
 *                                                                      *
 *     VERSION: 1.0                                                     *
 *     Programmers:                                                     *
 *             Christian Schafmeister                                   *
 *             David Rivkin                                             *
 *                                                                      *
 *     Principal Investigator: Peter A. Kollman                         *
 *                                                                      *
 ************************************************************************
 *
 *      Description:
 *              This object does very little except manage the
 *              reading and writing of several units into a 
 *              single DATABASE file.  It maintains an open
 *              DATABASE and a DICTIONARY of UNIT names.
 *              The object attatched to each UNIT name in the
 *              DICTIONARY is a dummy.
 */



#ifndef	LIBRARY_H
#define	LIBRARY_H

#ifndef	DICTIONARY_H
# include "dictionary.h"
#endif

#ifndef	DATABASE_H
# include "database.h"
#endif


typedef struct  {
	DICTIONARY	dLibrary;
	DATABASE	dbLibrary;
	DICTLOOP	dlContents;
} LIBRARYt;

typedef LIBRARYt	*LIBRARY;




extern LIBRARY	lLibraryOpen(char *sName, int iOpenMode);
extern OBJEKT	oLibraryLoad(LIBRARY ul, char *sName);
extern void	LibrarySave(LIBRARY ul, char *sName, OBJEKT oObj, PARMLIB plLibrary);
extern bool	bLibraryRemove(LIBRARY ul, char *sName);

extern void	LibraryLoop(LIBRARY ul);
extern char	*sLibraryNext(LIBRARY ul);
extern void	LibraryClose(LIBRARY *lPLibrary);

extern bool	bLibraryDBPrefix(LIBRARY lLib, char *cPResult, char *sStr );
#define	dbLibraryDatabase(l)	(l->dbLibrary)


#endif

