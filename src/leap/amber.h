/*
 *      File:   amber.h
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
 *              Read AMBER prep and parameter database files into
 *              UNITs and PARMSETs.
 *
 *              This module was written to maintain compatability with
 *              AMBER 3a.
 *
 */



#ifndef	AMBER_H
#define	AMBER_H

#ifndef	DICTIONARY_H
#include	"dictionary.h"
#endif

#define	AMBER_WILD_CARD	"X"


/*  amber.c  */

extern void		AmberAddAtomTypes( LIST lEntries );
extern DICTIONARY	dAmberReadPrepFile( char *sFilename );
extern PARMSET		psAmberReadParmSet( FILE *fIn, char *sFilename );
extern int		iAtomTypeHybridization( char *sType );

#endif
