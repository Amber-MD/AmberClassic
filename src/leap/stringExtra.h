/*
 *      File:   stringExtra.h
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
 *              Some extra routines to manipulate strings.
 */
 
#ifndef	STRINGEXTRA_H
#define	STRINGEXTRA_H



extern char	*sRemoveSpaces(char *sIn, char *sOut);
extern char	*sRemoveControlAndPadding(char *sRaw, char *sResult);
extern char	*sRemoveLeadingSpaces(char *sLine);
extern char	*sRemoveFirstString(char *sLine, char *sHead);
extern bool	bStringMatchPattern(char *sString, char *sPattern);
extern void	StringCopyMax(char *sDest, char *sSource, int iMax);

#endif
