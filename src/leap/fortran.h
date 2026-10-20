/*
 *      File:   fortran.h
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
 *
 *              This file contains routines for performing
 *              FORTRAN style formatted output.
 *              It contains a routine that allows the caller
 *              to specify the number of entries that should
 *              appear on each line and the format of each
 *              entry on the line.
 *              Then there is a routine that the user uses
 *              to write each subsequent piece of data.
 */

#ifndef	FORTRAN_H
# define FORTRAN_H

extern void	FortranFile(FILE *fOut);
extern void	FortranFormat(int iPerLine, char *sFormat);

extern void	FortranWriteInt(int iVal);
extern void	FortranWriteDouble(double dVal);
extern void	FortranWriteString(char *sVal);

extern void	FortranEndLine(void);

extern char	*sFortranReadString( char *sString );
extern char	*sFortranReadLabel( char *sString );

extern int	iFortranReadInt(void);
extern double	dFortranReadDouble(void);
extern void	FortranSkipLine(void);

extern void	FortranDebugOn(void);
extern void	FortranDebug(char *sStr);	

extern STRING	GsFortranDebug;


#endif
