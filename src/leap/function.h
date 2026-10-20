#ifndef FUNCTION_H
#define FUNCTION_H
/*
 *	File:	function.h
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
 *	Description:
 *		Store filename/start-stop line numbers/function name
 *		and index them by integer.
 *
 */

#ifndef	BASICS_H
#include	"basics.h"
#endif

#ifndef	VARARRAY_H
#include	"varArray.h"
#endif

typedef	struct	{
	STRING		sFilename;
	int		iStart;
	int		iStop;
	STRING		sFunction;
	bool		bPrint;
} FUNCTIONt;

#define	NO_FUNCTION	-1

extern	VARARRAY	GvaFunctions;
extern	VARARRAY	GvaFiles;

extern void	FunctionAdd(char *sFilename, int iStart, int iStop, 
			char *sFunction);
extern int	iFunctionFindWithFilenameLine(char *sFilename, int iLine);
extern int	iFunctionFindWithFunction(char *sFunction);


extern void	FunctionFileAdd(char *sName);
extern int	iFunctionFilenameFind(char *sFilename);



#define	sFunctionFilename(i)	(PVAI(GvaFunctions,FUNCTIONt,i)->sFilename)
#define	iFunctionStart(i)	(PVAI(GvaFunctions,FUNCTIONt,i)->iStart)
#define	iFunctionStop(i)	(PVAI(GvaFunctions,FUNCTIONt,i)->iStop)
#define	sFunctionFunction(i)	(i==NO_FUNCTION \
	? "(?)" : PVAI(GvaFunctions,FUNCTIONt,i)->sFunction)
#define	bFunctionPrint(i)	(i==NO_FUNCTION \
	? true : PVAI(GvaFunctions,FUNCTIONt,i)->bPrint)

#define	FunctionSetPrint(i,b)	(PVAI(GvaFunctions,FUNCTIONt,i)->bPrint=(b))

#define	iFunctionCount()	(iVarArrayElementCount(GvaFunctions))

#define	iFunctionFileCount()	(iVarArrayElementCount(GvaFiles))
#define	sFunctionFileFilename(i)	(PVAI(GvaFiles,STRING,i))


#endif
