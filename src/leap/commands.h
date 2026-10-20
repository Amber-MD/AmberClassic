/*
 *      File:   commands.h
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
 *      External declarations for ALL command handlers.
 *
 *      Each command handler has it's name prefixed with 'Cmd' and
 *      has two arguments, an integer representing the number
 *      of arguments and an array of OBJEKTS which are the arguments
 *
 */



#ifndef	COMMANDS_H
#define	COMMANDS_H

#ifndef	VARARRAY_H
#include	"varArray.h"
#endif

typedef OBJEKT  (*CMD_FUNCTION)(int iArgCount, ASSOC aaArgs[]);

typedef struct {
        STRING          sName;
        CMD_FUNCTION        fCallback;
} COMMANDt;

extern  COMMANDt        cCommands[];


typedef struct {
        STRING  sName;
        STRING  sCommand;
} ALIASt;

typedef ALIASt	*ALIAS;

extern	VARARRAY	GvaAlias;

OBJEKT
oCmd_charge( int, ASSOC* );

OBJEKT
oCmd_check( int, ASSOC* );
#endif
