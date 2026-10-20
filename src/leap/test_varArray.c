/*
 *	File:	test_varArray.c
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
 */

#include	"basics.h"

#include	"varArray.h"


int main(iint argc, char *argv[])
{
VARARRAY	vaTest;

int		i;

    BasicsInitialize();

    vaTest = vaVarArrayCreate(sizeof(int));
    i=1;
    VarArrayInsertBefore( vaTest, 0, &i );
    i=3;
    VarArrayInsertBefore( vaTest, 1, &i );
    i=2;
    VarArrayInsertBefore( vaTest, 1, &i );
    i = 0;
    VarArrayInsertBefore( vaTest, 0, &i );

    for ( i=0; i<iVarArrayElementCount(vaTest); i++ ) {
	printf( "varArray[%d] = %d\n", i, *PVAI(vaTest,int,i) );
    }
}

