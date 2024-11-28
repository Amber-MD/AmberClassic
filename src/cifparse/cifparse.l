%option prefix="cifp"
%option noyywrap

%{ 
/*

COPYRIGHT 1996 Rutgers - The State University of New Jersey

This software is provided WITHOUT WARRANTY OF MERCHANTABILITY OR
FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR
IMPLIED.  RUTGERS MAKE NO REPRESENTATION OR WARRANTY THAT THE
SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT OR OTHER
PROPRIETARY RIGHT.

The user of this software shall indemnify, hold harmless and defend
Rutgers, its governors, trustees, officers, employees, students,
agents and the authors against any and all claims, suits,
losses, liabilities, damages, costs, fees, and expenses including
reasonable attorneys' fees resulting from or arising out of the
use of this software.  This indemnification shall include, but is
not limited to, any and all claims alleging products liability.

*/

/* 
  FILE:         cifparse.l

  This file is part of the NDBQUERY application,
  a program component of the Nucleic Acid Database Project.

  H. M. Berman, W. K. Olson, D. L. Beveridge, J. K. Westbrook, A. Gelbin,
  T. Demeny, S. H. Shieh, A. R. Srinivasan, and B. Schneider. (1992).
  The Nucleic Acid Database: A Comprehensive Relational Database of 
  Three-Dimensional Structures of Nucleic Acids.  Biophys J., 63, 751-
  759.

  Questions about this software should be directed to:

              ndbadmin@ndbserver.rutgers.edu

  PURPOSE:    DDL 2.1 compliant CIF file lexer ...

  Modified to allow multiple lexers in single program, support flex, etc.
    by D.A. Case,  12/97
  Further modified to recognize some NMR_STAR syntax used for NEF-style 
    restraints -- DAC, 11/14
  Make use of name-prefix and file-prefix bison options, 

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cifparse.h"
#include "cifp.tab.h"

/* 
 *  Local prototypes
 */
static void alt_cifpmore();
static void delete_punctuation(char*, char*, int*);
static char tbuf[MAXVALUELENGTH+1];
/*
 *  Local static variables ...
 */

static int   isText = FALSE;
static int   i, j;
static int   buflen;

%}
%%
[\n];|[\n]    {                                       /* New Lines */
          lineNo++;
          if (cifptext[cifpleng-1] == ';') {
             if (isText == TRUE) {          /* end of text value */
   	        strncat(tbuf,cifptext,MAXVALUELENGTH);
		buflen = strlen(tbuf);
	        delete_punctuation(cifplval.TempBuffer,tbuf,&buflen);
                isText = FALSE;
                return(VALUE);
	     } else {  /* text value begins */
	       isText = TRUE;
	       strcpy(tbuf,"");
	     }
	  }
          if (isText == TRUE) alt_cifpmore();
       }

^;|^;([^;])+      {                      /* multi-line text  value */    
      if (cifptext[0] == ';' && cifptext[cifpleng-1] != '\n') {  /* text value begins*/
	 isText = TRUE;
	 strncpy(tbuf,cifptext,MAXVALUELENGTH);
      }
      else {  
          isText = FALSE;
	  strncat(tbuf,cifptext,MAXVALUELENGTH);
	  buflen = strlen(tbuf);
          for (i=0; i<buflen; i++)
	       if (tbuf[i] == '\n') lineNo++;

	  delete_punctuation(cifplval.TempBuffer,tbuf,&buflen); 
          return(VALUE);
        } 
    }

[\t ]+      {                                    /* white space */ 
      if (isText) alt_cifpmore();
       else {
	 for (i=0; i<cifpleng; i++)
	   if (cifptext[i] == '\n') lineNo++;
       }
    }

[Dd][Aa][Tt][Aa][_][^\n\t ]*   {                 /* data_ */ 
      if (isText) alt_cifpmore();
      else {
        strncpy(cifplval.TempBuffer,cifptext,MAXVALUELENGTH);
        return (DATABLOCK);
      }
    }

[Ss][Aa][Vv][Ee][_][^\n\t ]*   {                 /* save_ */ 
      if (isText) alt_cifpmore();
      else {
        strncpy(cifplval.TempBuffer,cifptext,MAXVALUELENGTH);
        return (DATABLOCK);
      }
    }

[Ss][Aa][Vv][Ee][_][\n\t ]*   {                 /* end of save_ */ 
      if (isText) alt_cifpmore();
      else lineNo++;
    }

[Ss][Tt][Oo][Pp][_][\n\t ]*   {                 /* end of loop_ */ 
      if (isText) alt_cifpmore();
      else lineNo++;
    }

[Ll][Oo][Oo][Pp][_]      {                        /* loop_ */
      if (isText) alt_cifpmore();
      else return (LOOP);
    }

[?]|[.]  {                                        /* Unknown */
      if (isText) alt_cifpmore();
      else return (UNKNOWN);
    }

[#][^\n]*[\n]       {                             /* COMMENT */
      if (isText) alt_cifpmore();
      else lineNo++;
    }

[_][^\n \t]+  {                /* Item names start with _ end with spaces*/
      if (isText) {
         alt_cifpmore();
      } else {
        /* If the beginning of text is in buffer cifptext */
         for (i=0; i<cifpleng; i++) 
           if (cifptext[i] == '_') break;
	 strncpy(cifplval.TempBuffer,cifptext,MAXVALUELENGTH);
	 return(ITEMNAME);
      }
   }


[^\n\t _][^\n \t]*    {  /* unquoted or mixed quote strings strings */
     if (!isText) {
        j=0;
        for (i=cifpleng-1; i >= 0; i--) {
            if ( cifptext[i] == '\'' || cifptext[i] == '\"') {
               cifptext[i]='\0';
               break;
            } else 
             break;
        }
        for (i=0; i < cifpleng; i++) {
             if (cifptext[i] == '\'' || cifptext[i] == '\"') {
                j++;
                break;
            } else
               break;
        }
        strncpy(cifplval.TempBuffer, &cifptext[j],MAXVALUELENGTH);
	return(VALUE);
     }
     else {
       alt_cifpmore();
     }
   }


[\'][^'\n]*[\'\n][ \t\n]   {  /* Single Quoted Strings */

     if (!isText) {

        j=0;
        for (i=cifpleng-1; i >= 0; i--) {
            if ( cifptext[i] == ' ' || cifptext[i] == '\t' ||  cifptext[i] == '\n')
                 cifptext[i]='\0';
            else if ( cifptext[i] == '\'') {
               cifptext[i]='\0';
               break;
            } else 
             break;
        }
        for (i=0; i < cifpleng; i++) {
            if ( cifptext[i] == ' ' || cifptext[i] == '\t') 
                j++;
             else if (cifptext[i] == '\'') {
                j++;
                break;
            } else
               break;
        }
        strncpy(cifplval.TempBuffer, &cifptext[j],MAXVALUELENGTH);
	return(VALUE);
     }
     else {

        alt_cifpmore();
        if (cifptext[cifpleng-1] == '\n') unput('\n');
     }
   }


[\"][^"\n]*[\"\n][ \t\n]    {  /* Double Quoted Strings */

     if (!isText) {
        j=0;
        for (i=cifpleng-1; i >= 0; i--) {
            if ( cifptext[i] == ' ' || cifptext[i] == '\t' ||  cifptext[i] == '\n')
                 cifptext[i]='\0';
            else if ( cifptext[i] == '\"') {
               cifptext[i]='\0';
               break;
            } else 
             break;
        }
        for (i=0; i < cifpleng; i++) {
            if ( cifptext[i] == ' ' || cifptext[i] == '\t') 
                j++;
             else if (cifptext[i] == '\"') {
                j++;
                break;
            } else
               break;
        }
        strncpy(cifplval.TempBuffer, &cifptext[j], MAXVALUELENGTH);
	return(VALUE);
     }
     else {
       alt_cifpmore();
        if (cifptext[cifpleng-1] == '\n') unput('\n');
     }
   }

%%
/*---------------------------------------------------------------------------
 *                     delete_punctuation
 * Purpose:                                                                  
 *         Deletes leading and trailing semicolons and quotes on text strings.
 * Input :                                                                   
 *         char *ibuf:   The string which cifplex finds                  
 *         char *obuf:   unquoted string
 *         int  *bufLen: The sting length of ibuf                   
 *---------------------------------------------------------------------------*/
static void delete_punctuation(char *obuf, char *ibuf, int *bufLen)
{
  int i;
  if (*bufLen == 0) return;
  if (*bufLen > MAXVALUELENGTH) {
    fprintf( stderr,"String value exceeds maximum length %d > %d",*bufLen,MAXVALUELENGTH);
    fprintf( stderr,"You must reset the CIFPARSE definition of MAXVALUELENGTH greater than %d to process this file",*bufLen);
    exit(1);
  }

  for (i=0; i<*bufLen; i++)
    if (ibuf[i] == ';') break;
  strcpy(obuf, &ibuf[i+1]);
  (*bufLen)-=i;
  
  for (i= *bufLen-1; i >=0; i--)
    if (obuf[i] == ';') break;
  if (obuf[i-1] == '\n') i--;			   
  (*bufLen)-= (*bufLen-i);
  obuf[*bufLen] = '\0';
}

static void alt_cifpmore()
{
   buflen = strlen(tbuf);
   if (buflen + strlen(cifptext) > MAXVALUELENGTH) {
     fprintf( stderr,"Exceeded maximum buffer size %d starting at line %d",MAXVALUELENGTH,lineNo);
     fprintf( stderr,"Increase the MAXVALUELENGTH to read this file");
   } else 
     strcat(tbuf, cifptext);
}

