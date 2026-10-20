/* 
 *      File:   help.h
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
 *              This object manages the database of help information.
 *              Help is stored as a DICTIONARY of keywords which
 *              are associated with very large strings which contain
 *              the help information for that keyword.
 *              There is only one HELP object per system.
 */


#ifndef	HELP_H
#define	HELP_H	

typedef struct  {
	char	*sUpSubject;
	char	*sSubject;
	char	*sText;
} HELPt;

typedef HELPt*  HELP;



/*
 *      Functions
 */

extern void	HelpInitialize(void);
extern void	HelpShutdown(void);
extern void     HTInit(void);
extern void	HelpLoop(void);
extern HELP	hHelpNext(void);

                /* sHelpText returns pointers to VERY VERY LONG STRINGS */

extern HELP	hHelp(char *sSubject);

#define sHelpSubject(h)         ((HELP)(h))->sSubject
#define sHelpUpSubject(h)       ((HELP)(h))->sUpSubject
#define sHelpText(h)            ((HELP)(h))->sText

                /* HelpAdd should not be called by the user */
                /* It should only be called from within helptext.c */
extern void	HelpAdd(char *sUpSubject, char *sText);


#endif
