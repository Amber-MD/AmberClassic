/*
 *      File: sysdepend.c
 *
 ************************************************************************
 *                            LEAP                                      *
 *                                                                      *
 *                   Copyright (c) 1992, 1995                           *
 *           Regents of the University of California                    *
 *                     All Rights Reserved.                             *
 *                                                                      *
 *  This software provided pursuant to a license agreement containing   *
 *  restrictions on its disclosure, duplication, and use. This software *
 *  contains confidential and proprietary information, and may not be   *
 *  extracted or distributed, in whole or in part, for any purpose      *
 *  whatsoever, without the express written permission of the authors.  *
 *  This notice, and the associated author list, must be attached to    *
 *  all copies, or extracts, of this software. Any additional           *
 *  restrictions set forth in the license agreement also apply to this  *
 *  software.                                                           *
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
 *              ALL system dependent code is in this file.
 */

	/*
	 * Define some macros that make determining the system 
	 */
	/*
	 * easier 
	 */

#include        <errno.h>
#include	<signal.h>

#include	<sys/types.h>

#include	<sys/stat.h>

#include	"basics.h"


#include        <stdarg.h>
#include        <stdlib.h>

/*
 *------------------------------------------------------------------
 *
 *      My own printf which checks verbosity levels
 *      and writes output to an optional log file.
 */

#define	MAXCHARSPERPRINTF	5000	/*
					 * 5000 characters max 
					 */
static char SsOutputBuffer[MAXCHARSPERPRINTF];
extern void myPrintString();



/*
 *      myPrintf
 *
 *      If the verbosity level in GiVerbosityLevel is greater than the
 *      value GiVerbosity then print the message to stdout.
 *      If the log file GFLog is defined then write the output to
 *      GFLog regardless of verbosity level, unless the verbosity level
 *      is -1.
 */

void 
myPrintf (char *fmt,...)
{
	va_list args;

	if (GfPrintStringCallback == NULL) {
		DFATAL (("Illegal print string callback!"));
	}

	va_start (args, fmt);

	vsprintf (SsOutputBuffer, fmt, args);
	myPrintString (SsOutputBuffer);

	/*
	 * Flush the file constantly so that we are sure 
	 * that if the program crashes it will contain 
	 * the session 
	 */

	if (GfLog != NULL)
		fflush (GfLog);

	va_end( args );
}

/*
 *    SysdependDirectoryList
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Get a list of names in the directory pointed to
 *      by cPPath.  Return an array of STRINGs that were MALLOCd
 */
#ifdef WIN32

void
SysdependDirectoryList( char *cPPath, STRING *saPNames[], int *iPNumber)
{
  // function derived from this example: https://msdn.microsoft.com/en-us/library/aa365200(VS.85).aspx
  // (I'd bet money that by the time anybody actually clicks the above link, it will be broken)

  // Prepare string for use with FindFile functions.  First, copy the
  // string to a buffer, then append '\*' to the directory name.
  TCHAR pathWithSuffix[MAX_PATH];


  strcpy(pathWithSuffix, cPPath);
  strcat(pathWithSuffix, "\\*");

  WIN32_FIND_DATA findData;
  memset(&findData, 0, sizeof(WIN32_FIND_DATA));

  HANDLE findHandle;

  if ((findHandle = FindFirstFile(pathWithSuffix, &findData)) == INVALID_HANDLE_VALUE) {
	VPFATALEXIT(("FindFirstFile(%s): %s\n", cPPath, strerror(errno) ));
    	*iPNumber = 0;
    	return;
  }

  int fileCount = 0;

  do
  {
	if (strcmp(findData.cFileName, ".") != 0) /*do not take into account . dir*/
	{
	  ++fileCount;
	}
  }
  while(FindNextFile(findHandle, &findData));

  // win32 has no "rewinddir" function, so we have to close and reopen the handle.
  FindClose(findHandle);
  findHandle = FindFirstFile(pathWithSuffix, &findData);

  MALLOC ((*saPNames), STRING *, sizeof (STRING) * fileCount);

  int index = 0;

  do
  {
	if (strcmp(findData.cFileName, ".") != 0) /*do not take into account . dir*/
	{
		(void) strcpy ((char*)((*saPNames)[index++]), findData.cFileName);
		++index;
	}
  }
  while(FindNextFile(findHandle, &findData));

  FindClose(findHandle);

  *iPNumber = fileCount;
}

#elif (defined SYSV || defined __i386__)

void 
SysdependDirectoryList( char *cPPath, STRING *saPNames[], int *iPNumber)
{
  struct dirent *namelist;
  DIR *dp;
  int iNumber, i;

  if ((dp = opendir(cPPath)) == NULL) {
	VPFATALEXIT(("opendir(%s): %s\n", cPPath, strerror(errno) ));
    	*iPNumber = 0;
    	return;
  }

  iNumber = 0;
  while ((namelist = readdir(dp)) != NULL) {
    if (strcmp(namelist->d_name, ".") == 0) /*do not take into account . dir*/
      continue;
    iNumber++;
  }

  rewinddir(dp);

  MALLOC ((*saPNames), STRING *, sizeof (STRING) * iNumber);

  i = 0;
  while ((namelist = readdir(dp)) != NULL) {
    if (strcmp(namelist->d_name, ".") == 0) /*do not take into account . dir*/
      continue;
    (void) strcpy ((char*)((*saPNames)[i++]), (char*)(namelist->d_name));
  }

  closedir(dp);
  *iPNumber = iNumber;
}

#else /* not SYSV == BSD */

void 
SysdependDirectoryList( char *cPPath, STRING *saPNames[], int *iPNumber)
{
  struct dirent **namelist;
  int iNumber, i;

  iNumber = scandir (cPPath, &namelist, NULL, NULL);
  if (iNumber == -1) {
	VPFATALEXIT(("scandir(%s): %s\n", cPPath, strerror(errno) ));
    	*iPNumber = 0;
    	return;
  }

  MALLOC ((*saPNames), STRING *, sizeof (STRING) * iNumber);
  for (i = 0; i < iNumber; i++) {
    (void) strcpy ((char*)((*saPNames)[i]), (char*)(namelist[i]->d_name));
    FREE(namelist[i]);
  }
  *iPNumber = iNumber;

  /*
   * Free the array that was returned 
   */
  FREE(namelist);
}

#endif


/*
 *    fsSysdependFileStatus
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the status of the file whose name is passed
 *      in (cPName).  Return the status in a FILESTATUSt 
 *      record.
 */
FILESTATUSt 
fsSysdependFileStatus (char *cPName)
{
  FILESTATUSt fsStatus;
  struct stat buf;

  fsStatus.fMode = FILEDOESNOTEXIST;
  fsStatus.iSize = -1;
  if (!stat (cPName, &buf)) {

    if (buf.st_mode & S_IFDIR)
      fsStatus.fMode |= FILEDIRECTORY;
    if (buf.st_mode & S_IFREG)
      fsStatus.fMode |= FILENORMAL;

    fsStatus.iSize = buf.st_size;
  }
  return(fsStatus);
}



/*
 *    SysdependCurrentWorkingDirectory
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the current working directory path
 *      in sPath.
 */
void 
SysdependCurrentWorkingDirectory (STRING sPath)
{
  char caPath[MAXPATHLEN];
  char *cPResult;

#ifdef	NeXT
  cPResult = (char *) getwd (caPath);
#elif defined(WIN32)
  cPResult = _getcwd (caPath, sizeof(caPath));
#else
  cPResult = getcwd (caPath, sizeof(caPath));
#endif

  if (cPResult == NULL) {
    DFATAL (("Could not get working path\n"));
  }
  if (strlen (caPath) > sizeof (STRING)) {
    DFATAL (("Path name %s is too long!", caPath));
  }
  strcpy (sPath, caPath);
  strcat (sPath, "/");
}




