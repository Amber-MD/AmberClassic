#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "Matrix.h"
#include "Topology.h"
#include "Constants.h"
#include "ParamRead.h"
#include "ConfigSamp.h"
#include "Peptide.h"
#include "Parse.h"
#include "Macros.h"
#include "mdgxVector.h"
#include "ptrajmask.h"

//-----------------------------------------------------------------------------
// CountWords: this function counts the number of words appearing on a line.
//             To store all words, sequentially and individually, use the
//             ParseWords function in this library.
//
// Arguments:
//   line:     the line of interest, containing spaces or tabs to delimit words
//-----------------------------------------------------------------------------
int CountWords(char* line)
{
  int i, nword, slen, onword, inquotes;

  nword = 0;
  slen = strlen(line);
  onword = 0;
  inquotes = 0;
  for (i = 0; i < slen; i++) {
    if (line[i] == 34 || line[i] == 39) {
      inquotes = 1 - inquotes;
    }
    if (line[i] != ' ' && line[i] != 9 && onword == 0 && line[i] != '\n') {
      onword = 1;
      nword++;
    }
    else if (inquotes == 0 && (line[i] == ' ' || line[i] == 9) &&
             onword == 1) {
      onword = 0;
    }
  }

  return nword;
}

//-----------------------------------------------------------------------------
// ToUpper: function for converting a character to upper case           
//
// Arguments:
//   c:     the character to convert
//-----------------------------------------------------------------------------
char ToUpper(char c)
{
  int buffer;

  if (c >= 'a' && c <= 'z') {
    buffer = c - 'a';
    return 'A' + buffer;
  }
  else {
    return c;
  }
}

//-----------------------------------------------------------------------------
// ToLower: function for converting a character to lower case           
//
// Arguments:
//   c:     the character to convert
//-----------------------------------------------------------------------------
char ToLower(char c)
{
  int buffer;

  if (c >= 'A' && c <= 'Z') {
    buffer = c - 'A';
    return 'a' + buffer;
  }
  else {
    return c;
  }
}

//-----------------------------------------------------------------------------
// StringToUpper: function for converting an entire string to upper case
//
// Arguments:
//   s:       the string to convert
//-----------------------------------------------------------------------------
void StringToUpper(char* s)
{
  int i, slen;

  slen = strlen(s);
  for (i = 0; i < slen; i++) {
    s[i] = ToUpper(s[i]);
  }
}

//-----------------------------------------------------------------------------
// StringToLower: function for converting an entire string to upper case
//
// Arguments:
//   s:       the string to convert
//-----------------------------------------------------------------------------
void StringToLower(char* s)
{
  int i, slen;

  slen = strlen(s);
  for (i = 0; i < slen; i++) {
    s[i] = ToLower(s[i]);
  }
}

//-----------------------------------------------------------------------------
// Type4Char: fill a type definition string up to 4 characters, padding it with
//            blank spaces as necessary.
//
// Arguments:
//   a:     the type string to write
//   t:     the type string input
//-----------------------------------------------------------------------------
void Type4Char(char* a, char* t)
{
  int i, pivot;

  strncpy(a, t, 4);
  pivot = 0;
  for (i = 0; i < 4; i++) {
    if (a[i] == '\0') {
      pivot = 1;
    }
    if (pivot == 1) {
      a[i] = ' ';
    }
  }
  a[4] = '\0';
}

//-----------------------------------------------------------------------------
// CaselessStrcmp: function that works like strcmp, but returns 0 for cases in
//                 which the two strings are identical if capitalization is
//                 not an issue.
//
// Arguments:
//   s1, s2:    the two strings to compare
//   nc:        the number of characters from each string to compare.  If set
//              to -1, all characters will be compared and the function will
//              return 1 (error) if the strings are of different lengths.  If
//              set zero, all characters up to the length of the shorter
//              string will be compared.  If set to any positive number, those
//              nc characters will be compared, assuming both strings are long
//              enough (as many characters as possible will be compared if nc
//              is larger than the length of either string).
//-----------------------------------------------------------------------------
int CaselessStrcmp(char* s1, char* s2, int nc)
{
  int i, ls1, ls2, minl, decision;
  char* s1work;
  char* s2work;

  ls1 = strlen(s1);
  ls2 = strlen(s2);
  s1work = (char*)malloc((ls1+1)*sizeof(char));
  s2work = (char*)malloc((ls2+1)*sizeof(char));
  strcpy(s1work, s1);
  strcpy(s2work, s2);
  StringToUpper(s1work);
  StringToUpper(s2work);
  decision = 0;
  if (nc == -1) {
    if (ls1 != ls2) {
      decision = 1;
    }
    else {
      for (i = 0; i < ls1; i++) {
        if (s1work[i] != s2work[i]) {
          decision = 1;
        }
      }
    }
  }
  else if (nc == 0) {
    minl = (ls1 < ls2) ? ls1 : ls2;
    for (i = 0; i < minl; i++) {
      if (s1work[i] != s2work[i]) {
        decision = 1;
      }
    }
  }
  else {
    minl = (ls1 < ls2) ? ls1 : ls2;
    if (nc > minl) {
      nc = minl;
    }
    for (i = 0; i < nc; i++) {
      if (s1work[i] != s2work[i]) {
        decision = 1;
      }
    }
  }

  // Free allocated memory
  free(s1work);
  free(s2work);

  return decision;
}

//-----------------------------------------------------------------------------
// RealXpYf: function to convert a character string, which is assumed to be a
//           number in decimal format, to a double-precision real.  Much faster
//           than scanf.
//
// Arguments:
//   word:    the word containing a formatted decimal number
//   X:       the number of characters expected in the number (including any
//            preceding whitespace)
//   Y:       the number of characters expected after the decimal
//-----------------------------------------------------------------------------
double RealXpYf(char* word, int X, int Y)
{
  int i, started, xym1;
  double mult, result, msign;
  char* buff;

  // One sanity check, try to salvage a number if it fails
  xym1 = X-Y-1;
  if (word[xym1] != '.') {
    buff = (char*)malloc(X+1*sizeof(char));
    strncpy(buff, word, X); 
    buff[X] = '\0';
    sscanf(buff, "%lf", &result);
    free(buff);
    return result;
  }

  // Loop over all digits
  mult = pow(10.0, X-Y-2);
  started = 0;
  msign = 1.0;
  result = 0.0;
  for (i = 0; i < X; i++) {
    if (word[i] >= '0' && word[i] <= '9') {
      result += mult*(word[i] - '0');
      started = 1;
    }
    else if (word[i] == '-') {
      if (started == 0) {
        msign = -1.0;
        started = 1;
      }
      else {
        printf("RealXpYf >> %s is a nonsensical number.\n", word);
        exit(1);
      }
    }
    else if (i == xym1) {
      continue;
    }
    else if (word[i] == ' ' && started == 1) {
      if (i < xym1) {
        printf("RealXpYf >> %s breaks the requested format.\n", word);
        exit(1);
      }
      break;
    }
    else if (word[i] != ' ') {
      buff = (char*)malloc(X*sizeof(char));
      strncpy(buff, word, X);
      sscanf(buff, "%lf", &result);
      free(buff);
      printf("RealXpYf >> %s is a nonsensical number.\n", word);
      exit(1);
    }
    mult *= 0.1;
  }

  return result*msign;
}

//-----------------------------------------------------------------------------
// WordIsInteger: function to determine whether a word really looks like an
//                integer, and then return 0 if false or 1 if true.
//
// Arguments:
//   word:       the word to analyze
//-----------------------------------------------------------------------------
int WordIsInteger(char* word)
{
  int i, wlen, onint, nnum;

  wlen = strlen(word);
  onint = 0;
  nnum = 0;
  for (i = 0; i < wlen; i++) {
    if ((word[i] < '0' || word[i] > '9') && word[i] != '-' && word[i] != '+' &&
	word[i] != ' ') {
      return 0;
    }
    if (onint == 0 && ((word[i] >= '0' && word[i] <= '9') ||
                       word[i] == '-' || word[i] == '+')) {
      nnum++;
      onint = 1;
      continue;
    }
    if (onint == 1 && word[i] == ' ') {
      onint = 0;
    }
  }
  if (nnum != 1) {
    return 0;
  }

  return 1;
}

//-----------------------------------------------------------------------------
// WordIsNumber: function to determine whether a word really looks like a
//               number, and then return 0 if false or 1 if true.     
//
// Arguments:
//   word:       the word to analyze
//-----------------------------------------------------------------------------
int WordIsNumber(char* word)
{
  int i, wlen, hitdot;

  if ((word[0] < '0' || word[0] > '9') &&
      word[0] != '-' && word[0] != '+' && word[0] != '.') {
    return 0;
  }
  wlen = strlen(word);
  hitdot = (word[0] == '.') ? 1 : 0;
  for (i = 1; i < wlen; i++) {
    if ((word[i] < '0' || word[i] > '9') && word[i] != '.' &&
        word[i] != 'E' && word[i] != 'e' && word[i] != 'D' &&
        word[i] != 'd') {
      return 0;
    }
    if (word[i] == 'D' || word[i] == 'd') {
      word[i] = 'E';
    }
    if (word[i] == 'E' || word[i] == 'e') {
      if (i < wlen-2 && (word[i+1] == '-' || word[i+1] == '+')) {
        i++;
        continue;
      }
      else {
        return 0;
      }
    }
    if (word[i] == '.') {
      if (hitdot == 1) {
        return 0;
      }
      else {
        hitdot = 1;
      }
    }
  }

  return 1;
}

//-----------------------------------------------------------------------------
// WordIsAtomType: function to determine whether a word really looks like an
//                 atom type, returning 0 if false or 1 if true.
//
// Arguments:
//   word:       the word to analyze
//-----------------------------------------------------------------------------
int WordIsAtomType(char* word)
{
  int wlen;

  wlen = strlen(word);
  if (wlen > 2 || wlen == 0) {
    return 0;
  }
  if (!(word[0] >= 'A' && word[0] <= 'Z') &&
      !(word[0] >= 'a' && word[0] <= 'z') &&
      !(word[0] >= '0' && word[0] <= '9')) {
    return 0;
  }
  if (!(word[1] >= 'A' && word[1] <= 'Z') && word[1] != '#' &&
      !(word[1] >= 'a' && word[1] <= 'z') && word[1] != '@' &&
      !(word[1] >= '0' && word[1] <= '9') && word[1] != '*' &&
      word[1] != ' ' && word[1] != '\0') {
    return 0;
  }

  return 1;
}

//-----------------------------------------------------------------------------
// ParseWords: this function stores all words appearing on a line.      
//
// Arguments:
//   line:       the line to decompose into words as delimited by whitespace
//               or tab (ascii text value 9) characters
//-----------------------------------------------------------------------------
cmat ParseWords(char* line)
{
  int i, j, nword, onword, inquotes, nc;
  cmat CW;

  nword = CountWords(line);
  if (nword == 0) {
    CW = CreateCmat(1, MAXNAME);
  }
  else {
    CW = CreateCmat(nword, MAXNAME);
  }
  i = 0;
  j = 0;
  while (i < nword) {

    // Advance to the next word
    while (line[j] == ' ' || line[j] == 9) {
      j++;
    }

    // Read this word
    onword = 1;
    if (line[j] == 34 || line[j] == 39) {
      inquotes = 1;
      j++;
    }
    else {
      inquotes = 0;
    }
    nc = 0;
    while (onword == 1) {
      if (inquotes == 0 && (line[j] == ' ' || line[j] == 9 ||
                            line[j] == '\n' || line[j] == '\0')) {
        onword = 0;
        CW.map[i][nc] = '\0';
      }
      else if (inquotes == 1 && (line[j] == 34 || line[j] == 39 || 
                                 line[j] == '\0')) {
        onword = 0;
        inquotes = 0;
        CW.map[i][nc] = '\0';
      }
      else {
        CW.map[i][nc] = line[j];
        nc++;
      }
      j++;
    }
    i++;
  }

  return CW;
}

//-----------------------------------------------------------------------------
// LocateWordInLine: find the first character of a word within a line, and
//                   verify that it is followed by white space.  Return the
//                   index of the character in the line / string, or -1 if it
//                   cannot be found.
//
//
// Arguments:
//   w:       the word to search for
//   s:       the line to search
//   minpos:  the starting place within the string to search
//-----------------------------------------------------------------------------
int LocateWordInLine(char* w, char* s, int minpos)
{
  int i, j, maxpos;
  
  if (minpos >= strlen(s)) {
    return -1;
  }
  if (minpos < 0) {
    minpos = 0;
  }
  maxpos = strlen(s) - strlen(w);
  for (i = minpos; i < maxpos; i++) {
    if (s[i] == w[0] && strncmp(w, &s[i], strlen(w)) == 0) {
      return i;
    }
  }
  
  return -1;
}


//-----------------------------------------------------------------------------
// RemoveWhiteSpace: function for removing white space at the beginning 
//                   (or end) of a string.                              
//                                                                      
// Arguments:                                                           
//   a:     the string                                                  
//   tail:  set to 0 to remove whitespace at the head, 1 at the tail    
//   ws:    define white space (could be ' ' or '0', for example)       
//-----------------------------------------------------------------------------
void RemoveWhiteSpace(char* a, int tail, char ws)
{
  int i, j, asize, num_rem;

  asize = strlen(a);
  if (tail == 0) {
    num_rem = 0;
    while (a[num_rem] == ws && num_rem < asize) {
      num_rem++;
    }
    j = 0;
    for (i = num_rem; i < asize; i++) {
      a[j] = a[i];
      j++;
    }
    a[j] = '\0';
  }
  else {
    i = strlen(a)-1;
    while (a[i] == ws) {
      a[i] = '\0';
      i--;
    }
  }
}

//-----------------------------------------------------------------------------
// EqualSpace: add white space before and after any equal symbols in a  
//             line UNLESS they occur within ('') or ("") quotes.       
//                                                                      
// Arguments:                                                           
//   line:    the line (from an input command file)                     
//-----------------------------------------------------------------------------
void EqualSpace(char* line)
{
  int i, j, slen, inquotes;
  slen = strlen(line);

  inquotes = 0;
  for (i = 0; i < slen; i++) {
    if (line[i] == '=' && inquotes == 0) {
      line[slen+2] = '\0';
      for (j = slen-1; j > i; j--) {
        line[j+2] = line[j];
      }
      line[i] = ' ';
      line[i+1] = '=';
      line[i+2] = ' ';
      slen += 2;
      i += 2;
    }

    // ASCII standard character set used to identify quotation symbols
    if (line[i] == 34 || line[i] == 39) {
      inquotes = 1 - inquotes;
    }
  }
}

//-----------------------------------------------------------------------------
// RemoveComments: this function removes commented text from a line.  Comments
//                 are denoted by the appearance of '%', '#', or '$' symbols.
//
// Arguments:
//   line:         the line to expunge comments from
//-----------------------------------------------------------------------------
void RemoveComments(char* line)
{
  int i, inquotes;
  const int slen = strlen(line);

  inquotes = 0;
  for (i = 0; i < slen; i++) {
    if (inquotes == 0 &&
        (line[i] == '$' || line[i] == '%' || line[i] == '#')) {
      line[i] = '\0';
      break;
    }

    // ASCII standard character set used to identify quotation symbols
    if (line[i] == 34 || line[i] == 39) {
      inquotes = 1 - inquotes;
    }
  }
}

//-----------------------------------------------------------------------------
// NixCommaCarriage: replaces all commas and carriage returns on a line 
//                   with white space.                                  
//
// Arguments:
//   line:         the line to expunge commas and carriage returns from
//-----------------------------------------------------------------------------
void NixCommaCarriage(char* line)
{
  int i, inquotes;
  const int slen = strlen(line);

  inquotes = 0;
  for (i = 0; i < slen; i++) {
    if (inquotes == 0 && (line[i] == ',' || line[i] == '\n')) {
      line[i] = ' ';
    }

    // ASCII standard character set used to identify quotation symbols
    if (line[i] == 34 || line[i] == 39) {
      inquotes = 1 - inquotes;
    }
  }
}

//-----------------------------------------------------------------------------
// AddToGlossary: add to a growing glossary of terms.  This function first
//                searches for how many terms are in the glossary, then adds
//                the specified number of new terms.  This function makes use
//                of a variable argument list.
//
// Arguments:
//   glossary:   the list of terms to build upon
//   nargs:      the number of additional arguments to take in
//   <args>:     nargs number of character strings
//-----------------------------------------------------------------------------
void AddToGlossary(cmat *glossary, int nargs, ...)
{
  int i, j, nterm, found;
  char buffer[32];
  va_list ap;

  // Count the number of words already in the glossary
  nterm = 0;
  for (i = 0; i < glossary->row; i++) {
    if (glossary->map[i][0] != '\0') {
      nterm++;
    }
  }
  if (nterm == glossary->row) {
    *glossary = ReallocCmat(glossary, glossary->row + 32, glossary->col);
  }

  va_start(ap, nargs);
  for (i = 0; i < nargs; i++) {
    strcpy(buffer, va_arg(ap, const char*));
    found = 0;
    for (j = 0; j < nterm; j++) {
      if (strcmp(buffer, glossary->map[j]) == 0) {
        found = 1;
      }
    }
    if (found == 0) {
      strcpy(glossary->map[nterm], buffer);
      nterm++;
      if (nterm == glossary->row) {
        *glossary = ReallocCmat(glossary, glossary->row + 32, glossary->col);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// WordInGlossary: function to determine whether a word is found in a glossary
//                 of many terms.  Returns 1 if true, zero if false.
//
// Arguments:
//   glossary:   the list of terms to match against
//   word:       the word to find in the glossary
//-----------------------------------------------------------------------------
int WordInGlossary(cmat *glossary, char* word)
{
  int i;

  for (i = 0; i < glossary->row; i++) {
    if (strcmp(word, glossary->map[i]) == 0) {
      return 1;
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// AdvanceToSegment: this function searches a command file from the current
//                   file pointer (or, from the beginning, if specified) and
//                   stops after finding the heading for the section of
//                   interest.                           
//                                                                      
// Arguments:                                                           
//   inp:     the input file                                            
//   segname: the name of the target segment                            
//   scan0:   if set to 1, causes the file to be rewound before beginning the
//            search; if set to 2, causes the segment name to be searched
//            directly, with no leading '&' symbol, and also prompts a file
//            rewind                                
//-----------------------------------------------------------------------------
int AdvanceToSegment(FILE *inp, char* segname, int scan0)
{
  int slen, collect;
  char line[MAXLINE];

  // Start by rewinding the file
  slen = strlen(segname);
  if (scan0 == 1 || scan0 == 2) {
    rewind(inp);
  }
  collect = 0;
  while (collect == 0) {
    if (fgets(line, MAXLINE, inp) == NULL) {
      break;
    }
    RemoveWhiteSpace(line, 0, ' ');
    if (scan0 < 2) {
      if (line[0] == '&' && strncmp(&line[1], segname, slen) == 0) {
        collect = 1;
      }
    }
    else if (scan0 == 2) {
      if (strncmp(&line[0], segname, slen) == 0) {
        collect = 1;
      }
    }
  }

  return collect;
}

//-----------------------------------------------------------------------------
// DetectNamelistEnd: detect the end of a namelist, signified by a line 
//                    containing "&end."                                
//                                                                      
// Arguments:                                                           
//   line:    the line (from an input command file)                     
//   errmsg:  the error message to display (the name of the function in which
//            which the call to DetectNamelistEnd originated)         
//-----------------------------------------------------------------------------
int DetectNamelistEnd(char* line, char* errmsg)
{
  if (line[0] == '&') {
    if (strncmp(&line[1], "end", 3) == 0) {
      return 0;
    }
    else {
      printf("%s >> Error.  New segment encountered before termination of "
             "%s >> current segment.\n", errmsg, errmsg);
      exit(1);
    }
  }
  else if (line[0] == '/' && (line[1] == ' ' || line[1] == '\0')) {
    return 0;
  }

  return 1;
}

//-----------------------------------------------------------------------------
// ReadNamelistLine: once a namelist has been detected, this function will read
//                   one line from the input file, determine whether a namelist
//                   end has been encountered, and return the processed line.
//                                                                      
// Arguments:                                                           
//   line:      character string, allocated buffer                      
//   lwords:    pointer to character matrix, allocated in this function 
//   callfunc:  string indicating the function that called this reader  
//   inp:       the input file                                          
//-----------------------------------------------------------------------------
int ReadNamelistLine(char* line, cmat *lwords, char* callfunc, FILE *inp)
{
  // Read the next line
  fgets(line, MAXLINE, inp);
  RemoveWhiteSpace(line, 0, ' ');
  RemoveComments(line);

  // Break if the line is "&end"
  if (DetectNamelistEnd(line, callfunc) == 0) {
    return 0;
  }

  // Eliminate and add spaces between special characters
  // "=", "\n", and ","                                 
  NixCommaCarriage(line);
  EqualSpace(line);
  *lwords = ParseWords(line);

  return 1;
}

//-----------------------------------------------------------------------------
// SeekString: find a tag for a particular string.                      
//                                                                      
// Arguments:                                                           
//   L:       the list of words on this line of the input command file  
//   val:     the value to set if the string flag is found              
//   sname:   the string flag being sought                              
//   salias:  alias for the string flag being sought                    
//-----------------------------------------------------------------------------
void SeekString(cmat L, char* val, char* sname, char* salias)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(val, L.map[i+1]);
      }
      else if (i+2 < L.row) {
        strcpy(val, L.map[i+2]);
      }
      else {
        printf("SeekString >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekMultiString: after finding a keyword (tag), continue to read strings
//                  that follow until another keyword is found.
//
// Arguments:
//   L:        the list of words on this line of the input command file  
//   val:      the values to set if the string flag is found (this matrix will
//             be expanded if necessary by this function)
//   sname:    the string flag being sought                              
//   salias:   alias for the string flag being sought                    
//   glossary: a list of keywords to match against, to keep this open-ended
//             declaration from bleeding into other declarations within the
//             same namelist
//   append:   flag to have this multi-string append to a growing list of
//             values rather than just creating a new one from the last such
//             instance of the keyword it found
//-----------------------------------------------------------------------------
void SeekMultiString(cmat L, cmat *val, char* sname, char* salias,
                     cmat *glossary, int append)
{
  int i, j, nw;

  if (append == 1) {
    nw = 0;
    for (i = 0; i < val->row; i++) {
      if (val->map[i][0] != '\0') {
        nw++;
      }
    }
    *val = ReallocCmat(val, nw+32, val->col);
  }
  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+1 == L.row || WordInGlossary(glossary, L.map[i+1]) == 1) {
        printf("SeekMultiString >> Error.  No value specified for identifier "
               "%s/%s.\n", sname, salias);
        exit(1);
      }
      j = (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) ? i+1 : i+2;
      if (append == 0) {
        nw = 0;
      }
      while (j < L.row && WordInGlossary(glossary, L.map[j]) == 0) {
        strcpy(val->map[nw], L.map[j]);
        nw++;
        if (nw == val->row) {
          *val = ReallocCmat(val, val->row+32, val->col);
        }
        j++;
      }
      *val = ReallocCmat(val, nw, val->col);
    }
  }
  if (append == 1) {
    nw = 0;
    for (i = 0; i < val->row; i++) {
      if (val->map[i][0] != '\0') {
        nw++;
      }
    }
    if (nw == 0) {
      nw = 1;
    }
    *val = ReallocCmat(val, nw, val->col);
  }
}

//-----------------------------------------------------------------------------
// SeekString: find a tag for a particular string.                      
//                                                                      
// Arguments:                                                           
//   L:       the list of words on this line of the input command file  
//   val:     the growing list of values set by finding the appropriate flag
//   sname:   the string flag being sought                              
//   salias:  alias for the string flag being sought                    
//   counter: the counter variable (gets incremented)                   
//-----------------------------------------------------------------------------
void SeekStringInc(cmat L, cmat *val, char* sname, char* salias, int *counter)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(val->map[*counter], L.map[i+1]);
        *counter += 1;
      }
      else if (i+2 < L.row) {
        strcpy(val->map[*counter], L.map[i+2]);
        *counter += 1;
      }
      else {
        printf("SeekStringInc >> Error.  No value specified for identifier "
               "%s/%s.\n", sname, salias);
      }
      if (*counter == val->row) {
        *val = ReallocCmat(val, 2*val->row, val->col);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekSSR: find a tag for a particular string triplet, and increment a counter
//          variable each time the tag is found.                
//
// Arguments:                                                           
//   L:       the list of words on this line of the input command file  
//   val1:    the first value to set if the string flag is found        
//   val2:    the second value to set if the string flag is found       
//   val3:    the third (numerical) value to set if the flag is found   
//   sname:   the string flag being sought                              
//   salias:  alias for the string flag being sought                    
//   counter: the counter variable (gets incremented)                   
//-----------------------------------------------------------------------------
void SeekSSR(cmat L, char* val1, char* val2, double *val3, char* sname,
             char* salias, int *counter)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+3 < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(val1, L.map[i+1]);
        strcpy(val2, L.map[i+2]);
        *val3 = atof(L.map[i+3]);
        *counter += 1;
      }
      else if (i+4 < L.row) {
        strcpy(val1, L.map[i+2]);
        strcpy(val2, L.map[i+3]);
        *val3 = atof(L.map[i+4]);
        *counter += 1;
      }
      else {
        printf("SeekSSR >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekS3R: find a tag for a particular string triplet, and increment a counter
//          variable each time the tag is found.                
//                                                                      
// Arguments:                                                           
//   L:       the list of words on this line of the input command file  
//   val1:    the first value to set if the string flag is found        
//   val2:    the second value to set if the string flag is found       
//   val3:    the third value to set if the string flag is found        
//   val4:    the fourth (numerical) value to set if the flag is found  
//   sname:   the string flag being sought                              
//   salias:  alias for the string flag being sought                    
//   counter: the counter variable (gets incremented)                   
//-----------------------------------------------------------------------------
void SeekS3R(cmat L, char* val1, char* val2, char* val3, double *val4,
             char* sname, char* salias, int *counter)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+4 < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(val1, L.map[i+1]);
        strcpy(val2, L.map[i+2]);
        strcpy(val3, L.map[i+3]);
        *val4 = atof(L.map[i+4]);
        *counter += 1;
      }
      else if (i+5 < L.row) {
        strcpy(val1, L.map[i+2]);
        strcpy(val2, L.map[i+3]);
        strcpy(val3, L.map[i+4]);
        *val4 = atof(L.map[i+5]);
        *counter += 1;
      }
      else {
        printf("SeekS3R >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// TouchupAtomName: fill out an atom name so that it has precisely four 
//                  characters by adding whitespace where necessary.    
//                                                                      
// Arguments:                                                           
//   atn:     the atom name, or atom type name                          
//-----------------------------------------------------------------------------
void TouchupAtomName(char* atn)
{
  int i, slen;

  slen = strlen(atn);
  for (i = slen; i < 4; i++) {
    atn[i] = ' ';
  }
  atn[4] = '\0';
}

//-----------------------------------------------------------------------------
// SeekSinglePoint: seek a single point calculation or a series of single point
//                  calculations.  The format of this command may be simply
//                  "<keyword> <topology> <coordinates> <energy>", where
//                  keyword is either System or sys, and energy is either a
//                  real number or the name of yet another file containing a
//                  list of real numbers.  The format also accepts its own
//                  list of keywords for more detailed input.
//
// Arguments:
//   L:         the list of words on this line of the input command file
//   mp:        the master parameter development set (this routine does not
//              simply return a value like other Seek(...) functions)
//   sname:     the string flag being sought
//   salias:    alias for the string flag being sought
//   maxconf:   the maximum number of conformations (mmsys structures) that
//              have been allocated (not stored in mp)
//   glossary:  glossary of terms in the overarching namelist
//-----------------------------------------------------------------------------
void SeekSinglePoint(cmat L, prmset *mp, char* sname, char* salias,
                     int *maxconf, cmat *glossary)
{
  int i, j;

  for (i = 0; i < L.row-3; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {

      // This looks like a single point conformation.  Allocate space
      // to store the relevant information and then read it.
      mp->conf[mp->nconf].crdsrc = (char*)calloc(MAXNAME, sizeof(char));
      mp->conf[mp->nconf].tpsrc = (char*)calloc(MAXNAME, sizeof(char));
      mp->conf[mp->nconf].esrc = (char*)calloc(MAXNAME, sizeof(char));
      j = i+1;
      while (j < L.row && WordInGlossary(glossary, L.map[j]) == 0) {
        if (j < L.row-1) {
          if (strcmp(L.map[j], "-p") == 0) {
            strcpy(mp->conf[mp->nconf].tpsrc, L.map[j+1]);
          }
          if (strcmp(L.map[j], "-c") == 0) {
            strcpy(mp->conf[mp->nconf].crdsrc, L.map[j+1]);
          }
          if (strcmp(L.map[j], "-e") == 0) {
            strcpy(mp->conf[mp->nconf].esrc, L.map[j+1]);
          }
        }
        j++;
      }
      if (mp->conf[mp->nconf].tpsrc[0] == '\0' ||
          mp->conf[mp->nconf].crdsrc[0] == '\0' ||
          mp->conf[mp->nconf].esrc[0] == '\0') {

        // Not all of the critical pieces were found.  Try reading the
        // information in the simplest way possible.
        strcpy(mp->conf[mp->nconf].tpsrc, L.map[i+1]);
        strcpy(mp->conf[mp->nconf].crdsrc, L.map[i+2]);
        strcpy(mp->conf[mp->nconf].esrc, L.map[i+3]);
      }
      mp->nconf += 1;
      if (mp->nconf == *maxconf) {
        *maxconf *= 2;
        mp->conf = (mmsys*)realloc(mp->conf, (*maxconf)*sizeof(mmsys));
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekBondTermID: detect a bond, angle, or torsion identifier for
//                 optimization.  GetParamNamelist in the Command.c library
//                 uses this function, but it is placed here to keep Command.c
//                 tidy.                              
//                                                                      
// Arguments:                                                           
//   L:       the list of words on this line of the input command file  
//   mp:      the main parameter set, storing a growing list of bond,   
//            angle, and torsion IDs                                    
//   sname:   the string flag being sought                              
//   salias:  alias for the string flag being sought                    
//   maxadj:  the maximum number of adjustable terms that can be held without
//            reallocating the array (this number will be incremented and the
//            array will be expanded if necessary)  
//   order:   the order of the bonded term (as always, bond=2, angle=3, and
//            torsion=4)                                            
//-----------------------------------------------------------------------------
void SeekBondTermID(cmat L, prmset *mp, char* sname, char* salias,
                    int *maxadj, int order)
{
  int i, j, alen, blen, clen, dlen, nfield;
  xbonddef *xb;
  xangldef *xa;
  torterm *xt;
  xcmapdef *xm;
  bondrst *restb;
  anglrst *resta;
  torrst *resth;

  // Need one extra field if this is seeking a restraint weight
  if (order == 7) {
    nfield = 1;
  }
  else if (order == -7) {
    nfield = 2;
  }
  else {
    nfield = (order > 0) ? order : order + 1;
  }

  // Loop over all words
  for (i = 0; i < L.row; i++) {
    if (order == 2)       xb = &mp->badj[mp->nbadj];
    else if (order == 3)  xa = &mp->aadj[mp->naadj];
    else if (order == 4)  xt = &mp->hadj[mp->nhadj];
    else if (order == 7)  xm = &mp->madj[mp->nmadj];
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i + nfield < L.row && strcmp(L.map[i+1], "=") != 0) {
        if (order == 2) {
          Type4Char(xb->atype, L.map[i+1]);
          Type4Char(xb->btype, L.map[i+2]);
          xb->comment = (char*)malloc(MAXNAME*sizeof(char));
          mp->nbadj += 1;
        }
        else if (order == 3) {
          Type4Char(xa->atype, L.map[i+1]);
          Type4Char(xa->btype, L.map[i+2]);
          Type4Char(xa->ctype, L.map[i+3]);
          xa->comment = (char*)malloc(MAXNAME*sizeof(char));
          mp->naadj += 1;
        }
        else if (order == 4) {
          Type4Char(xt->atype, L.map[i+1]);
          Type4Char(xt->btype, L.map[i+2]);
          Type4Char(xt->ctype, L.map[i+3]);
          Type4Char(xt->dtype, L.map[i+4]);
          xt->comment = (char*)malloc(MAXNAME*sizeof(char));
          mp->nhadj += 1;
        }
	else if (order == 7) {
	  xm->resnames = CreateCmat(1, 5);
	  strncpy(xm->resnames.map[0], L.map[i+1], 4);
	  mp->nmadj += 1;
	}
      }
      else if (i + nfield + 1 < L.row) {
        if (order == 2) {
          Type4Char(xb->atype, L.map[i+2]);
          Type4Char(xb->btype, L.map[i+3]);
          xb->comment = (char*)malloc(MAXNAME*sizeof(char));
          mp->nbadj += 1;
        }
        else if (order == 3) {
          Type4Char(xa->atype, L.map[i+2]);
          Type4Char(xa->btype, L.map[i+3]);
          Type4Char(xa->ctype, L.map[i+4]);
          xa->comment = (char*)malloc(MAXNAME*sizeof(char));
          mp->naadj += 1;
        }
        else if (order == 4) {
          Type4Char(xt->atype, L.map[i+2]);
          Type4Char(xt->btype, L.map[i+3]);
          Type4Char(xt->ctype, L.map[i+4]);
          Type4Char(xt->dtype, L.map[i+5]);
          xt->comment = (char*)malloc(MAXNAME*sizeof(char));
          mp->nhadj += 1;
        }
	else if (order == 7) {
	  xm->resnames = CreateCmat(1, 5);
	  strncpy(xm->resnames.map[0], L.map[i+2], 4);
	  mp->nmadj += 1;
	}
      }
      else {
        printf("SeekTorsionID >> Error.  Invalid inputs for identifier "
               "%s/%s.\n", sname, salias);
        printf("SeekTorsionID >> This identifier requires four consecutive "
               "atom type names.\n");
        exit(1);
      }

      // Touch up the atom names
      if (order == 2) {
        TouchupAtomName(xb->atype);
        TouchupAtomName(xb->btype);
      }
      else if (order == 3) {
        TouchupAtomName(xa->atype);
        TouchupAtomName(xa->btype);
        TouchupAtomName(xa->ctype);
      }
      else if (order == 4) {
        TouchupAtomName(xt->atype);
        TouchupAtomName(xt->btype);
        TouchupAtomName(xt->ctype);
        TouchupAtomName(xt->dtype);
      }
      else if (order == 7) {
        TouchupAtomName(xm->resnames.map[0]);
        xm->nres = 1;
      }

      // Extend the adjustable terms array if needed
      if (order == 2 && mp->nbadj == *maxadj) {
        *maxadj += 32;
        mp->badj = (xbonddef*)realloc(mp->badj, (*maxadj)*sizeof(xbonddef));
      }
      else if (order == 3 && mp->naadj == *maxadj) {
        *maxadj += 32;
        mp->aadj = (xangldef*)realloc(mp->aadj, (*maxadj)*sizeof(xangldef));
      }
      else if (order == 4 && mp->nhadj == *maxadj) {
        *maxadj += 32;
        mp->hadj = (torterm*)realloc(mp->hadj, (*maxadj)*sizeof(torterm));
      }
      else if (order == 7 && mp->nmadj == *maxadj) {
	*maxadj += 32;
	mp->madj = (xcmapdef*)realloc(mp->madj, (*maxadj)*sizeof(xcmapdef));
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekRecast: detect directives to change the name of an atom type, perhaps
//             only in specific instances.                      
//                                                                      
// Arguments:                                                           
//   L:        the list of words on this line of the input command file 
//   mp:       the main parameter set                                   
//   sname:    the string flag being sought                             
//   salias:   alias for the string flag being sought                   
//   maxhold:  the maximum number of type recasts that can be held without
//             reallocating the array (this number will be incremented and the
//             array will be expanded if necessary) 
//   specinst: flag to indicate that this atom recasting takes place only in
//             specific instances (1) or if it is general (0)   
//-----------------------------------------------------------------------------
void SeekRecast(cmat L, prmset *mp, char* sname, char* salias, int *maxhold,
                int specinst)
{
  int i, j, ilmin, olen, nlen;
  char *orig, *new;
  typeswitch *xs;
  typebranch *xb;

  for (i = 0; i < L.row; i++) {
    if (specinst == 1) {
      xb = &mp->cleave[mp->ncleave];
      orig = xb->orig;
      new = xb->pnew;
    }
    else {
      xs = &mp->recast[mp->nrecast];
      orig = xs->orig;
      new = xs->pnew;
    }
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      ilmin = i + 2 + specinst;
      if (ilmin < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(orig, L.map[i+1]);
        strcpy(new, L.map[i+2]);
        if (specinst == 1) {
          strcpy(mp->cleave[mp->ncleave].instances, L.map[i+3]);
        }
      }
      else if (ilmin + 1 < L.row) {
        strcpy(orig, L.map[i+2]);
        strcpy(new, L.map[i+3]);
        if (specinst == 1) {
          strcpy(mp->cleave[mp->ncleave].instances, L.map[i+4]);
        }
      }
      else {
        printf("SeekRecast >> Error.  Invalid inputs for identifier "
               "%s/%s.\n", sname, salias);
        if (specinst == 0) {
          printf("SeekRecast >> This identifier requires two "
                 "consecutive atom type names.\n");
        }
        else {
          printf("SeekRecast >> This identifier requires two "
                 "consecutive atom type names.\nSeekRecast >> and an ambmask "
                 "string describing instances of the new type.\n");
        }
        exit(1);
      }

      // Touch up the atom type names
      TouchupAtomName(orig);
      TouchupAtomName(new);

      // Extend the array if needed
      if (specinst == 0) {
        mp->nrecast += 1;
        if (mp->nrecast == *maxhold) {
          *maxhold += 32;
          mp->recast = (typeswitch*)realloc(mp->recast,
                                            (*maxhold)*sizeof(typeswitch));
        }
      }
      else {
        mp->ncleave += 1;
        if (mp->ncleave == *maxhold) {
          *maxhold += 32;
          mp->cleave = (typebranch*)realloc(mp->recast,
                                            (*maxhold)*sizeof(typebranch));
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekSpecParmRest: seek a specific parameter restraint.  This replaces
//                   functionality originally incorporated into SeekBondTermID,
//                   building flexibility in the number of details that can be
//                   added and also the order in which they are specified.
//                                                                      
// Arguments:                                                           
//   L:         the list of words on this line of the input command file 
//   mp:        the main parameter set                                   
//   sname:     the string flag being sought                             
//   salias:    alias for the string flag being sought                   
//   order:     the order of the restraint (2=bond, 3=angle, 4=torsion)  
//   maxsr:     the maximum number of restraints that can be held        
//   glossary:  a list of keywords to match against, to keep this open-ended
//              declaration from bleeding into other declarations within the
//              same namelist
//-----------------------------------------------------------------------------
void SeekSpecParmRest(cmat L, prmset *mp, char* sname, char* salias, int order,
                      int *maxsr, cmat *glossary)
{
  int i, j, nt;
  bondrst *restb;
  anglrst *resta;
  torrst *resth;
  nmroper *restr;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) != 0 && strcmp(L.map[i], salias) != 0) {
      continue;
    }

    // This is a request for a specific parameter restraint
    if (order == 2) {
      restb = &mp->userrstB[mp->nuserrstB];
      restb->has_wK = 0;
      restb->has_wl0 = 0;
      restb->has_tK = 0;
      restb->has_tl0 = 0;
      restb->rstwK = 0.0;
      restb->rstwl0 = 0.0;
      restb->targK = 0.0;
      restb->targl0 = 0.0;
    }
    else if (order == 3) {
      resta = &mp->userrstA[mp->nuserrstA];
      resta->has_wK = 0;
      resta->has_wTh0 = 0;
      resta->has_tK = 0;
      resta->has_tTh0 = 0;
      resta->rstwK = 0.0;
      resta->rstwTh0 = 0.0;
      resta->targK = 0.0;
      resta->targTh0 = 0.0;
    }
    else if (order == 4) {
      resth = &mp->userrstH[mp->nuserrstH];
      resth->rstw = 0.0;
      resth->pn = 0.0;
      resth->target = 0.0;
    }
    else if (order == 6) {
      restr = &mp->userrstR[mp->nuserrstR];
      restr->rstw = 0.0;
      restr->target = 1.0;
    }

    // Loop over all words in the line and see what we can find
    nt = 0;
    for (j = i+1; j < L.row; j++) {
      if (strcmp(L.map[j], "=") == 0) {
        continue;
      }

      // This is relevant information
      if (strcmp(L.map[j], "Keq") == 0 && j < L.row-1) {
        if (order == 2) {
          restb->has_tK = 1;
          restb->targK = atof(L.map[j+1]);
          j++;
        }
        else if (order == 3) {
          resta->has_tK = 1;
          resta->targK = atof(L.map[j+1]);
          j++;
        }
      }
      else if (strcmp(L.map[j], "wtKeq") == 0 && j < L.row-1) {
        if (order == 2) {
          restb->has_wK = 1;
          restb->rstwK = atof(L.map[j+1]);
          j++;
        }
        else if (order == 3) {
          resta->has_wK = 1;
          resta->rstwK = atof(L.map[j+1]);
          j++;
        }
      }
      else if (strcmp(L.map[j], "Leq") == 0 && order == 2) {
        restb->has_tl0 = 1;
        restb->targl0 = atof(L.map[j+1]);
        j++;
      }
      else if (strcmp(L.map[j], "Theq") == 0 && order == 3) {
        resta->has_tTh0 = 1;
        resta->targTh0 = atof(L.map[j+1])*PI/180.0;
        j++;
      }
      else if (strcmp(L.map[j], "wtLeq") == 0 && order == 2) {
        restb->has_wl0 = 1;
        restb->rstwl0 = atof(L.map[j+1]);
        j++;
      }
      else if (strcmp(L.map[j], "wtTheq") == 0 && order == 3) {
        resta->has_wTh0 = 1;
        resta->rstwTh0 = atof(L.map[j+1]);
        j++;
      }
      else if ((strcmp(L.map[j], "period") == 0 ||
                strcmp(L.map[j], "per") == 0) && j < L.row-1 && order == 4) {
        resth->pn = atof(L.map[j+1]);
        j++;
      }
      else if ((strcmp(L.map[j], "weight") == 0 ||
                strcmp(L.map[j], "rwt") == 0) && j < L.row-1 && order == 4) {
        resth->rstw = atof(L.map[j+1]);
        j++;
      }
      else if ((strcmp(L.map[j], "target") == 0 ||
                strcmp(L.map[j], "trg") == 0) && j < L.row-1) {
        if (order == 4) {
          resth->target = atof(L.map[j+1]);
          j++;
        }
        else if (order == 6) {
          restr->target = atof(L.map[j+1]);
          j++;
        }
      }
      else if (WordIsNumber(L.map[j]) == 1 && order == 4) {
        resth->rstw = atof(L.map[j]);
      }
      else if (strcmp(L.map[j], "break") == 0 ||
               WordInGlossary(glossary, L.map[j]) == 1) {
        break;
      }
      else if (nt < order && nt < 4 && strlen(L.map[j]) <= 4) {
        if (nt == 0) {
          if (order == 2) strcpy(restb->atype, L.map[j]);
          else if (order == 3) strcpy(resta->atype, L.map[j]);
          else if (order == 4) strcpy(resth->atype, L.map[j]);
          else if (order == 6) strcpy(restr->amask.map[0], L.map[j]);
          nt++;
        }
        else if (nt == 1) {
          if (order == 2) strcpy(restb->btype, L.map[j]);
          else if (order == 3) strcpy(resta->btype, L.map[j]);
          else if (order == 4) strcpy(resth->btype, L.map[j]);
          else if (order == 6) strcpy(restr->bmask.map[0], L.map[j]);
          nt++;
        }
        else if (nt == 2) {
          if (order == 3) strcpy(resta->ctype, L.map[j]);
          else if (order == 4) strcpy(resth->ctype, L.map[j]);
          else if (order == 6) strcpy(restr->cmask.map[0], L.map[j]);
          nt++;
        }
        else if (nt == 3) {
          if (order == 4) {
            strcpy(resth->dtype, L.map[j]);
            nt++;
          }
          else if (order == 6) {
            strcpy(restr->dmask.map[0], L.map[j]);
            nt++;
          }
        }
      }
    }
    i = j;

    // Fill out atom type names (if this is not an NMR operation restraint)
    // and expand the array of restraints if necessary.
    if (order == 2) {
      TouchupAtomName(restb->atype);
      TouchupAtomName(restb->btype);
      mp->nuserrstB += 1;
      if (mp->nuserrstB == *maxsr) {
        *maxsr += 32;
        mp->userrstB = (bondrst*)realloc(mp->userrstB,
                                         (*maxsr)*sizeof(bondrst));
      }
    }
    else if (order == 3) {
      TouchupAtomName(resta->atype);
      TouchupAtomName(resta->btype);
      TouchupAtomName(resta->ctype);
      mp->nuserrstA += 1;
      if (mp->nuserrstA == *maxsr) {
        *maxsr += 32;
        mp->userrstA = (anglrst*)realloc(mp->userrstA,
                                         (*maxsr)*sizeof(anglrst));
      }
    }
    else if (order == 4) {
      TouchupAtomName(resth->atype);
      TouchupAtomName(resth->btype);
      TouchupAtomName(resth->ctype);
      TouchupAtomName(resth->dtype);
      mp->nuserrstH += 1;
      if (mp->nuserrstH == *maxsr) {
        *maxsr += 32;
        mp->userrstH = (torrst*)realloc(mp->userrstH, (*maxsr)*sizeof(torrst));
      }
    }
    else if (order == 6) {
      mp->nuserrstR += 1;
      if (mp->nuserrstR == *maxsr) {
        *maxsr += 32;
        mp->userrstR = (nmroper*)realloc(mp->userrstR,
                                         (*maxsr)*sizeof(nmroper));
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekGeomRest: seek a geometric constraint.  This applies to angles only:
//               some groups of angles need to add up to specific values (i.e.
//               360 degrees).  The format of a geometry restraint is 'angl XX
//               XX XX' repeated as many times as necessary, where 'XX' defines
//               one atom type, and also 'target <float>' for the target value
//               of the sum of all angles.
//
// Arguments:                                                           
//   L:        the list of words on this line of the input command file 
//   mp:       the main parameter set                                   
//   sname:    the string flag being sought                             
//   salias:   alias for the string flag being sought                   
//   maxsr:    the maximum number of geometry requests that can be held 
//-----------------------------------------------------------------------------
void SeekGeomRest(cmat L, prmset *mp, char* sname, char* salias, int *maxgeom)
{
  int i, j, k, nangl, natmtype;
  cmat gtype;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) != 0 && strcmp(L.map[i], salias) != 0) {
      continue;
    }

    // This is a geometry request.  Count the number of angles.
    gtype = CreateCmat(24, 5);
    nangl = 0;
    for (j = i+1; j < L.row-3; j++) {
      if (strcmp(L.map[j], "angl") == 0) {
        natmtype = 0;
        for (k = j+1; k < j+4; k++) {
          if (strlen(L.map[k]) <= 4) {
            natmtype++;
          }
        }
        if (natmtype == 3 && nangl < 8) {
          for (k = j+1; k < j+4; k++) {
            strcpy(gtype.map[3*nangl+k-j-1], L.map[k]);
            TouchupAtomName(gtype.map[3*nangl+k-j-1]);
          }
          nangl++;
        }
      }
    }
    if (nangl < 2) {
      continue;
    }

    // If we are still here, there are angles in the request.
    mp->anglsum[mp->nanglsum].atmtype = ReallocCmat(&gtype, 3*nangl, 5);
    mp->anglsum[mp->nanglsum].nvar = nangl;
    mp->anglsum[mp->nanglsum].fitcol = (int*)malloc(2*nangl*sizeof(int));
    mp->anglsum[mp->nanglsum].rstrow = (int*)malloc(nangl*sizeof(int));
    mp->anglsum[mp->nanglsum].basis = (double*)malloc(2*nangl*sizeof(double));
    SetIVec(mp->anglsum[mp->nanglsum].fitcol, 2*nangl, -1);
    mp->anglsum[mp->nanglsum].target = 2.0*PI;
    mp->anglsum[mp->nanglsum].srange = 5.0*PI/180.0;
    for (j = i+1; j < L.row-1; j++) {
      if (strcmp(L.map[j], "target") == 0 || strcmp(L.map[j], "trg") == 0) {
        mp->anglsum[mp->nanglsum].target = atof(L.map[j+1]) * PI/180.0;
      }
    }
    for (j = i+1; j < L.row-1; j++) {
      if (strcmp(L.map[j], "range") == 0 || strcmp(L.map[j], "lim") == 0) {
        mp->anglsum[mp->nanglsum].srange = atof(L.map[j+1]) * PI/180.0;
      }
    }
    mp->nanglsum += 1;
    if (mp->nanglsum >= *maxgeom) {
      *maxgeom += 32;
      mp->anglsum = (geomrst*)realloc(mp->anglsum, (*maxgeom)*sizeof(geomrst));
    }
  }
}

//-----------------------------------------------------------------------------
// SeekSpecReq: seek a spectral resampling request.  Adjustable dihedrals
//              within the fitting problem which matches the criteria given in
//              one or more of this array's spectrum requests will be flagged,
//              isolated in a sub-problem, and then held tightly to each of a
//              range of values enumerated on a one- or two-dimensional grid
//              as all other dihedrals flagged for spectrum optimization or    
//              support of it are also re-optimized given the restrained
//              variables.  The set of solutions to these sub-problems is the
//              spectrum of sub-optimal solutions which can then be evaluated
//              for promixity to the global minimum in the figure of merit
//              (energy RMSE from training data) and proximity to each other
//              in terms of the solution vectors.  One data set can thereby
//              yield many distinct solutions of comparable quality and
//              possibly different behavior.
//
// Arguments:                                                           
//   L:         the list of words on this line of the input command file 
//   mp:        the main parameter set                                   
//   sname:     the string flag being sought                             
//   salias:    alias for the string flag being sought                   
//   maxsr:     the maximum number of spectrum requests that can be held 
//   glossary:  a list of keywords to match against, to keep this open-ended
//              declaration from bleeding into other declarations within the
//              same namelist
//-----------------------------------------------------------------------------
void SeekSpecReq(cmat L, prmset *mp, char* sname, char* salias, int *maxsr,
                 cmat *glossary)
{
  int i, j, k, nt, slen;
  specreq *sr;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) != 0 && strcmp(L.map[i], salias) != 0) {
      continue;
    }

    // This is a spectrum request
    mp->spectrum = 1;
    sr = &mp->spectralH[mp->nspectralH];
    sr->level = 1;
    sr->order = 4;
    sr->minval = -2.0;
    sr->maxval = 2.0;
    sr->atmtypes = CreateCmat(4, 8);
    nt = 0;
    for (j = i+1; j < L.row; j++) {
      if (strcmp(L.map[j], "=") == 0) {
        continue;
      }

      // This is relevant information
      if (strcmp(L.map[j], "retain") == 0) {
        sr->level = 1;
      }
      else if (strcmp(L.map[j], "sample") == 0) {
        sr->level = 2;
      }
      else if (strcmp(L.map[j], "order") == 0 && j < L.row-1) {
        sr->order = atoi(L.map[j+1]);
        j++;
      }
      else if (strcmp(L.map[j], "min") == 0 && j < L.row-1) {
        sr->minval = atof(L.map[j+1]);
        j++;
      }
      else if (strcmp(L.map[j], "max") == 0 && j < L.row-1) {
        sr->maxval = atof(L.map[j+1]);
        j++;
      }
      else if (strcmp(L.map[j], "spc") == 0 && j < L.row-1) {
        sr->gspc = atof(L.map[j+1]);
        j++;
      }
      else if (strlen(L.map[j]) <= 4 && nt < 4) {
        strcpy(sr->atmtypes.map[nt], L.map[j]);
        nt++;
      }
      else if (strcmp(L.map[j], "break") == 0 ||
               WordInGlossary(glossary, L.map[j]) == 1) {
        break;
      }
    }
    i = j;

    // Compact the atom type matrix and properly name its types
    if (nt < 4) {
      sr->atmtypes = ReallocCmat(&sr->atmtypes, nt, 8);
    }
    for (j = 0; j < nt; j++) {
      TouchupAtomName(sr->atmtypes.map[j]);
    }

    // Expand the array of spectrum requests if necessary
    mp->nspectralH += 1;
    if (mp->nspectralH == *maxsr) {
      j = *maxsr + 32;
      mp->spectralH = (specreq*)realloc(mp->spectralH, j*sizeof(specreq));
    }
  }
}

//-----------------------------------------------------------------------------
// SeekStringPlusVal: seek a string plus a real number variable, then increment
//                    a counter.                              
//
// Arguments:                                                           
//   L:        the list of words on this line of the input command file 
//   val1:     character string storing the first input value           
//   val2:     double-precision pointer storing the second input value  
//   sname:    the string flag being sought                             
//   salias:   alias for the string flag being sought                   
//   counter:  counter variable                                         
//-----------------------------------------------------------------------------
void SeekStringPlusVal(cmat L, char* val1, double *val2, char* sname,
                       char* salias, int *counter)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+2 < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(val1, L.map[i+1]);
        *val2 = atof(L.map[i+2]);
        *counter += 1;
      }
      else if (i+3 < L.row) {
        strcpy(val1, L.map[i+2]);
        *val2 = atof(L.map[i+3]);
        *counter += 1;
      }
      else {
        printf("SeekStringPlusVal >> Error.  No value specified for "
               "identifier %s/%s.\n", sname, salias);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekRecord: seek a string input, but increment a counter each time an input
//             match can be found.  Increment a counter and store the input in
//             a character matrix according to the counter.
//
// Arguments:                                                           
//   L:        the list of words on this line of the input command file 
//   C:        character matrix storing the input values (records)      
//   sname:    the string flag being sought                             
//   salias:   alias for the string flag being sought                   
//   counter:  counter variable                                         
//-----------------------------------------------------------------------------
void SeekRecord(cmat L, cmat *C, char* sname, char* salias, int *counter)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(C->map[*counter], L.map[i+1]);
        *counter += 1;
      }
      else if (i+2 < L.row) {
        strcpy(C->map[*counter], L.map[i+2]);
        *counter += 1;
      }
      else {
        printf("SeekRecord >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
      }
    }
  }

  // Extend matrix if needed
  if (*counter == C->row) {
    *C = ReallocCmat(C, C->row+32, C->col);
  }
}

//-----------------------------------------------------------------------------
// ReadNumericalSuffix: read the numerical suffix to a control variable.
//                                                                      
// Arguments:                                                           
//   term:     the control variable                                     
//   nstart:   the anticipated start of the numerical part of term      
//   notanum:  flag to indicate that the numerical part contained non-integer
//             characters                                   
//-----------------------------------------------------------------------------
static int ReadNumericalSuffix(char* term, int nstart, int *notanum)
{
  int i, vpos, lflag;

  // Determine any numerical suffix to the flag
  lflag = strlen(term);
  *notanum = 0;
  if (lflag == nstart) {
    *notanum = 1;
    return -1;
  }
  for (i = nstart; i < lflag; i++) {
    if (term[i] < '0' || term[i] > '9') {
      *notanum = 1;
      return -1;
    }
  }
  vpos = atoi(&term[nstart])-1;

  return vpos;
}

//-----------------------------------------------------------------------------
// SeekNString: find a tag for a particular string and then check it for
//              numerical extensions denoting that the string denotes an input
//              for an alternate molecular system.                
//                                                                      
// Arguments:                                                           
//   Similar to SeekString above, but val is now a cmat struct with more than
//   one row.  The integer array fspec denotes whether any rows of val had
//   already been specified before this function was called, and therefore
//   should not be reassigned.                                
//-----------------------------------------------------------------------------
cmat SeekNString(cmat L, cmat *val, int* fspec, char* sname, char* salias)
{
  int i, nstart, notanumber, vpos, sfound, afound;

  const int lsname = strlen(sname);
  const int lsalias = strlen(salias);
  for (i = 0; i < L.row; i++) {
    sfound = 0;
    afound = 0;
    if (strncmp(L.map[i], sname, lsname) == 0) {
      sfound = 1;
      nstart = lsname;
    }
    else if (strncmp(L.map[i], salias, lsalias) == 0) {
      afound = 1;
      nstart = lsalias;
    }
    if (sfound == 1 || afound == 1) {
      vpos = ReadNumericalSuffix(L.map[i], nstart, &notanumber);
      if (notanumber == 1) {
        continue;
      }
      if (vpos < 0 || vpos > MAXSYS) {
        printf("SeekNString >> Error.  Invalid numerical value %d specified "
               "for identifier\nSeekNString >> %s/%s.\n", vpos+1, sname,
               salias);
        exit(1);
      }
      if (vpos >= val->row) {
        *val = ReallocCmat(val, vpos+1, MAXNAME);
      }
      else if (fspec[vpos] == 1) {
        continue;
      }

      // This string value may be set
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(val->map[vpos], L.map[i+1]);
      }
      else if (i+2 < L.row) {
        strcpy(val->map[vpos], L.map[i+2]);
      }
      else {
        printf("SeekNString >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
        exit(1);
      }
    }
  }

  return *val;
}

//-----------------------------------------------------------------------------
// SeekFixQ:  finds a tag for a charge fix restraint                    
//                                                                      
// Arguments:                                                           
//   L:       the list of words on this line of the input command file  
//   val1:    first destination of ambmask string                       
//   val2:    solvated charge target                                    
//   val3:    second destination of ambmask string                      
//   val4:    vacuum charge target, defaults to solvated charge target  
//   sname:   the string flag being sought                              
//   salias:  alias for the string flag being sought                    
//   counter: the counter variable (gets incremented)                   
//-----------------------------------------------------------------------------
void SeekFixQ(cmat L, char* val1, double *val2, char* val3, double *val4,
              char* sname, char* salias, int *counter)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+2 < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(val1, L.map[i+1]);
        strcpy(val3, L.map[i+1]);
        *val2 = atof(L.map[i+2]);
        if (i+3 < L.row) {
          *val4 = atof(L.map[i+3]);
        }
        else {
          *val4 = atof(L.map[i+2]);
        }
        *counter += 1;
      }
      else if (i+3 < L.row) {
        strcpy(val1, L.map[i+2]);
        strcpy(val3, L.map[i+2]);
        *val2 = atof(L.map[i+3]);
        if (i+4< L.row) {
          *val4 = atof(L.map[i+4]);
        }
        else {
          *val4 = atof(L.map[i+3]);
        }
        *counter += 1;
      }
      else {
        printf("SeekFixQ >> Error.  Invalid inputs for identifier "
               "%s/%s.]n", sname, salias);
        printf("SeekFixQ >> This identifier requires an amber mask "
               "a solution-phase atomic charge to be fixed for the mask, "
               "and optionally a vacuum-phase atomic charge to be fixed for "
               "the mask.\n");
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekAmberLibrary: find input to write an Amber-format library based on an
//                   existing template file.
//
// Arguments:
//   L:        the list of words on this line of the input command file
//-----------------------------------------------------------------------------
void SeekAmberLibrary(cmat L, fset *myfit, cmat *glossary)
{
  int i, j, hasout, hastmp, hasstl;
  char lbuff[MAXNAME];
  
  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], "amblib") == 0 ||
        strcmp(L.map[i], "AmberLibrary") == 0) {

      // Increment the number of libraries to process
      if (myfit->namblib >= myfit->maxamblib) {
        myfit->maxamblib += 8;
        if (myfit->namblib == 0) {
          myfit->qlibout = CreateCmat(myfit->maxamblib, MAXNAME);
          myfit->qlibtmp = CreateCmat(myfit->maxamblib, MAXNAME);
          myfit->qlibstyle = CreateImat(1, myfit->maxamblib);
        }
        else {
          myfit->qlibout = ReallocCmat(&myfit->qlibout, myfit->maxamblib,
                                       MAXNAME);
          myfit->qlibtmp = ReallocCmat(&myfit->qlibtmp, myfit->maxamblib,
                                       MAXNAME);
          myfit->qlibstyle = ReallocImat(&myfit->qlibstyle, 1,
                                         myfit->maxamblib);
        }
      }
      hasout = 0;
      hastmp = 0;
      hasstl = 0;
      for (j = i+1; j < L.row - 1; j++) {
        if (WordInGlossary(glossary, L.map[j]) == 1) {
          break;
        }
        else if (strcmp(L.map[j], "out") == 0) {
          strcpy(myfit->qlibout.map[myfit->namblib], L.map[j+1]);
          hasout = 1;
        }
        else if (strcmp(L.map[j], "template") == 0) {
          strcpy(myfit->qlibtmp.map[myfit->namblib], L.map[j+1]);
          hastmp = 1;
        }
        else if (strcmp(L.map[j], "style") == 0) {
          strcpy(lbuff, L.map[j+1]);
          StringToLower(lbuff);
          if (strcmp(lbuff, "resp") == 0) {
            myfit->qlibstyle.map[0][myfit->namblib] = 0;
            hasstl = 1;
          }
          else if (strcmp(lbuff, "vacuum") == 0) {
            myfit->qlibstyle.map[0][myfit->namblib] = 1;
            hasstl = 1;
          }
          else if (strcmp(lbuff, "condensed") == 0) {
            myfit->qlibstyle.map[0][myfit->namblib] = 2;
            hasstl = 1;
          }
          else if (strcmp(lbuff, "hyperpol") == 0) {
            myfit->qlibstyle.map[0][myfit->namblib] = 3;
            hasstl = 1;
          }
        }
      }
      if (hasout == 1 && hastmp == 1 && hasstl == 1) {
        myfit->namblib += 1;
      }
      else {
        if (hasout == 0) {
          printf("SeekAmberLibrary >> Error.  No ouput amber library was "
                 "specified.\n");
        }
        if (hastmp == 0) {
          printf("SeekAmberLibrary >> Error.  No template amber library was "
                 "specified.\n");
        }
        if (hasstl == 0) {
          printf("SeekAmberLibrary >> Error.  No partial charge style was "
                 "specified.\n");
        }
        exit(1);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekManipulator: find input commands for an operation to perform on a
//                  molecule that will alter its configuration / conformation.
//                  This function can get called with many different aliases,
//                  but unlike other functions like SeekInt this one detects
//                  the input control parameter name or alias and does
//                  different things based on what those terms are. 
//
// Arguments:
//   L:        the list of words on this line of the input command file
//   cfsinp:   the configuration generator input control data
//   maxops:   the maximum number of operations that cfsinp can hold (the data
//             structure itself stores the exact number, but the maximum number
//             is probably a local variable in the function calling this one so
//             that it doesn't have to be stored forever).
//   sname:    the string flag being sought
//   salias:   alias for the string flag being sought
//   glossary: a glossary of all keywords in the namelist where Manipulators
//             may be found (&configs).  This is needed for open-ended commands
//             like this one so that the contents of other commands don't get
//             tossed in by mistake.
//-----------------------------------------------------------------------------
void SeekManipulator(cmat L, configs *cfsinp, int *maxops, char* sname,
                     char* salias, cmat *glossary)
{
  int i, j, k, minfound, maxfound, absrange, scatter, strategy;
  int natom, FtopGiven, quadwinGiven, labelGiven;
  double Ftop;
  cfigop *myop;

  // Figure out what flavor of manipulation this is
  absrange = 0;
  scatter = 0;
  if (strcmp(salias, "random") == 0) {
    absrange = 1;
    scatter = 1;
    strategy = 0;
  }
  if (strcmp(salias, "uniform") == 0) {
    absrange = 1;
    scatter = 0;
    strategy = 1;
  }
  if (strcmp(salias, "rpert") == 0) {
    absrange = 0;
    scatter = 1;
    strategy = 2;
  }
  if (strcmp(salias, "gpert") == 0) {
    absrange = 0;
    scatter = 0;
    strategy = 3;
  }
  if (strcmp(salias, "set") == 0) {
    absrange = 1;
    scatter = 0;
    strategy = 4;
  }

  // Start recording new operations
  myop = &cfsinp->ops[cfsinp->nops];
  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {

      // We've found what appears to be an operation.  Start reading the
      // input.  The number of atoms involved will be determined by the
      // first recognizable keyword (up until then, all words will be
      // assumed to be atom names).
      natom = 0;
      minfound = 0;
      maxfound = 0;
      myop->Krst = 64.0;
      myop->distance = 0.0;
      myop->atommasks = CreateCmat(4, 64);
      myop->fbhw = 0.0;
      myop->quadwin = 0.5;

      // Keep track of what the user has said--this will be important when
      // interpreting what the user wants the restraint to look like.
      FtopGiven = 0;
      quadwinGiven = 0;
      labelGiven = 0;

      // Look at what the user has told mdgx
      for (j = i+1; j < L.row; j++) {

        // Detect various range specifications or a stiffness constant
        if (j < L.row-3 &&
            strcmp(L.map[j], "{") == 0 && strcmp(L.map[j+3], "}") == 0 &&
            WordIsNumber(L.map[j+1]) == 1 && WordIsNumber(L.map[j+2])) {
          myop->minval = atof(L.map[j+1]);
          myop->maxval = atof(L.map[j+2]);
          minfound = 1;
          maxfound = 1;
          j += 3;
        }
        else if (j < L.row-8 && (strcmp(L.map[j], "box") == 0 ||
                                 strcmp(L.map[j], "line") == 0) &&
                 strcmp(L.map[j+1], "{") == 0 &&
                 strcmp(L.map[j+8], "}") == 0 &&
                 WordIsNumber(L.map[j+2]) == 1 && WordIsNumber(L.map[j+3]) &&
                 WordIsNumber(L.map[j+4]) == 1 && WordIsNumber(L.map[j+5]) &&
                 WordIsNumber(L.map[j+6]) == 1 && WordIsNumber(L.map[j+7])) {
          if (strcmp(L.map[j], "box") == 0) {
            myop->pegtype = 0;
          }
          else {
            myop->pegtype = 1;
          }
          myop->minvx = atof(L.map[j+2]);
          myop->minvy = atof(L.map[j+3]);
          myop->minvz = atof(L.map[j+4]);
          myop->maxvx = atof(L.map[j+5]);
          myop->maxvy = atof(L.map[j+6]);
          myop->maxvz = atof(L.map[j+7]);
          minfound = 1;
          maxfound = 1;
          j += 8;
        }
        else if (j < L.row-5 && strcmp(L.map[j], "peg") == 0 &&
                 strcmp(L.map[j+1], "{") == 0 &&
                 strcmp(L.map[j+5], "}") == 0 &&
                 WordIsNumber(L.map[j+2]) == 1 && WordIsNumber(L.map[j+3]) &&
                 WordIsNumber(L.map[j+4]) == 1) {
          myop->minvx = atof(L.map[j+2]);
          myop->minvy = atof(L.map[j+3]);
          myop->minvz = atof(L.map[j+4]);
          myop->maxvx = myop->minvx;
          myop->maxvy = myop->minvy;
          myop->maxvz = myop->minvz;
          minfound = 1;
          maxfound = 1;
          j += 5;
        }
        else if (j < L.row-1 && WordIsNumber(L.map[j+1]) == 1 &&
                 (strcmp(L.map[j], "min") == 0 ||
                  strcmp(L.map[j], "max") == 0 ||
                  strcmp(L.map[j], "value") == 0 ||
                  strcmp(L.map[j], "Krst") == 0 ||
                  strcmp(L.map[j], "distance") == 0 ||
                  strcmp(L.map[j], "fbhw") == 0 ||
                  strcmp(L.map[j], "quadratic") == 0 ||
                  strcmp(L.map[j], "Ftop") == 0)) {
          if (strcmp(L.map[j], "min") == 0 && strategy < 4) {
            myop->minval = atof(L.map[j+1]);
            minfound = 1;
          }
          if (strcmp(L.map[j], "max") == 0 && strategy < 4) {
            myop->maxval = atof(L.map[j+1]);
            maxfound = 1;
          }
          if (strcmp(L.map[j], "value") == 0 && strategy == 4) {
            myop->minval = atof(L.map[j+1]);
            myop->maxval = myop->minval; 
            minfound = 1;
            maxfound = 1;
          }
          if (strcmp(L.map[j], "Krst") == 0) {
            myop->Krst = atof(L.map[j+1]);
          }
          if (strcmp(L.map[j], "distance") == 0) {
            myop->distance = atof(L.map[j+1]);
          }
          if (strcmp(L.map[j], "fbhw") == 0) {
            myop->fbhw = atof(L.map[j+1]);
          }
          if (strcmp(L.map[j], "Ftop") == 0) {
            Ftop = atof(L.map[j+1]);
            FtopGiven = 1;
          }
          if (strcmp(L.map[j], "quadratic") == 0) {
            myop->quadwin = atof(L.map[j+1]);
            quadwinGiven = 1;
          }
          j++;
        }
        else if (j < L.row-1 && strcmp(L.map[j], "label") == 0) {
          if (strlen(L.map[j+1]) >= 32) {
            L.map[j+1][31] = '\0';
          }
          strcpy(myop->label, L.map[j+1]);
          labelGiven = 1;
          j++;
        }

        // Detect a break
        else if (strcmp(L.map[j], "break") == 0 ||
                 WordInGlossary(glossary, L.map[j]) == 1) {
          break;
        }

        // If we're still here this must be an ambmask string for an atom
        else {

          // Check on the number of atoms
          if (natom == 4) {
            printf("SeekManipulator >> Error.  Manipulator command contains "
                   "too many atom names:\n");
            for (k = 0; k < L.row; k++) {
              printf("%s ", L.map[k]);
            }
            printf("\n");
            exit(1);
          }
          strcpy(myop->atommasks.map[natom], L.map[j]);
          natom++;
        }
      }

      // Add to the list of operations if this looks to be in order, or
      // clean up the operation in memory if we are not making it official.
      if (natom >= 1 && minfound == 1 && maxfound == 1 &&
          (natom == 1 || absrange == 0 || myop->minval <= myop->maxval)) {
        myop->order = natom;
        if (myop->order >= 3) {
          myop->minval *= PI/180.0;
          myop->maxval *= PI/180.0;
          myop->fbhw *= PI/180.0;
        }
        myop->absrange = absrange;
        myop->scatter = scatter;
        if (quadwinGiven == 1 && natom >= 3) {
          myop->quadwin *= PI/180.0;
        }
        if (quadwinGiven == 0 && FtopGiven == 1) {
          myop->quadwin = fabs(Ftop)/fabs(2.0*myop->Krst);
        }
        if (labelGiven == 0) {
          sprintf(myop->label, "%d\n", cfsinp->nops+1);
        }
        cfsinp->nops += 1;
        if (cfsinp->nops >= *maxops) {
          *maxops += 32;
          cfsinp->ops = (cfigop*)realloc(cfsinp->ops,
                                         (*maxops)*sizeof(cfigop));
        }
        myop = &cfsinp->ops[cfsinp->nops];
      }
      else {
        DestroyCmat(&myop->atommasks);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekPeptideSystem: find a directive to launch a peptide system in implicit
//                    solvent.
//
// Arguments:
//
//-----------------------------------------------------------------------------
void SeekPeptideSystem(cmat L, pepcon *ppctrl, int *maxsys, char* sname,
		       char* salias, cmat *glossary)
{
  int i, j, nsys;

  nsys = ppctrl->nsys;
  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {

      // This is a peptide system.  Initialize and then log its contents.
      ppctrl->Tranges[nsys].x = -298.0;
      ppctrl->Tranges[nsys].y = -298.0;
      ppctrl->Treplicas[nsys] = 1;
      ppctrl->Preplicas[nsys] = 1;
      ppctrl->Sreplicas[nsys] = 1;
      ppctrl->tpinames.map[nsys][0] = '\0';
      ppctrl->tpfnames.map[nsys][0] = '\0';
      for (j = i+1; j < L.row; j++) {

	// Detect various range specifications or a stiffness constant
	if (j < L.row-1 && strcmp(L.map[j], "Tmin") == 0 &&
	    WordIsNumber(L.map[j+1])) {
	  ppctrl->Tranges[nsys].x = atof(L.map[j+1]);
	  j++;
	}
	else if (j < L.row-1 && strcmp(L.map[j], "Tmax") == 0 &&
		 WordIsNumber(L.map[j+1])) {
	  ppctrl->Tranges[nsys].y = atof(L.map[j+1]);
	  j++;
	}
	else if (j < L.row-1 && strcmp(L.map[j], "temp0") == 0 &&
		 WordIsNumber(L.map[j+1])) {
	  ppctrl->Tranges[nsys].x = atof(L.map[j+1]);
	  ppctrl->Tranges[nsys].y = atof(L.map[j+1]);
	  j++;
	}
	else if (j < L.row-1 && strcmp(L.map[j], "T-rep") == 0 &&
		 WordIsInteger(L.map[j+1])) {
	  ppctrl->Treplicas[nsys] = atoi(L.map[j+1]);
	  j++;
	}
        else if (j < L.row-1 && strcmp(L.map[j], "P-rep") == 0 &&
                 WordIsInteger(L.map[j+1])) {
	  ppctrl->Preplicas[nsys] = atoi(L.map[j+1]);
	  j++;
	}
	else if (j < L.row-1 && strcmp(L.map[j], "N-rep") == 0 &&
		 WordIsInteger(L.map[j+1])) {
	  ppctrl->Sreplicas[nsys] = atoi(L.map[j+1]);
	  j++;
	}
	else if (j < L.row-1 && strcmp(L.map[j], "-p") == 0) {
	  strcpy(ppctrl->tpinames.map[nsys], L.map[j+1]);
	  j++;
	}
	else if (j < L.row-1 && strcmp(L.map[j], "-pi") == 0) {
	  strcpy(ppctrl->tpinames.map[nsys], L.map[j+1]);
	  j++;
	}
	else if (j < L.row-1 && strcmp(L.map[j], "-pf") == 0) {
	  strcpy(ppctrl->tpfnames.map[nsys], L.map[j+1]);
	  j++;
	}
	else if (j < L.row-1 && strcmp(L.map[j], "-c") == 0) {
	  strcpy(ppctrl->crdnames.map[nsys], L.map[j+1]);
	  j++;
	}
	else if (j < L.row-1 && strcmp(L.map[j], "-x") == 0) {
	  strcpy(ppctrl->trjbases.map[nsys], L.map[j+1]);
	  j++;
	}
	else if (j < L.row-1 && strcmp(L.map[j], "-o") == 0) {
	  strcpy(ppctrl->outbases.map[nsys], L.map[j+1]);
	  j++;
	}
	else if (j < L.row-1 && strcmp(L.map[j], "-r") == 0) {
	  strcpy(ppctrl->rstrtbases.map[nsys], L.map[j+1]);
	  j++;
	}

	// Detect a break
	else if (strcmp(L.map[j], "break") == 0 ||
		 WordInGlossary(glossary, L.map[j]) == 1) {
	  break;
	}
      }

      // Check input
      if (ppctrl->tpinames.map[nsys][0] == '\0') {
	if (ppctrl->tpfnames.map[nsys][0] != '\0') {
          printf("SeekPeptideSystem >> Error.  A base topology must be "
		 "specified for system\nSeekPeptideSystem >> %d (perturbed "
		 "topology %s)\n", nsys, ppctrl->tpfnames.map[nsys]);
	}
	else {
          printf("SeekPeptideSystem >> Error.  A base topology must be "
		 "specified for system\nSeekPeptideSystem >> %d.\n", nsys);
	}
      }
      
      // Increment the number of peptide systems,
      // and allocate new memory if necessary.
      ppctrl->nsys += 1;
      if (ppctrl->nsys >= *maxsys) {
	*maxsys *= 2;
	ppctrl->tpinames = ReallocCmat(&ppctrl->tpinames, *maxsys, MAXNAME);
	ppctrl->tpfnames = ReallocCmat(&ppctrl->tpfnames, *maxsys, MAXNAME);
	ppctrl->outbases = ReallocCmat(&ppctrl->outbases, *maxsys, MAXNAME);
	ppctrl->crdnames = ReallocCmat(&ppctrl->crdnames, *maxsys, MAXNAME);
	ppctrl->trjbases = ReallocCmat(&ppctrl->trjbases, *maxsys, MAXNAME);
	ppctrl->rstrtbases = ReallocCmat(&ppctrl->rstrtbases, *maxsys,
					 MAXNAME);
	ppctrl->Tranges = (double2*)realloc(ppctrl->Tranges,
					    (*maxsys) * sizeof(double2));
	ppctrl->Treplicas = (int*)realloc(ppctrl->Treplicas,
					  (*maxsys) * sizeof(int));
	ppctrl->Preplicas = (int*)realloc(ppctrl->Preplicas,
					  (*maxsys) * sizeof(int));
	ppctrl->Sreplicas = (int*)realloc(ppctrl->Sreplicas,
					  (*maxsys) * sizeof(int));
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekOpsCombo: find a directive to combine two or three configuration
//               sampling operations.  This will instate two or three
//               dimensional grid sampling.
//
// Arguments:
//   Same as for SeekManipulator above, except that maxops is now maxcombos,
//   still storing the maximum number of these directives that can currently
//   be stored.
//-----------------------------------------------------------------------------
void SeekOpsCombo(cmat L, configs *cfsinp, int *maxcombos, char* sname,
                  char* salias, cmat *glossary)
{
  int i, j;
  coupler *mycombo;

  mycombo = &cfsinp->combo[cfsinp->ncombo];
  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {

      // This appears to be a combination.  Start reading.
      j = i+1;
      mycombo->npcs = 0;
      while (j < L.row && strcmp(L.map[j], "break") != 0 &&
             WordInGlossary(glossary, L.map[j]) == 0 && mycombo->npcs < 3) {
        strcpy(mycombo->labels[mycombo->npcs], L.map[j]);
        mycombo->npcs += 1;
        j++;
      }
      cfsinp->ncombo += 1;
      if (cfsinp->ncombo >= *maxcombos) {
        *maxcombos += 32;
        cfsinp->combo = (coupler*)realloc(cfsinp->combo,
                                          (*maxcombos)*sizeof(coupler));
      }
      mycombo = &cfsinp->combo[cfsinp->ncombo];
    }
  }
}

//-----------------------------------------------------------------------------
// SeekSPItem: find a piece of single point data.
//
// Arguments:
//   L:        the list of words taken from the line
//   spelist:  the data structure storing the single point energy items
//   maxitems: the maximum number of operations that spelist is ready to hold
//             (the data structure itself stores the exact number, but the
//             maximum number is a local variable in the function calling this
//             one so that it doesn't have to be stored forever).
//   sname:    the string flag being sought
//   salias:   alias for the string flag being sought
//   glossary: a glossary of all keywords in the namelist where Manipulators
//             may be found (&configs).  This is needed for open-ended commands
//             like this one so that the contents of other commands don't get
//             tossed in by mistake.
//-----------------------------------------------------------------------------
void SeekSPItem(cmat L, spdata *spelist, int *maxitems, char* sname,
                char* salias, cmat *glossary)
{
  int i, j, hascrd, hase;
  sptask *myitem;

  for (i = 0; i < L.row-1; i++) {
    if (strcmp(L.map[i], sname) != 0 && strcmp(L.map[i], salias) != 0) {
      continue;
    }

    // This hit the keyword for a single point data item.
    // Now find out what's here.
    j = i+1;
    myitem = &spelist->items[spelist->nitems];
    myitem->takeqmcrd = 1;
    myitem->crdsrc = (char*)malloc(MAXNAME*sizeof(char));
    myitem->esrc = (char*)malloc(MAXNAME*sizeof(char));
    hase = 0;
    hascrd = 0;
    while (j < L.row && WordInGlossary(glossary, L.map[j]) != 1) {
      if (strcmp(L.map[j], "-nrg") == 0 && j < L.row-1) {
        strcpy(myitem->esrc, L.map[j+1]);
        hase = 1;
      }
      else if (strcmp(L.map[j], "-coord") == 0 && j < L.row-1) {
        strcpy(myitem->crdsrc, L.map[j+1]);
        myitem->takeqmcrd = 0;
        hascrd = 1;
      }
      else if (hase == 0 && j < L.row) {
        strcpy(myitem->esrc, L.map[j]);
        hase = 1;
      }
      else if (hascrd == 0 && j < L.row) {
        strcpy(myitem->crdsrc, L.map[j]);
        hascrd = 1;
      }
      j++;
    }
    if (hase == 1) {
      spelist->nitems += 1;
      if (spelist->nitems >= *maxitems) {
        *maxitems *= 2;
        spelist->items = (sptask*)realloc(spelist->items,
                                          (*maxitems)*sizeof(sptask));
      }
    }
    else {

      // Free allocated memory if this failed to find a complete item
      free(myitem->crdsrc);
      free(myitem->esrc);
    }
  }
}

//-----------------------------------------------------------------------------
// SeekNumberedName: find a name string and then a tag within it which delimits
//                   a prefix, place to put a number for multiple different 
//                   instances of the information (i.e. a molecular
//                   configuration), and a suffix.
//
// Arguments:
//   L:       the list of words on this line of the input command file
//   base:    the base name (this pre-allocated character array will be filled
//            as a means of returning a value)
//   suffix:  similar to base but for the file suffix
//   sname:   the string flag being sought
//   salias:  alias for the string flag being sought
//-----------------------------------------------------------------------------
void SeekNumberedName(cmat L, char* base, char* suffix, char* sname,
                      char* salias)
{
  int i, found, slen, pivot;

  // Strategy: store the value found entirely in base, then parse from there.
  found = 0;
  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        strcpy(base, L.map[i+1]);
        found = 1;
      }
      else if (i+2 < L.row) {
        strcpy(base, L.map[i+2]);
        found = 1;
      }
      else {
        printf("SeekNumberedName >> Error.  No value specified for identifier "
               "%s/%s.\n", sname, salias);
        exit(1);
      }
    }
  }

  // Now parse the name and find the "[N]" delimiter that indicates
  // 'number goes here.'
  if (found == 1) {

    // Find the location of the delimiter
    slen = strlen(base);
    pivot = -1;
    for (i = 0; i <= slen-3; i++) {
      if (strncmp(&base[i], "[N]", 3) == 0) {
        pivot = i;
      }
    }

    // If there was no delimiter, the file name is all base.
    if (pivot == -1 || pivot == slen-3) {
      suffix[0] = '\0';
    }
    else {
      for (i = pivot+3; i < slen; i++) {
        suffix[i-pivot-3] = base[i];
        base[i] = '\0';
      }
      suffix[slen-pivot-3] = '\0';
      base[pivot] = '\0';
      base[pivot+1] = '\0';
      base[pivot+2] = '\0';
    }
  }
}

//-----------------------------------------------------------------------------
// SeekReal: find a tag for a particular real value.                    
//
// Arguments:
//   L:       the list of words on this line of the input command file
//   val:     string containing the possible real number
//   sname:   the string flag being sought
//   salias:  alias for the string flag being sought
//-----------------------------------------------------------------------------
void SeekReal(cmat L, double *val, char* sname, char* salias)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        *val = atof(L.map[i+1]);
      }
      else if (i+2 < L.row) {
        *val = atof(L.map[i+2]);
      }
      else {
        printf("SeekReal >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
        exit(1);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekNReal: find a tag for a particular real value, one of a series.  Unlike
//            SeekNString, this function does not check for the
//            pre-specification of values and does not limit the series number
//            to within [ 0, MAXSYS ].                    
//                                                                      
// Arguments:                                                           
//   Same as SeekReal above, except in this case val is an array not a pointer.
//   The val array must be pre-allocated.  maxidx indicates the maximum
//   permissible value ID.                                  
//-----------------------------------------------------------------------------
void SeekNReal(cmat L, double* val, char* sname, char* salias, int maxidx)
{
  int i, lsname, lsalias, vidx, match, nstart, notanumber;

  lsname = strlen(sname);
  lsalias = strlen(salias);
  for (i = 0; i < L.row; i++) {
    match = 0;
    if (strncmp(L.map[i], sname, lsname) == 0) {
      match = 1;
      nstart = lsname;
    }
    else if (strncmp(L.map[i], salias, lsalias) == 0) {
      match = 1;
      nstart = lsalias;
    }
    if (match == 1) {
      vidx = ReadNumericalSuffix(L.map[i], nstart, &notanumber);
      if (notanumber == 1) {
        continue;
      }
      if (vidx < 0 || vidx >= maxidx) {
        printf("SeekNReal >> Error.  Invalid numerical value %d specified "
               "for identifier\nSeekNReal >> %s/%s.\n", vidx+1, sname,
               salias);
        exit(1);
      }
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        val[vidx] = atof(L.map[i+1]);
      }
      else if (i+2 < L.row) {
        val[vidx] = atof(L.map[i+2]);
      }
      else {
        printf("SeekNReal >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
        exit(1);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekInt: find a tag for a particular integer value.                  
//
// Arguments:
//   L:       the list of words on this line of the input command file
//   val:     string containing the possible integer
//   sname:   the string flag being sought
//   salias:  alias for the string flag being sought
//-----------------------------------------------------------------------------
void SeekInt(cmat L, int *val, char* sname, char* salias)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        *val = strtol(L.map[i+1], NULL, 10);
      }
      else if (i+2 < L.row) {
        *val = strtol(L.map[i+2], NULL, 10);
      }
      else {
        printf("SeekInt >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
        exit(1);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// SeekLLInt: find a tag for a particular long long integer value.      
//
// Arguments:
//   L:       the list of words on this line of the input command file
//   val:     string containing the possible long long int
//   sname:   the string flag being sought
//   salias:  alias for the string flag being sought
//-----------------------------------------------------------------------------
void SeekLLInt(cmat L, long long int *val, char* sname, char* salias)
{
  int i;

  for (i = 0; i < L.row; i++) {
    if (strcmp(L.map[i], sname) == 0 || strcmp(L.map[i], salias) == 0) {
      if (i+1 < L.row && strcmp(L.map[i+1], "=") != 0) {
        sscanf(L.map[i+1], "%lld", val);
      }
      else if (i+2 < L.row) {
        sscanf(L.map[i+2], "%lld", val);
      }
      else {
        printf("SeekLLInt >> Error.  No value specified for identifier %s/%s"
               ".\n", sname, salias);
        exit(1);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// ParseAmbmask: parse an ambmask string.  Uses ptrajmask.c             
//                                                                      
// Arguments:                                                           
//   maskstr:   the mask string                                         
//   tp:        the topology                                            
//   crd:       the coordinates                                         
//-----------------------------------------------------------------------------
int* ParseAmbMask(char* maskstr, prmtop *tp, coord *crd)
{
  int i, j;
  int* imask;
  char* cmask;
  NAME* atomNames;
  NAME* residueNames;
  NAME* atomTypes;

  // Reshape the atom and residue name arrays
  atomNames = (NAME*)malloc(tp->natom*sizeof(NAME));
  atomTypes = (NAME*)malloc(tp->natom*sizeof(NAME));
  residueNames = (NAME*)malloc(tp->nres*sizeof(NAME));
  for (i = 0; i < tp->natom; i++) {
    for (j = 0; j < 4; j++) {
      atomNames[i][j] = tp->AtomNames[4*i+j];
      atomTypes[i][j] = tp->AtomTypes[4*i+j];
    }
    atomNames[i][4] = '\0';
    atomTypes[i][4] = '\0';
  }
  for (i = 0; i < tp->nres; i++) {
    for (j = 0; j < 4; j++) {
      residueNames[i][j] = tp->ResNames[4*i+j];
    }
    residueNames[i][4] = '\0';
  }

  // Call to Dan Roe's C-implementation parser
  cmask = parseMaskString(maskstr, tp->natom, tp->nres, atomNames,
                          residueNames, tp->ResLims, crd->loc, atomTypes, 0);

  // Conversion to mdgx mask format
  imask = (int*)malloc(tp->natom*sizeof(int));
  if (cmask == NULL) {
    SetIVec(imask, tp->natom, 0);
  }
  else {
    for (i = 0; i < tp->natom; i++) {
      imask[i] = (cmask[i] == 'T') ? 1 : 0;
    }
  }

  // Free allocated memory
  free(cmask);
  free(atomNames);
  free(residueNames);
  free(atomTypes);

  return imask;
}

//-----------------------------------------------------------------------------
// FOpenSafe: open a new file, if and only if it does not already exist when
//            file overwriting is not permitted.                     
//                                                                      
// Arguments:                                                           
//   fname:  the name of the file to open                               
//   ovrwrt: flag to authorize overwriting (1 permits, 0 restricts)     
//-----------------------------------------------------------------------------
FILE* FOpenSafe(char* fname, int ovrwrt)
{
  FILE* outp;

  if (ovrwrt == 1 || (ovrwrt == 0 && (outp = fopen(fname, "r")) == NULL)) {
    outp = fopen(fname, "w");
  }
  else {
    printf("FOpenSafe >> Error.  File %s already exists.\n", fname);
    exit(1);
  }

  return outp;
}

//-----------------------------------------------------------------------------
// ReadNumericalShorthand: converts a string, which may contain terms such as
//                         "GB", "MB", "gb", or "KB" into a long long int. 
//                         "GB" or "G" or "gb" or "g" equals gigabyte, "MB"
//                         megabtye, etc.  Numbers parsed by this routine
//                        cannot be negative.          
//
// Arguments:                                                           
//   numstr:      the string that will be converted into a number       
//-----------------------------------------------------------------------------
long long int ReadNumericalShorthand(char* numstr)
{
  int i, slen;
  long long int product, tensplace, multiplier;

  // Parse the string for special characters
  slen = strlen(numstr);
  multiplier = 1;
  for (i = 0; i < slen; i++) {
    if (numstr[i] < '0' || numstr[i] > '9') {
      if (numstr[i] == 'G' || numstr[i] == 'g') {
        multiplier = 1073741824;
      }
      else if (numstr[i] == 'M' || numstr[i] == 'm') {
        multiplier = 1048576;
      }
      else if (numstr[i] == 'K' || numstr[i] == 'k') {
        multiplier = 1024;
      }
      else {
        printf("ReadNummericalShorthand >> Error.  Unable to parse %s into "
               "digits.\n", numstr);
        exit(1);
      }
      if (i < slen-2 ||
          (i < slen-1 && !(numstr[i+1] == 'B' || numstr[i+1] == 'b'))) {
        printf("ReadNummericalShorthand >> Error.  Unable to parse %s into "
               "digits.\n", numstr);
        exit(1);
      }
      slen = i;
    }
  }

  // Compose the number
  tensplace = 1;
  product = 0;
  for (i = slen-1; i >= 0; i--) {
    product += (numstr[i] - '0')*tensplace;
    tensplace *= 10;
  }
  product *= multiplier;

  return product;
}

//-----------------------------------------------------------------------------
// DetectBinaryFile: runs a simple test to verify whether a file is ascii text
//                   or instead can only be parsed as a binary file.  Returns
//                   0 for ascii, 1 for binary.
//
// Arguments:
//   fname:          the name of the file
//-----------------------------------------------------------------------------
int DetectBinaryFile(char* fname)
{
  int isbinary;
  long long int il, fsize;
  int testchar;
  FILE *finp;

  // Detect binary files as those having any byte whose value is outside the
  // range [0, 127].
  if ((finp = fopen(fname, "r")) == NULL) {
    printf("DetectBinaryFile >> Error.  File %s does not exist.\n", fname);
    exit(1);
  }
  isbinary = 0;
  while ((testchar = fgetc(finp)) != EOF && isbinary == 0) {
    if (testchar < 0) {
      isbinary = 1;
      break;
    }
  }
  fclose(finp);

  return isbinary;
}

//-----------------------------------------------------------------------------
// FindFileType: function to determine whether a string names a regular file
//               or, instead, a directory or regular expression. 
//               Returns 0 if it's a regular file, 1 if it's a directory, and
//               2 otherwise.
//
// Arguments:
//   path:      the name of the file or directory to test
//-----------------------------------------------------------------------------
int FindFileType(const char* path)
{
  struct stat path_stat;

  if (stat(path, &path_stat) == 0) {
    if (S_ISREG(path_stat.st_mode) == 1) {
      return 0;
    }
    else if (S_ISDIR(path_stat.st_mode) == 1) {
      return 1;
    }
  }

  return 2;
}

//-----------------------------------------------------------------------------
// Alphabetical: function for comparing two strings.  Each must end in '\0'.
//
// Arguments:
//   a, b:       the objects of type WordNdex (see ParseDS.h) to compare
//-----------------------------------------------------------------------------
int Alphabetical(const void *a, const void *b)
{
  int i, aend, bend;
  char *worda, *wordb;

  worda = ((WordNdex*)a)[0].word;
  wordb = ((WordNdex*)b)[0].word;

  i = 0;
  aend = 0;
  bend = 0;
  while (aend == 0 || bend == 0) {
    if (worda[i] < wordb[i]) {
      return -1;
    }
    if (worda[i] > wordb[i]) {
      return 1;
    }
    if (worda[i] == '\0') {
      aend = 1;
    }
    if (wordb[i] == '\0') {
      bend = 1;
    }
    i++;
  }

  return 0;
}

//-----------------------------------------------------------------------------
// NumberComponent: function for comparing the numerical parts of two strings,
//                  which must already have been extracted.
//
// Arguments:
//   a, b:       the objects of type WordNdex (see ParseDS.h) to compare
//-----------------------------------------------------------------------------
int NumberComponent(const void *a, const void *b)
{
  int i, numa, numb;
  char *worda, *wordb;

  numa = ((WordNdex*)a)[0].numeric;
  numb = ((WordNdex*)b)[0].numeric;

  if (numa < numb) {
    return -1;
  }
  if (numa > numb) {
    return 1;
  }

  return 0;
}

//-----------------------------------------------------------------------------
// ExtractDigits: extract from a string characters that are not digits (0-9).
//                The result of either all digits within a string or all non-
//                digits within a string (depending on the flag takedigits)
//                will be stored in another pre-allocated string.
//
// Arguments: 
//   nns:         character array to hold the result
//   mystr:       string to operate on
//   takedigits:  all non-digits will be extracted if this flag is set to 0.
//                All digits will be extracted otherwise.
//-----------------------------------------------------------------------------
static void ExtractDigits(char* nns, char* mystr, int takedigits)
{
  int i, j, slen;

  slen = strlen(mystr);
  j = 0;
  if (takedigits == 0) {
    for (i = 0; i < slen; i++) {
      if (mystr[i] < '0' || mystr[i] > '9') {
        nns[j] = mystr[i];
        j++;
      }
    }
  }
  else {
    for (i = 0; i < slen; i++) {
      if (mystr[i] >= '0' && mystr[i] <= '9') {
        nns[j] = mystr[i];
        j++;
      }
    }
  }
  nns[j] = '\0';
}

//-----------------------------------------------------------------------------
// AlphabetizeCmat: return an array of integers that, when read in order,
//                  gives the alphabetical order of strings in the rows of a
//                  character matrix.  This routine will try to alphabetize
//                  things as cleanly as possible: first organizing them in a
//                  standard alphabetical order, then seeking out blocks that
//                  might differ only by a numerical component and sorting
//                  those separately.
//
// Arguments:
//   C:        the character matrix to alphabetize
//-----------------------------------------------------------------------------
int* AlphabetizeCmat(cmat *C)
{
  int i, j, k, sdiff;
  int* order;
  char* nni;
  char* nnj;
  WordNdex* pods;

  // Create a scratch array with pointers into the character matrix and
  // numerical tags to go along with 'em
  pods = (WordNdex*)malloc(C->row*sizeof(WordNdex));

  // Load up the pods array, sort, and read back the results
  for (i = 0; i < C->row; i++) {
    pods[i].pos = i;
    pods[i].word = C->map[i];
  }
  qsort(pods, C->row, sizeof(WordNdex), Alphabetical);

  // Sub-sort based on common non-numerical components
  nni = (char*)malloc(C->col*sizeof(char));
  nnj = (char*)malloc(C->col*sizeof(char));
  i = 0;
  while (i < C->row) {
    ExtractDigits(nni, pods[i].word, 0);
    j = i;
    sdiff = 0;
    while (j < C->row && sdiff == 0) {
      ExtractDigits(nnj, pods[j].word, 0);
      sdiff = strcmp(nni, nnj);
      j++;
    }
    if (j > i+1) {
      for (k = i; k < j; k++) {
        ExtractDigits(nnj, pods[k].word, 1);
        pods[k].numeric = (strlen(nnj) > 0) ? atoi(nnj) : 0;
      }
      qsort(&pods[i], j-i, sizeof(WordNdex), NumberComponent);
    }
    i = j;
  }

  // Load up the ordering for return
  order = (int*)malloc(C->row*sizeof(int));
  for (i = 0; i < C->row; i++) {
    order[pods[i].pos] = i;
  }

  // Free allocated memory
  free(pods);
  free(nni);
  free(nnj);

  return order;
}

//-----------------------------------------------------------------------------
// InterpretRegExp: break down a regular expression into its components
//
// Arguments:
//   query:     the regular expression to break down
//-----------------------------------------------------------------------------
static regexp InterpretRegExp(char* query)
{
  int i, j, k, qlen, partcon, pos, onpart, nchoice, maxchoice;
  regexp RE;

  // First, figure out whether the start and the end are wildcards
  qlen = strlen(query);
  RE.matchhead = (query[0] == '*') ? 0 : 1;
  RE.matchtail = (query[qlen-1] == '*') ? 0 : 1;

  // Count the number of parts
  onpart = 0;
  RE.nparts = 0;
  for (i = 0; i < qlen; i++) {
    if (query[i] != '*' && onpart == 0) {
      onpart = 1;
      RE.nparts += 1;
    }
    if (query[i] == '*' && onpart == 1) {
      onpart = 0;
    }
  }
  RE.parts = (cmat*)malloc(RE.nparts*sizeof(cmat));
  for (i = 0; i < RE.nparts; i++) {
    RE.parts[i] = CreateCmat(qlen, 128);
  }

  // Scan over all parts and log them
  onpart = 0;
  partcon = 0;
  for (i = 0; i < qlen; i++) {
    if (query[i] != '*') {
      pos = 0;
      j = i;
      maxchoice = 1;
      while (j < qlen && query[j] != '*') {

        // If we hit a '[', seek the corresponding ']'.  If no such thing
        // can be found, take '[' as just another character.  If the closing
        // bracket can be found, take all the characters in between as
        // choices 
        if (query[j] == '[') {
          k = j;
          while (k < qlen && query[k] != ']' && query[k] != '*') {
            k++;
          }
          if (k < qlen) {
            nchoice = k-j-1;
            if (nchoice > 0) {
              for (k = 0; k < nchoice; k++) {
                RE.parts[partcon].map[pos][k] = query[j+1+k];
              }
              pos++;
              if (nchoice > maxchoice) {
                maxchoice = nchoice;
              }
            }
            j += nchoice + 1;
          }
          else {

            // This is an odd case in which '[' really does appear to be
            // part of the word.
            RE.parts[partcon].map[pos][0] = query[j];
            pos++;
          }
        }
        else {
          RE.parts[partcon].map[pos][0] = query[j];
          pos++;
        }
        j++;
      }
      RE.parts[partcon] = ReallocCmat(&RE.parts[partcon], pos, maxchoice+1);
      partcon++;
      i = j;
    }
  }

  return RE;
}

//-----------------------------------------------------------------------------
// MatchSubstringRegExp: match a substring with a regular expression part.
//                       Returns 1 if the substring can be matched, 0 if not.
//
// Arguments:
//   S:       the substring to match
//   slen:    the length of S
//   C:       the character matrix to use in matching
//-----------------------------------------------------------------------------
static int MatchSubstringRegExp(char* S, cmat *C)
{
  int i, j, match;

  for (i = 0; i < C->row; i++) {
    match = 0;
    for (j = 0; j < C->col; j++) {
      if (C->map[i][j] == '?' || S[i] == C->map[i][j]) {
        match = 1;
      }
    }
    if (match == 0) {
      return 0;
    }
  }

  return 1;
}

//-----------------------------------------------------------------------------
// JoinPath: joint two pieces of a path, similar to the python function
//           os.path.join.
//
// Arguments:
//   newpath:    the (pre-allocated) array to hold the joined path
//   path1:      path of the parent directory (before the slash)
//   slash:      slash character to invoke
//   path2:      path of the file or continuation directory (after the slash)
//-----------------------------------------------------------------------------
void JoinPath(char* newpath, char* path1, char slash, char* path2)
{
  int i;

  if (path1[strlen(path1)-1] == slash) {
    sprintf(newpath, "%s%s", path1, path2);
  }
  else {
    sprintf(newpath, "%s%c%s", path1, slash, path2);
  }
}

//-----------------------------------------------------------------------------
// DirMatchRegExp: list all files (including directories) within a directory
//                 matching the target regular expression.
//
// Arguments:
//   query:     the regular expression to match
//   rootdir:   the root directory in which to search
//   namelim:   the maximum number of names for which space has been allocated
//   nameL:     the list of names into which hits matching the query will go
//-----------------------------------------------------------------------------
static void DirMatchRegExp(char* rootdir, regexp* levelkeys, int mylevel,
                           int nlevel, char slash, cmat *nameL, int *nnames)
{
  int i, j, pass, nmlen;
  char* newrootdir;
  regexp *RE;
  DIR *dir;
  struct dirent *ent;

  // Allocate for a growing root directory name
  newrootdir = (char*)malloc(MAXNAME*sizeof(char));
  RE = &levelkeys[mylevel];

  // Return immediately if the root directory cannot be found
  if ((dir = opendir(rootdir)) == NULL) {
    return;
  }
  while ((ent = readdir(dir)) != NULL) {

    // Skip the "right here" and "go back" entries
    if (strcmp(ent->d_name, "..") == 0 || strcmp(ent->d_name, ".") == 0) {
      continue;
    }

    // There could be a case where the regular expression was '*',
    // meaning that there are no parts to match and anything goes.
    if (RE->nparts == 0) {

      // Have we reached the end?  If so, check that this is really a file,
      // not some other directory.  If not, make a recursive call.
      if (mylevel == nlevel-1) {
        JoinPath(nameL->map[*nnames], rootdir, slash, ent->d_name);
        if (FindFileType(nameL->map[*nnames]) == 0) {
          *nnames += 1;
          if (*nnames == nameL->row) {
            *nameL = ReallocCmat(nameL, 2*nameL->row, nameL->col); 
          }
        }
      }
      else {
        JoinPath(newrootdir, rootdir, slash, ent->d_name);
        DirMatchRegExp(newrootdir, levelkeys, mylevel+1, nlevel, slash, nameL,
                       nnames);
      }
      continue;
    }

    // Match all names against the regular expression at the right level.
    // Step through all of the regular expression's parts and increment
    // the counter npos to make sure everything gets matched in order.
    nmlen = strlen(ent->d_name) - RE->parts[RE->nparts-1].row;

    // If the tail has to match, then the final part
    // of the name should be tested now.
    if (RE->matchtail == 1) {
      pass = MatchSubstringRegExp(&ent->d_name[nmlen],
                                  &RE->parts[RE->nparts-1]);
      if (pass == 0) {
        continue;
      }
    }
    else {
      pass = 1;
    }

    // Step through the file name and verify that it matches
    i = 0;
    j = 0;
    while (i < nmlen && j < RE->nparts - RE->matchtail && pass == 1) {
      if (MatchSubstringRegExp(&ent->d_name[i], &RE->parts[j]) == 1) {
        i += RE->parts[j].row;
        j++;
      }
      else {

        // If the head had to match, we have to have advanced by now.
        if (RE->matchhead == 1 && i == 0) {
          pass = 0;
        }
        i++;
      }
    }

    // If not all parts got matched, this didn't work.
    if (j < RE->nparts - RE->matchtail) {
      continue;
    }

    // If all parts did get matched, descend further
    // into the tree or log the new name.
    if (mylevel == nlevel-1) {
      JoinPath(nameL->map[*nnames], rootdir, slash, ent->d_name);
      if (FindFileType(nameL->map[*nnames]) == 0) {
        *nnames += 1;
        if (*nnames == nameL->row) {
          *nameL = ReallocCmat(nameL, 2*nameL->row, nameL->col); 
        }
      }
    }
    else {
      JoinPath(newrootdir, rootdir, slash, ent->d_name);
      DirMatchRegExp(newrootdir, levelkeys, mylevel+1, nlevel, slash, nameL,
                     nnames);
    }
  }

  // Free allocated memory
  free(newrootdir);

  // Close the root directory
  closedir(dir);
}

//-----------------------------------------------------------------------------
// RegExpFileSearch: search for all regular files matching an input string.
//                   The list of files will be returned in an alphabetically
//                   ordered character matrix.
//
// Arguments:
//   query:          the regular expression search query
//-----------------------------------------------------------------------------
cmat RegExpFileSearch(char* query)
{
  int h, i, j, istart, iend, nlevel, nsearch, qlen, maxlevel, leadingslash;
  int nchar, nfi;
  char slash;
  char* rootdir;
  cmat levels, allfi;
  regexp* levelkeys;

  // Platform-specific directory demarcations
#ifdef linux
  slash = '/';
#endif
#ifdef unix
  slash = '/';
#endif
#ifdef _WIN32
  slash = '\\';
#endif
#ifdef _WIN64
  slash = '\\';
#endif

  // First, formulate a list of directory levels in the regular expression
  // by taking whatever appears between slash characters.  Detect leading
  // slashes and log them appropriately.
  nlevel = 0;
  nchar = 0;
  istart = 0;
  leadingslash = 0;
  qlen = strlen(query);
  if (query[0] == slash) {
    leadingslash = 1;
    while (istart < qlen && query[istart] == slash) {
      istart++;
    }
  }
  iend = qlen-1;
  while (query[iend] == slash) {
    iend--;
  }
  maxlevel = 32;
  levels = CreateCmat(maxlevel, MAXNAME);
  for (i = istart; i <= iend; i++) {
    if (query[i] == slash) {
      j = i;
      while (j < iend && query[j] == slash) {
        j++;
      }
      i = j-1;
      levels.map[nlevel][nchar] = '\0';
      nlevel++;
      if (nlevel == maxlevel) {
        maxlevel += 32;
        levels = ReallocCmat(&levels, maxlevel, MAXNAME);
      }
      nchar = 0;
    }
    else {
      levels.map[nlevel][nchar] = query[i];
      nchar++;
    }
  }
  if (nchar > 0) {
    nlevel++;
  }

  // Find the regular expression wildcards '*', '[...]', and '?'.
  levelkeys = (regexp*)malloc(nlevel*sizeof(regexp));
  for (i = 0; i < nlevel; i++) {
    levelkeys[i] = InterpretRegExp(levels.map[i]);
  }

  // Now that the number of searches is known, we must allocate
  // memory for each of them and perform the tree search.
  rootdir = (char*)malloc(MAXNAME*sizeof(char));
  allfi = CreateCmat(32, MAXNAME);
  if (leadingslash == 1) {
    sprintf(rootdir, "%c", slash);
  }
  else {
    sprintf(rootdir, ".%c", slash);
  }
  nfi = 0;
  DirMatchRegExp(rootdir, levelkeys, 0, nlevel, slash, &allfi, &nfi);

  // Free allocated memory 
  free(rootdir);

  // Check for blank lines in allfi to get the true number of entries.
  i = 0;
  while (allfi.map[i][0] != '\0') {
    i++;
  }
  allfi = ReallocCmat(&allfi, i, MAXNAME);

  return allfi;
}

//-----------------------------------------------------------------------------
// DirectoryFileSearch: search a directory for all regular files (do not
//                      recursively search sub-directories).
//
// Arguments:
//   mydir:     the name of the directory to search
//-----------------------------------------------------------------------------
cmat DirectoryFileSearch(char* mydir)
{
  int i, namelim;
  char slash;
  cmat allfi;
  DIR *dir;
  struct dirent *ent;

  // OS-dependent directory demarcation
#ifdef linux
  slash = '/';
#endif
#ifdef unix
  slash = '/';
#endif
#ifdef _WIN32
  slash = '\\';
#endif
#ifdef _WIN64
  slash = '\\';
#endif

  if ((dir = opendir(mydir)) == NULL) {
    allfi = CreateCmat(1, MAXNAME);
    allfi.row = 0;
    return allfi;
  }
  namelim = 32;
  allfi = CreateCmat(namelim, MAXNAME);
  i = 0;
  while ((ent = readdir(dir)) != NULL) {
    if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) {
      continue;
    }
    JoinPath(allfi.map[i], mydir, slash, ent->d_name);
    i++;
    if (i == namelim) {
      namelim *= 2;
      allfi = ReallocCmat(&allfi, namelim, MAXNAME);
    }
  }
  closedir(dir);
  allfi = ReallocCmat(&allfi, i, MAXNAME);

  return allfi;
}

//-----------------------------------------------------------------------------
// DoubleFromLine: grab a double-precision number from a character string.
//                 Returns 0 on success, 1 on error.
//
// Arguments:
//   line:    the character string to search
//   widx:    the index of the word to grab and turn into a number
//-----------------------------------------------------------------------------
int DoubleFromLine(char* line, int widx, double *val)
{
  int errcode;
  cmat Lwords;

  Lwords = ParseWords(line);
  if (Lwords.row >= widx && WordIsNumber(Lwords.map[widx]) == 1) {
    *val = atof(Lwords.map[widx]);
    errcode = 0;
  }
  else {
    errcode = 1;
    *val = 0.0;
  }
  DestroyCmat(&Lwords);

  return errcode;
}

//-----------------------------------------------------------------------------
// ValidRealData: function for determining whether a block of text is valid
//                real data, containing only comments behind the accepted
//                punctuation, commas, or numbers.
//
// Arguments:
//   fname:     the name of the file to read
//-----------------------------------------------------------------------------
int ValidRealData(char* fname)
{
  int i, j, valid;
  cmat cfi, Lwords;

  // Read the file into memory
  cfi = Ascii2Mem(fname, 256, 8, "File not found");

  // Search for anything that's not a number
  valid = 1;
  i = 0;
  while (i < cfi.row && valid == 1) {
    Lwords = ParseWords(cfi.map[i]);
    for (j = 0; j < Lwords.row; j++) {
      if (WordIsNumber(Lwords.map[j]) == 0) {
        valid = 0;
      }
    }
    DestroyCmat(&Lwords);
    i++;
  }

  // Free allocated memory
  DestroyCmat(&cfi);

  return valid;
}

//-----------------------------------------------------------------------------
// ReadRealData: function for taking in a series of doubles from a file.  The
//               result is returned as a double-precision matrix, a row vector.
//
// Arguments:
//   fname:     the name of the file to read
//-----------------------------------------------------------------------------
dmat ReadRealData(char* fname)
{
  int i, j, npt;
  dmat vals;
  cmat cfi, Lwords;

  // Read the file into memory
  cfi = Ascii2Mem(fname, 256, 8, "File not found");
  npt = 0;
  for (i = 0; i < cfi.row; i++) {
    Lwords = ParseWords(cfi.map[i]);
    npt += Lwords.row;
    DestroyCmat(&Lwords);
  }
  vals = CreateDmat(1, npt, 0);
  npt = 0;
  for (i = 0; i < cfi.row; i++) {
    Lwords = ParseWords(cfi.map[i]);
    for (j = 0; j < Lwords.row; j++) {
      vals.data[npt] = atof(Lwords.map[j]);
      npt++;
    }
  }

  return vals;
}

//-----------------------------------------------------------------------------
// ReadNamelistsFromFile: function for extrating namelists from a file as an
//                        array of cmat structs.  This is how I should have
//                        done things from the beginning, so future patches
//                        should incorporate this strategy to make it possible
//                        to read commands that span multiple lines.
//
// Arguments:
//   fname:     the name of the file
//-----------------------------------------------------------------------------
nmlgroup ReadNamelistsFromFile(char* fname, int nargs, ...)
{
  int i, j, k, onlist, nw, wlen, maxlen;
  int* nmlen;
  char buffer[64];
  cmat cfi, Lwords, allwords;
  nmlgroup mynml;
  va_list ap;

  // Read the file verbatim into memory, first as a
  // list of lines, then as a list of individual words
  cfi = Ascii2Mem(fname, 512, 2, "Namelist file not found.");
  nw = 0;
  maxlen = 1;
  for (i = 0; i < cfi.row; i++) {
    RemoveWhiteSpace(cfi.map[i], 0, ' ');
    RemoveComments(cfi.map[i]);
    NixCommaCarriage(cfi.map[i]);
    EqualSpace(cfi.map[i]);
    Lwords = ParseWords(cfi.map[i]);
    if (Lwords.map[i][0] == '\0') {
      continue;
    }
    for (j = 0; j < Lwords.row; j++) {
      wlen = strlen(Lwords.map[j]);
      nw++;
      if (wlen > maxlen) {
        maxlen = wlen;
      }
    }
  }
  allwords = CreateCmat(nw, maxlen+1);
  nw = 0;
  for (i = 0; i < cfi.row; i++) {
    Lwords = ParseWords(cfi.map[i]);
    if (Lwords.map[i][0] == '\0') {
      continue;
    }
    for (j = 0; j < Lwords.row; j++) {
      strcpy(allwords.map[nw], Lwords.map[j]);
      nw++;
    }
  }

  // Parse the list of words and extract namelists
  mynml.count = 0;
  j = 0;
  for (i = 0; i < allwords.row; i++) {
    if (allwords.map[i][0] == '&' && strcmp(allwords.map[i], "&end") != 0) {
      j++;
    }
  }
  nmlen = (int*)malloc(2*j*sizeof(int));
  if (nargs > 0) {
    va_start(ap, nargs);
    for (i = 0; i < nargs; i++) {
      strcpy(buffer, va_arg(ap, const char*));
      onlist = 0;
      for (j = 0; j < allwords.row; j++) {
        if (onlist == 0 && allwords.map[j][0] == '&' &&
            strcmp(&allwords.map[j][1], buffer) == 0) {
          nmlen[2*mynml.count] = j+1;
          onlist = 1;
        }
        if (onlist == 1 && (strcmp(allwords.map[j], "&end") == 0 ||
                            strcmp(allwords.map[j], "//") == 0)) {
          nmlen[2*mynml.count+1] = j;
          mynml.count += 1;
          onlist = 0;
        }
      }
    }
  }
  else {
    onlist = 0;
    for (j = 0; j < allwords.row; j++) { 
      if (onlist == 0 && allwords.map[j][0] == '&' &&
          strcmp(allwords.map[j], "&end") == 0) {
        nmlen[2*mynml.count] = j+1;
        onlist = 1;
      }
      if (onlist == 1 && (strcmp(allwords.map[j], "&end") == 0 ||
                          strcmp(allwords.map[j], "//") == 0)) {
        nmlen[2*mynml.count+1] = j;
        mynml.count += 1;
        onlist = 0;
      }
    }
  }

  // Allocate space for each namelist
  mynml.titles = CreateCmat(mynml.count, 64);
  mynml.nml = (cmat*)malloc(mynml.count*sizeof(cmat));
  for (i = 0; i < mynml.count; i++) {
    strcpy(mynml.titles.map[i], &allwords.map[nmlen[2*i]-1][1]);
    mynml.nml[i] = CreateCmat(nmlen[2*i+1]-nmlen[2*i], maxlen+1);
    k = 0;
    for (j = nmlen[2*i]; j < nmlen[2*i+1]; j++) {
      strcpy(mynml.nml[i].map[k], allwords.map[j]);
      k++;
    }
  }

  // Free allocated memory
  free(nmlen);

  return mynml;
}
