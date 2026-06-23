/* LINTLIBRARY */

#include       <stdio.h>
#include       <ctype.h>
#include       <stdarg.h>
#include       <stdlib.h>
#include       "hybrid36.h"

/*
 *      pdb_sscanf performs similarly to sscanf, execept that fields are of
 *      fixed length and a complete line is always consumed.  The field
 *      width defaults to one.  If the line is shorter than expected then
 *      the default is returned.
 *
 *      No Length modifier characters allowed e.g. %lf (%f = double)
 *
 *              d       get an integer.  Default:  0.
 *              x       get a hexadecimal integer.  Default:  0.
 *              h       get a hybrid36 encoded integer.  Default:  0.
 *              f       get a floating point number (C double).  Default:  0.0.
 *              (space) ignore characters within field
 *              s       get a C string, trailing spaces are stripped (leading
 *                      spaces retained); the field width is used as a limit on
 *                      the string length, the null character is appended
 *                      to the end of the string.  Default:  empty string.
 *              c       get a character(s); no stripping of spaces, nor is
 *                      a null character appended.  Default:  space(s).
 */


void 
pdb_sscanf(char *buffer, char *fmt, ...)
{
        va_list         ap;
        register int    i, field_width;
        char            *s, *t;
        char            c;

        va_start(ap, fmt);
        for (; *fmt != '\0'; fmt++)
                if (*fmt != '%') {
                        if (*buffer != '\0')
                                buffer++;
                } else {

                        /* calculate field_width */
                        field_width = 0;
                        for (++fmt; isdigit(*fmt); fmt++)
                                field_width = field_width * 10 + *fmt - '0';
                        if (field_width == 0)
                                field_width = 1;        /* default */

                        switch (*fmt) {

                        case 'd':                       /* integer */
                                s = buffer;
                                for (i = 0; i < field_width; i++, buffer++) {
                                        if (*buffer == '\0')
                                                break;
                                }
                                c = *buffer;
                                *buffer = '\0';
                                *(va_arg(ap, int *)) = atoi(s);
                                *buffer = c;
                                break;

                        case 'x':                       /* hexadecimal integer */
                                s = buffer;
                                for (i = 0; i < field_width; i++, buffer++) {
                                        if (*buffer == '\0')
                                                break;
                                }
                                c = *buffer;
                                *buffer = '\0';
                                *(va_arg(ap, int *)) = (int)strtol(s, NULL, 16);
                                *buffer = c;
                                break;
                       case 'h':                       /* hybrid36 */
                                s = buffer;
                                for (i = 0; i < field_width; i++, buffer++) {
                                        if (*buffer == '\0') {
                                                field_width = i;
                                                break;
                                        }
                                }
                                int *result = va_arg(ap, int *);
                                if (hy36decode(field_width,s,result)>0) *result=0;
                                break;
                        case 'f':                       /* floating point */
                                s = buffer;
                                for (i = 0; i < field_width; i++, buffer++) {
                                        if (*buffer == '\0')
                                                break;
                                }
                                c = *buffer;
                                *buffer = '\0';
                                *(va_arg(ap, double *)) = atof(s);
                                *buffer = c;
                                break;

                        case 's':                       /* string */
                                s = t = va_arg(ap, char *);
                                for (i = 0; i < field_width; i++) {
                                        if (*buffer == '\0')
                                                break;
                                        *s++ = *buffer++;
                                }
                                *s = '\0';
                                /* remove trailing spaces */
                                while (s > t && isspace(*--s))
                                        *s = '\0';
                                break;

                        case 'c':                       /* character(s) */
                                s = va_arg(ap, char *);
                                for (i = 0; i < field_width; i++)
                                        s[i] = ' ';     /* default */

                                for (i = 0; i < field_width; i++) {
                                        if (*buffer == '\0')
                                                break;
                                        *s++ = *buffer++;
                                }
                                break;

                        case ' ':                       /* space (ignore) */
                                for (i = 0; i < field_width; i++, buffer++)
                                        if (*buffer == '\0')
                                                break;
                                break;

                        default:
                                /* error message ? */
                                break;
                        }
                }
        va_end(ap);
}
