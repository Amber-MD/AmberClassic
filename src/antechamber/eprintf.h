// Error reporting facilities and wrappers for standard library functions.
// All function names start with 'e' and for std lib functions follow with
// the std lib function names.

#ifndef EPRINTF_H
#define EPRINTF_H

#include <stddef.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif

extern void *ecalloc(const size_t num, const size_t size);
extern FILE *efopen(const char *filename, const char *mode);
extern char *egetenv(const char *name);
extern void *emalloc(const size_t n);
extern void enonfatal(const char *format, ...);
extern void eprintf(const char *format, ...);
extern char *eprogramname();
extern void *erealloc(void *p, const size_t n);
extern char *estrdup(const char *from);
extern void esetprogramname(const char *name);
extern int esystem(const char *command);

#ifdef __cplusplus
}
#endif

#endif                          // EPRINTF_H
