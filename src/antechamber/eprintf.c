// Error reporting facilities and wrappers for standard library functions.
// All function names start with 'e' and for std lib functions follow with
// the std lib function names.
// Based on Kernighan and Pike.

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "eprintf.h"


static char *programname = NULL;


// Allocate and clear memory, report if error.

void *ecalloc(const size_t num, const size_t size)
{
    void *p = calloc(num, size);
    if (p == NULL)
        eprintf("Memory allocation and clearing of (%u) bytes failed.", num*size);
    return p;
}


// Open a file, report if error.

FILE *efopen(const char *filename, const char *mode)
{
    FILE *p = fopen(filename, mode);
    if (p == NULL)
        eprintf("Cannot open file (%s) with mode (%s).", filename, mode);
    return p;
}


// Free the memory for the copied name of the program.

static void efreeprogramname()
{
    free(programname);
}


// Get an environment variable, report if error.

char *egetenv(const char *name)
{
    char *p = getenv(name);
    if (p == NULL)
        eprintf("Environment variable (%s) not found.", name);
    return p;
}


// Allocate memory, report if error.

void *emalloc(const size_t size)
{
    void *p = malloc(size);
    if (p == NULL)
        eprintf("Memory allocation of (%u) bytes failed.", size);
    return p;
}


// Print error message and soldier on.  E.g.:
// antechamber: Non Fatal Error!
// Bla, bla as in argument format.
// If the caller's errno is nonzero, emit the std library errno message.

void enonfatal(const char *format, ...)
{
    int callerserrno = errno;
    va_list args;

    fflush(stdout);
    if (eprogramname() != NULL)
        fprintf(stderr, "%s: Non Fatal Error!\n", eprogramname());

    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);

    // if the caller's errno is nonzero, emit the std library errno message.
    if (callerserrno != 0) {
        fprintf(stderr, "\n%s", strerror(errno));
    }
    fprintf(stderr, "\n");

    errno = callerserrno;
    return;
}


// Print error message and exit.  E.g.:
// antechamber: Fatal Error!
// Bla, bla as in argument format.
// If the caller's errno is nonzero, emit the std library errno message.

void eprintf(const char *format, ...)
{
    int callerserrno = errno;
    va_list args;

    fflush(stdout);
    if (eprogramname() != NULL)
        fprintf(stderr, "%s: Fatal Error!\n", eprogramname());

    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);

    // if the caller's errno is nonzero, emit the std library errno message.
    if (callerserrno != 0) {
        fprintf(stderr, "\n%s", strerror(errno));
    }
    fprintf(stderr, "\n");

    efreeprogramname();
    errno = callerserrno;
    exit(EXIT_FAILURE);
}


// Reallocate memory, report if error.

void *erealloc(void *p, const size_t size)
{
    p = realloc(p, size);
    if (p == NULL)
        eprintf("Memory reallocation of (%u) bytes failed.", size);
    return p;
}


// Duplicate a string, report if error.

char *estrdup(const char *from)
{
    char *to;
    to = (char *) malloc(strlen(from) + 1);
    if (to == NULL)
        eprintf("Memory allocation (%.20s) failed.", from);
    strcpy(to, from);
    return to;
}


// Return the name of the program.

char *eprogramname()
{
    return programname;
}


// Set the name of the program to be used in error messages.

void esetprogramname(const char *name)
{
    programname = estrdup(name);
}


// Execute a system command carefully.
// In particular, flush the standard output streams, use errono, and report errors.

int esystem(const char *command)
{
    int callerserrno = errno;
    fflush(stdout);
    fflush(stderr);
    errno = 0;
    int status = system(command);
    if (status != 0) {
        eprintf("Cannot properly run \"%s\".", command);
    }
    errno = callerserrno;
    return status;
}


