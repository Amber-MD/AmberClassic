#include <stdio.h>
#include <time.h>

/* NAB interface to date and time. Note, repeated calls will leak memory! */

char	*date( void )
{

	static char string[11];
	size_t smax=11;
	time_t now;

	now = time( NULL );
	strftime( string, smax, "%m/%d/%Y", localtime( &now ) );

	return( string );

}

char	*timeofday( void )
{

	static char string[9];
	size_t smax=9;
	time_t now;

	now = time( NULL );
	strftime( string, smax, "%H:%M:%S", localtime( &now ) );

	return( string );

}

#if !(defined(WIN32) || defined(WIN64))
char	*ftime( char *fmt )
{

/*   NAB interface to system routine strftime   */

	static char string[50];
	size_t smax=50;
	time_t now;

	now = time( NULL );
	strftime( string, smax, fmt, localtime( &now ) );
	string[49] = '\0';  /* in case of overflow, no explicit checks here */

	return( string );

}
#endif

