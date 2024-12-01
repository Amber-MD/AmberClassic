#include <stdio.h>

#include "errormsg.h"

static	int	errs = 0;

void	errormsg_s( int fatal, char fmt[], char str[] )
{
	char	e_msg[ 256 ];

	sprintf( e_msg, fmt, str );
	errormsg( fatal, e_msg );
}
void	errormsg_2s( int fatal, char fmt[], char str1[], char str2[] )
{
	char	e_msg[ 256 ];

	sprintf( e_msg, fmt, str1, str2 );
	errormsg( fatal, e_msg );
}
void	errormsg_d( int fatal, char fmt[], int i )
{
	char	e_msg[ 256 ];

	sprintf( e_msg, fmt, i );
	errormsg( fatal, e_msg );
}
int	errors( void )
{

	return( errs );
}
