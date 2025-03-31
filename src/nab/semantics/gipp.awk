nawk 'BEGIN{
	DEBUG = 0
	nt = nr = 0
}
/OP_/ { nt++ ; nr = 0; next }
$0 ~ /[0-9]+ T_/ {
	nr++
	gsub( "\\|   ","| 0 ")
	if ( DEBUG ) print
	data[ nt FS nr ] = $0
	next
}
END{
# scan one element to get types
	for ( r = 1; r <= nr; r++ ) {
		ne = split( data[ 1 FS r ], a, "|" )
		split( a[1], b )
		type[b[1]+0] = b[2]
	}
	if ( DEBUG )
	for ( r in type )
		print r, type[r]
		
	printf "static int crud[][%d][%d] = {\n", nr, nr
	for ( t = 1; t <= nt; t++ ) {
		print "{"
		for ( r = 1; r <= nr; r++ ) {
			ne = split( data[ t FS r ], a, "|" ) - 1
			printf "{"
			for ( e = 2; e <= ne; e++ ) {
				q = a[e]+0
				if ( !( q in type ) )
					printf " %s", "T_SPECIAL"
				else
					printf " %s", type[a[e]+0]
				if ( e != ne ) printf ","
			}
			printf "}"
			if ( r != nr ) printf ","
			printf "\n"
		}
		printf "}"
		if ( t != nt ) printf ","
		printf "\n"
	}
	print "};"
}' $*

