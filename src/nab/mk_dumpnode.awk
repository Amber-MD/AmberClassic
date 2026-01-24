BEGIN	{ 	printf( "#include <stdio.h>\n" );
		printf( "#include \"nab.h\"\n" );
		printf( "#include \"y.tab.h\"\n" );
		printf( "\n" );
		printf( "#define\tINCR\t2\n" );
		printf( "#define\tDUMPFILE\t\"dump.file\"\n" );
		printf( "\n" );
		printf( "static\tFILE\t*dfp = NULL;\n" );
		printf( "void\tdumpnode();\n" );
		printf( "\n" );
		printf( "void\tdumpexpr( fp, np, indent )\n" );
		printf( "FILE	*fp;\n" );
		printf( "NODE_T\t*np;\n" );
		printf( "int\tindent;\n" );
		printf( "{\n" );
		printf( "\n" );
		printf( "\tif( !fp ){\n" );
		printf( "\t\tif( !dfp ){\n" );
		printf( "\t\t\tif( ( dfp = fopen( DUMPFILE, \"w\" ) ) == NULL )\n" );
		printf( "\t\t\treturn;\n" );
		printf( "\t\t}\n" );
		printf( "\t\tfp = dfp;\n" );
		printf( "\t}\n" );
		printf( "\tif( np ){\n" );
		printf( "\t\tdumpnode( fp, np, indent );\n" );
		printf( "\t\tdumpexpr( fp, np->n_left, indent + INCR );\n" );
		printf( "\t\tdumpexpr( fp, np->n_right, indent + INCR );\n" );
		printf( "\t}\n" );
		printf( "}\n" );
		printf( "\n" );
		printf( "void\tdumpnode( fp, np, indent )\n" );
		printf( "FILE\t*fp;\n" );
		printf( "NODE_T\t*np;\n" );
		printf( "int\tindent;\n" );
		printf( "{\n" );
		printf( "\tint\ti;\n" );
		printf( "\n" );
		printf( "\tfprintf( fp, \"%%*s\", indent, \"\" );\n" );
		printf( "\tfprintf( fp, \"%%7d: lf = %%7d rt = %%7d tp,cl,kn = %%2d,%%2d,%%2d \",\n" );
		printf( "\t\tnp, np->n_left, np->n_right, np->n_type, np->n_class, np->n_kind );\n" );
		printf( "\tswitch( np->n_sym ){\n" );
	}
	{ if( $1 ~/^%token/ ){
		printf( "\tcase %s :\n", $2 );
		printf( "\t\tfprintf( fp, \"%s", $2 );
		if( $2 == "SYM_IDENT" || $2 == "SYM_ATTRIBUTE" )
			printf( " = %%s\\n\", np->n_val.v_value.v_cval );\n" );
		else if( $2 == "SYM_INT_LIT" )
			printf( " = %%d\\n\", np->n_val.v_value.v_ival );\n" );
		else if( $2 == "SYM_FLOAT_LIT" )
			printf( " = %%8.3f\\n\", np->n_val.v_value.v_fval );\n" );
		else if( $2 == "SYM_STRING_LIT" )
			printf( " = \\\"%%s\\\"\\n\", np->n_val.v_value.v_cval );\n" );
		else if( $2 == "SYM_TYPE" )
			printf( " = %%d\\n\", np->n_val.v_value.v_ival );\n" );
		else if( $2 == "SYM_INDEX" )
			printf( " = %%d\\n\", np->n_val.v_value.v_ival );\n" );
		else
			printf( "\\n\" );\n" );
		printf( "\t\tbreak;\n" );
	  }
	}
END	{	printf( "\tdefault :\n" );
		printf( "\t\tfprintf( fp, \"dumpnode: unknown symbol %%d\\n\"," );
		printf( " np->n_sym );\n" );
		printf( "\t\tbreak;\n" );
		printf( "\t}\n" );
		printf( "}\n" ); }
