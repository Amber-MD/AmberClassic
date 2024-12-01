function update( t, w )
{
	var	i, j, k, v, f, s, c;
	var	found, load, url;

	load = false;
	if( t == "d" ){
		f = document.f_dop;
		c = f.d_current;
		s = f.s_dop;
	}else{
		f = document.f_uop;
		c = f.u_current;
		s = f.s_uop;
	}

	if( w == "s" ){
		i = s.selectedIndex;
		if( i == -1 )
			c.value = " -- Not Set -- ";
		else{
			v = s.options[ i ].value;
			c.value = v;
			if( t == "d" ){
				load = true;
				url = v + ".html";
			}
		}
	}else if( w == "p" ){
		i = s.selectedIndex;
		if( i == -1 )
			c.value = " -- Not Set -- ";
		else if( i == 0 ){
			v = s.options[ i ].value;
			alert( v + " is 1st value" );
		}else{
			i--;
			s.selectedIndex = i;
			v = s.options[ i ].value;
			c.value = v;
			if( t == "d" ){
				load = true;
				url = v + ".html";
			}
		}
		
	}else if( w == "n" ){
		i = s.selectedIndex;
		if( i == s.options.length - 1 ){
			v = s.options[ i ].value;
			alert( v + " is last value" );
		}else{
			i++;
			s.selectedIndex = i;
			v = s.options[ i ].value;
			c.value = v;
			if( t == "d" ){
				load = true;
				url = v + ".html";
			}
		}
	}else if( w == "f" ){
		i = 0; j = s.length - 1;
		for( found = false; i <= j; ){
			k = Math.floor( ( i + j ) / 2 );
			if( s.options[k].value == c.value ){
				found = true; 
				break;
			}else if( s.options[k].value < c.value )
				i = k + 1;
			else
				j = k - 1;
		}
		if( !found )
			alert( "Not Found" );
		else{
			s.selectedIndex = k;
			v = s.options[ k ].value;
			c.value = v;
			if( t == "d" ){
				load = true;
				url = v + ".html";
			}
		}
	}else
		alert( "update: illegal value of w: '" + w + "'" );

	if( load )
		parent.RuleFrame.location = url;
}
