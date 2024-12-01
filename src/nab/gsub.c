#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "nab.h"

/*
 *	sub() uses POSIX rules.  POSIX rules differ from gensub rules only
 *	for the case \q, where q is not a recognized escape.  POSIX simplye
 *	copies the \q, but gensub() drops the \ and copies only the q.
 *
 *
 *	Auth	You type	sub() sees	sub() generates
 *
 *	p,g	      &		   &		the matched text
 *	p,g	    \\&		  \&		a literal '&'
 *	p,g	   \\\\		  \\		a literal '\\'
 *	p,g	  \\\\&		 \\&		a literal '\\', then the matched text
 *	p,g	\\\\\\&		\\\&		a literal "\\&"
 *	p	    \\q		  \q		a literal "\\q"
 *	g	    \\q		  \q		a literal 'q'
*/

#define	EXPBUF_SIZE	256
#define	MAX_MATCH	100

extern	char	*loc1, *loc2;

static	int
cnt_amps(const char *);

static	char	*
mk_substr(char *, int, char *, char *);

int	NAB_gsub(int all, char **r, char **s, char **t)
{
	char	expbuf[EXPBUF_SIZE];
	char	*sp, *tp, *ntp;
	int	circ;
	char	*s_match[MAX_MATCH];
	char	*e_match[MAX_MATCH];
	int	m, n_match = 0;
	int	a, n_amp;
	size_t	s_len, t_len;
	char	*new_s = NULL;
	char	*new_t = NULL;
	size_t	nti, ntimax;

	if(*r == NULL){
		fprintf(stderr, "gsub: ERROR: regexp r is NULL\n");
		return 0;
	}else if(**r == '\0'){
		fprintf(stderr, "gsub: ERROR: regexp r is \"\"\n");
		return 0;
	}
	if(*s == NULL){
		fprintf(stderr, "gsub: ERROR: substitute string s is NULL\n");
		return 0;
	}
	if(*t == NULL){
		fprintf(stderr, "gsub: ERROR: target string t is NULL\n");
		return 0;
	}

	compile(*r, expbuf, &expbuf[EXPBUF_SIZE], '\0');

	n_amp = cnt_amps(*s);
	s_len = strlen(*s);
	t_len = *t ? strlen(*t) : 0;
	circ = **r == '^';
	for(n_match = 0, tp = *t; /* *tp &&*/ step(tp, expbuf); ){

		if(n_match < MAX_MATCH){
			s_match[n_match] = loc1;
			e_match[n_match] = loc2;
			if(n_amp != 0)
				t_len += s_len - n_amp + n_amp * (loc2 - loc1);
			else
				t_len += s_len - (loc2 - loc1);
		}
		n_match++;
		if((circ || !all) && n_match == 1)
			break;
		if(!*tp)
			break;
		tp = loc2 != loc1 ? loc2 : loc2 + 1;
	}
	if(n_match > MAX_MATCH){
		fprintf(stderr, "gsub: too many matches %d: only the first %d substitutions performed\n", n_match, MAX_MATCH);
		n_match = MAX_MATCH;
	}
	// make room for the \0
	t_len += 1;
	new_t = (char *)malloc(t_len * sizeof(char));
	if(new_t == NULL){
		fprintf(stderr, "gsub: can't allocate nt\n");
		n_match = 0;
	}else{
		for(ntp = new_t, tp = *t, m = 0; m < n_match; m++){
			ntimax = s_match[m] - tp;
			for(nti = 0; nti < ntimax; nti++)
				*ntp++ = *tp++;
			new_s = mk_substr(*s, n_amp, s_match[m], e_match[m]);
			if(new_s == NULL){
				n_match = 0;
				goto CLEAN_UP;
			}
			strcpy(ntp, new_s);
			ntp += strlen(new_s);
			free(new_s);
			new_s = NULL;
			tp = e_match[m];
		}
		for( ; *tp; )
			*ntp++ = *tp++;
		*ntp = '\0';
		free(*t);
		*t = new_t;
	}

CLEAN_UP : ;

	if(new_s != NULL)
		free(new_s);

	return n_match;
}

static	int
cnt_amps(const char *s)
{
	const char	*sp;
	int	n_amp;

	for(n_amp = 0, sp = s; *sp; ){
		if(*sp == '\\'){	// look for \&, which is a literal &
			sp++;
			if(*sp == '&')
				sp++;
		}else if(*sp == '&'){	// unescaped & -> what was matched
			n_amp++;
			sp++;
		}else
			sp++;
	}

	return n_amp;
}

static	char	*
mk_substr(char *s, int n_amp, char *s_match, char *e_match)
{
	char	*new_s = NULL;
	size_t	s_new_s;
	int	a, i;
	char	*sp, *nsp, *mp;

	s_new_s = n_amp == 0 ? strlen(s) + 1 : strlen(s) - n_amp + n_amp * (e_match - s_match) + 1;
	new_s = (char *)malloc(s_new_s * sizeof(char));
	if(new_s == NULL)
		fprintf(stderr, "mk_substr: can't allocate new_s\n");
	else{
		for(nsp = new_s, sp = s; *sp; ){
			if(*sp == '\\'){
				sp++;
				if(*sp == '&')			// \& -> insert '&'
					*nsp++ = '&';
				else if(*sp == '\\')		// \\ -> insert '\'
					*nsp++ = '\\';
				else{				// \q -> insert  "\q"
					*nsp++ = '\\';
					*nsp++ = *sp;
				}
			}else if(*sp == '&'){			// insert the match
				for(mp = s_match; mp < e_match; mp++)
					*nsp++ = *mp;
			}else					// insert the current char
				*nsp++ = *sp;
			sp++;
		}
		*nsp = '\0';
	}

	return new_s;
}
