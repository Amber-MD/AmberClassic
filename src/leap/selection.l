%{
#include "selection.h"
#include "selection.tab.h"
#include <stdlib.h>
#include <string.h>

int selparse(void);
%}

%option prefix="sel"
%option noyywrap noinput nounput

%x TILDE_MODE

WS         [ \t\r\n]+
REALNUM    [0-9]+[.][0-9]+
NUM        [0-9]+

  /* BRACKET: fnmatch character class or literal escape within a name.
   * [CA] = C or A (character class, passed to fnmatch as-is)
   * [*]  = literal asterisk, [?] = literal question mark, [-] = literal dash.
   * Full POSIX bracket expression syntax is supported via fnmatch().
   * Note: [[:upper:]] style collating classes work but are unlikely to be
   * needed in practice — [A-Z] covers the same ground.                   */
BRACKET    \[[^\]]*\]

  /* A name starts with a letter, single-quote (nucleic acid atoms: O3', C5'),
   * a glob wildcard, or a bracket expression e.g. [ABC]CA.
   * Leading digits are NOT allowed — bare integers lex as INTEGER for index
   * selection. Legacy leading-digit atom names like 1HG1 require ~: ~1HG1.
   * Continuation allows letters, digits, and common punctuation found in
   * real molecular names. # added for force field atom types.            */
NAMESTART  [A-Za-z'*?]
NAMECONT   [A-Za-z0-9_+\-*.'?#]
NAME       ({BRACKET}|{NAMESTART})({NAMECONT}|{BRACKET})*

  /* Tilde literal: like NAME but - is just another char (not DASH token),
   * and the only terminators are structural delimiters that start a new
   * selection spec. Bracket escapes work identically to NAME.
   * + not * so that empty tilde literal is caught by the error rule.    */
TILDECHAR  [^ \t\r\n,()&|!@:^]
TILDENAME  ({TILDECHAR}|{BRACKET})+

%%

{WS}                    ;

"("                     return LPAREN;
")"                     return RPAREN;
","                     return COMMA;
"::"                    return COLON_COLON;      /* chain ID               */
":%"                    return COLON_PCT;         /* residue type           */
":;"                    return COLON_SEMICOLON;   /* RES_PDBSEQ: always PDB resSeq */
":#"                    return COLON_HASH;        /* RES_INDEX: always flat sequential index */
":"                     return COLON;             /* PDB resSeq (default)   */
"^"                     return CARAT;             /* molecule index         */
"@%"                    return AT_PCT;            /* atom type              */
"@/"                    return AT_SLASH;          /* atom element           */
"@"                     return AT;               /* atom name/index        */
"<@"                    return LESS_AT;
"<:"                    return LESS_COLON;
"<;"                    return LESS_SEMICOLON;
"<^"                    return LESS_CARAT;
">@"                    return GREATER_AT;
">:"                    return GREATER_COLON;
">;"                    return GREATER_SEMICOLON;
">^"                    return GREATER_CARAT;
"!"                     return NOT;
"&"                     return AND;
"|"                     return OR;
"-"                     return DASH;

"~"                     { BEGIN(TILDE_MODE); return TILDE; }

{REALNUM}               {
                            sellval.rnum = atof(yytext);
                            return REAL;
                        }

{NUM}                   {
                            sellval.inum = atol(yytext);
                            return INTEGER;
                        }

{NAME}                  {
                            sellval.text = strdup(yytext);
                            if (!sellval.text) abort();
                            return NAME;
                        }

  /* ------------------------------------------------------------------ */
  /* TILDE_MODE: after ~, consume one forced-literal name token.         */
  /*                                                                     */
  /* ~ exists because:                                                   */
  /*   (1) names like "017" must not be parsed as integer index 17 —    */
  /*       leading zeros are significant (e.g. modified nucleotide 017) */
  /*   (2) '-' can appear in residue names so cannot be a terminator    */
  /*   (3) single-quote cannot serve as a quote char since it is valid  */
  /*       in atom names (O3', C5' in nucleic acids)                    */
  /*   (4) double-quote adds complexity for no real gain                */
  /*                                                                     */
  /* The forced flag propagated to mk_text_node() disables glob matching */
  /* on the node — ~ means "literal string, no wildcards".              */
  /*                                                                     */
  /* Note: the lexer delivers the raw string including any [..] brackets */
  /* to the parser unchanged. Bracket interpretation (character classes  */
  /* vs literal escapes) is the pattern evaluator's responsibility.     */
  /* ------------------------------------------------------------------ */
<TILDE_MODE>{WS}        ;

<TILDE_MODE>{TILDENAME} {
                            sellval.text = strdup(yytext);
                            if (!sellval.text) abort();
                            BEGIN(INITIAL);
                            return NAME;
                        }

<TILDE_MODE>.           {
                            /* Empty or invalid tilde literal — structural
                             * delimiter immediately followed ~            */
                            fprintf(stderr, "selection lex error: empty ~ literal\n");
                            BEGIN(INITIAL);
                            return 0;
                        }

.                       { return yytext[0]; }

%%

extern SELNODE selection_root;
SELNODE selParseAtomMask(const char *mask) {
    YY_BUFFER_STATE b;
    SELNODE root;
    b = sel_scan_string(mask);
    selection_root = NULL;
    if (selparse() != 0) {
        sel_delete_buffer(b);
        return NULL;
    }
    root = selection_root;
    sel_delete_buffer(b);
    return root;
}
