%{
//----------------------------------------------------------------------
//   AtomMask selection syntax:
//
//  expression  := operand
//               | '!' expression
//               | expression '&' expression
//               | expression '|' expression
//               | '(' expression ')'
//               | expression distance_op
//
//  operand     := selection_spec_list | '*'
//
//  selection_spec_list := selection_spec
//               | res_selection  atom_selection      -- implied AND (narrow down)
//               | mol_selection  atom_selection      -- implied AND
//               | mol_selection  res_selection       -- implied AND
//               | mol_res_selection atom_selection   -- implied AND (3-level chain)
//               | selection_spec_list selection_spec -- implied OR (same/other level)
//
//  res_selection   := ':' value_list          -- RESNUM: default (PDB resSeq or flat in compat mode)
//                   | ':;' value_list         -- RES_PDBSEQ: always PDB resSeq
//                   | ':#' value_list         -- RES_INDEX: always flat sequential index
//                   | ':%' name               -- RESTYPE: residue type (single name, always literal)
//                   | '::' value_list         -- chain ID
//                   | '::' value_list ':' value_list  -- chain ID + PDB resSeq (: forced to PDB resSeq)
//                   | ':(' expr ')'           -- residues containing expr
//
//  atom_selection  := '@' value_list          -- atom number or name
//                   | '@%' value_list         -- atom type (always literal)
//                   | '@/' value_list         -- element
//
//  mol_selection   := '^' value_list          -- molecule number
//                   | '^(' expr ')'           -- molecules containing expr
//
//  value_list  := value (',' value)*
//
//  value       := integer
//               | integer '-' integer         -- numeric range
//               | name                        -- glob active, fnmatch semantics
//               | '~' name                   -- forced literal, no glob
//               -- note: * and ? are valid name chars, lex as NAME via NAMESTART
//               -- bare * at expression level means select-all
//
//  distance_op := ('<@' | '>@' |             -- by atom
//                  '<:' | '>:' |             -- by residue
//                  '<;' | '>;' |             -- by residue center
//                  '<^' | '>^')              -- by molecule
//                 number                     -- REAL or INTEGER (promoted to double)
//
//  Residue numbering modes:
//    RESNUM    ':'   default — PDB resSeq normally, flat index in compat mode
//    RES_PDBSEQ':;'  always PDB resSeq (cpptraj called this "original resnum")
//    RES_INDEX ':#'  always flat sequential index
//    RESTYPE   ':%'  residue type — single literal name, not a number
//
//  Implied AND fires when moving DOWN the hierarchy: mol -> res -> atom
//  Implied OR  fires for same-level or upward juxtaposition
//
//  Distance ops are postfix on a primary expression. The left operand is
//  the reference selection; atoms within/beyond the cutoff of any reference
//  atom are selected. The nonbond grid is built lazily on first eval and
//  cached on the AST node for reuse across trajectory frames.
//  dist is stored as-is; the grid setup squares it internally.
//------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "select_mask.h"

static SELNODE mk_node(SELNODEKINDt k, SELNODE l, SELNODE r);
static SELNODE mk_text_node(SELNODEKINDt k, char *text, int forced);
static SELNODE mk_int_node(SELNODEKINDt k, long v);
static SELNODE mk_range(SELNODEKINDt k, long a, long b);
static SELNODE mk_dist_node(SELNODEKINDt k, SELNODE ref, double dist);

extern int sellex(void);
extern void selerror(const char *s);

extern char sel_error_token[];  /* last token seen, set by lexer via YY_USER_ACTION */

SELNODE selection_root = NULL;
%}

%union {
    char   *text;
    long    inum;
    double  rnum;
    SELNODE node;
}

%token <text> NAME           /* string */
%token <inum> INTEGER        /* integer */
%token <rnum> REAL           /* double, used for distance cutoffs */
%token TILDE                 /* ~ force literal string for next token */
%token COLON                 /* :  residue RESNUM (default) */
%token COLON_COLON           /* :: chain ID */
%token COLON_SEMICOLON       /* :; RES_PDBSEQ — always PDB resSeq */
%token COLON_HASH            /* :# RES_INDEX  — always flat sequential index */
%token COLON_PCT             /* :% RESTYPE    — residue type, single literal name */
%token AT                    /* @  atom name or global index (cpptraj compat) */
%token AT_HASH               /* @# explicit global flat atom index */
%token AT_SEMICOLON          /* @; atom index within residue (PDB order, 1-based) */
%token AT_PCT                /* @% atom type (always literal) */
%token AT_SLASH              /* @/ atom element */
%token LESS_AT               /* <@ distance within, by atom */
%token LESS_COLON            /* <: distance within, by residue */
%token LESS_SEMICOLON        /* <; distance within, by residue center */
%token LESS_CARAT            /* <^ distance within, by molecule */
%token GREATER_AT            /* >@ distance beyond, by atom */
%token GREATER_COLON         /* >: distance beyond, by residue */
%token GREATER_SEMICOLON     /* >; distance beyond, by residue center */
%token GREATER_CARAT         /* >^ distance beyond, by molecule */
%token CARAT                 /* ^  molecule number or containment */
%token COMMA                 /* , */
%token AND                   /* & explicit AND */
%token OR                    /* | explicit OR */
%token NOT                   /* ! */
%token LPAREN                /* ( */
%token RPAREN                /* ) */
%token DASH                  /* - */
/* NOTE: no STAR token — * and ? are valid NAMESTART chars and lex as NAME.
 * fnmatch handles glob matching. Bare * at expression level is select-all,
 * handled by checking NAME == "*" in primary.                             */

%type <node> expr term factor primary
%type <node> selection_spec_list selection_spec
%type <node> res_selection atom_selection mol_selection mol_res_selection
%type <node> chain_res_selection
%type <node> res_list_default res_item_default
%type <node> res_list_pdbseq  res_item_pdbseq
%type <node> res_list_index   res_item_index
%type <node> atom_list atom_item
%type <node> atom_list_index   atom_item_index
%type <node> atom_list_residx  atom_item_residx
%type <node> chain_list chain_item
%type <node> elem_list elem_item
%type <node> type_list type_item
%type <node> mol_list mol_item
%type <rnum> number
%type <text> name forced_name any_name

/* Implied operator precedence: AND binds tighter than OR */
%left  IMPLIED_OR
%right IMPLIED_AND

/* Explicit operator precedence (low to high) */
%left OR
%left AND
%right NOT

%%

input
    : expr                          { selection_root = $1; }
    | error                         { selerror("invalid selection expression"); YYERROR; }
    ;

expr
    : expr OR term                  { $$ = mk_node(SEL_NODE_OR,  $1, $3); }
    | term                          { $$ = $1; }
    ;

term
    : term AND factor               { $$ = mk_node(SEL_NODE_AND, $1, $3); }
    | factor                        { $$ = $1; }
    ;

factor
    : NOT factor                    { $$ = mk_node(SEL_NODE_NOT, $2, NULL); }
    | primary                       { $$ = $1; }
    ;

/* --- distance ops are postfix on any primary expression ---
 *
 * expr <@ 5.0  selects all atoms within 5.0 Å of any atom matched by expr.
 * The reference set is the left child; the nonbond grid is built from it
 * lazily on first eval and cached on the node. dist is stored as-is;
 * neighbor_grid_setup() squares it internally.
 * Both REAL and INTEGER are accepted for distance — INTEGER is promoted
 * to double so <@5 and <@5.0 are equivalent.                            */
primary
    : LPAREN expr RPAREN             { $$ = $2; }
    | LPAREN error RPAREN            { selerror("invalid expression inside parentheses");
                                       YYERROR; }
    | NAME                           { if (strcmp($1, "*") != 0) {
                                           selerror("bare name requires : @ ^ prefix — "
                                                    "did you mean :name or @name? "
                                                    "Use * alone to select all atoms.");
                                           free($1); YYERROR;
                                       }
                                       free($1);
                                       $$ = mk_node(SEL_NODE_ALL, NULL, NULL); }
    | selection_spec_list            { $$ = $1; }
    | primary LESS_AT        number  { $$ = mk_dist_node(SEL_NODE_DIST_WITHIN_ATOM,    $1, $3); }
    | primary LESS_COLON     number  { $$ = mk_dist_node(SEL_NODE_DIST_WITHIN_RES,     $1, $3); }
    | primary LESS_SEMICOLON number  { $$ = mk_dist_node(SEL_NODE_DIST_WITHIN_RESCEN,  $1, $3); }
    | primary LESS_CARAT     number  { $$ = mk_dist_node(SEL_NODE_DIST_WITHIN_MOL,     $1, $3); }
    | primary GREATER_AT     number  { $$ = mk_dist_node(SEL_NODE_DIST_BEYOND_ATOM,    $1, $3); }
    | primary GREATER_COLON  number  { $$ = mk_dist_node(SEL_NODE_DIST_BEYOND_RES,     $1, $3); }
    | primary GREATER_SEMICOLON number { $$ = mk_dist_node(SEL_NODE_DIST_BEYOND_RESCEN,$1, $3); }
    | primary GREATER_CARAT  number  { $$ = mk_dist_node(SEL_NODE_DIST_BEYOND_MOL,     $1, $3); }
    | primary LESS_AT        NAME    { selerror("distance requires a number e.g. <@5.0");
                                       free($3); YYERROR; }
    | primary GREATER_AT     NAME    { selerror("distance requires a number e.g. >@5.0");
                                       free($3); YYERROR; }
    ;

/* number: accept REAL or INTEGER, promote INTEGER to double.
 * Allows <@5 and <@5.0 to be equivalent.                     */
number
    : REAL                          { $$ = $1; }
    | INTEGER                       { $$ = (double)$1; }
    ;

selection_spec_list
    : selection_spec
        { $$ = $1; }
    | res_selection  atom_selection      %prec IMPLIED_AND
        { $$ = mk_node(SEL_NODE_AND, $1, $2); }
    | mol_selection  atom_selection      %prec IMPLIED_AND
        { $$ = mk_node(SEL_NODE_AND, $1, $2); }
    | mol_res_selection                  %prec IMPLIED_AND
        { $$ = $1; }
    | mol_res_selection atom_selection   %prec IMPLIED_AND
        { $$ = mk_node(SEL_NODE_AND, $1, $2); }
    | selection_spec_list selection_spec %prec IMPLIED_OR
        { $$ = mk_node(SEL_NODE_OR,  $1, $2); }
    ;

/* three-level chain: ^1:ALA@CA — mol+res reduces here first via its own
 * %prec so bison doesn't confuse it with selection_spec_list mol res.
 * The combined node can then accept an atom_selection for the third level,
 * or stand alone as a two-level mol+res selection.                        */
mol_res_selection
    : mol_selection res_selection        %prec IMPLIED_AND
        { $$ = mk_node(SEL_NODE_AND, $1, $2); }
    ;

/* --- typed selection specs --- */

res_selection
    : COLON           res_list_default  { $$ = $2; }
    | COLON_SEMICOLON res_list_pdbseq   { $$ = $2; }
    | COLON_HASH      res_list_index    { $$ = $2; }
    | COLON_PCT       name              { $$ = mk_text_node(SEL_NODE_RESTYPE, $2, 1); }
    | COLON_PCT       INTEGER           { selerror(":% residue type requires a name, not a number — "
                                                   "did you mean :# for residue index?");
                                          YYERROR; }
    | COLON_COLON     chain_list        { $$ = $2; }
    | chain_res_selection               { $$ = $1; }
    | COLON LPAREN expr RPAREN          { $$ = mk_node(SEL_NODE_RES_CONTAINS, $3, NULL); }
    ;

/* '::' chain ':' resnum — ':' is always PDB resSeq here regardless of
 * cpptraj_compat, since an explicit chain ID implies PDB coordinate context.
 * ::A:#17 and ::A:%name are handled via implied AND since :# and :% are
 * already unambiguous regardless of context.                              */
chain_res_selection
    : COLON_COLON chain_list COLON res_list_pdbseq
        { $$ = mk_node(SEL_NODE_AND, $2, $4); }
    ;

atom_selection
    : AT           atom_list           { $$ = $2; }
    | AT_HASH      atom_list_index     { $$ = $2; }   /* explicit global flat index */
    | AT_SEMICOLON atom_list_residx    { $$ = $2; }   /* index within residue, PDB order */
    | AT_PCT       type_list           { $$ = $2; }
    | AT_SLASH     elem_list           { $$ = $2; }
    ;

/* --- atom index lists --- */

atom_list_index
    : atom_item_index                           { $$ = $1; }
    | atom_list_index COMMA atom_item_index     { $$ = mk_node(SEL_NODE_OR, $1, $3); }
    ;

atom_item_index
    : INTEGER               { $$ = mk_int_node(SEL_NODE_ATOM_INDEX,       $1); }
    | INTEGER DASH INTEGER  { $$ = mk_range(   SEL_NODE_RANGE_ATOM_INDEX, $1, $3); }
    | NAME                  { selerror("@# requires integer — use @ for atom name");
                              free($1); YYERROR; }
    ;

/* --- atom index within residue (PDB order, 1-based) --- */

atom_list_residx
    : atom_item_residx                          { $$ = $1; }
    | atom_list_residx COMMA atom_item_residx   { $$ = mk_node(SEL_NODE_OR, $1, $3); }
    ;

atom_item_residx
    : INTEGER               { $$ = mk_int_node(SEL_NODE_ATOM_RESIDX,       $1); }
    | INTEGER DASH INTEGER  { $$ = mk_range(   SEL_NODE_RANGE_ATOM_RESIDX, $1, $3); }
    | NAME                  { selerror("@; requires integer — atom index within residue");
                              free($1); YYERROR; }
    ;

mol_selection
    : CARAT mol_list                    { $$ = $2; }
    | CARAT LPAREN expr RPAREN          { $$ = mk_node(SEL_NODE_MOL_CONTAINS, $3, NULL); }
    ;

selection_spec
    : res_selection                     { $$ = $1; }
    | atom_selection                    { $$ = $1; }
    | mol_selection                     { $$ = $1; }
    ;

/* --- residue lists — three separate sets to carry numbering scheme ---
 *
 * Each scheme gets its own list, item, and range node kinds so the
 * evaluator knows which numbering to use without inspecting context.
 * default (RESNUM) allows names since : takes both names and numbers.
 * pdbseq and index are number-only contexts.                          */

res_list_default
    : res_item_default                           { $$ = $1; }
    | res_list_default COMMA res_item_default    { $$ = mk_node(SEL_NODE_OR, $1, $3); }
    ;

res_item_default
    : name                  { $$ = mk_text_node(SEL_NODE_RESNAME,      $1, 0); }
    | forced_name           { $$ = mk_text_node(SEL_NODE_RESNAME,      $1, 1); }
    | INTEGER               { $$ = mk_int_node( SEL_NODE_RESNUM,       $1); }
    | INTEGER DASH INTEGER  { $$ = mk_range(    SEL_NODE_RANGE_RESNUM, $1, $3); }
    ;

res_list_pdbseq
    : res_item_pdbseq                          { $$ = $1; }
    | res_list_pdbseq COMMA res_item_pdbseq    { $$ = mk_node(SEL_NODE_OR, $1, $3); }
    ;

res_item_pdbseq
    : INTEGER               { $$ = mk_int_node( SEL_NODE_RES_PDBSEQ,       $1); }
    | INTEGER DASH INTEGER  { $$ = mk_range(    SEL_NODE_RANGE_RES_PDBSEQ, $1, $3); }
    | NAME                  { selerror(":; PDB resSeq requires integer — "
                                       "use : for residue name selection");
                              free($1); YYERROR; }
    ;

res_list_index
    : res_item_index                         { $$ = $1; }
    | res_list_index COMMA res_item_index    { $$ = mk_node(SEL_NODE_OR, $1, $3); }
    ;

res_item_index
    : INTEGER               { $$ = mk_int_node( SEL_NODE_RES_INDEX,       $1); }
    | INTEGER DASH INTEGER  { $$ = mk_range(    SEL_NODE_RANGE_RES_INDEX, $1, $3); }
    | NAME                  { selerror(":# residue index requires integer — "
                                       "use : for residue name selection");
                              free($1); YYERROR; }
    ;

/* --- atom list --- */

atom_list
    : atom_item                         { $$ = $1; }
    | atom_list COMMA atom_item         { $$ = mk_node(SEL_NODE_OR, $1, $3); }
    ;

atom_item
    : name                  { $$ = mk_text_node(SEL_NODE_ATOMNAME,   $1, 0); }
    | forced_name           { $$ = mk_text_node(SEL_NODE_ATOMNAME,   $1, 1); }
    | INTEGER               { $$ = mk_int_node( SEL_NODE_INDEX,      $1); }
    | INTEGER DASH INTEGER  { $$ = mk_range(    SEL_NODE_RANGE_ATOM, $1, $3); }
    ;

/* --- chain list --- */

chain_list
    : chain_item                        { $$ = $1; }
    | chain_list COMMA chain_item       { $$ = mk_node(SEL_NODE_OR, $1, $3); }
    ;

chain_item
    : any_name              { $$ = mk_text_node(SEL_NODE_CHAINID, $1, 0); }
    | INTEGER               { /* Chain IDs have no integer index meaning so promote
                               * to string. Only context where INTEGER -> string is
                               * unambiguous. Note: atol() loses leading zeros so
                               * ::01 still requires ::~01.                        */
                              char buf[32];
                              snprintf(buf, sizeof(buf), "%ld", $1);
                              $$ = mk_text_node(SEL_NODE_CHAINID, strdup(buf), 1); }
    ;

/* --- element list --- */

elem_list
    : elem_item                         { $$ = $1; }
    | elem_list COMMA elem_item         { $$ = mk_node(SEL_NODE_OR, $1, $3); }
    ;

elem_item
    : any_name              { $$ = mk_text_node(SEL_NODE_ELEMENT,    $1, 0); }
    | INTEGER               { $$ = mk_int_node( SEL_NODE_ELEMENT_NUM,$1); }
    | INTEGER DASH INTEGER  { $$ = mk_range(    SEL_NODE_RANGE_ELEM, $1, $3); }
    ;

/* --- atom type list (always forced literal — no glob in type names) --- */

type_list
    : type_item                         { $$ = $1; }
    | type_list COMMA type_item         { $$ = mk_node(SEL_NODE_OR, $1, $3); }
    ;

type_item
    : any_name              { $$ = mk_text_node(SEL_NODE_ATOMTYPE, $1, 1); }
    | INTEGER               { selerror("@% atom type requires a name, not a number");
                              YYERROR; }
    ;

/* --- molecule list (numbers and ranges only — no string names) --- */

mol_list
    : mol_item                          { $$ = $1; }
    | mol_list COMMA mol_item           { $$ = mk_node(SEL_NODE_OR, $1, $3); }
    ;

mol_item
    : INTEGER               { $$ = mk_int_node(SEL_NODE_MOLNUM,    $1); }
    | INTEGER DASH INTEGER  { $$ = mk_range(   SEL_NODE_RANGE_MOL, $1, $3); }
    | NAME                  { selerror("^ molecule selector requires integer — "
                                       "molecule names are not supported");
                              free($1); YYERROR; }
    ;

/* --- name nonterminals --- */

name
    : NAME                  { $$ = $1; }
    ;

forced_name
    /* The lexer handles ~ as a unit: ~017 is emitted as TILDE NAME("017"),
     * never as TILDE INTEGER. This preserves leading zeros so ~017 != ~17.
     * The parser never sees TILDE INTEGER — adding that rule would corrupt
     * numeric names by routing through atol() and losing leading zeros.
     * The 'forced' flag in mk_text_node() disables glob pattern matching:
     * ~ means "literal string, no wildcards, no fnmatch".                */
    : TILDE NAME            { $$ = $2; }
    | TILDE INTEGER         { selerror("~ followed by bare integer — "
                                       "leading zeros are preserved so ~017 != ~17, "
                                       "but the lexer should have emitted this as NAME. "
                                       "This is an internal error.");
                              YYERROR; }
    ;

any_name
    : name                  { $$ = $1; }
    | forced_name           { $$ = $1; }
    ;

%%

void selerror(const char *s)
{
    VPFATAL("selection parse error near '%s': %s\n", sel_error_token, s);
}

static SELNODE alloc_node(SELNODEKINDt k)
{
    SELNODE n = (SELNODE)calloc(1, sizeof(SELNODEt));
    if (!n) abort();
    n->kind = k;
    n->cache_object = SEL_CACHE_INVALID;
    return n;
}

static SELNODE mk_node(SELNODEKINDt k, SELNODE l, SELNODE r)
{
    SELNODE n = alloc_node(k);
    n->left = l;
    n->right = r;
    return n;
}

static SELNODE mk_text_node(SELNODEKINDt k, char *text, int forced)
{
    SELNODE n = alloc_node(k);
    n->text = text;
    n->forced_string = forced;
    n->has_glob = !forced && (strchr(text, '*') || strchr(text, '?') ||
                              strchr(text, '['));
    return n;
}

static SELNODE mk_int_node(SELNODEKINDt k, long v)
{
    SELNODE n = alloc_node(k);
    n->a = v;
    return n;
}

static SELNODE mk_range(SELNODEKINDt k, long a, long b)
{
    SELNODE n = alloc_node(k);
    n->a = a;
    n->b = b;
    return n;
}

static SELNODE mk_dist_node(SELNODEKINDt k, SELNODE ref, double dist)
{
    /* dist stored as-is — neighbor_grid_setup() squares it internally.
     * ngGrid is left NULL by calloc — built lazily on first eval and
     * cached here for reuse across trajectory frames. Caller must call
     * sel_invalidate_coords(root) between frames if coordinates change. */
    SELNODE n = alloc_node(k);
    n->left = ref;
    n->dist = dist;
    return n;
}
