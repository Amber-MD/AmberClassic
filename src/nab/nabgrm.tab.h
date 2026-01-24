/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

#ifndef YY_YY_NABGRM_TAB_H_INCLUDED
# define YY_YY_NABGRM_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token kinds.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    YYEMPTY = -2,
    YYEOF = 0,                     /* "end of file"  */
    YYerror = 256,                 /* error  */
    YYUNDEF = 257,                 /* "invalid token"  */
    SYM_ADDRESS = 258,             /* SYM_ADDRESS  */
    SYM_ALLOCATE = 259,            /* SYM_ALLOCATE  */
    SYM_AND = 260,                 /* SYM_AND  */
    SYM_ASSERT = 261,              /* SYM_ASSERT  */
    SYM_ASSIGN = 262,              /* SYM_ASSIGN  */
    SYM_ATOM = 263,                /* SYM_ATOM  */
    SYM_ATSIGN = 264,              /* SYM_ATSIGN  */
    SYM_ATTRIBUTE = 265,           /* SYM_ATTRIBUTE  */
    SYM_BOUNDS = 266,              /* SYM_BOUNDS  */
    SYM_BREAK = 267,               /* SYM_BREAK  */
    SYM_CALL = 268,                /* SYM_CALL  */
    SYM_COMMA = 269,               /* SYM_COMMA  */
    SYM_CONTINUE = 270,            /* SYM_CONTINUE  */
    SYM_DEALLOCATE = 271,          /* SYM_DEALLOCATE  */
    SYM_DEBUG = 272,               /* SYM_DEBUG  */
    SYM_DECL = 273,                /* SYM_DECL  */
    SYM_DELETE = 274,              /* SYM_DELETE  */
    SYM_DONT_MATCH = 275,          /* SYM_DONT_MATCH  */
    SYM_DYNAMIC = 276,             /* SYM_DYNAMIC  */
    SYM_ELSE = 277,                /* SYM_ELSE  */
    SYM_EQUAL = 278,               /* SYM_EQUAL  */
    SYM_ERROR = 279,               /* SYM_ERROR  */
    SYM_FILE = 280,                /* SYM_FILE  */
    SYM_FLOAT = 281,               /* SYM_FLOAT  */
    SYM_FLOAT_LIT = 282,           /* SYM_FLOAT_LIT  */
    SYM_FOR = 283,                 /* SYM_FOR  */
    SYM_FOREACH = 284,             /* SYM_FOREACH  */
    SYM_GREATER = 285,             /* SYM_GREATER  */
    SYM_GREATER_EQUAL = 286,       /* SYM_GREATER_EQUAL  */
    SYM_HASHED = 287,              /* SYM_HASHED  */
    SYM_IDENT = 288,               /* SYM_IDENT  */
    SYM_IF = 289,                  /* SYM_IF  */
    SYM_IN = 290,                  /* SYM_IN  */
    SYM_INDEX = 291,               /* SYM_INDEX  */
    SYM_INDIRECT = 292,            /* SYM_INDIRECT  */
    SYM_INT = 293,                 /* SYM_INT  */
    SYM_INT_LIT = 294,             /* SYM_INT_LIT  */
    SYM_LBRACE = 295,              /* SYM_LBRACE  */
    SYM_LBRACK = 296,              /* SYM_LBRACK  */
    SYM_LESS = 297,                /* SYM_LESS  */
    SYM_LESS_EQUAL = 298,          /* SYM_LESS_EQUAL  */
    SYM_LIST = 299,                /* SYM_LIST  */
    SYM_LPAREN = 300,              /* SYM_LPAREN  */
    SYM_MATCH = 301,               /* SYM_MATCH  */
    SYM_MATRIX = 302,              /* SYM_MATRIX  */
    SYM_MINUS = 303,               /* SYM_MINUS  */
    SYM_MINUS_ASSIGN = 304,        /* SYM_MINUS_ASSIGN  */
    SYM_MINUS_MINUS = 305,         /* SYM_MINUS_MINUS  */
    SYM_MODULUS = 306,             /* SYM_MODULUS  */
    SYM_MODULUS_ASSIGN = 307,      /* SYM_MODULUS_ASSIGN  */
    SYM_MOLECULE = 308,            /* SYM_MOLECULE  */
    SYM_NEGATE = 309,              /* SYM_NEGATE  */
    SYM_NOT = 310,                 /* SYM_NOT  */
    SYM_NOT_EQUAL = 311,           /* SYM_NOT_EQUAL  */
    SYM_OR = 312,                  /* SYM_OR  */
    SYM_PARM = 313,                /* SYM_PARM  */
    SYM_PERIOD = 314,              /* SYM_PERIOD  */
    SYM_PLUS = 315,                /* SYM_PLUS  */
    SYM_PLUS_ASSIGN = 316,         /* SYM_PLUS_ASSIGN  */
    SYM_PLUS_PLUS = 317,           /* SYM_PLUS_PLUS  */
    SYM_POINT = 318,               /* SYM_POINT  */
    SYM_POINTS_TO = 319,           /* SYM_POINTS_TO  */
    SYM_RBRACE = 320,              /* SYM_RBRACE  */
    SYM_RBRACK = 321,              /* SYM_RBRACK  */
    SYM_RESIDUE = 322,             /* SYM_RESIDUE  */
    SYM_RETURN = 323,              /* SYM_RETURN  */
    SYM_RPAREN = 324,              /* SYM_RPAREN  */
    SYM_SEMICOLON = 325,           /* SYM_SEMICOLON  */
    SYM_SIZE_T = 326,              /* SYM_SIZE_T  */
    SYM_SLASH = 327,               /* SYM_SLASH  */
    SYM_SLASH_ASSIGN = 328,        /* SYM_SLASH_ASSIGN  */
    SYM_STAR = 329,                /* SYM_STAR  */
    SYM_STAR_ASSIGN = 330,         /* SYM_STAR_ASSIGN  */
    SYM_STMTLIST = 331,            /* SYM_STMTLIST  */
    SYM_STRING = 332,              /* SYM_STRING  */
    SYM_STRING_LIT = 333,          /* SYM_STRING_LIT  */
    SYM_STRUCT = 334,              /* SYM_STRUCT  */
    SYM_TEST = 335,                /* SYM_TEST  */
    SYM_TYPE = 336,                /* SYM_TYPE  */
    SYM_UPARROW = 337,             /* SYM_UPARROW  */
    SYM_UPARROW_ASSIGN = 338,      /* SYM_UPARROW_ASSIGN  */
    SYM_WHILE = 339                /* SYM_WHILE  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif

/* Value type.  */


extern YYSTYPE yylval;


int yyparse (void);


#endif /* !YY_YY_NABGRM_TAB_H_INCLUDED  */
