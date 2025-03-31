/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison implementation for Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30802

/* Bison version string.  */
#define YYBISON_VERSION "3.8.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* First part of user prologue.  */
#line 2 "nabgrm.y"

#include <stdio.h>
#include "nab.h"
#include "cgen.h"
#include "errormsg.h"

extern	VALUE_T	val;
static	VALUE_T	v_type;

int yyerror();

# define YYSTYPE_IS_DECLARED 1
# define YYDEBUG 1


#line 87 "nabgrm.tab.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

#include "nabgrm.tab.h"
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_SYM_ADDRESS = 3,                /* SYM_ADDRESS  */
  YYSYMBOL_SYM_ALLOCATE = 4,               /* SYM_ALLOCATE  */
  YYSYMBOL_SYM_AND = 5,                    /* SYM_AND  */
  YYSYMBOL_SYM_ASSERT = 6,                 /* SYM_ASSERT  */
  YYSYMBOL_SYM_ASSIGN = 7,                 /* SYM_ASSIGN  */
  YYSYMBOL_SYM_ATOM = 8,                   /* SYM_ATOM  */
  YYSYMBOL_SYM_ATSIGN = 9,                 /* SYM_ATSIGN  */
  YYSYMBOL_SYM_ATTRIBUTE = 10,             /* SYM_ATTRIBUTE  */
  YYSYMBOL_SYM_BOUNDS = 11,                /* SYM_BOUNDS  */
  YYSYMBOL_SYM_BREAK = 12,                 /* SYM_BREAK  */
  YYSYMBOL_SYM_CALL = 13,                  /* SYM_CALL  */
  YYSYMBOL_SYM_COMMA = 14,                 /* SYM_COMMA  */
  YYSYMBOL_SYM_CONTINUE = 15,              /* SYM_CONTINUE  */
  YYSYMBOL_SYM_DEALLOCATE = 16,            /* SYM_DEALLOCATE  */
  YYSYMBOL_SYM_DEBUG = 17,                 /* SYM_DEBUG  */
  YYSYMBOL_SYM_DECL = 18,                  /* SYM_DECL  */
  YYSYMBOL_SYM_DELETE = 19,                /* SYM_DELETE  */
  YYSYMBOL_SYM_DONT_MATCH = 20,            /* SYM_DONT_MATCH  */
  YYSYMBOL_SYM_DYNAMIC = 21,               /* SYM_DYNAMIC  */
  YYSYMBOL_SYM_ELSE = 22,                  /* SYM_ELSE  */
  YYSYMBOL_SYM_EQUAL = 23,                 /* SYM_EQUAL  */
  YYSYMBOL_SYM_ERROR = 24,                 /* SYM_ERROR  */
  YYSYMBOL_SYM_FILE = 25,                  /* SYM_FILE  */
  YYSYMBOL_SYM_FLOAT = 26,                 /* SYM_FLOAT  */
  YYSYMBOL_SYM_FLOAT_LIT = 27,             /* SYM_FLOAT_LIT  */
  YYSYMBOL_SYM_FOR = 28,                   /* SYM_FOR  */
  YYSYMBOL_SYM_FOREACH = 29,               /* SYM_FOREACH  */
  YYSYMBOL_SYM_GREATER = 30,               /* SYM_GREATER  */
  YYSYMBOL_SYM_GREATER_EQUAL = 31,         /* SYM_GREATER_EQUAL  */
  YYSYMBOL_SYM_HASHED = 32,                /* SYM_HASHED  */
  YYSYMBOL_SYM_IDENT = 33,                 /* SYM_IDENT  */
  YYSYMBOL_SYM_IF = 34,                    /* SYM_IF  */
  YYSYMBOL_SYM_IN = 35,                    /* SYM_IN  */
  YYSYMBOL_SYM_INDEX = 36,                 /* SYM_INDEX  */
  YYSYMBOL_SYM_INDIRECT = 37,              /* SYM_INDIRECT  */
  YYSYMBOL_SYM_INT = 38,                   /* SYM_INT  */
  YYSYMBOL_SYM_INT_LIT = 39,               /* SYM_INT_LIT  */
  YYSYMBOL_SYM_LBRACE = 40,                /* SYM_LBRACE  */
  YYSYMBOL_SYM_LBRACK = 41,                /* SYM_LBRACK  */
  YYSYMBOL_SYM_LESS = 42,                  /* SYM_LESS  */
  YYSYMBOL_SYM_LESS_EQUAL = 43,            /* SYM_LESS_EQUAL  */
  YYSYMBOL_SYM_LIST = 44,                  /* SYM_LIST  */
  YYSYMBOL_SYM_LPAREN = 45,                /* SYM_LPAREN  */
  YYSYMBOL_SYM_MATCH = 46,                 /* SYM_MATCH  */
  YYSYMBOL_SYM_MATRIX = 47,                /* SYM_MATRIX  */
  YYSYMBOL_SYM_MINUS = 48,                 /* SYM_MINUS  */
  YYSYMBOL_SYM_MINUS_ASSIGN = 49,          /* SYM_MINUS_ASSIGN  */
  YYSYMBOL_SYM_MINUS_MINUS = 50,           /* SYM_MINUS_MINUS  */
  YYSYMBOL_SYM_MODULUS = 51,               /* SYM_MODULUS  */
  YYSYMBOL_SYM_MODULUS_ASSIGN = 52,        /* SYM_MODULUS_ASSIGN  */
  YYSYMBOL_SYM_MOLECULE = 53,              /* SYM_MOLECULE  */
  YYSYMBOL_SYM_NEGATE = 54,                /* SYM_NEGATE  */
  YYSYMBOL_SYM_NOT = 55,                   /* SYM_NOT  */
  YYSYMBOL_SYM_NOT_EQUAL = 56,             /* SYM_NOT_EQUAL  */
  YYSYMBOL_SYM_OR = 57,                    /* SYM_OR  */
  YYSYMBOL_SYM_PARM = 58,                  /* SYM_PARM  */
  YYSYMBOL_SYM_PERIOD = 59,                /* SYM_PERIOD  */
  YYSYMBOL_SYM_PLUS = 60,                  /* SYM_PLUS  */
  YYSYMBOL_SYM_PLUS_ASSIGN = 61,           /* SYM_PLUS_ASSIGN  */
  YYSYMBOL_SYM_PLUS_PLUS = 62,             /* SYM_PLUS_PLUS  */
  YYSYMBOL_SYM_POINT = 63,                 /* SYM_POINT  */
  YYSYMBOL_SYM_POINTS_TO = 64,             /* SYM_POINTS_TO  */
  YYSYMBOL_SYM_RBRACE = 65,                /* SYM_RBRACE  */
  YYSYMBOL_SYM_RBRACK = 66,                /* SYM_RBRACK  */
  YYSYMBOL_SYM_RESIDUE = 67,               /* SYM_RESIDUE  */
  YYSYMBOL_SYM_RETURN = 68,                /* SYM_RETURN  */
  YYSYMBOL_SYM_RPAREN = 69,                /* SYM_RPAREN  */
  YYSYMBOL_SYM_SEMICOLON = 70,             /* SYM_SEMICOLON  */
  YYSYMBOL_SYM_SIZE_T = 71,                /* SYM_SIZE_T  */
  YYSYMBOL_SYM_SLASH = 72,                 /* SYM_SLASH  */
  YYSYMBOL_SYM_SLASH_ASSIGN = 73,          /* SYM_SLASH_ASSIGN  */
  YYSYMBOL_SYM_STAR = 74,                  /* SYM_STAR  */
  YYSYMBOL_SYM_STAR_ASSIGN = 75,           /* SYM_STAR_ASSIGN  */
  YYSYMBOL_SYM_STMTLIST = 76,              /* SYM_STMTLIST  */
  YYSYMBOL_SYM_STRING = 77,                /* SYM_STRING  */
  YYSYMBOL_SYM_STRING_LIT = 78,            /* SYM_STRING_LIT  */
  YYSYMBOL_SYM_STRUCT = 79,                /* SYM_STRUCT  */
  YYSYMBOL_SYM_TEST = 80,                  /* SYM_TEST  */
  YYSYMBOL_SYM_TYPE = 81,                  /* SYM_TYPE  */
  YYSYMBOL_SYM_UPARROW = 82,               /* SYM_UPARROW  */
  YYSYMBOL_SYM_UPARROW_ASSIGN = 83,        /* SYM_UPARROW_ASSIGN  */
  YYSYMBOL_SYM_WHILE = 84,                 /* SYM_WHILE  */
  YYSYMBOL_YYACCEPT = 85,                  /* $accept  */
  YYSYMBOL_program = 86,                   /* program  */
  YYSYMBOL_defpart = 87,                   /* defpart  */
  YYSYMBOL_defs = 88,                      /* defs  */
  YYSYMBOL_stmtpart = 89,                  /* stmtpart  */
  YYSYMBOL_stmts = 90,                     /* stmts  */
  YYSYMBOL_def = 91,                       /* def  */
  YYSYMBOL_type_decl = 92,                 /* type_decl  */
  YYSYMBOL_var_decl = 93,                  /* var_decl  */
  YYSYMBOL_type = 94,                      /* type  */
  YYSYMBOL_simple_type = 95,               /* simple_type  */
  YYSYMBOL_struct_type = 96,               /* struct_type  */
  YYSYMBOL_field_list = 97,                /* field_list  */
  YYSYMBOL_field = 98,                     /* field  */
  YYSYMBOL_id_list = 99,                   /* id_list  */
  YYSYMBOL_var_list = 100,                 /* var_list  */
  YYSYMBOL_var = 101,                      /* var  */
  YYSYMBOL_aspec = 102,                    /* aspec  */
  YYSYMBOL_as_list = 103,                  /* as_list  */
  YYSYMBOL_asize = 104,                    /* asize  */
  YYSYMBOL_func_decl = 105,                /* func_decl  */
  YYSYMBOL_func_def = 106,                 /* func_def  */
  YYSYMBOL_107_1 = 107,                    /* $@1  */
  YYSYMBOL_func_hdr = 108,                 /* func_hdr  */
  YYSYMBOL_109_2 = 109,                    /* $@2  */
  YYSYMBOL_formals = 110,                  /* formals  */
  YYSYMBOL_fp_list = 111,                  /* fp_list  */
  YYSYMBOL_f_parm = 112,                   /* f_parm  */
  YYSYMBOL_func_body = 113,                /* func_body  */
  YYSYMBOL_114_3 = 114,                    /* $@3  */
  YYSYMBOL_115_4 = 115,                    /* $@4  */
  YYSYMBOL_116_5 = 116,                    /* $@5  */
  YYSYMBOL_f_defpart = 117,                /* f_defpart  */
  YYSYMBOL_lv_decls = 118,                 /* lv_decls  */
  YYSYMBOL_f_stmtpart = 119,               /* f_stmtpart  */
  YYSYMBOL_stmt = 120,                     /* stmt  */
  YYSYMBOL_alloc_stmt = 121,               /* alloc_stmt  */
  YYSYMBOL_assert_stmt = 122,              /* assert_stmt  */
  YYSYMBOL_break_stmt = 123,               /* break_stmt  */
  YYSYMBOL_cmpd_stmt = 124,                /* cmpd_stmt  */
  YYSYMBOL_125_6 = 125,                    /* $@6  */
  YYSYMBOL_continue_stmt = 126,            /* continue_stmt  */
  YYSYMBOL_dealloc_stmt = 127,             /* dealloc_stmt  */
  YYSYMBOL_debug_stmt = 128,               /* debug_stmt  */
  YYSYMBOL_delete_stmt = 129,              /* delete_stmt  */
  YYSYMBOL_expr_stmt = 130,                /* expr_stmt  */
  YYSYMBOL_if_stmt = 131,                  /* if_stmt  */
  YYSYMBOL_132_7 = 132,                    /* $@7  */
  YYSYMBOL_for_stmt = 133,                 /* for_stmt  */
  YYSYMBOL_return_stmt = 134,              /* return_stmt  */
  YYSYMBOL_135_8 = 135,                    /* $@8  */
  YYSYMBOL_136_9 = 136,                    /* $@9  */
  YYSYMBOL_while_stmt = 137,               /* while_stmt  */
  YYSYMBOL_if_hdr = 138,                   /* if_hdr  */
  YYSYMBOL_139_10 = 139,                   /* $@10  */
  YYSYMBOL_140_11 = 140,                   /* $@11  */
  YYSYMBOL_141_12 = 141,                   /* @12  */
  YYSYMBOL_for_hdr = 142,                  /* for_hdr  */
  YYSYMBOL_143_13 = 143,                   /* $@13  */
  YYSYMBOL_144_14 = 144,                   /* $@14  */
  YYSYMBOL_for_ctrl = 145,                 /* for_ctrl  */
  YYSYMBOL_for_in = 146,                   /* for_in  */
  YYSYMBOL_for_count = 147,                /* for_count  */
  YYSYMBOL_148_15 = 148,                   /* $@15  */
  YYSYMBOL_149_16 = 149,                   /* $@16  */
  YYSYMBOL_for_expr = 150,                 /* for_expr  */
  YYSYMBOL_for_test_expr = 151,            /* for_test_expr  */
  YYSYMBOL_while_hdr = 152,                /* while_hdr  */
  YYSYMBOL_153_17 = 153,                   /* $@17  */
  YYSYMBOL_154_18 = 154,                   /* $@18  */
  YYSYMBOL_155_19 = 155,                   /* @19  */
  YYSYMBOL_dbg_list = 156,                 /* dbg_list  */
  YYSYMBOL_e_list = 157,                   /* e_list  */
  YYSYMBOL_expr = 158,                     /* expr  */
  YYSYMBOL_lval = 159,                     /* lval  */
  YYSYMBOL_ar_lval = 160,                  /* ar_lval  */
  YYSYMBOL_at_lval = 161,                  /* at_lval  */
  YYSYMBOL_rval = 162,                     /* rval  */
  YYSYMBOL_disj = 163,                     /* disj  */
  YYSYMBOL_conj = 164,                     /* conj  */
  YYSYMBOL_a_expr = 165,                   /* a_expr  */
  YYSYMBOL_term = 166,                     /* term  */
  YYSYMBOL_factor = 167,                   /* factor  */
  YYSYMBOL_primary = 168,                  /* primary  */
  YYSYMBOL_incr = 169,                     /* incr  */
  YYSYMBOL_actuals = 170,                  /* actuals  */
  YYSYMBOL_ap_list = 171,                  /* ap_list  */
  YYSYMBOL_a_parm = 172,                   /* a_parm  */
  YYSYMBOL_i_list = 173,                   /* i_list  */
  YYSYMBOL_assignop = 174,                 /* assignop  */
  YYSYMBOL_relop = 175,                    /* relop  */
  YYSYMBOL_addop = 176,                    /* addop  */
  YYSYMBOL_mulop = 177,                    /* mulop  */
  YYSYMBOL_incrop = 178,                   /* incrop  */
  YYSYMBOL_unop = 179,                     /* unop  */
  YYSYMBOL_id = 180,                       /* id  */
  YYSYMBOL_num = 181,                      /* num  */
  YYSYMBOL_string = 182                    /* string  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

/* Work around bug in HP-UX 11.23, which defines these macros
   incorrectly for preprocessor constants.  This workaround can likely
   be removed in 2023, as HPE has promised support for HP-UX 11.23
   (aka HP-UX 11i v2) only through the end of 2022; see Table 2 of
   <https://h20195.www2.hpe.com/V2/getpdf.aspx/4AA4-7673ENW.pdf>.  */
#ifdef __hpux
# undef UINT_LEAST8_MAX
# undef UINT_LEAST16_MAX
# define UINT_LEAST8_MAX 255
# define UINT_LEAST16_MAX 65535
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))


/* Stored state numbers (used for stacks). */
typedef yytype_int16 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if !defined yyoverflow

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* !defined yyoverflow */

#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  27
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   352

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  85
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  98
/* YYNRULES -- Number of rules.  */
#define YYNRULES  188
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  268

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   339


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   181,   181,   182,   183,   184,   185,   186,   187,   188,
     189,   191,   192,   193,   194,   196,   199,   202,   203,   204,
     208,   212,   216,   220,   224,   228,   232,   236,   240,   244,
     248,   252,   256,   257,   259,   261,   262,   264,   265,   266,
     267,   269,   270,   271,   272,   274,   275,   277,   280,   284,
     284,   286,   286,   289,   290,   291,   292,   294,   297,   299,
     301,   297,   304,   305,   306,   307,   308,   309,   311,   312,
     313,   314,   315,   316,   317,   318,   319,   320,   321,   322,
     323,   325,   330,   334,   338,   338,   342,   346,   351,   355,
     360,   363,   365,   364,   366,   367,   370,   367,   373,   375,
     376,   377,   375,   381,   382,   381,   385,   386,   387,   390,
     392,   389,   393,   394,   395,   397,   399,   400,   401,   399,
     405,   407,   408,   409,   412,   413,   415,   416,   417,   418,
     420,   423,   424,   426,   427,   429,   430,   432,   433,   435,
     436,   438,   439,   441,   442,   443,   444,   445,   446,   448,
     450,   451,   452,   453,   454,   455,   457,   458,   459,   461,
     462,   464,   466,   468,   470,   472,   474,   475,   477,   478,
     479,   481,   482,   483,   485,   486,   487,   488,   489,   490,
     491,   492,   493,   495,   496,   498,   499,   500,   501
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "SYM_ADDRESS",
  "SYM_ALLOCATE", "SYM_AND", "SYM_ASSERT", "SYM_ASSIGN", "SYM_ATOM",
  "SYM_ATSIGN", "SYM_ATTRIBUTE", "SYM_BOUNDS", "SYM_BREAK", "SYM_CALL",
  "SYM_COMMA", "SYM_CONTINUE", "SYM_DEALLOCATE", "SYM_DEBUG", "SYM_DECL",
  "SYM_DELETE", "SYM_DONT_MATCH", "SYM_DYNAMIC", "SYM_ELSE", "SYM_EQUAL",
  "SYM_ERROR", "SYM_FILE", "SYM_FLOAT", "SYM_FLOAT_LIT", "SYM_FOR",
  "SYM_FOREACH", "SYM_GREATER", "SYM_GREATER_EQUAL", "SYM_HASHED",
  "SYM_IDENT", "SYM_IF", "SYM_IN", "SYM_INDEX", "SYM_INDIRECT", "SYM_INT",
  "SYM_INT_LIT", "SYM_LBRACE", "SYM_LBRACK", "SYM_LESS", "SYM_LESS_EQUAL",
  "SYM_LIST", "SYM_LPAREN", "SYM_MATCH", "SYM_MATRIX", "SYM_MINUS",
  "SYM_MINUS_ASSIGN", "SYM_MINUS_MINUS", "SYM_MODULUS",
  "SYM_MODULUS_ASSIGN", "SYM_MOLECULE", "SYM_NEGATE", "SYM_NOT",
  "SYM_NOT_EQUAL", "SYM_OR", "SYM_PARM", "SYM_PERIOD", "SYM_PLUS",
  "SYM_PLUS_ASSIGN", "SYM_PLUS_PLUS", "SYM_POINT", "SYM_POINTS_TO",
  "SYM_RBRACE", "SYM_RBRACK", "SYM_RESIDUE", "SYM_RETURN", "SYM_RPAREN",
  "SYM_SEMICOLON", "SYM_SIZE_T", "SYM_SLASH", "SYM_SLASH_ASSIGN",
  "SYM_STAR", "SYM_STAR_ASSIGN", "SYM_STMTLIST", "SYM_STRING",
  "SYM_STRING_LIT", "SYM_STRUCT", "SYM_TEST", "SYM_TYPE", "SYM_UPARROW",
  "SYM_UPARROW_ASSIGN", "SYM_WHILE", "$accept", "program", "defpart",
  "defs", "stmtpart", "stmts", "def", "type_decl", "var_decl", "type",
  "simple_type", "struct_type", "field_list", "field", "id_list",
  "var_list", "var", "aspec", "as_list", "asize", "func_decl", "func_def",
  "$@1", "func_hdr", "$@2", "formals", "fp_list", "f_parm", "func_body",
  "$@3", "$@4", "$@5", "f_defpart", "lv_decls", "f_stmtpart", "stmt",
  "alloc_stmt", "assert_stmt", "break_stmt", "cmpd_stmt", "$@6",
  "continue_stmt", "dealloc_stmt", "debug_stmt", "delete_stmt",
  "expr_stmt", "if_stmt", "$@7", "for_stmt", "return_stmt", "$@8", "$@9",
  "while_stmt", "if_hdr", "$@10", "$@11", "@12", "for_hdr", "$@13", "$@14",
  "for_ctrl", "for_in", "for_count", "$@15", "$@16", "for_expr",
  "for_test_expr", "while_hdr", "$@17", "$@18", "@19", "dbg_list",
  "e_list", "expr", "lval", "ar_lval", "at_lval", "rval", "disj", "conj",
  "a_expr", "term", "factor", "primary", "incr", "actuals", "ap_list",
  "a_parm", "i_list", "assignop", "relop", "addop", "mulop", "incrop",
  "unop", "id", "num", "string", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-184)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-52)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     209,  -184,  -184,  -184,  -184,  -184,  -184,  -184,  -184,  -184,
    -184,  -184,    -2,    37,   182,  -184,   209,  -184,  -184,    -2,
    -184,   -25,  -184,  -184,   -20,  -184,    19,  -184,   250,   250,
      -8,     2,   250,   274,   250,  -184,  -184,  -184,  -184,  -184,
     250,  -184,  -184,  -184,  -184,  -184,  -184,  -184,  -184,   182,
    -184,  -184,  -184,  -184,  -184,  -184,  -184,  -184,  -184,  -184,
    -184,  -184,  -184,  -184,   182,   182,   182,     3,    44,  -184,
    -184,  -184,     7,    55,   295,    -6,    26,     0,  -184,    -2,
     250,    31,  -184,  -184,  -184,    14,    72,    -7,  -184,  -184,
      49,    20,   243,    29,    34,  -184,  -184,    38,   250,    39,
    -184,    78,    40,    62,    66,   182,    45,   250,    68,  -184,
      94,  -184,  -184,  -184,  -184,   250,  -184,  -184,    85,  -184,
    -184,  -184,  -184,   250,  -184,   250,   250,  -184,  -184,  -184,
    -184,  -184,  -184,  -184,  -184,  -184,   250,  -184,  -184,   250,
    -184,  -184,  -184,  -184,   250,   250,   -15,  -184,   -23,  -184,
     250,  -184,    -2,   105,    79,  -184,  -184,  -184,    -2,    58,
     243,  -184,  -184,  -184,    59,    -4,  -184,   250,  -184,  -184,
    -184,    13,  -184,  -184,  -184,  -184,   116,    65,  -184,  -184,
    -184,  -184,  -184,  -184,  -184,  -184,  -184,    64,  -184,   121,
    -184,    95,  -184,  -184,    75,  -184,   131,  -184,   209,   209,
      77,   135,  -184,  -184,  -184,  -184,   250,   250,  -184,    81,
     250,   182,   250,  -184,  -184,   250,  -184,   113,    -2,  -184,
      87,  -184,   140,  -184,    -2,  -184,   209,  -184,    -2,    88,
    -184,  -184,    89,  -184,    22,  -184,  -184,  -184,  -184,  -184,
    -184,  -184,  -184,  -184,   209,   182,  -184,  -184,  -184,  -184,
      -2,    93,    96,  -184,   182,    99,   250,  -184,  -184,  -184,
    -184,   100,  -184,   101,  -184,  -184,   250,  -184
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       4,    27,    25,    23,    21,    19,    26,    29,    24,    28,
      20,    22,     0,     0,     8,     3,     5,    11,    12,     0,
      17,    18,    13,    14,    49,   185,    30,     1,     0,     0,
       0,     0,     0,     0,     0,   187,   103,    99,   186,    84,
       0,   183,   182,   184,   181,    95,   188,   116,     2,     7,
       9,    69,    70,    71,    73,    72,    74,    75,    76,    68,
      78,    77,    79,    80,     0,     0,     0,     0,   143,   127,
     128,   124,   131,   133,   135,   137,   139,   141,   146,     0,
       0,   126,   144,   145,     6,     0,    37,    39,    15,    47,
       0,     0,     0,     0,     0,    83,    86,     0,     0,     0,
     121,   122,     0,     0,     0,     0,     0,     0,     0,    10,
      91,    94,    98,    90,   159,     0,   161,   164,     0,   160,
     163,   162,   165,     0,   151,     0,     0,   173,   168,   171,
     170,   174,   166,   167,   172,   169,     0,   176,   175,     0,
     180,   179,   178,   177,     0,     0,   150,   126,   143,   147,
     153,    16,     0,     0,     0,    58,    50,    48,     0,     0,
      32,    81,    82,    87,     0,     0,    88,     0,    89,   104,
     100,     0,   149,    96,   117,    92,   157,     0,   130,   125,
     132,   134,   136,   138,   140,   142,   156,     0,   152,   154,
      38,    39,    46,    42,     0,    41,    43,    45,    54,    63,
       0,    35,    31,    33,   120,   123,   113,     0,    85,     0,
       0,     0,     0,   129,   148,     0,    40,     0,     0,    18,
       0,    53,    55,    64,     0,    59,    62,    34,     0,     0,
     106,   107,     0,   112,   126,   101,    97,   118,    93,   158,
     155,    44,    57,    52,     0,    67,    65,    36,   105,   109,
       0,     0,     0,    56,    66,     0,   115,   108,   102,   119,
      60,     0,   114,     0,   110,    61,   113,   111
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -184,  -184,  -184,   150,  -184,  -102,  -184,  -184,  -183,  -178,
     -77,    33,    12,  -184,   -55,    24,   -44,  -184,   -40,  -184,
    -184,  -184,  -184,  -184,  -184,  -184,   -64,  -184,  -184,  -184,
    -184,  -184,  -184,  -184,  -184,   -42,  -184,  -184,  -184,  -184,
    -184,  -184,  -184,  -184,  -184,  -184,  -184,  -184,  -184,  -184,
    -184,  -184,  -184,  -184,  -184,  -184,  -184,  -184,  -184,  -184,
    -184,  -184,  -184,  -184,  -184,   -85,  -184,  -184,  -184,  -184,
    -184,  -184,   -87,   -28,   -24,  -184,  -184,    60,    67,  -184,
     -65,    46,    47,   115,  -184,  -184,   -13,  -184,   -12,  -184,
    -184,  -184,  -184,   -60,  -184,   -10,  -184,  -184
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
       0,    13,    14,    15,    48,    49,    16,    17,    18,    19,
      20,   219,   159,   160,   200,    85,    86,   194,   195,   196,
      22,    23,    90,    24,   154,   220,   221,   222,   156,   199,
     245,   263,   225,   226,   255,    50,    51,    52,    53,    54,
     105,    55,    56,    57,    58,    59,    60,   211,    61,    62,
     107,   209,    63,    64,   104,   207,   251,    65,   103,   206,
     229,   230,   231,   256,   266,   232,   261,    66,   108,   210,
     252,    99,   100,    67,    68,    69,    70,    71,    72,    73,
      74,    75,    76,    77,    78,   187,   188,   189,   177,   123,
     136,   139,   144,    79,    80,    81,    82,    83
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      93,    94,    26,   171,    97,   101,   102,   109,   124,    87,
     167,   164,   106,    25,    91,   158,   223,    28,   115,    29,
     218,   224,   110,   111,   112,    30,   115,    42,    31,    32,
      33,    25,    34,    21,   153,   140,   118,    27,   -51,    44,
      35,    36,   137,   246,   118,    88,    25,    37,   224,    21,
      89,   114,    38,    39,   138,   146,   148,   250,    40,    92,
     126,    41,    95,    42,   125,   172,   218,   150,    43,   147,
     165,   182,    96,   113,   183,    44,   150,   141,   208,   173,
     205,    45,   145,   158,   151,   115,   152,   176,   124,   155,
     157,    46,   167,   116,    42,   179,   117,    47,   142,   161,
     143,   148,   148,   118,   162,   119,    44,   169,   163,   166,
     168,   170,   148,   174,   172,   148,   175,   120,   178,   121,
     148,   148,   186,   202,   198,   197,   192,   122,   204,   109,
     212,   213,    35,   214,   192,   215,   153,   193,    25,   101,
      35,   216,   191,   254,    38,   217,    25,   227,   201,   228,
      40,   236,    38,    41,   244,    42,   243,   248,    40,   249,
      43,    41,   258,    42,   260,   259,    84,    44,    43,   238,
     264,   265,   203,   247,   242,    44,   190,   241,   233,   235,
     253,   267,   237,    46,   176,   180,    28,   186,    29,   197,
     184,    46,   185,   181,    30,   149,   234,    31,    32,    33,
     239,    34,   240,     0,     0,     0,     0,     0,   191,    35,
      36,     0,   109,     0,   191,    25,    37,     1,   201,     0,
       2,    38,    39,     0,     0,     0,     0,    40,   262,     0,
      41,     0,    42,     0,     3,     4,     0,    43,   233,     0,
     257,     0,     0,     0,    44,     0,     0,     5,     0,     0,
      45,     1,     0,     0,     2,     0,     6,     0,     0,     0,
      46,     0,     7,     0,     0,     0,    47,     0,     3,     4,
       0,     0,     8,     0,     0,     0,     9,    35,     0,     0,
      10,     5,     0,    25,     0,     0,    11,     0,    12,    38,
       6,     0,     0,     0,     0,    40,     7,     0,    41,     0,
      42,    35,     0,     0,     0,    43,     8,    25,     0,     0,
       9,     0,    44,    38,    10,   127,     0,     0,   128,    98,
      11,     0,    41,     0,    42,   129,   130,     0,    46,    43,
     131,     0,     0,     0,     0,     0,    44,   132,   133,     0,
       0,   134,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   135,    46
};

static const yytype_int16 yycheck[] =
{
      28,    29,    12,   105,    32,    33,    34,    49,    68,    19,
      14,    98,    40,    33,    24,    92,   199,     4,    41,     6,
     198,   199,    64,    65,    66,    12,    41,    50,    15,    16,
      17,    33,    19,     0,    41,     9,    59,     0,    45,    62,
      27,    28,    48,   226,    59,    70,    33,    34,   226,    16,
      70,     7,    39,    40,    60,    79,    80,    35,    45,    40,
       5,    48,    70,    50,    57,    69,   244,    45,    55,    79,
      98,   136,    70,    70,   139,    62,    45,    51,    65,   107,
     167,    68,    82,   160,    70,    41,    14,   115,   148,    40,
      70,    78,    14,    49,    50,   123,    52,    84,    72,    70,
      74,   125,   126,    59,    70,    61,    62,    45,    70,    70,
      70,    45,   136,    45,    69,   139,    22,    73,    33,    75,
     144,   145,   150,    65,    45,   153,    21,    83,    69,   171,
      14,    66,    27,    69,    21,    14,    41,    32,    33,   167,
      27,    66,   152,   245,    39,    14,    33,    70,   158,    14,
      45,    70,    39,    48,    14,    50,    69,    69,    45,    70,
      55,    48,    69,    50,    65,    69,    16,    62,    55,   211,
      70,    70,   160,   228,   218,    62,   152,   217,   206,   207,
     244,   266,   210,    78,   212,   125,     4,   215,     6,   217,
     144,    78,   145,   126,    12,    80,   206,    15,    16,    17,
     212,    19,   215,    -1,    -1,    -1,    -1,    -1,   218,    27,
      28,    -1,   254,    -1,   224,    33,    34,     8,   228,    -1,
      11,    39,    40,    -1,    -1,    -1,    -1,    45,   256,    -1,
      48,    -1,    50,    -1,    25,    26,    -1,    55,   266,    -1,
     250,    -1,    -1,    -1,    62,    -1,    -1,    38,    -1,    -1,
      68,     8,    -1,    -1,    11,    -1,    47,    -1,    -1,    -1,
      78,    -1,    53,    -1,    -1,    -1,    84,    -1,    25,    26,
      -1,    -1,    63,    -1,    -1,    -1,    67,    27,    -1,    -1,
      71,    38,    -1,    33,    -1,    -1,    77,    -1,    79,    39,
      47,    -1,    -1,    -1,    -1,    45,    53,    -1,    48,    -1,
      50,    27,    -1,    -1,    -1,    55,    63,    33,    -1,    -1,
      67,    -1,    62,    39,    71,    20,    -1,    -1,    23,    45,
      77,    -1,    48,    -1,    50,    30,    31,    -1,    78,    55,
      35,    -1,    -1,    -1,    -1,    -1,    62,    42,    43,    -1,
      -1,    46,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    56,    78
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     8,    11,    25,    26,    38,    47,    53,    63,    67,
      71,    77,    79,    86,    87,    88,    91,    92,    93,    94,
      95,    96,   105,   106,   108,    33,   180,     0,     4,     6,
      12,    15,    16,    17,    19,    27,    28,    34,    39,    40,
      45,    48,    50,    55,    62,    68,    78,    84,    89,    90,
     120,   121,   122,   123,   124,   126,   127,   128,   129,   130,
     131,   133,   134,   137,   138,   142,   152,   158,   159,   160,
     161,   162,   163,   164,   165,   166,   167,   168,   169,   178,
     179,   180,   181,   182,    88,   100,   101,   180,    70,    70,
     107,   180,    40,   158,   158,    70,    70,   158,    45,   156,
     157,   158,   158,   143,   139,   125,   158,   135,   153,   120,
     120,   120,   120,    70,     7,    41,    49,    52,    59,    61,
      73,    75,    83,   174,   178,    57,     5,    20,    23,    30,
      31,    35,    42,    43,    46,    56,   175,    48,    60,   176,
       9,    51,    72,    74,   177,    82,   159,   180,   159,   168,
      45,    70,    14,    41,   109,    40,   113,    70,    95,    97,
      98,    70,    70,    70,   157,   158,    70,    14,    70,    45,
      45,    90,    69,   158,    45,    22,   158,   173,    33,   158,
     162,   163,   165,   165,   166,   167,   158,   170,   171,   172,
     100,   180,    21,    32,   102,   103,   104,   158,    45,   114,
      99,   180,    65,    97,    69,   157,   144,   140,    65,   136,
     154,   132,    14,    66,    69,    14,    66,    14,    94,    96,
     110,   111,   112,    93,    94,   117,   118,    70,    14,   145,
     146,   147,   150,   158,   180,   158,    70,   158,   120,   173,
     171,   103,   101,    69,    14,   115,    93,    99,    69,    70,
      35,   141,   155,   111,    90,   119,   148,   180,    69,    69,
      65,   151,   158,   116,    70,    70,   149,   150
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_uint8 yyr1[] =
{
       0,    85,    86,    87,    87,    88,    88,    89,    89,    90,
      90,    91,    91,    91,    91,    92,    93,    94,    94,    95,
      95,    95,    95,    95,    95,    95,    95,    95,    95,    95,
      96,    96,    97,    97,    98,    99,    99,   100,   100,   101,
     101,   102,   102,   103,   103,   104,   104,   105,   105,   107,
     106,   109,   108,   110,   110,   111,   111,   112,   114,   115,
     116,   113,   117,   117,   118,   118,   119,   119,   120,   120,
     120,   120,   120,   120,   120,   120,   120,   120,   120,   120,
     120,   121,   122,   123,   125,   124,   126,   127,   128,   129,
     130,   131,   132,   131,   133,   135,   136,   134,   137,   139,
     140,   141,   138,   143,   144,   142,   145,   145,   146,   148,
     149,   147,   150,   150,   151,   151,   153,   154,   155,   152,
     156,   156,   157,   157,   158,   158,   159,   159,   159,   160,
     161,   162,   162,   163,   163,   164,   164,   165,   165,   166,
     166,   167,   167,   168,   168,   168,   168,   168,   168,   168,
     169,   169,   170,   170,   171,   171,   172,   173,   173,   174,
     174,   174,   174,   174,   174,   174,   175,   175,   175,   175,
     175,   175,   175,   175,   175,   176,   176,   177,   177,   177,
     177,   178,   178,   179,   179,   180,   181,   181,   182
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     2,     1,     0,     1,     2,     1,     0,     1,
       2,     1,     1,     1,     1,     2,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       2,     5,     1,     2,     3,     1,     3,     1,     3,     1,
       4,     1,     1,     1,     3,     1,     1,     2,     3,     0,
       3,     0,     6,     1,     0,     1,     3,     2,     0,     0,
       0,     8,     1,     0,     1,     2,     1,     0,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     3,     3,     2,     0,     4,     2,     3,     3,     3,
       2,     2,     0,     5,     2,     0,     0,     5,     2,     0,
       0,     0,     7,     0,     0,     6,     1,     1,     3,     0,
       0,     7,     1,     0,     1,     0,     0,     0,     0,     7,
       3,     1,     1,     3,     1,     3,     1,     1,     1,     4,
       3,     1,     3,     1,     3,     1,     3,     1,     3,     1,
       3,     1,     3,     1,     1,     1,     1,     2,     4,     3,
       2,     2,     1,     0,     1,     3,     1,     1,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYNOMEM         goto yyexhaustedlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use YYerror or YYUNDEF. */
#define YYERRCODE YYUNDEF


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)




# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Kind, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  if (!yyvaluep)
    return;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp,
                 int yyrule)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)]);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep)
{
  YY_USE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/* Lookahead token kind.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;




/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yychar = YYEMPTY; /* Cause a token to be read.  */

  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    YYNOMEM;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        YYNOMEM;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          YYNOMEM;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */


  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token\n"));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else if (yychar == YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      goto yyerrlab1;
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 2: /* program: defpart stmtpart  */
#line 181 "nabgrm.y"
                                   { CG_genend(); }
#line 1502 "nabgrm.tab.c"
    break;

  case 4: /* defpart: %empty  */
#line 183 "nabgrm.y"
                                { (yyval.npval) = NULL; }
#line 1508 "nabgrm.tab.c"
    break;

  case 8: /* stmtpart: %empty  */
#line 187 "nabgrm.y"
                                { (yyval.npval) = NULL; }
#line 1514 "nabgrm.tab.c"
    break;

  case 15: /* type_decl: struct_type SYM_SEMICOLON  */
#line 197 "nabgrm.y"
                                { (yyval.npval) = node( SYM_DECL, 0, (yyvsp[-1].npval), 0 );
				  CG_genvardecl( (yyval.npval), 0, 0, 0 ); }
#line 1521 "nabgrm.tab.c"
    break;

  case 16: /* var_decl: type var_list SYM_SEMICOLON  */
#line 200 "nabgrm.y"
                                { (yyval.npval) = node( SYM_DECL, 0, (yyvsp[-2].npval), (yyvsp[-1].npval) );
				  CG_genvardecl( (yyval.npval), 0, 0, 0 ); }
#line 1528 "nabgrm.tab.c"
    break;

  case 17: /* type: simple_type  */
#line 202 "nabgrm.y"
                        { (yyval.npval) = (yyvsp[0].npval); }
#line 1534 "nabgrm.tab.c"
    break;

  case 18: /* type: struct_type  */
#line 203 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 1540 "nabgrm.tab.c"
    break;

  case 19: /* simple_type: SYM_INT  */
#line 204 "nabgrm.y"
                                        { 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_INT;
			(yyval.npval) = node( SYM_TYPE, &v_type, 0, 0 ); }
#line 1549 "nabgrm.tab.c"
    break;

  case 20: /* simple_type: SYM_SIZE_T  */
#line 208 "nabgrm.y"
                                {
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_SIZE_T;
			(yyval.npval) = node( SYM_TYPE, &v_type, 0, 0 ); }
#line 1558 "nabgrm.tab.c"
    break;

  case 21: /* simple_type: SYM_FLOAT  */
#line 212 "nabgrm.y"
                                { 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_FLOAT;
			(yyval.npval) = node( SYM_TYPE, &v_type, 0, 0 ); }
#line 1567 "nabgrm.tab.c"
    break;

  case 22: /* simple_type: SYM_STRING  */
#line 216 "nabgrm.y"
                                { 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_STRING;
			(yyval.npval) = node( SYM_TYPE, &v_type, 0, 0 ); }
#line 1576 "nabgrm.tab.c"
    break;

  case 23: /* simple_type: SYM_FILE  */
#line 220 "nabgrm.y"
                                { 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_FILE;
			(yyval.npval) = node( SYM_TYPE, &v_type, 0, 0 ); }
#line 1585 "nabgrm.tab.c"
    break;

  case 24: /* simple_type: SYM_POINT  */
#line 224 "nabgrm.y"
                                { 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_POINT;
			(yyval.npval) = node( SYM_TYPE, &v_type, 0, 0 ); }
#line 1594 "nabgrm.tab.c"
    break;

  case 25: /* simple_type: SYM_BOUNDS  */
#line 228 "nabgrm.y"
                                { 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_BOUNDS;
			(yyval.npval) = node( SYM_TYPE, &v_type, 0, 0 ); }
#line 1603 "nabgrm.tab.c"
    break;

  case 26: /* simple_type: SYM_MATRIX  */
#line 232 "nabgrm.y"
                                { 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_MATRIX;
			(yyval.npval) = node( SYM_TYPE, &v_type, 0, 0 ); }
#line 1612 "nabgrm.tab.c"
    break;

  case 27: /* simple_type: SYM_ATOM  */
#line 236 "nabgrm.y"
                                { 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_ATOM;
			(yyval.npval) = node( SYM_TYPE, &v_type, 0, 0 ); }
#line 1621 "nabgrm.tab.c"
    break;

  case 28: /* simple_type: SYM_RESIDUE  */
#line 240 "nabgrm.y"
                                { 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_RESIDUE;
			(yyval.npval) = node( SYM_TYPE, &v_type, 0, 0 ); }
#line 1630 "nabgrm.tab.c"
    break;

  case 29: /* simple_type: SYM_MOLECULE  */
#line 244 "nabgrm.y"
                                { 
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_MOLECULE;
			(yyval.npval) = node( SYM_TYPE, &v_type, 0, 0 ); }
#line 1639 "nabgrm.tab.c"
    break;

  case 30: /* struct_type: SYM_STRUCT id  */
#line 248 "nabgrm.y"
                                {
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_USER;
			(yyval.npval) = node( SYM_TYPE, &v_type, node( SYM_STRUCT, 0, (yyvsp[0].npval), 0 ), 0 ); }
#line 1648 "nabgrm.tab.c"
    break;

  case 31: /* struct_type: SYM_STRUCT id SYM_LBRACE field_list SYM_RBRACE  */
#line 252 "nabgrm.y"
                                                                 {
			v_type.v_type = T_INT;
			v_type.v_value.v_ival = T_USER;
			(yyval.npval) = node( SYM_TYPE, &v_type, node( SYM_STRUCT, 0, (yyvsp[-3].npval), (yyvsp[-1].npval) ), 0 ); }
#line 1657 "nabgrm.tab.c"
    break;

  case 32: /* field_list: field  */
#line 256 "nabgrm.y"
                        { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[0].npval), 0 ); }
#line 1663 "nabgrm.tab.c"
    break;

  case 33: /* field_list: field field_list  */
#line 258 "nabgrm.y"
                                        { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[-1].npval), (yyvsp[0].npval) ); }
#line 1669 "nabgrm.tab.c"
    break;

  case 34: /* field: simple_type id_list SYM_SEMICOLON  */
#line 260 "nabgrm.y"
                                        { (yyval.npval) = node( SYM_DECL, 0, (yyvsp[-2].npval), (yyvsp[-1].npval) ); }
#line 1675 "nabgrm.tab.c"
    break;

  case 35: /* id_list: id  */
#line 261 "nabgrm.y"
                        { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[0].npval), 0 ); }
#line 1681 "nabgrm.tab.c"
    break;

  case 36: /* id_list: id SYM_COMMA id_list  */
#line 263 "nabgrm.y"
                                        { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 1687 "nabgrm.tab.c"
    break;

  case 37: /* var_list: var  */
#line 264 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[0].npval), 0 ); }
#line 1693 "nabgrm.tab.c"
    break;

  case 38: /* var_list: var SYM_COMMA var_list  */
#line 265 "nabgrm.y"
                                         { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 1699 "nabgrm.tab.c"
    break;

  case 39: /* var: id  */
#line 266 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 1705 "nabgrm.tab.c"
    break;

  case 40: /* var: id SYM_LBRACK aspec SYM_RBRACK  */
#line 268 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LBRACK, 0, (yyvsp[-3].npval), (yyvsp[-1].npval) ); }
#line 1711 "nabgrm.tab.c"
    break;

  case 41: /* aspec: as_list  */
#line 269 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 1717 "nabgrm.tab.c"
    break;

  case 42: /* aspec: SYM_HASHED  */
#line 270 "nabgrm.y"
                                { (yyval.npval) = node( SYM_HASHED, 0, 0, 0 ); }
#line 1723 "nabgrm.tab.c"
    break;

  case 43: /* as_list: asize  */
#line 271 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[0].npval), 0 ); }
#line 1729 "nabgrm.tab.c"
    break;

  case 44: /* as_list: asize SYM_COMMA as_list  */
#line 273 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 1735 "nabgrm.tab.c"
    break;

  case 45: /* asize: expr  */
#line 274 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 1741 "nabgrm.tab.c"
    break;

  case 46: /* asize: SYM_DYNAMIC  */
#line 275 "nabgrm.y"
                                { (yyval.npval) = node( SYM_DYNAMIC, 0, 0, 0 ); }
#line 1747 "nabgrm.tab.c"
    break;

  case 47: /* func_decl: func_hdr SYM_SEMICOLON  */
#line 278 "nabgrm.y"
                                { CG_genop( NULL, SYM_SEMICOLON);
				  CG_genfend( NULL ); }
#line 1754 "nabgrm.tab.c"
    break;

  case 48: /* func_decl: func_hdr id SYM_SEMICOLON  */
#line 281 "nabgrm.y"
                                { CG_genop( NULL, SYM_SEMICOLON );
				  CG_genfend( (yyvsp[-1].npval) ); }
#line 1761 "nabgrm.tab.c"
    break;

  case 49: /* $@1: %empty  */
#line 284 "nabgrm.y"
                                { CG_genfstart(); }
#line 1767 "nabgrm.tab.c"
    break;

  case 50: /* func_def: func_hdr $@1 func_body  */
#line 285 "nabgrm.y"
                                { CG_genfend( NULL ); }
#line 1773 "nabgrm.tab.c"
    break;

  case 51: /* $@2: %empty  */
#line 286 "nabgrm.y"
                                { CG_genfhdr( (yyvsp[-1].npval), (yyvsp[0].npval) ); }
#line 1779 "nabgrm.tab.c"
    break;

  case 52: /* func_hdr: type id $@2 SYM_LPAREN formals SYM_RPAREN  */
#line 288 "nabgrm.y"
                                { CG_genplist( (yyvsp[-1].npval) ); }
#line 1785 "nabgrm.tab.c"
    break;

  case 53: /* formals: fp_list  */
#line 289 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 1791 "nabgrm.tab.c"
    break;

  case 54: /* formals: %empty  */
#line 290 "nabgrm.y"
                                { (yyval.npval) = NULL; }
#line 1797 "nabgrm.tab.c"
    break;

  case 55: /* fp_list: f_parm  */
#line 291 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[0].npval), 0 ); }
#line 1803 "nabgrm.tab.c"
    break;

  case 56: /* fp_list: f_parm SYM_COMMA fp_list  */
#line 293 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 1809 "nabgrm.tab.c"
    break;

  case 57: /* f_parm: type var  */
#line 295 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[0].npval), NULL ); 
			  	  (yyval.npval) = node( SYM_DECL, 0, (yyvsp[-1].npval), (yyval.npval) ); }
#line 1816 "nabgrm.tab.c"
    break;

  case 58: /* $@3: %empty  */
#line 297 "nabgrm.y"
                                { CG_genpdecls();
				  CG_genop( NULL, SYM_LBRACE ); }
#line 1823 "nabgrm.tab.c"
    break;

  case 59: /* $@4: %empty  */
#line 299 "nabgrm.y"
                                { CG_genedefs( TRUE ); }
#line 1829 "nabgrm.tab.c"
    break;

  case 60: /* $@5: %empty  */
#line 301 "nabgrm.y"
                                { CG_genestmts( TRUE );
				  CG_genop( NULL, SYM_RBRACE ); }
#line 1836 "nabgrm.tab.c"
    break;

  case 61: /* func_body: SYM_LBRACE $@3 f_defpart $@4 f_stmtpart SYM_RBRACE $@5 SYM_SEMICOLON  */
#line 303 "nabgrm.y"
                                { (yyval.npval)=NULL; }
#line 1842 "nabgrm.tab.c"
    break;

  case 63: /* f_defpart: %empty  */
#line 305 "nabgrm.y"
                                { (yyval.npval) = NULL; }
#line 1848 "nabgrm.tab.c"
    break;

  case 67: /* f_stmtpart: %empty  */
#line 309 "nabgrm.y"
                                { (yyval.npval) = NULL; }
#line 1854 "nabgrm.tab.c"
    break;

  case 81: /* alloc_stmt: SYM_ALLOCATE expr SYM_SEMICOLON  */
#line 326 "nabgrm.y"
                                { CG_genmain();
				  (yyval.npval) = node( SYM_ALLOCATE, 0, 0, (yyvsp[-1].npval) );
				  CG_genexpr( (yyval.npval) );
				  CG_genop( NULL, SYM_SEMICOLON ); }
#line 1863 "nabgrm.tab.c"
    break;

  case 82: /* assert_stmt: SYM_ASSERT expr SYM_SEMICOLON  */
#line 331 "nabgrm.y"
                                { CG_genmain();
				  (yyval.npval) = node( SYM_ASSERT, 0, 0, (yyvsp[-1].npval) );
				  CG_genassert( (yyval.npval) ); }
#line 1871 "nabgrm.tab.c"
    break;

  case 83: /* break_stmt: SYM_BREAK SYM_SEMICOLON  */
#line 335 "nabgrm.y"
                                { CG_genmain();
				  CG_genrword( SYM_BREAK );
				  CG_genop( NULL, SYM_SEMICOLON ); }
#line 1879 "nabgrm.tab.c"
    break;

  case 84: /* $@6: %empty  */
#line 338 "nabgrm.y"
                                { CG_genmain();
				  CG_genop( NULL, SYM_LBRACE ); }
#line 1886 "nabgrm.tab.c"
    break;

  case 85: /* cmpd_stmt: SYM_LBRACE $@6 stmts SYM_RBRACE  */
#line 341 "nabgrm.y"
                                { CG_genop( NULL, SYM_RBRACE ); }
#line 1892 "nabgrm.tab.c"
    break;

  case 86: /* continue_stmt: SYM_CONTINUE SYM_SEMICOLON  */
#line 343 "nabgrm.y"
                                { CG_genmain();
				  CG_genrword( SYM_CONTINUE );
				  CG_genop( NULL, SYM_SEMICOLON ); }
#line 1900 "nabgrm.tab.c"
    break;

  case 87: /* dealloc_stmt: SYM_DEALLOCATE expr SYM_SEMICOLON  */
#line 347 "nabgrm.y"
                                { CG_genmain();
				  (yyval.npval) = node( SYM_DEALLOCATE, 0, 0, (yyvsp[-1].npval) );
			 	  CG_genexpr((yyval.npval));
				  CG_genop( NULL, SYM_SEMICOLON ); }
#line 1909 "nabgrm.tab.c"
    break;

  case 88: /* debug_stmt: SYM_DEBUG dbg_list SYM_SEMICOLON  */
#line 352 "nabgrm.y"
                                { CG_genmain();
				  (yyval.npval) = node( SYM_DEBUG, 0, 0, (yyvsp[-1].npval) );
				  CG_gendebug( (yyval.npval) ); }
#line 1917 "nabgrm.tab.c"
    break;

  case 89: /* delete_stmt: SYM_DELETE expr SYM_SEMICOLON  */
#line 356 "nabgrm.y"
                                { CG_genmain();
				  (yyval.npval) = node( SYM_DELETE, 0, 0, (yyvsp[-1].npval) );
			 	  CG_genexpr((yyval.npval));
				  CG_genop( NULL, SYM_SEMICOLON ); }
#line 1926 "nabgrm.tab.c"
    break;

  case 90: /* expr_stmt: expr SYM_SEMICOLON  */
#line 361 "nabgrm.y"
                                { CG_genmain(); CG_genexpr( (yyvsp[-1].npval) );
				  CG_genop( NULL, SYM_SEMICOLON ); }
#line 1933 "nabgrm.tab.c"
    break;

  case 92: /* $@7: %empty  */
#line 365 "nabgrm.y"
                                { CG_genrword( SYM_ELSE ); }
#line 1939 "nabgrm.tab.c"
    break;

  case 95: /* $@8: %empty  */
#line 367 "nabgrm.y"
                                { CG_genmain();
				  CG_genrword(SYM_RETURN);
				  CG_genop( NULL, SYM_LPAREN ); }
#line 1947 "nabgrm.tab.c"
    break;

  case 96: /* $@9: %empty  */
#line 370 "nabgrm.y"
                                { CG_genexpr( (yyvsp[0].npval) ); }
#line 1953 "nabgrm.tab.c"
    break;

  case 97: /* return_stmt: SYM_RETURN $@8 expr $@9 SYM_SEMICOLON  */
#line 371 "nabgrm.y"
                                { CG_genop( NULL, SYM_RPAREN );
				  CG_genop( NULL, SYM_SEMICOLON ); }
#line 1960 "nabgrm.tab.c"
    break;

  case 99: /* $@10: %empty  */
#line 375 "nabgrm.y"
                                { CG_genmain(); CG_genrword( SYM_IF ); }
#line 1966 "nabgrm.tab.c"
    break;

  case 100: /* $@11: %empty  */
#line 376 "nabgrm.y"
                                { CG_genop( NULL, SYM_LPAREN ); }
#line 1972 "nabgrm.tab.c"
    break;

  case 101: /* @12: %empty  */
#line 377 "nabgrm.y"
                                { (yyval.npval) = node( SYM_TEST, 0, 0, (yyvsp[0].npval) );
				  CG_genexpr( (yyval.npval) ); }
#line 1979 "nabgrm.tab.c"
    break;

  case 102: /* if_hdr: SYM_IF $@10 SYM_LPAREN $@11 expr @12 SYM_RPAREN  */
#line 379 "nabgrm.y"
                                { CG_genop( NULL, SYM_RPAREN ); }
#line 1985 "nabgrm.tab.c"
    break;

  case 103: /* $@13: %empty  */
#line 381 "nabgrm.y"
                                { CG_genmain(); CG_genrword( SYM_FOR ); }
#line 1991 "nabgrm.tab.c"
    break;

  case 104: /* $@14: %empty  */
#line 382 "nabgrm.y"
                                { CG_genop( NULL, SYM_LPAREN ); }
#line 1997 "nabgrm.tab.c"
    break;

  case 105: /* for_hdr: SYM_FOR $@13 SYM_LPAREN $@14 for_ctrl SYM_RPAREN  */
#line 384 "nabgrm.y"
                                { CG_genop( NULL, SYM_RPAREN ); }
#line 2003 "nabgrm.tab.c"
    break;

  case 108: /* for_in: id SYM_IN id  */
#line 387 "nabgrm.y"
                                { (yyval.npval) = node( SYM_FOREACH, 0, (yyvsp[-2].npval), (yyvsp[0].npval) );
				  CG_genexpr( (yyval.npval) ); }
#line 2010 "nabgrm.tab.c"
    break;

  case 109: /* $@15: %empty  */
#line 390 "nabgrm.y"
                                { CG_genop( NULL, SYM_SEMICOLON ); }
#line 2016 "nabgrm.tab.c"
    break;

  case 110: /* $@16: %empty  */
#line 392 "nabgrm.y"
                                { CG_genop( NULL, SYM_SEMICOLON ); }
#line 2022 "nabgrm.tab.c"
    break;

  case 112: /* for_expr: expr  */
#line 393 "nabgrm.y"
                                { CG_genexpr( (yyvsp[0].npval) ); }
#line 2028 "nabgrm.tab.c"
    break;

  case 113: /* for_expr: %empty  */
#line 394 "nabgrm.y"
                                { (yyval.npval) = NULL; }
#line 2034 "nabgrm.tab.c"
    break;

  case 114: /* for_test_expr: expr  */
#line 395 "nabgrm.y"
                                { (yyval.npval) = node( SYM_TEST, 0, 0, (yyvsp[0].npval) );
				  CG_genexpr( (yyval.npval) ); }
#line 2041 "nabgrm.tab.c"
    break;

  case 115: /* for_test_expr: %empty  */
#line 397 "nabgrm.y"
                                { (yyval.npval) = NULL; }
#line 2047 "nabgrm.tab.c"
    break;

  case 116: /* $@17: %empty  */
#line 399 "nabgrm.y"
                                { CG_genmain(); CG_genrword( SYM_WHILE ); }
#line 2053 "nabgrm.tab.c"
    break;

  case 117: /* $@18: %empty  */
#line 400 "nabgrm.y"
                                { CG_genop( NULL, SYM_LPAREN ); }
#line 2059 "nabgrm.tab.c"
    break;

  case 118: /* @19: %empty  */
#line 401 "nabgrm.y"
                                { (yyval.npval) = node( SYM_TEST, 0, 0, (yyvsp[0].npval) );
				  CG_genexpr( (yyval.npval) ); }
#line 2066 "nabgrm.tab.c"
    break;

  case 119: /* while_hdr: SYM_WHILE $@17 SYM_LPAREN $@18 expr @19 SYM_RPAREN  */
#line 403 "nabgrm.y"
                                { CG_genop( NULL, SYM_RPAREN ); }
#line 2072 "nabgrm.tab.c"
    break;

  case 120: /* dbg_list: SYM_LPAREN e_list SYM_RPAREN  */
#line 406 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[-1].npval); }
#line 2078 "nabgrm.tab.c"
    break;

  case 121: /* dbg_list: e_list  */
#line 407 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2084 "nabgrm.tab.c"
    break;

  case 122: /* e_list: expr  */
#line 408 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[0].npval), 0 ); }
#line 2090 "nabgrm.tab.c"
    break;

  case 123: /* e_list: expr SYM_COMMA e_list  */
#line 410 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 2096 "nabgrm.tab.c"
    break;

  case 124: /* expr: rval  */
#line 412 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2102 "nabgrm.tab.c"
    break;

  case 125: /* expr: lval assignop expr  */
#line 414 "nabgrm.y"
                                { (yyval.npval) = node( (yyvsp[-1].ival), 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 2108 "nabgrm.tab.c"
    break;

  case 126: /* lval: id  */
#line 415 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2114 "nabgrm.tab.c"
    break;

  case 127: /* lval: ar_lval  */
#line 416 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2120 "nabgrm.tab.c"
    break;

  case 128: /* lval: at_lval  */
#line 417 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2126 "nabgrm.tab.c"
    break;

  case 129: /* ar_lval: lval SYM_LBRACK i_list SYM_RBRACK  */
#line 419 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LBRACK, 0, (yyvsp[-3].npval), (yyvsp[-1].npval) ); }
#line 2132 "nabgrm.tab.c"
    break;

  case 130: /* at_lval: lval SYM_PERIOD SYM_IDENT  */
#line 421 "nabgrm.y"
                                { (yyval.npval) = node( SYM_PERIOD, 0, (yyvsp[-2].npval),
				  node( SYM_ATTRIBUTE, &val, 0, 0 ) ); }
#line 2139 "nabgrm.tab.c"
    break;

  case 131: /* rval: disj  */
#line 423 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2145 "nabgrm.tab.c"
    break;

  case 132: /* rval: disj SYM_OR rval  */
#line 425 "nabgrm.y"
                                { (yyval.npval) = node( SYM_OR, 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 2151 "nabgrm.tab.c"
    break;

  case 133: /* disj: conj  */
#line 426 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2157 "nabgrm.tab.c"
    break;

  case 134: /* disj: conj SYM_AND disj  */
#line 428 "nabgrm.y"
                                { (yyval.npval) = node( SYM_AND, 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 2163 "nabgrm.tab.c"
    break;

  case 135: /* conj: a_expr  */
#line 429 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2169 "nabgrm.tab.c"
    break;

  case 136: /* conj: a_expr relop a_expr  */
#line 431 "nabgrm.y"
                                { (yyval.npval) = node( (yyvsp[-1].ival), 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 2175 "nabgrm.tab.c"
    break;

  case 137: /* a_expr: term  */
#line 432 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2181 "nabgrm.tab.c"
    break;

  case 138: /* a_expr: term addop a_expr  */
#line 434 "nabgrm.y"
                                { (yyval.npval) = node( (yyvsp[-1].ival), 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 2187 "nabgrm.tab.c"
    break;

  case 139: /* term: factor  */
#line 435 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2193 "nabgrm.tab.c"
    break;

  case 140: /* term: factor mulop term  */
#line 437 "nabgrm.y"
                                { (yyval.npval) = node( (yyvsp[-1].ival), 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 2199 "nabgrm.tab.c"
    break;

  case 141: /* factor: primary  */
#line 438 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2205 "nabgrm.tab.c"
    break;

  case 142: /* factor: primary SYM_UPARROW factor  */
#line 440 "nabgrm.y"
                                { (yyval.npval) = node( SYM_UPARROW, 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 2211 "nabgrm.tab.c"
    break;

  case 143: /* primary: lval  */
#line 441 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2217 "nabgrm.tab.c"
    break;

  case 144: /* primary: num  */
#line 442 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2223 "nabgrm.tab.c"
    break;

  case 145: /* primary: string  */
#line 443 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2229 "nabgrm.tab.c"
    break;

  case 146: /* primary: incr  */
#line 444 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2235 "nabgrm.tab.c"
    break;

  case 147: /* primary: unop primary  */
#line 445 "nabgrm.y"
                                { (yyval.npval) = node( (yyvsp[-1].ival), 0, NULL, (yyvsp[0].npval) ); }
#line 2241 "nabgrm.tab.c"
    break;

  case 148: /* primary: id SYM_LPAREN actuals SYM_RPAREN  */
#line 447 "nabgrm.y"
                                { (yyval.npval) = node( SYM_CALL, 0, (yyvsp[-3].npval), (yyvsp[-1].npval) ); }
#line 2247 "nabgrm.tab.c"
    break;

  case 149: /* primary: SYM_LPAREN expr SYM_RPAREN  */
#line 449 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LPAREN, 0, NULL, (yyvsp[-1].npval) ); }
#line 2253 "nabgrm.tab.c"
    break;

  case 150: /* incr: incrop lval  */
#line 450 "nabgrm.y"
                                { (yyval.npval) = node( (yyvsp[-1].ival), 0, 0, (yyvsp[0].npval) ); }
#line 2259 "nabgrm.tab.c"
    break;

  case 151: /* incr: lval incrop  */
#line 451 "nabgrm.y"
                                { (yyval.npval) = node( (yyvsp[0].ival), 0, (yyvsp[-1].npval), 0 ); }
#line 2265 "nabgrm.tab.c"
    break;

  case 152: /* actuals: ap_list  */
#line 452 "nabgrm.y"
                                { (yyval.npval) = (yyvsp[0].npval); }
#line 2271 "nabgrm.tab.c"
    break;

  case 153: /* actuals: %empty  */
#line 453 "nabgrm.y"
                                { (yyval.npval) = NULL; }
#line 2277 "nabgrm.tab.c"
    break;

  case 154: /* ap_list: a_parm  */
#line 454 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[0].npval), 0 ); }
#line 2283 "nabgrm.tab.c"
    break;

  case 155: /* ap_list: a_parm SYM_COMMA ap_list  */
#line 456 "nabgrm.y"
                                { (yyval.npval) = node( SYM_LIST, 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 2289 "nabgrm.tab.c"
    break;

  case 156: /* a_parm: expr  */
#line 457 "nabgrm.y"
                                { (yyval.npval) = node( SYM_PARM, 0, 0, (yyvsp[0].npval) ); }
#line 2295 "nabgrm.tab.c"
    break;

  case 157: /* i_list: expr  */
#line 458 "nabgrm.y"
                                { (yyval.npval) = node( SYM_INDEX, 0, (yyvsp[0].npval), 0 ); }
#line 2301 "nabgrm.tab.c"
    break;

  case 158: /* i_list: expr SYM_COMMA i_list  */
#line 460 "nabgrm.y"
                                { (yyval.npval) = node( SYM_INDEX, 0, (yyvsp[-2].npval), (yyvsp[0].npval) ); }
#line 2307 "nabgrm.tab.c"
    break;

  case 159: /* assignop: SYM_ASSIGN  */
#line 461 "nabgrm.y"
                                { (yyval.ival) = SYM_ASSIGN; }
#line 2313 "nabgrm.tab.c"
    break;

  case 160: /* assignop: SYM_PLUS_ASSIGN  */
#line 463 "nabgrm.y"
                                { (yyval.ival) = SYM_PLUS_ASSIGN; }
#line 2319 "nabgrm.tab.c"
    break;

  case 161: /* assignop: SYM_MINUS_ASSIGN  */
#line 465 "nabgrm.y"
                                { (yyval.ival) = SYM_MINUS_ASSIGN; }
#line 2325 "nabgrm.tab.c"
    break;

  case 162: /* assignop: SYM_STAR_ASSIGN  */
#line 467 "nabgrm.y"
                                { (yyval.ival) = SYM_STAR_ASSIGN; }
#line 2331 "nabgrm.tab.c"
    break;

  case 163: /* assignop: SYM_SLASH_ASSIGN  */
#line 469 "nabgrm.y"
                                { (yyval.ival) = SYM_SLASH_ASSIGN; }
#line 2337 "nabgrm.tab.c"
    break;

  case 164: /* assignop: SYM_MODULUS_ASSIGN  */
#line 471 "nabgrm.y"
                                { (yyval.ival) = SYM_MODULUS_ASSIGN; }
#line 2343 "nabgrm.tab.c"
    break;

  case 165: /* assignop: SYM_UPARROW_ASSIGN  */
#line 473 "nabgrm.y"
                                { (yyval.ival) = SYM_UPARROW_ASSIGN; }
#line 2349 "nabgrm.tab.c"
    break;

  case 166: /* relop: SYM_LESS  */
#line 474 "nabgrm.y"
                                { (yyval.ival) = SYM_LESS; }
#line 2355 "nabgrm.tab.c"
    break;

  case 167: /* relop: SYM_LESS_EQUAL  */
#line 476 "nabgrm.y"
                                { (yyval.ival) = SYM_LESS_EQUAL; }
#line 2361 "nabgrm.tab.c"
    break;

  case 168: /* relop: SYM_EQUAL  */
#line 477 "nabgrm.y"
                                { (yyval.ival) = SYM_EQUAL; }
#line 2367 "nabgrm.tab.c"
    break;

  case 169: /* relop: SYM_NOT_EQUAL  */
#line 478 "nabgrm.y"
                                { (yyval.ival) = SYM_NOT_EQUAL; }
#line 2373 "nabgrm.tab.c"
    break;

  case 170: /* relop: SYM_GREATER_EQUAL  */
#line 480 "nabgrm.y"
                                { (yyval.ival) = SYM_GREATER_EQUAL; }
#line 2379 "nabgrm.tab.c"
    break;

  case 171: /* relop: SYM_GREATER  */
#line 481 "nabgrm.y"
                                { (yyval.ival) = SYM_GREATER; }
#line 2385 "nabgrm.tab.c"
    break;

  case 172: /* relop: SYM_MATCH  */
#line 482 "nabgrm.y"
                                { (yyval.ival) = SYM_MATCH; }
#line 2391 "nabgrm.tab.c"
    break;

  case 173: /* relop: SYM_DONT_MATCH  */
#line 484 "nabgrm.y"
                                { (yyval.ival) = SYM_DONT_MATCH; }
#line 2397 "nabgrm.tab.c"
    break;

  case 174: /* relop: SYM_IN  */
#line 485 "nabgrm.y"
                                { (yyval.ival) = SYM_IN; }
#line 2403 "nabgrm.tab.c"
    break;

  case 175: /* addop: SYM_PLUS  */
#line 486 "nabgrm.y"
                                { (yyval.ival) = SYM_PLUS; }
#line 2409 "nabgrm.tab.c"
    break;

  case 176: /* addop: SYM_MINUS  */
#line 487 "nabgrm.y"
                                { (yyval.ival) = SYM_MINUS; }
#line 2415 "nabgrm.tab.c"
    break;

  case 177: /* mulop: SYM_STAR  */
#line 488 "nabgrm.y"
                                { (yyval.ival) = SYM_STAR; }
#line 2421 "nabgrm.tab.c"
    break;

  case 178: /* mulop: SYM_SLASH  */
#line 489 "nabgrm.y"
                                { (yyval.ival) = SYM_SLASH; }
#line 2427 "nabgrm.tab.c"
    break;

  case 179: /* mulop: SYM_MODULUS  */
#line 490 "nabgrm.y"
                                { (yyval.ival) = SYM_MODULUS; }
#line 2433 "nabgrm.tab.c"
    break;

  case 180: /* mulop: SYM_ATSIGN  */
#line 491 "nabgrm.y"
                                { (yyval.ival) = SYM_ATSIGN; }
#line 2439 "nabgrm.tab.c"
    break;

  case 181: /* incrop: SYM_PLUS_PLUS  */
#line 492 "nabgrm.y"
                                { (yyval.ival) = SYM_PLUS_PLUS; }
#line 2445 "nabgrm.tab.c"
    break;

  case 182: /* incrop: SYM_MINUS_MINUS  */
#line 494 "nabgrm.y"
                                { (yyval.ival) = SYM_MINUS_MINUS; }
#line 2451 "nabgrm.tab.c"
    break;

  case 183: /* unop: SYM_MINUS  */
#line 495 "nabgrm.y"
                                { (yyval.ival) = SYM_NEGATE; }
#line 2457 "nabgrm.tab.c"
    break;

  case 184: /* unop: SYM_NOT  */
#line 496 "nabgrm.y"
                                { (yyval.ival) = SYM_NOT; }
#line 2463 "nabgrm.tab.c"
    break;

  case 185: /* id: SYM_IDENT  */
#line 498 "nabgrm.y"
                                { (yyval.npval) = node( SYM_IDENT, &val, 0, 0 ); }
#line 2469 "nabgrm.tab.c"
    break;

  case 186: /* num: SYM_INT_LIT  */
#line 499 "nabgrm.y"
                                { (yyval.npval) = node( SYM_INT_LIT, &val, 0, 0 ); }
#line 2475 "nabgrm.tab.c"
    break;

  case 187: /* num: SYM_FLOAT_LIT  */
#line 500 "nabgrm.y"
                                { (yyval.npval) = node( SYM_FLOAT_LIT, &val, 0, 0 ); }
#line 2481 "nabgrm.tab.c"
    break;

  case 188: /* string: SYM_STRING_LIT  */
#line 502 "nabgrm.y"
                                { (yyval.npval) = node( SYM_STRING_LIT, &val, 0, 0 ); }
#line 2487 "nabgrm.tab.c"
    break;


#line 2491 "nabgrm.tab.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (YY_("syntax error"));
    }

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;
  ++yynerrs;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  YY_ACCESSING_SYMBOL (yystate), yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturnlab;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturnlab;


/*-----------------------------------------------------------.
| yyexhaustedlab -- YYNOMEM (memory exhaustion) comes here.  |
`-----------------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturnlab;


/*----------------------------------------------------------.
| yyreturnlab -- parsing is finished, clean up and return.  |
`----------------------------------------------------------*/
yyreturnlab:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 503 "nabgrm.y"


#include "lex.yy.c"
