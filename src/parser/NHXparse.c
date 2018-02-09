/* A Bison parser, made by GNU Bison 2.5.  */

/* Bison implementation for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2011 Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

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

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.5"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse         yytree_parse
#define yylex           yytree_lex
#define yyerror         yytree_error
#define yylval          yytree_lval
#define yychar          yytree_char
#define yydebug         yytree_debug
#define yynerrs         yytree_nerrs


/* Copy the first part of user declarations.  */

/* Line 268 of yacc.c  */
#line 1 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "NHXnode.h"
#include "NHXtree.h"
#include "NHXannotation.h"


/* 
   Defined this to get more detailed error messages from parser. Only
   for debugging really. To messy for a user. 
*/
/* #define YYERROR_VERBOSE */ 

/* Here comes some C declarations */

extern FILE *yytree_in;
extern char *yytree_text;
extern unsigned int lineno; /* Current line number in input file */
extern unsigned int n_left_parens; /* Number of matched left parens '(' */
extern unsigned int n_right_parens; /* Number of matched right parens ')' */
extern unsigned int n_leaves;	/* Number of found tree leaves*/

int err_flag = 0; /* arve's workaround for the non-working YYRECOVERING() macro */

/* The list of found trees */
struct NHXtree *input_trees;

/* For better error messages, we store the name of the current file */
char *current_filename;

/* When we see an annotation in the form "TAG = 4711", an 
   annotation structure is created already after the equal sign.
   That way we have some space to put the value in. OK, bad 
   explanation, but look at the code... (rules for 'annotation'
   and 'value'.)
*/
struct NHXannotation *current_annotation;

/* Some prototypes */
void set_str_annotation(char *str);
void set_int_annotation(int i);
void set_float_annotation(float f);
void set_int_list_annotation(struct int_list *il);
typedef enum {
  string_type=1, 
  int_type=2, 
  float_type=4, 
  number_type=6,	/* == int_type | float_type */
  int_list_type=8
} type;
type get_annotation_type();
void check_annotation_type(type actual_type);

void err_msg(char *s);
void yyerror(char *s);


/* Line 268 of yacc.c  */
#line 139 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.c"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     LEFT_PAREN = 258,
     RIGHT_PAREN = 259,
     EQUAL = 260,
     COLON = 261,
     SEMICOLON = 262,
     COMMA = 263,
     SEPARATOR = 264,
     APOSTROPHE = 265,
     TAG = 266,
     STRING = 267,
     FLOAT = 268,
     INTEGER = 269,
     NHX_ANNOTATION_START = 270,
     BEEP_ANNOTATION_START = 271,
     ANNOTATION_END = 272,
     SPECIES_NAME = 273,
     IS_DUPLICATION = 274,
     NODE_ID = 275,
     UNKNOWN_ANNOTATION = 276
   };
#endif
/* Tokens.  */
#define LEFT_PAREN 258
#define RIGHT_PAREN 259
#define EQUAL 260
#define COLON 261
#define SEMICOLON 262
#define COMMA 263
#define SEPARATOR 264
#define APOSTROPHE 265
#define TAG 266
#define STRING 267
#define FLOAT 268
#define INTEGER 269
#define NHX_ANNOTATION_START 270
#define BEEP_ANNOTATION_START 271
#define ANNOTATION_END 272
#define SPECIES_NAME 273
#define IS_DUPLICATION 274
#define NODE_ID 275
#define UNKNOWN_ANNOTATION 276




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 293 of yacc.c  */
#line 61 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"

  struct NHXtree *t;		/* For returning full trees */
  struct NHXnode *v;		/* For returning tree nodes */
  struct NHXannotation *a;	/* For handling node annotations */
  float branch_time;
  char *str;         /* Dealing with leaf names */
  unsigned integer;
  struct int_list *il;



/* Line 293 of yacc.c  */
#line 229 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 343 of yacc.c  */
#line 241 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

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
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
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
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  13
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   63

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  22
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  19
/* YYNRULES -- Number of rules.  */
#define YYNRULES  45
/* YYNRULES -- Number of states.  */
#define YYNSTATES  66

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   276

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
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
      15,    16,    17,    18,    19,    20,    21
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     4,     6,     8,    11,    13,    15,    19,
      23,    30,    33,    35,    39,    40,    42,    44,    46,    50,
      51,    54,    56,    58,    60,    61,    63,    65,    68,    73,
      78,    81,    84,    86,    90,    91,    93,    94,    99,   101,
     103,   105,   107,   111,   113,   115
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      23,     0,    -1,    -1,    24,    -1,    25,    -1,    24,    25,
      -1,    27,    -1,    27,    -1,    26,     8,    27,    -1,    28,
      30,    32,    -1,     3,    26,     4,    29,    30,    32,    -1,
       3,     1,    -1,    12,    -1,    10,    12,    10,    -1,    -1,
      14,    -1,    13,    -1,    12,    -1,    10,    12,    10,    -1,
      -1,     6,    31,    -1,     1,    -1,    13,    -1,    14,    -1,
      -1,    33,    -1,    34,    -1,    33,    34,    -1,    15,    36,
      35,    17,    -1,    16,    36,    35,    17,    -1,    15,     1,
      -1,    16,     1,    -1,    37,    -1,    37,    36,    35,    -1,
      -1,     9,    -1,    -1,    12,     5,    38,    39,    -1,     1,
      -1,    12,    -1,    14,    -1,    13,    -1,     3,    40,     4,
      -1,     1,    -1,    14,    -1,    40,    14,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,   109,   109,   110,   113,   114,   117,   121,   122,   130,
     136,   142,   145,   146,   151,   152,   153,   154,   155,   159,
     160,   161,   164,   165,   168,   169,   172,   173,   176,   177,
     178,   179,   182,   183,   186,   187,   191,   190,   193,   196,
     197,   198,   199,   200,   203,   204
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "LEFT_PAREN", "RIGHT_PAREN", "EQUAL",
  "COLON", "SEMICOLON", "COMMA", "SEPARATOR", "APOSTROPHE", "TAG",
  "STRING", "FLOAT", "INTEGER", "NHX_ANNOTATION_START",
  "BEEP_ANNOTATION_START", "ANNOTATION_END", "SPECIES_NAME",
  "IS_DUPLICATION", "NODE_ID", "UNKNOWN_ANNOTATION", "$accept",
  "tree_file", "tree_list", "tree", "subtree_list", "subtree", "leaf",
  "label", "newick_weight", "number", "possible_x_annotation",
  "ext_annotations", "ext_annotation", "annotation_list",
  "possible_separator", "annotation", "$@1", "value", "int_list", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    22,    23,    23,    24,    24,    25,    26,    26,    27,
      27,    27,    28,    28,    29,    29,    29,    29,    29,    30,
      30,    30,    31,    31,    32,    32,    33,    33,    34,    34,
      34,    34,    35,    35,    36,    36,    38,    37,    37,    39,
      39,    39,    39,    39,    40,    40
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     1,     1,     2,     1,     1,     3,     3,
       6,     2,     1,     3,     0,     1,     1,     1,     3,     0,
       2,     1,     1,     1,     0,     1,     1,     2,     4,     4,
       2,     2,     1,     3,     0,     1,     0,     4,     1,     1,
       1,     1,     3,     1,     1,     2
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     0,    12,     0,     3,     4,     6,     0,    11,
       0,     7,     0,     1,     5,    21,     0,    24,    14,     0,
      13,    22,    23,    20,     0,     0,     9,    25,    26,     0,
      17,    16,    15,     0,     8,    30,    35,     0,    31,     0,
      27,     0,    24,    38,     0,     0,    34,     0,    18,    10,
      36,    28,     0,    29,     0,    33,    43,     0,    39,    41,
      40,    37,    44,     0,    42,    45
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     4,     5,     6,    10,     7,     8,    33,    17,    23,
      26,    27,    28,    45,    37,    46,    54,    61,    63
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -33
static const yytype_int8 yypact[] =
{
      32,    24,    -3,   -33,    51,    32,   -33,   -33,     2,   -33,
      44,   -33,    18,   -33,   -33,   -33,    19,    34,    33,    32,
     -33,   -33,   -33,   -33,    28,    29,   -33,    34,   -33,    41,
     -33,   -33,   -33,     2,   -33,   -33,   -33,     3,   -33,     3,
     -33,    45,    34,   -33,    49,    39,    22,    40,   -33,   -33,
     -33,   -33,     3,   -33,    10,   -33,   -33,    46,   -33,   -33,
     -33,   -33,   -33,    12,   -33,   -33
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -33,   -33,   -33,    53,   -33,     0,   -33,   -33,    26,   -33,
      20,   -33,    36,   -32,   -25,   -33,   -33,   -33,   -33
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -35
static const yytype_int8 yytable[] =
{
      39,    11,   -19,    15,    43,   -19,   -19,    47,    16,    12,
     -19,    56,   -19,    57,   -19,    44,    64,   -19,   -19,    34,
      55,    52,    58,    59,    60,     9,    65,     1,    20,    35,
      38,    36,    21,    22,     2,     1,     3,    36,    36,   -32,
     -34,   -34,     2,    29,     3,    30,    31,    32,    18,    24,
      25,    13,    19,    41,    50,    48,    51,    53,    14,    42,
      62,     0,    49,    40
};

#define yypact_value_is_default(yystate) \
  ((yystate) == (-33))

#define yytable_value_is_error(yytable_value) \
  YYID (0)

static const yytype_int8 yycheck[] =
{
      25,     1,     0,     1,     1,     3,     4,    39,     6,    12,
       8,     1,    10,     3,    12,    12,     4,    15,    16,    19,
      52,    46,    12,    13,    14,     1,    14,     3,    10,     1,
       1,     9,    13,    14,    10,     3,    12,     9,     9,    17,
      12,    12,    10,    10,    12,    12,    13,    14,     4,    15,
      16,     0,     8,    12,     5,    10,    17,    17,     5,    33,
      14,    -1,    42,    27
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,    10,    12,    23,    24,    25,    27,    28,     1,
      26,    27,    12,     0,    25,     1,     6,    30,     4,     8,
      10,    13,    14,    31,    15,    16,    32,    33,    34,    10,
      12,    13,    14,    29,    27,     1,     9,    36,     1,    36,
      34,    12,    30,     1,    12,    35,    37,    35,    10,    32,
       5,    17,    36,    17,    38,    35,     1,     3,    12,    13,
      14,    39,    14,    40,     4,    14
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* This macro is provided for backward compatibility. */

#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
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


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (0, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  YYSIZE_T yysize1;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = 0;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                yysize1 = yysize + yytnamerr (0, yytname[yyx]);
                if (! (yysize <= yysize1
                       && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                  return 2;
                yysize = yysize1;
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  yysize1 = yysize + yystrlen (yyformat);
  if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
    return 2;
  yysize = yysize1;

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

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

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
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

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

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
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:

/* Line 1806 of yacc.c  */
#line 109 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { err_msg("No input tree!"); }
    break;

  case 3:

/* Line 1806 of yacc.c  */
#line 110 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { input_trees = (yyvsp[(1) - (1)].t); }
    break;

  case 4:

/* Line 1806 of yacc.c  */
#line 113 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.t) = new_tree((yyvsp[(1) - (1)].v), NULL); }
    break;

  case 5:

/* Line 1806 of yacc.c  */
#line 114 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.t) = new_tree((yyvsp[(2) - (2)].v), (yyvsp[(1) - (2)].t));   }
    break;

  case 8:

/* Line 1806 of yacc.c  */
#line 122 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.v) = new_node(NULL); 
					    (yyvsp[(1) - (3)].v)->parent = (yyval.v);
					    (yyvsp[(3) - (3)].v)->parent = (yyval.v);
					    (yyval.v)->left = (yyvsp[(1) - (3)].v);
					    (yyval.v)->right = (yyvsp[(3) - (3)].v);
					  }
    break;

  case 9:

/* Line 1806 of yacc.c  */
#line 132 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.v) = (yyvsp[(1) - (3)].v);
				  annotate_node((yyval.v), append_annotations((yyvsp[(2) - (3)].a), (yyvsp[(3) - (3)].a)));
				}
    break;

  case 10:

/* Line 1806 of yacc.c  */
#line 138 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyvsp[(2) - (6)].v)->name = (yyvsp[(4) - (6)].str); /* label */
				  annotate_node((yyvsp[(2) - (6)].v), append_annotations((yyvsp[(5) - (6)].a), (yyvsp[(6) - (6)].a)));
				  (yyval.v) = (yyvsp[(2) - (6)].v);
		                }
    break;

  case 11:

/* Line 1806 of yacc.c  */
#line 142 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.v) = NULL; err_msg("Could not parse subtree"); }
    break;

  case 12:

/* Line 1806 of yacc.c  */
#line 145 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.v) = new_node((yyvsp[(1) - (1)].str)); n_leaves++;}
    break;

  case 13:

/* Line 1806 of yacc.c  */
#line 146 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.v) = new_node((yyvsp[(2) - (3)].str)); n_leaves++; }
    break;

  case 14:

/* Line 1806 of yacc.c  */
#line 151 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.str) = NULL; }
    break;

  case 15:

/* Line 1806 of yacc.c  */
#line 152 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { char str[10]; sprintf(str, "%d", (yyvsp[(1) - (1)].integer)); (yyval.str) = strdup(str);}
    break;

  case 16:

/* Line 1806 of yacc.c  */
#line 153 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { char str[10]; sprintf(str, "%f", (yyvsp[(1) - (1)].branch_time)); (yyval.str) = strdup(str);}
    break;

  case 17:

/* Line 1806 of yacc.c  */
#line 154 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.str) = (yyvsp[(1) - (1)].str);}
    break;

  case 18:

/* Line 1806 of yacc.c  */
#line 155 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.str) = (yyvsp[(2) - (3)].str);}
    break;

  case 19:

/* Line 1806 of yacc.c  */
#line 159 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    {(yyval.a) = NULL; }
    break;

  case 20:

/* Line 1806 of yacc.c  */
#line 160 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.a) = new_newick_weight((yyvsp[(2) - (2)].branch_time), NULL); }
    break;

  case 21:

/* Line 1806 of yacc.c  */
#line 161 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.a) = new_newick_weight(0.0, NULL); err_msg("Expected a branchlength");}
    break;

  case 22:

/* Line 1806 of yacc.c  */
#line 164 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.branch_time) = (yyvsp[(1) - (1)].branch_time);}
    break;

  case 23:

/* Line 1806 of yacc.c  */
#line 165 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.branch_time) = (float) (yyvsp[(1) - (1)].integer); }
    break;

  case 24:

/* Line 1806 of yacc.c  */
#line 168 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.a) = NULL;}
    break;

  case 27:

/* Line 1806 of yacc.c  */
#line 173 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.a) = append_annotations((yyvsp[(1) - (2)].a), (yyvsp[(2) - (2)].a)); }
    break;

  case 28:

/* Line 1806 of yacc.c  */
#line 176 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.a) = (yyvsp[(3) - (4)].a); }
    break;

  case 29:

/* Line 1806 of yacc.c  */
#line 177 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.a) = (yyvsp[(3) - (4)].a); }
    break;

  case 30:

/* Line 1806 of yacc.c  */
#line 178 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { err_msg("Syntax error in extended annotations"); (yyval.a) = NULL; }
    break;

  case 31:

/* Line 1806 of yacc.c  */
#line 179 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { err_msg("Syntax error in extended annotations"); (yyval.a) = NULL; }
    break;

  case 33:

/* Line 1806 of yacc.c  */
#line 183 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.a) = append_annotations((yyvsp[(1) - (3)].a), (yyvsp[(3) - (3)].a)); }
    break;

  case 36:

/* Line 1806 of yacc.c  */
#line 191 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { current_annotation = new_annotation((yyvsp[(1) - (2)].str), NULL);}
    break;

  case 37:

/* Line 1806 of yacc.c  */
#line 192 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.a) = current_annotation; }
    break;

  case 38:

/* Line 1806 of yacc.c  */
#line 193 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { err_msg("Syntax error in extended annotations");}
    break;

  case 39:

/* Line 1806 of yacc.c  */
#line 196 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { set_str_annotation((yyvsp[(1) - (1)].str));}
    break;

  case 40:

/* Line 1806 of yacc.c  */
#line 197 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { set_int_annotation((yyvsp[(1) - (1)].integer)); }
    break;

  case 41:

/* Line 1806 of yacc.c  */
#line 198 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { set_float_annotation((yyvsp[(1) - (1)].branch_time)); }
    break;

  case 42:

/* Line 1806 of yacc.c  */
#line 199 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { set_int_list_annotation((yyvsp[(2) - (3)].il)); }
    break;

  case 43:

/* Line 1806 of yacc.c  */
#line 200 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { err_msg("Wrong value type"); }
    break;

  case 44:

/* Line 1806 of yacc.c  */
#line 203 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.il) = new_int_list((yyvsp[(1) - (1)].integer), NULL); }
    break;

  case 45:

/* Line 1806 of yacc.c  */
#line 204 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"
    { (yyval.il) = new_int_list((yyvsp[(2) - (2)].integer), (yyvsp[(1) - (2)].il)); }
    break;



/* Line 1806 of yacc.c  */
#line 1789 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.c"
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
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
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

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
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
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
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
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 2067 of yacc.c  */
#line 207 "/home/auwnm/newPrimeTV/PrimeTV2/src/parser/NHXparse.y"


void
yyerror (char *s) {
  fprintf (stderr, "%s:line %d: %s\n", current_filename, lineno, s);
  if (n_right_parens > n_left_parens) {
     fprintf(stderr, "\tUnbalanced parenthesis!\n");
  }
  if (n_left_parens == 1) {	/* Style is everything! */
    fprintf (stderr, "\tAfter 1 leaf, %d '(' and %d ')' parens: '%s'\n", n_left_parens, n_right_parens, yytree_text);
  } else {
    fprintf (stderr, "\tAfter %d leaves, %d '(' and %d ')' parens: '%s'\n", n_leaves, n_left_parens, n_right_parens, yytree_text);
  }
  err_flag = 1;
}

void
err_msg(char *s) {
  if (err_flag) {	     
     fprintf(stderr, "%s\n", s);
     err_flag = 0;
  }
}	     

/*
  The following functions update the annotation with the actual values, 
  but also performs some type checking with the known tags.
*/

void
set_str_annotation(char *str) {
  check_annotation_type(string_type);
  current_annotation->arg.str = str;
}
void
set_int_annotation(int i) {
  type t = get_annotation_type();
  if (t == int_type) {
    current_annotation->arg.i = i;
  } else if (t == float_type) {
    current_annotation->arg.t = (float) i;
  }
}
void
set_int_list_annotation(struct int_list *il) {
  check_annotation_type(int_list_type);
  current_annotation->arg.il = il;
}
void
set_float_annotation(float f) {
  check_annotation_type(float_type);
  current_annotation->arg.t = f;
}



char *arb_tags[] = {"S",         "AC",          "ID",     "NT",       "BL",       "ET",       "NW",      "EX",      "D",     "TT",   NULL};
type arb_types[] = {string_type, int_list_type, int_type, float_type, float_type, float_type, float_type,int_type,int_type,float_type};


type
get_annotation_type() {
  int i;
  
  for (i=0; arb_tags[i] != NULL; i++) {
    if (strcmp(current_annotation->anno_type, arb_tags[i]) == 0) {
      	return arb_types[i];
    }
  } 
  
  fprintf(stderr, "%s:%d: Error, tag without known type: %s\n", 
	  current_filename, lineno, current_annotation->anno_type);
  exit(18);
}

/*
  A function to verify the type of an annotation given the tag
*/
void 
check_annotation_type(type actual_type) {
  int i;
/*   char errbuf[1024];  */
  
  for (i=0; arb_tags[i] != NULL; i++) {
    if (strcmp(current_annotation->anno_type, arb_tags[i]) == 0) {
      if (arb_types[i] & actual_type) {
	return;
      } else {
	fprintf(stderr, "%s:%d:  Error, wrong type for tag %s!\n", 
		current_filename,
		lineno,
		current_annotation->anno_type);
	exit(17);
      }
    }
  }
}


void
set_globals(const char *filename) {
   /* For error messages */
  current_filename = strdup(filename);
  n_left_parens = 0;
  n_right_parens = 0;
  n_leaves = 0;
  lineno = 1;
  	
#if YYDEBUG == 1
   yydebug = 1;
#endif

}	


/*
  read_tree
  
  Instruct the parser to read a tree from 'filename' or, if filename
  is NULL, STDIN. Return the read tree on success, or NULL.
*/
struct NHXtree *
read_tree(const char *filename) {
  FILE *f = NULL;
  int ret_val;

  if (filename == NULL) {
    yytree_in = stdin;
    set_globals("STDIN");	/* For better error messages */
  } else {
    f = fopen(filename, "r");
    set_globals(filename);	/* For better error messages */
    if (!f) {
      fprintf(stderr, "Could not open tree file '%s' for reading.\n", filename);
      return NULL;
    } else {
      yytree_in = f;
    }
  } 
  
  ret_val = yyparse();

  /* Cleanup */
  if (f != NULL) {
    close(f);
    yytree_in = stdin;
  }

  if (ret_val == 1)
    return NULL;
  else
    return input_trees;
}



/*
  read_tree_from_file_stream
  
  Instruct the parser to read a tree from a file stream.
*/

struct NHXtree *
read_tree_from_file_stream( FILE * f ) {
  int ret_val;
  set_globals("");
  yytree_in = f;
  ret_val = yyparse();
  if (ret_val == 1)
    return NULL;
  else
    return input_trees;
}



/*
  read_tree_string
  
  Instruct the parser to read from a string rather than a
  file. Otherwise behaves like read_tree.
*/
struct NHXtree *
read_tree_string(const char *str) {
  int ret_val;

  if (str == NULL) {
    fprintf(stderr, "Warning: Tried to read a tree from a NULL string.\n");
    return NULL;
  } 
    
  set_globals("<input string>");	/* For better error messages */
  read_from_string(str);

  ret_val = yyparse();

  /* Cleanup */
  close_string_buffer();

  if (ret_val == 1)
    return NULL;
  else
    return input_trees;
}

