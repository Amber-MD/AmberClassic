/*
 *   MxNameLen applies to block, category, item and keyword names...
 *   
 */
#define MxNameLen    		 256 
/*
 *   This is the size of the buffer length allocated to individual
 *   values.
 */

#define CifDefaultSpace           2

#undef   YYLMAX
#define  YYLMAX       1024
#if 0
#undef   YYMAXDEPTH  
#define  YYMAXDEPTH  20000
#undef   YYINITDEPTH 
#define  YYINITDEPTH  1000
#define  YYPRINT(file, type, value) cifpprint(file, type, value)
#endif

#define MAXVALUELENGTH   8182

#define TRUE                   1
#define FALSE                  0

typedef struct NdbCifRowFormat {
  char **columns;
} NdbCifRowFormat;

typedef struct NdbCifCategoryFormat {
  int numCol;
  int allCol;
  int curCol;
  int allRow;
  int numRow;
  int curRow;
  char categoryName[MxNameLen];
  char **colNames;
  NdbCifRowFormat *rows;
} NdbCifCategoryFormat;

typedef struct NdbCifDatablockFormat {
  int numCategory; /* Number of categories in this datablock */
  int allCategory; /* Allocated space */
  int curCategory; /* index of the current category */
  char datablockName[MxNameLen]; 
  NdbCifCategoryFormat *categories;
} NdbCifDatablockFormat;  

typedef struct NdbCifDatablocksFormat {
  int numDatablock; /* Number of datablocks in this structure */
  int curDatablock; /* index of the current datablock */
  int allDatablock; /* Allocated space */
  NdbCifDatablockFormat *datablocks;
} NdbCifDatablocksFormat;


int cifpparse();
// Parser routines:
void ndb_cif_process_item_name_list();
void ndb_cif_process_value_list();
void ndb_cif_process_item_name_value_pair();
void ndb_cif_process_loop_declaration();

char ndb_cif_set_null_char(char new_null);
int ndb_cif_read_file(FILE *fp);
int ndb_cif_write_file(FILE *fp);
int ndb_nef_write_file(FILE *fp);
void ndb_cif_print_item_name(FILE *fp, char *itemName, int *linePos);
void ndb_cif_print_item_value(FILE *fp, char *itemValue, int *linePos);
void ndb_cif_write_category(FILE *fp);
int ndb_cif_init(void);
int ndb_cif_close(void);
int ndb_cif_new_datablock(char *datablockName);
int  ndb_cif_put_datablock_name(char *datablockName);
int ndb_cif_new_category(char *categoryName);
int ndb_cif_new_row(void);
int ndb_cif_insert_new_row(int rowNo);
int ndb_cif_rewind_datablock(void);
int ndb_cif_rewind_category(void);
int ndb_cif_rewind_row(void);
int ndb_cif_rewind_column(void);
int ndb_cif_reset_datablocks(void);
int ndb_cif_reset_datablock(void);
int ndb_cif_reset_datablock_by_id(int datablockId);
int ndb_cif_reset_category(void);
int ndb_cif_reset_category_by_id(int datablockId, int categoryId);
int ndb_cif_remove_datablock(void);
int ndb_cif_remove_datablock_by_id(int datablockId);
int ndb_cif_remove_datablock_by_name(char *datablockName);
int ndb_cif_remove_category(void);
int ndb_cif_remove_category_by_id(int datablockId, int categoryId);
int ndb_cif_remove_category_by_name(char *datablockName, char *categoryName);
int ndb_cif_remove_row(void);
int ndb_cif_remove_row_by_id(int datablockId, int categoryId, int rowId);
int ndb_cif_next_datablock(void);
int ndb_cif_next_category(void);
int ndb_cif_next_row(void);
int ndb_cif_next_column(void);
int ndb_cif_move_datablock(char *datablockName);
int ndb_cif_move_datablock_by_id(int datablockId);
int ndb_cif_move_category_by_name(char *datablockName, char *categoryName);
int ndb_cif_move_category_by_id(int datablockId, int categoryId);
int ndb_cif_move_category(int categoryId);
int ndb_cif_move_row_by_name(char *datablockName, char *categoryName, int rowId);
int ndb_cif_move_row_by_id(int datablockId, int categoryId, int rowId);
int ndb_cif_move_row(int rowId);
int ndb_cif_count_datablock(void);
int ndb_cif_count_category(void);
int ndb_cif_count_row(void);
int ndb_cif_count_column(void);
int ndb_cif_current_datablock(void);
int ndb_cif_current_datablock_name(char *datablockName);
int ndb_cif_current_category(void);
int ndb_cif_current_category_name(char *categoryName);
int ndb_cif_current_row(void);
int ndb_cif_current_col(void);
int ndb_cif_put_item_keyword(char *itemKeyword);
int ndb_cif_remove_item_keyword(char *itemKeyword);
int ndb_cif_get_item_name(int colId, char *itemName);
int ndb_cif_get_item_shname(int colId, char *itemName);
int ndb_cif_get_item_value(int colId, char *fieldValue, int maxFieldLen);
char *ndb_cif_copy_item_value(int colId);
int ndb_cif_output_item(FILE *fp, int colId);
int ndb_cif_item_row_1_key(int colId, char *fieldValue);
int ndb_cif_put_item_value(int colId, char *fieldValue);
int ndb_cif_get_category_name_from_item_name(char *categoryName, char *itemName);
int ndb_cif_get_item_keyword_from_item_name(char *itemKeyword, char *itemName);
int ndb_cif_get_datablock_id(char *datablockName);
int ndb_cif_get_category_id(char *datablockName, char *categoryName);
int ndb_cif_get_column_id(char *datablockName, char *categoryName, char *itemKeyword);
void ndb_cif_print_datablock(FILE *fp);
void ndb_cif_pretty_print_datablock(FILE *fp);
void ndb_cif_print_category(FILE *fp, char *category);
void ndb_cif_print_datablocks(FILE *fp);
int get_column_index(int blockIndex, int categoryIndex, char *columnName);
int get_category_index(int blockIndex, char *categoryName);


#ifdef CIF_GLOBAL
	extern FILE *cifpin;
	char TempKeyword[MxNameLen+1], TempValue[MAXVALUELENGTH+1];
	NdbCifDatablocksFormat cifFiles;
	int  lineNo;      
#else
	extern char TempKeyword[MxNameLen+1], TempValue[MAXVALUELENGTH+1];
	extern FILE *cifpin;
	extern int  lineNo;
	extern NdbCifDatablocksFormat cifFiles;
#endif
