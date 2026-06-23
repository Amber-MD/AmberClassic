
/*
 *  'inner' avl stuff
 */
/* way3.h */

typedef int way3; /* -1, 0, 1 */

#define way3stop  ((way3)0)
#define way3left ((way3)-1)
#define way3right ((way3)1)

way3 makeway3(void);

#define way3sum(x,y) ((x)+(y))
/* assume x!=y */

#define way3opp(x) (-(x))

way3 way3opp2(void);

way3 way3random(void);

/* node.h */

typedef struct _node {
   struct _node *ptr[2]; /* left, right */
   way3 balance:2, trace:2;
   rectype data;
} node;

#define stepway(n,x) (((n)->ptr)[way3ix(x)])
#define stepopp(n,x) (((n)->ptr)[way3ix(way3opp(x))])

node *allocnode(void);
void freenode(void);
node *swapptr(void);
int way3ix(void);


/* tree.h */

#define SRF_FINDEQUAL 1
#define SRF_FINDLESS  2
#define SRF_FINDGREAT 4
#define SRF_SETMARK   8
#define SRF_FROMMARK 16

rectype *avltree_search(void);
rectype *avltree_insert(void);
rectype *avltree_delete(void);
void avltree_first(void);
void avltree_last(void);
long avltree_clear(void);

#define avltree_init(x) (*(x)=NULL)

int compkey(void);
void copydata(void);

/* 'PLUS' interface */

#define IX_OK   1
#define IX_FAIL 0
#define EOIX (-2)

int create_index(void);
int destroy_index(void);
int find_key(void);
int locate_key(void);
int add_key(void);
int delete_key(void);
int first_key(void);
int last_key(void);
int next_key(void);
int prev_key(void);
int find_exact(void);

