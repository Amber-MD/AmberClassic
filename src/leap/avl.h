#ifndef AVL_H
#define AVL_H
/*
 *  'application' avl stuff (BPLUS-style interface)
 */

/*
 * Duplicate key modes (IX_DESC.dup_keys):
 *
 *   IX_NODUPKEY  (0) -- keys must be unique. Adding a duplicate key is an error.
 *
 *   IX_DUPKEY    (1) -- the same key may appear with different recptr values.
 *                       Internally, key+recptr together form the composite tree key,
 *                       so each unique key+recptr pair occupies its own node.
 *
 *   IX_DUPKEYREC (2) -- identical key+recptr pairs are allowed; a repeat count is
 *                       maintained per node. Adding an existing key+recptr increments
 *                       the count; deleting decrements it and only removes the node
 *                       when the count reaches zero.
 *
 * Note: in all modes, compkey() uses key+recptr as the full comparison when
 * dup_keys is nonzero, so the tree structure is always unambiguous.
 */

/*
 * find_key(IX_REC *pe, IX_DESC *pix, int direction)
 *
 * Searches for pe->key and positions the traversal cursor for subsequent
 * next_key()/prev_key() calls. In dup modes, direction determines which
 * record of the matching key group the cursor lands on: direction > 0 lands
 * on the first record, direction < 0 lands on the last record.
 *
 * To traverse all records sharing a key in IX_DUPKEY mode, call find_key()
 * to land on the first (or last) record of the group, then call next_key()
 * (or prev_key()) repeatedly until pe->key changes.
 *
 * The direction argument widens the search when no exact match exists:
 *   direction == 0  -- exact match only
 *   direction <  0  -- if not found, accept the nearest lesser key
 *   direction >  0  -- if not found, accept the nearest greater key
 * This affects where the cursor lands for subsequent iteration, but the
 * return value always reflects whether an exact key match was found.
 *
 * On success, pe->recptr and pe->count are updated but pe->key is not copied.
 *
 * Returns IX_OK if an exact match was found, IX_FAIL otherwise (including
 * when only a neighbor was found via direction).
 */

/*
 * locate_key(IX_REC *pe, IX_DESC *pix, int direction)
 *
 * Same search logic and cursor positioning as find_key(), but intended for
 * range/proximity queries where the caller wants the actual tree record nearest
 * to the search key, not just a hit/miss result.
 *
 * Unlike find_key(), always copies the full found record (key+recptr+count)
 * back into pe, even when the result is only a neighbor.
 *
 * Returns:
 *   IX_OK   -- exact key match found
 *   IX_FAIL -- only a neighbor was found (pe holds that neighbor's data)
 *   IX_END  -- nothing found in the requested direction (tree boundary reached)
 *
 * The IX_END vs IX_FAIL distinction is important for scanning: IX_END means
 * the iteration has walked off the edge of the tree, while IX_FAIL means the
 * tree has content but not an exact match.
 */

/*
 * has_key(IX_REC *pe, IX_DESC *pix)
 *
 * Tests whether pe->key exists in the tree without disturbing the traversal
 * cursor. In IX_DUPKEY mode, pe->recptr may be set to test for a specific
 * key+recptr pair; in other modes pe->recptr is ignored.
 *
 * Returns IX_OK if found, IX_FAIL if not.
 */

/* NOTE: IX_REC with IX_LEN_CSTRING means it is in reality variable-length
 * C struct, but we just use a fixed size good enough for the code base
 */
#define IX_DEFAULTKEYLEN 256     /*  size of default key  */
#define IX_LEN_CSTRING 0

typedef void *IX_RECPOS;

typedef struct {
	IX_RECPOS recptr;           /* User data pointer */
	unsigned int count;      /* duplication count (dupkeys==2), decrement with delete */
	char key[IX_DEFAULTKEYLEN]; /* actually can be of any length */
} rectype;

typedef rectype	IX_REC;

typedef struct {
	void *root;
	int keylength; /* zero for null-terminated strings */
	int dup_keys;
		/*
		 *  0 -- repeated key causes an error message;
		 *  1 -- repeated key & rec cause an error message;
		 *  2 -- complete repetitions allowed, use repetition count.
		 */
} IX_DESC;
#define IX_NODUPKEY  0
#define IX_DUPKEY    1
#define IX_DUPKEYREC 2

#define IX_OK   1
#define IX_FAIL 0
#define IX_END (-2)   // End of index -- for locate_key(), next_key(), prev_key()

// Create & Destroy AVL tree dictionary
extern int	create_index(IX_DESC *pix, int dup, int keylength);
extern int	destroy_index(IX_DESC *pix);

extern int	find_key(IX_REC *pe, IX_DESC *pix, int direction);
extern int	locate_key(IX_REC *pe, IX_DESC *pix, int direction);
extern int      has_key(IX_REC *pe, IX_DESC *pix);

extern int	add_key(IX_REC *pe, IX_DESC *pix);
extern int	delete_key(IX_REC *pe, IX_DESC *pix);

// Move index state pointer to first or last key
extern int	first_key(IX_DESC *pix);
extern int	last_key(IX_DESC *pix);
// Use this to retreive subsequent key
// NOTE: tree is re-searched for current key, not cached, so state change won't segfault.
// IX_END means no more keys (not last key)
extern int	next_key(IX_REC *pe, IX_DESC *pix);
extern int	prev_key(IX_REC *pe, IX_DESC *pix);

#endif // AVL_H
