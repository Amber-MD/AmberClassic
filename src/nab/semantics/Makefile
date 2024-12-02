include ../../../config.h

RULES= \
SYM_ALLOCATE.rules \
SYM_AND.rules \
SYM_ASSERT.rules \
SYM_ASSIGN.rules \
SYM_ATSIGN.rules \
SYM_CALL.rules \
SYM_DEALLOCATE.rules \
SYM_DELETE.rules \
SYM_DONT_MATCH.rules \
SYM_EQUAL.rules \
SYM_FOREACH.rules \
SYM_GREATER.rules \
SYM_GREATER_EQUAL.rules \
SYM_IN.rules \
SYM_INDEX.rules \
SYM_LBRACK.rules \
SYM_LESS.rules \
SYM_LESS_EQUAL.rules \
SYM_LPAREN.rules \
SYM_MATCH.rules \
SYM_MINUS.rules \
SYM_MINUS_ASSIGN.rules \
SYM_MINUS_MINUS.rules \
SYM_MODULUS.rules \
SYM_MODULUS_ASSIGN.rules \
SYM_NEGATE.rules \
SYM_NOT.rules \
SYM_NOT_EQUAL.rules \
SYM_OR.rules \
SYM_PARM.rules \
SYM_PLUS.rules \
SYM_PLUS_ASSIGN.rules \
SYM_PLUS_PLUS.rules \
SYM_SLASH.rules \
SYM_SLASH_ASSIGN.rules \
SYM_STAR.rules \
SYM_STAR_ASSIGN.rules \
SYM_TEST.rules \
SYM_UPARROW.rules \
SYM_UPARROW_ASSIGN.rules

checkexpr.c: $(RULES) attribute.tab parse2 ../nab.h ../nabgrm.y
	@echo "[NAB] PARSE2 RULES"
	$(VB)./parse2 --y-tab-h ../nabgrm.tab.h --attributes-tab \
       ./attribute.tab --nab-h-path ../nab.h -o $@ *.rules

parse2: parse2.c
	@echo "[NAB] CC $<"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

clean::
	/bin/rm -f parse2 *.html
