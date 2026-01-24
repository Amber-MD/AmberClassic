# Makefile for ucpp
#
# (c) Thomas Pornin 1999 - 2002
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 4. The name of the authors may not be used to endorse or promote
#    products derived from this software without specific prior written
#    permission.
#
# THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR 
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
# OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

include ../../config.h

.SUFFIXES:
.SUFFIXES:  .c .o

.c.o:
	@echo "[UCPP]  CC $<"
	$(VB)$(CC) -c $(CFLAGS) -DSTAND_ALONE $<

COBJ = mem.o nhash.o cpp.o lexer.o assert.o macro.o eval.o

$(BINDIR)/ucpp: $(COBJ)
	@echo "[UCPP]  CC $@"
	$(VB)$(CC) $(LDFLAGS) -o $@ $(COBJ) $(LIBS)

assert.o: tune.h ucppi.h cpp.h nhash.h mem.h
cpp.o: tune.h ucppi.h cpp.h nhash.h mem.h
eval.o: tune.h ucppi.h cpp.h nhash.h mem.h arith.c arith.h
lexer.o: tune.h ucppi.h cpp.h nhash.h mem.h
macro.o: tune.h ucppi.h cpp.h nhash.h mem.h
mem.o: mem.h
nhash.o: nhash.h mem.h

ucpp: $(BINDIR)/ucpp

uninstall: clean
	-rm -f $(BINDIR)/ucpp

skip:
	echo Skipping installation of the ucpp preprocessor

clean:
	rm -f *.o ucpp core

