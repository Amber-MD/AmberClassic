include ../../config.h

#
# reduce top level Makefile
#

reduce:
	@(cd toolclasses && $(MAKE) toolclasses)
	@(cd libpdb && $(MAKE) libpdb)
	@(cd reduce_src && $(MAKE) install_reduce)

clean:
	@(cd toolclasses && $(MAKE) clean)
	@(cd libpdb && $(MAKE) clean)
	@(cd reduce_src && $(MAKE) clean)
