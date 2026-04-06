#  just redirect things to lower-level Makefiles

include config.h

install::
	cd src && $(MAKE) install

test::
	cd test && $(MAKE) test.$(INSTALLTYPE)

clean::
	cd src && $(MAKE) clean

uninstall::
	touch config.h
	cd src && $(MAKE) uninstall

distclean::
	touch config.h
	cd src && $(MAKE) distclean

