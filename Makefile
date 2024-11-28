#  just redirect things to lower-level Makefiles

include config.h

install::
	cd src && $(MAKE) install

test::
	cd test && $(MAKE) test.$(INSTALLTYPE)

clean::
	cd src && $(MAKE) clean

uninstall::
	cd src && $(MAKE) uninstall

distclean::
	cd src && $(MAKE) distclean
