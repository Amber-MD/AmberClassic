#  just redirect things to lower-level Makefiles

install::
	cd src && $(MAKE) install

test::
	cd test && $(MAKE) test.$(INSTALLTYPE)

clean::
	cd src && $(MAKE) clean

uninstall:: clean
	touch config.h
	cd src && $(MAKE) uninstall
	/bin/rm -f config.h

distclean:: clean uninstall
	touch config.h
	cd src && $(MAKE) distclean
	/bin/rm -f config.h

