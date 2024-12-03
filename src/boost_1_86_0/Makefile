include ../../config.h

install_boost:
	@echo "[BOOST] bootstrap"
	@./bootstrap.sh --with-libraries=iostreams,program_options \
          --prefix=$(BASEDIR)  > boost_build.log
	@echo "[BOOST] b2 install"
	@./b2 install >> boost_build.log

clean:
	/bin/rm -f b2 bjam project-config.jam boost_config.log boost_build.log
	/bin/rm -rf bin.v2
