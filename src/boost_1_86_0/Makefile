include ../../config.h

install_boost:
	./bootstrap.sh --with-libraries=atomic,chrono,date_time,exception,filesystem,graph,iostreams,mpi,predef,program_options,python,regex,serialization,system,test,thread,timer --prefix=$(BASEDIR)
	./b2 install

clean:
	/bin/rm -f b2 bjam project-config.jam boost_config.log
	/bin/rm -rf bin.v2
