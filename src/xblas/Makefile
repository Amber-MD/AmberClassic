include ../../config.h

lib-amb: $(LIBDIR)/libxblas-amb.a

$(LIBDIR)/libxblas-amb.a:
	cd src/common && $(MAKE) lib-amb
	cd src/dot && $(MAKE) lib-amb
	cd src/axpby && $(MAKE) lib-amb
	cd src/waxpby && $(MAKE) lib-amb
	cd src/gemv && $(MAKE) lib-amb
	cd src/gemm && $(MAKE) lib-amb
	cd src/gemv2 && $(MAKE) lib-amb

clean:
	cd src/common && $(MAKE) clean
	cd src/dot && $(MAKE) clean
	cd src/axpby && $(MAKE) clean
	cd src/waxpby && $(MAKE) clean
	cd src/gemv && $(MAKE) clean
	cd src/gemm && $(MAKE) clean
	cd src/gemv2 && $(MAKE) clean
