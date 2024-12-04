include ../../config.h

HEADERS= \
   metaFFT/ConvCorr.hpp          \
   metaFFT/fftwmain.hpp          \
   metaFFT/fftwwrap.hpp          \
   metatwist/Density.hpp         \
   metatwist/boosthead.hpp       \
   metatwist/lutilities.hpp      \
   metatwist/readccp4.hpp        \
   metatwist/readdx.hpp          \
   metatwist/scattering.hpp      \
   utilities/boosthead.hpp       \
   utilities/interpolation.hpp   \
   utilities/lapack.hpp          \
   utilities/rotations.hpp       \
   utilities/statistics.hpp      \
   utilities/transrot.hpp        \
   utilities/utilities.hpp      

skip:
	@echo "Skipping build of metatwist"

moft: $(BINDIR)/metatwist

metatwist.o:  metatwist.cpp $(HEADERS)
	@echo "[MoFT]  CXX $<"
	$(VB)$(CXX) $(COPTFLAGS) $(CFLAGS) -c -std=c++14 -I. -o $@ $<

$(BINDIR)/metatwist:  metatwist.o
	@echo "[MoFT]  CXX $@"
	$(VB)$(CXX) -std=c++14 -o $@ $< -L$(LIBDIR) $(BOOSTLIB) -lfftw3_omp -lfftw3

clean:
	/bin/rm -f *.o
