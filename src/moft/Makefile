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

# BOOSTLIB=-lboost_program_options -lboost_iostreams -lboost_regex -lboost_system -lboost_timer -lboost_chrono -lboost_filesystem -lboost_graph

moft: $(BINDIR)/metatwist

metatwist.o:  metatwist.cpp $(HEADERS)
	@echo "[MoFT]  CXX $<"
	$(VB)$(CXX) $(COPTFLAGS) $(CFLAGS) -c -std=c++14 -I. -o $@ $<

$(BINDIR)/metatwist:  metatwist.o
	@echo "[MoFT]  CXX $@"
	$(VB)$(CXX) -std=c++14 -o $@ $< -L$(LIBDIR) $(BOOSTLIB) $(FLIBSF)

clean:
	/bin/rm -f *.o
