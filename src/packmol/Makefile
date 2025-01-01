include ../../config.h

oall = cenmass.o \
       gencan.o \
       pgencan.o \
       initial.o \
       title.o \
       setsizes.o \
       exit_codes.o \
       getinp.o \
       strlength.o \
       output.o \
       checkpoint.o \
       writesuccess.o \
       fparc.o \
       gparc.o \
       gwalls.o \
       comprest.o \
       comparegrad.o \
       packmol.o \
       polartocart.o \
       resetboxes.o \
       tobar.o \
       setijk.o \
       setibox.o \
       restmol.o \
       swaptype.o \
       swaptypemod.o \
       ahestetic.o \
       heuristics.o \
       flashsort.o \
       jacobi.o \
       random.o \
       sizes.o \
       usegencan.o \
       compute_data.o \
       flashmod.o \
       computef.o \
       computeg.o \
       input.o \
	   gencan_ieee_signal_routines.o

install_packmol: $(BINDIR)/packmol

$(BINDIR)/packmol: $(oall)
	@echo "[PACKMOL]  FC $@"
	$(VB)$(FC) -o $@ $(oall)

modules = exit_codes.o sizes.o compute_data.o usegencan.o input.o flashmod.o swaptypemod.o ahestetic.o

exit_codes.o : exit_codes.f90
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c exit_codes.f90
sizes.o : sizes.f90 
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c sizes.f90
compute_data.o : compute_data.f90 sizes.o
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c compute_data.f90
input.o : input.f90 sizes.o 
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c input.f90
flashmod.o : flashmod.f90 sizes.o 
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c flashmod.f90
usegencan.o : usegencan.f90 sizes.o
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c usegencan.f90
swaptypemod.o : swaptypemod.f90 
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c swaptypemod.f90
ahestetic.o : ahestetic.f90 
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c ahestetic.f90
#
# Code compiled only for all versions
#
cenmass.o : cenmass.f90 $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c cenmass.f90
initial.o : initial.f90 $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c initial.f90
title.o : title.f90 $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c title.f90
setsizes.o : setsizes.f90 $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c setsizes.f90
getinp.o : getinp.f90  $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c getinp.f90
strlength.o : strlength.f90  
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c strlength.f90
output.o : output.f90  $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c output.f90
checkpoint.o : checkpoint.f90  $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c checkpoint.f90
writesuccess.o : writesuccess.f90  $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c writesuccess.f90
fparc.o : fparc.f90 $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c fparc.f90
gparc.o : gparc.f90 $(modules)   
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c gparc.f90
gwalls.o : gwalls.f90 $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c gwalls.f90
comprest.o : comprest.f90 $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c comprest.f90
comparegrad.o : comparegrad.f90 $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c comparegrad.f90
packmol.o : packmol.f90 $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c packmol.f90
polartocart.o : polartocart.f90 $(modules)   
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c polartocart.f90
resetboxes.o : resetboxes.f90 $(modules)   
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c resetboxes.f90
tobar.o : tobar.f90 $(modules)   
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c tobar.f90
setijk.o : setijk.f90 $(modules)   
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c setijk.f90
setibox.o : setibox.f90 $(modules)   
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c setibox.f90
restmol.o : restmol.f90 $(modules)   
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c restmol.f90
swaptype.o : swaptype.f90 $(modules)   
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c swaptype.f90
heuristics.o : heuristics.f90 $(modules)   
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c heuristics.f90
flashsort.o : flashsort.f90 $(modules)   
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c flashsort.f90
jacobi.o : jacobi.f90
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c jacobi.f90
pgencan.o : pgencan.f90 $(modules)
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c pgencan.f90
random.o : random.f90 
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c random.f90
computef.o : computef.f90 $(modules)   
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c computef.f90
computeg.o : computeg.f90 $(modules)   
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(FFLAGS) -c computeg.f90
gencan_ieee_signal_routines.o : gencan_ieee_signal_routines.f90
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(IEEE_SIGNAL_FFLAGS) -c gencan_ieee_signal_routines.f90
gencan.o : gencan.f
	@echo "[PACKMOL]  FC $<"
	$(VB)$(FC) $(GENCANFFLAGS) -c gencan.f 

clean: 
	@\rm -f ./*.o ./*.mod 
