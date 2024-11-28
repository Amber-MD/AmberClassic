include ../../config.h
CFLAGS+=$(COPTFLAGS)

SOURCE=calc_r_squared.c create_input.c defaults.c \
	elements.c error_messages.c eval_amber_std.c \
	file_io.c fitting_control.c help_functions.c \
	mem_alloc.c misc_utils.c options_summary.c \
	parameter_optimiser.c param_summary.c print_program_info.c \
	process_command_line.c process_job_control_setting.c \
	process_prmtop.c read_prmtop.c simplex.c write_input.c \
	genetic_algorithm.c bounds_check.c eval_amber_forces.c \
	dihedral_fitting.c wizard.c read_mdcrd.c read_energy.c

OBJECTS=calc_r_squared.o create_input.o defaults.o \
        elements.o error_messages.o eval_amber_std.o \
        file_io.o fitting_control.o help_functions.o \
        mem_alloc.o misc_utils.o options_summary.o \
        parameter_optimiser.o param_summary.o print_program_info.o \
        process_command_line.o process_job_control_setting.o \
        process_prmtop.o read_prmtop.o simplex.o write_input.o \
	genetic_algorithm.o bounds_check.o eval_amber_forces.o \
	dihedral_fitting.o wizard.o read_mdcrd.o read_energy.o

HEADERS=constants.h function_def.h forces.h prmtop_params.h

calc_r_squared.o: $(HEADERS)

paramfit: $(BINDIR)/paramfit

$(BINDIR)/paramfit: 	$(OBJECTS) 
	@echo "[PARAMFIT]  CC $<"
	$(VB)$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(LDFLAGS) $(LM) 

.c.o:
	@echo "[PARAMFIT]  CC $<"
	$(VB)$(CC) -c $(CFLAGS) $(AMBERCFLAGS) -o $@ $<
	
#

clean:
	/bin/rm -f *.o

uninstall: clean
	/bin/rm -f $(BINDIR)/paramfit

