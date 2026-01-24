#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Manual.h"
#include "Constants.h"

//-----------------------------------------------------------------------------
// HorizontalRule: prints a horizontal bar to whatever file (including stdout)
//                 specified.                                   
//                                                                      
// Arguments:                                                           
//   outp:   the target file                                            
//   int n:  the number of carriage returns to add after the rule       
//-----------------------------------------------------------------------------
void HorizontalRule(FILE *outp, int n)
{
  int i;

  fprintf(outp, "<++>---------------------------------------------------------"
         "--------------<++>\n");
  for (i = 0; i < n; i++) {
    fprintf(outp, "\n");
  }
}

//-----------------------------------------------------------------------------
// PrintSplash: print the splash lines for mdgx.  This is where to put  
//              information about the program's primary authorship and  
//              copyrights.                                             
//                                                                      
// Arguments:                                                           
//   outp:  the target file                                             
//-----------------------------------------------------------------------------
void PrintSplash(FILE *outp)
{
  HorizontalRule(outp, 0);
  fprintf(outp,
          "<++> mdgx: A molecular dynamics engine in the AMBER suite of "
          "programs      <++>\n"
          "<++>                                                             "
          "          <++>\n"
          "<++> Written by David S. Cerutti, Case Group (2009 - 2020)       "
          "          <++>\n");
  HorizontalRule(outp, 1);
}

//-----------------------------------------------------------------------------
// PrintVADesc: this function prints a variable name, alias, and description
//              using the specified formatting.             
//                                                                      
// Arguments:                                                           
//   leadspace: the number of leading white space characters            
//   vname:     the variable name                                       
//   vnlen:     the amount of space to give the variable name           
//   valias:    the variable alias                                      
//   valen:     the amount of space to give the variable alias          
//   vdesc:     the variable description                                
//   vdlen:     the amount of space to give the variable description (if the
//              description is longer than the amount of space alotted,   
//              additional lines will be used)               
//   vdindent:  the amount of indentation to apply to extra lines of the
//              variable description                                  
//   outp:      the output file                                         
//-----------------------------------------------------------------------------
void PrintVADesc(int leadspace, char* vname, int vnlen, char* valias,
                 int valen, char* vdesc, int vdlen, int vdindent, FILE *outp)
{
  int i, j, endfound, vpivot;
  char scratch[4096];

  // Print leading spaces
  for (i = 0; i < leadspace; i++) {
    fprintf(outp, " ");
  }

  // Extend the variable name, then print it
  j = strlen(vname);
  strcpy(scratch, vname);
  for (i = j; i < vnlen; i++) {
    scratch[i] = ' ';
  }
  scratch[vnlen] = '\0';
  fprintf(outp, "%s", scratch);

  // Extend the alias name, then print it
  j = strlen(valias);
  strcpy(scratch, valias);
  for (i = j; i < valen; i++) {
    scratch[i] = ' ';
  }
  scratch[valen] = '\0';
  fprintf(outp, "%s", scratch);

  // Parse the description
  endfound = 0;
  j = 0;
  while (endfound == 0) {
    for (i = 0; i < vdlen; i++) {
      if (vdesc[j+i] == ' ') {
        i++;
        vpivot = j+i;
      }
      if (vdesc[j+i] == '\0') {
        vpivot = j+i;
        endfound = 1;
        break;
      }
    }
    strncpy(scratch, &vdesc[j], vpivot-j);
    scratch[vpivot-j] = '\0';
    fprintf(outp, "%s\n", scratch);
    if (j == 0 && endfound == 0) {
      vdlen -= vdindent;
    }
    j = vpivot;
    while (vdesc[j] == ' ') {
      j++;
    }
    if (vdesc[j] == '\0') {
      endfound = 1;
    }
    if (endfound == 0) {
      for (i = 0; i < leadspace + vnlen + valen + vdindent; i++) {
        fprintf(outp, " ");
      }
    }
  }
}

//-----------------------------------------------------------------------------
// PrintParagraph: prints a paragraph within a specified width.  No leading
//                 white space or other columns are provided.  However, each  
//                 paragraph is terminated by printing one additional carriage
//                 return (so that paragaphs are separated by a blank line).
//
// Arguments:                                                           
//   vpar:   the paragraph string (carriage returns should not be included in
//           the string)                                    
//   width:  the width of the text to print (the text will be broken up at
//           whitespace characters to prevent lines from exceeding this width)
//   leader: leading words to put at the start of each line (if NULL, it will
//           be ignored)
//   outp:   the target file for this paragraph                         
//-----------------------------------------------------------------------------
void PrintParagraph(char *vpar, int width, char* leader, FILE *outp)
{
  int i, j, k, vpivot, endfound, leadlen;
  char* scratch;

  // Determine if there is a leading string to print at the start of each line
  if (leader != NULL) {
    width -= strlen(leader)+1;
  }
  if (width < 0) {
    printf("PrintParagraph >> Error.  Leading text exceeds the stated line "
           "width.\n");
    exit(1);
  }

  // Allocate space for staging the line to be written
  scratch = (char*)malloc(MAXLINE*sizeof(char));
  endfound = 0;
  j = 0;
  while (endfound == 0) {
    vpivot = j + width;
    for (i = 0; i < width; i++) {
      if (vpar[j+i] == ' ') {
        i++;
        vpivot = j+i;
      }
      if (vpar[j+i] == '\0') {
        vpivot = j+i;
        endfound = 1;
        break;
      }
    }

    // If we over-ran the length of the printable text,
    // insert an elipsis and then a break
    if (endfound == 0 && i == width && vpivot == j+width) {
      vpivot = j+i-3;
      strncpy(scratch, &vpar[j], vpivot-j);
      scratch[vpivot-j] = '.';
      scratch[vpivot-j+1] = '.';
      scratch[vpivot-j+2] = '.';
      scratch[vpivot-j+3] = '\0';
    }
    else {
      strncpy(scratch, &vpar[j], vpivot-j);
      scratch[vpivot-j] = '\0';
    }
    if (leader != NULL) {
      fprintf(outp, "%s %s\n", leader, scratch);
    }
    else {
      fprintf(outp, "%s\n", scratch);
    }
    j = vpivot;
  }
  fprintf(outp, "\n");

  // Free allocated memory
  free(scratch);
}

//-----------------------------------------------------------------------------
// PrintUsage: print a brief set of usage instructions, to help users get
//             started.                                             
//-----------------------------------------------------------------------------
void PrintUsage()
{
  printf("Usage:\n"
         ">> mdgx -i <input file>         (simplest operation)\n"
         ">> mdgx <arguments>             (operation similar to SANDER)\n"
         ">> mdgx <documentation name>    (prints documentation)\n\n");
  PrintParagraph("Command-line information may be entered using the same "
                 "arguments as the SANDER and PMEMD programs.  Alternatively, "
                 "all input may be provided in a single file using the -i "
                 "argument.  Any of the following arguments may also be used "
                 "to obtain additional documentation:", 79, NULL, stdout);
  printf("  -INPUT:    print all command line input options\n"
         "  -IFILE:    documentation on input file format\n"
         "  -FILES:    print descriptions of &files namelist variables "
         "(these may also\n"
         "             be entered as command line input)\n"
         "  -CNTRL:    print descriptions of &cntrl namelist variables (most "
         "are similar\n"
         "             to SANDER variables, but some are unique to mdgx and "
         "some SANDER\n"
         "             variables are not supported)\n"
         "  -EWALD:    print &ewald namelist variables\n"
         "  -FORCE:    print &force namelist variables\n"
         "  -FITQ:     print &fitq (charge fitting) namelist variables\n"
         "  -PARAM:    print &param (bonded term fitting) namelist "
	 "variables\n"
         "  -IPOLQ:    print &ipolq (Implicitly Polarized Charge) namelist "
         "variables\n"
         "  -CONFIGS:  print &configs (small molecule conformation "
	 "generation) keywords\n"
         "  -PPTD:     print &pptd (small oligomer molecular dynamics) "
	 "keywords\n"
	 "  -RULE:     print &rule keywords for specifying virtual sites, "
	 "extra points\n"
	 "  -PARMEDIT: print topology modification keywords, essential for "
	 "committing\n"
	 "             added virtual sites to files for later use\n"
         "  -ATTR:     attributions of certain aspects of the code, with "
         "references\n\n");
}

//-----------------------------------------------------------------------------
// PrintCommandLineInputOptions: this function essentially reproduces what the
//                               AMBER manual already has for SANDER and PMEMD 
//                               command line input, but since mdgx has some
//                               new features it is necessary to have
//                               independent documentation.
//-----------------------------------------------------------------------------
void PrintCommandLineInputOptions()
{
  PrintSplash(stdout);
  PrintVADesc(2, "-O", 5, " ", 2, "Overwrite output files if they exist "
              "(appending files with the -A option found in SANDER is "
              "currently not supported)", 71, 0, stdout);
  PrintVADesc(2, "-i", 5, " ", 2, "(input) control data for an energy "
              "minimization / molecular dynamics run", 71, 0, stdout);
  PrintVADesc(2, "-o", 5, " ", 2, "(output) human-readable state information "
              "and diagnostics", 71, 0, stdout);
  PrintVADesc(2, "-p", 5, " ", 2, "(input) molecular topology file (AMBER "
              "prmtop format)", 71, 0, stdout);
  PrintVADesc(2, "-p2", 5, " ", 2, "(input) molecular topology file; if "
              "thermodynamic integration is active, this topology describes "
              "the final state while the original topology describes the "
              "initial state", 71, 0, stdout);
  PrintVADesc(2, "-xpt", 5, " ", 2, "(input) extra points rule file directing "
              "mdgx to add extra points to the topology at run time", 71, 0,
              stdout);
  PrintVADesc(2, "-xpt2", 5, " ", 2, "(input) extra points rule file for the "
              "topology specified by the -p2 flag", 71, 0, stdout);
  PrintVADesc(2, "-c", 5, " ", 2, "(input) initial coordinates (and, "
              "optionally, velocities) and periodic box dimensions", 71, 0,
              stdout);
  PrintVADesc(2, "-c2", 5, " ", 2, "(input) initial coordinates; if "
              "thermodynamic integration is active, this second set of input "
              "coordinates pertains to the initial coordinates of the system "
              "in its final state as the mixing parameter lambda goes to 1",
              71, 0, stdout);
  PrintVADesc(2, "-d", 5, " ", 2, "(output) comprehensive force / energy "
              "report file", 71, 0, stdout);
  PrintVADesc(2, "-x", 5, " ", 2, "(output) coordinate trajectory file", 71,
              0, stdout);
  PrintVADesc(2, "-x2", 5, " ", 2, "(output) coordinate trajectory file; only "
              "used when thermodynamic integration is active", 71, 0, stdout);
  PrintVADesc(2, "-v", 5, " ", 2, "(output) velocity trajectory file", 71, 0,
              stdout);
  PrintVADesc(2, "-v2", 5, " ", 2, "(output) velocity trajectory file; only "
              "used when thermodynamic integration is active", 71, 0, stdout);
  PrintVADesc(2, "-e", 5, " ", 2, "(output) energy data over trajectory", 71,
              0, stdout);
  PrintVADesc(2, "-r", 5, " ", 2, "(output) checkpoint (and final) "
              "coordinates and periodic unit cell dimensions from energy "
              "minimization runs, plus final velocities from molecular "
              "dynamics runs\n", 71, 0, stdout);
  PrintVADesc(2, "-r2", 5, " ", 2, "(output) checkpoint file; only used when "
              "thermodynamic integration is active\n", 71, 0, stdout);
}

//-----------------------------------------------------------------------------
// PrintInputFormat: helpful documentation on the format of mdgx input files.
//-----------------------------------------------------------------------------
void PrintInputFormat()
{
  PrintSplash(stdout);
  PrintParagraph("The typical mdgx input file is designed to look very much "
                 "like a typical SANDER input file.  However, there are some "
                 "key differences implemented to make the mdgx input file "
                 "format more flexible and the control data more intuitive.",
                 79, NULL, stdout);
  PrintParagraph("The only significant restriction introduced to the mdgx "
                 "input file format is that different segments of the input "
                 "file must begin with the namelist identifier (e.g. &cntrl, "
                 "&ewald) on its own line, and be terminated with the "
                 "identifier &end, again on its own line.", 79, NULL, stdout);
  printf("Here is an example input file:\n\n"
         "&files\n"
         "  -p       Tests/ions.top\n"
         "  -c       Tests/ions.min\n"
         "  -rst     Tests/ionsMDGX\n"
         "  -rstsuff .rst\n"
         "&end\n\n"
         "&cntrl\n"
         "  DoRATTLE = 1,   LJ14Fac = 2.0,   Elec14Fac = 1.2,\n"
         "  ElecCut = 9.0,  vdw_cutoff = 15.0,\n"
         "  dt = 0.001,   nstlim 500000,  nfistep = 1000,\n"
         "  ntpr = 1000,   ntwr 1000,  ntwf = 1000,\n"
         "  Temperature = 0.0,\n"
         "  SplnSpc = 0.015625,\n"
         "&end\n\n"
         "&ewald\n"
         "  ordr1 = 4,\n"
         "  ordr2 = 4,\n"
         "  ordr3 = 4,\n"
         "  nfft1 = 64,\n"
         "  nfft2 = 64,\n"
         "  nfft3 = 64;\n"
         "&end\n\n");
  PrintParagraph("Note the presence of the familiar <namelist identifier> "
                 "<arguments> <&end> format, carried over from SANDER. "
                 "However, mdgx includes new namelists such as &files (which "
                 "allows the bulk of the command line information to be given "
                 "in the input file) and new variables such as file suffixes "
                 "(for specifying multiple files in a single run).  There are "
                 "also aliases for the familiar (though sometimes "
                 "unintelligible) SANDER namelist variables, to provide a "
                 "format that will accept SANDER input while also permitting "
                 "users to write less cryptic command files.", 79, NULL,
                 stdout);
  PrintParagraph("Another important change to the mdgx file format is that = "
                 "signs are no longer required between a variable name and "
                 "the desired value.  In fact, commas are no longer required "
                 "either, though they are useful for separating different "
                 "attributes in each namelist.", 79, NULL, stdout);
}

//-----------------------------------------------------------------------------
// PrintFilesNamelistVariables: this function provides documentation on 
//                              &files namelist variables, with whatever
//                              aliases may be available.               
//-----------------------------------------------------------------------------
void PrintFilesNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("The files namelist is a means of specifying many input and "
                 "output files, which would otherwise be given on the command "
                 "line, in the input file along with other namelists. "
                 "However, if files are specified on the command line, they "
                 "will take precedence over any data in the &files namelist.",
                 79, NULL, stdout);
  PrintParagraph("The &files namelist can also be used to specify suffixes "
                 "for output files, in the event that multiple sequential "
                 "output files are to be generated in a single run.  This can "
                 "be useful for runs on managed resources that do not allow "
                 "multiple initializations of a program in a single job "
                 "submission, or for managing many segments of a very long "
                 "trajectory.  The suffixes only come into use if the "
                 "variable nfistep (alias FileStepCount) is set in the &cntrl "
                 "namelist.  Otherwise, only one file of each output type "
                 "will be written and only the base file names will be used. "
                 "File overwriting has slightly different behavior when "
                 "multiple sequential output files are specified.  If "
                 "overwriting is activated in such a case, mdgx will search "
                 "for the lowest file number such that a complete state "
                 "information file (specified by the -o variable) or a "
                 "complete restart file (specified by the -r variable) are "
                 "unavailable.  The dynamics will begin, or resume, at that "
                 "point.", 79, NULL, stdout);
  printf("  Name      Alias      Description\n"
         " ------ ------------- ---------------------------------------------"
         "-------------\n");
  PrintVADesc(1, "-p", 7, "Topology", 14, "(input) molecular topology file "
              "(AMBER prmtop format)", 58, 2, stdout);
  PrintVADesc(1, "-p2", 7, "Topology2", 14, "(input) molecular topology file; "
              "if thermodynamic integration is active, this topology "
              "describes the final state while the original topology "
              "describes the initial state", 58, 2, stdout);
  PrintVADesc(1, "-xpt", 7, "EPRules", 14, "(input) extra points rule file "
              "directing mdgx to add extra points to the topology at run time",
              58, 2, stdout);
  PrintVADesc(1, "-xpt2", 7, "EPRules2", 14, "(input) extra points rule file "
              "for the topology specified by the -p2 flag", 58, 2, stdout);
  PrintVADesc(1, "-c", 7, "StartCrd", 14, "(input) initial coordinates (and, "
              "optionally, velocities) and periodic unit cell size (mdgx does "
              "not run with non-periodic unit cells)", 58, 2, stdout);
  PrintVADesc(1, "-d", 7, "ForceDump", 14, "(output) a comprehensive force "
              "and energy report; this file is analogous to forcedump.dat as "
              "produced by SANDER", 58, 2, stdout);
  PrintVADesc(1, "-rrp", 7, "ResReport", 14, "(output) a complete description "
              "of the various residue types in the system; does not include "
              "information on connections between residues", 58, 2, stdout);
  PrintVADesc(1, "-o", 7, "OutputBase", 14, "(output) human-readable state "
              "information and diagnostics", 58, 2, stdout);
  PrintVADesc(1, "-e", 7, "EnergyBase", 14, "(output) energy data over "
              "trajectory", 58, 2, stdout);
  PrintVADesc(1, "-x", 7, "CrdTrajBase", 14, "(output) coordinate trajectory "
              "file; coordinate sets saved at specified time intervals", 58,
              2, stdout);
  PrintVADesc(1, "-v", 7, "VelTrajBase", 14, "(output) velocity trajectory "
              "file; velocity sets saved at specified time intervals", 58, 2,
              stdout);
  PrintVADesc(1, "-f", 7, "FrcTrajBase", 14, "(output) force trajectory "
              "file; forces on all particles saved at specified time"
              "intervals", 58, 2, stdout);
  PrintVADesc(1, "-r", 7, "RestartBase", 14, "(output) checkpoint (and final) "
              "coordinates and periodic unit cell dimensions from energy "
              "minimization runs, plus final velocities from molecular"
              "dynamics runs", 58, 2, stdout);
  PrintVADesc(1, "-osf", 7, "OutputSuff", 14, "output state information data "
              "file suffix (if this or other suffixes are not specified, the "
              "base name is taken to be the complete file name)", 58, 2,
              stdout);
  PrintVADesc(1, "-esf", 7, "EnergySuff", 14, "Energy data file suffix", 58,
              2, stdout);
  PrintVADesc(1, "-xsf", 7, "CrdTrajSuff", 14, "Coordinate trajectory file "
              "suffix", 58, 2, stdout);
  PrintVADesc(1, "-vsf", 7, "VelTrajSuff", 14, "Velocity trajectory file "
              "suffix", 58, 2, stdout);
  PrintVADesc(1, "-fsf", 7, "FrcTrajSuff", 14, "Force trajectory file suffix",
              58, 2, stdout);
  PrintVADesc(1, "-rsf", 7, "RestartSuff", 14, "Restart file suffix\n", 58, 2,
              stdout);
}

//-----------------------------------------------------------------------------
// PrintCntrlNamelistVariables: this function provides documentation on &cntrl
//                              namelist variables, with whatever aliases may
//                              be available.               
//-----------------------------------------------------------------------------
void PrintCntrlNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("The &cntrl namelist is required for any SANDER run, and is "
                 "present with little modification in mdgx.  Most SANDER "
                 "input files can therefore be read directly by mdgx or "
                 "adapted without much effort.  However, there are some "
                 "variables that have either been replaced or discarded.  The "
                 "following list describes all variables available in the "
                 "mdgx &cntrl namelist.", 79, NULL, stdout);
  PrintParagraph("Molecular Dynamics Variables:", 79, NULL, stdout);
  printf("    Name        Alias      Description\n"
         " ---------- ------------- -----------------------------------------"
         "-------------\n");
  PrintVADesc(1, "imin", 11, "RunMode", 14, "The run mode (0 = molecular "
              "dynamics, 1 = energy minimization, 2 = force computation for "
              "the input coordinates)", 54, 2, stdout);
  PrintVADesc(1, "irest", 11, "RestartMD", 14, "Set to 1 to request that "
              "molecular dynamics be restarted using velocities given in the "
              "input coordiantes file; set to 0 to request that initial " 
              "velocities be assigned to random values in a Maxwell "
              "distribution", 54, 2, stdout);
  PrintVADesc(1, "ioutfm", 11, "CoordFormat", 14, "The format of trajectory "
              "output coordinates, 0 (default) being ascii format with three "
              "decimal places and 1 being binary NetCDF format, specifying "
              "all coordinates to single floating point precision", 54, 2,
              stdout);
  PrintVADesc(1, "nstlim", 11, "StepCount", 14, "Number of MD steps to be "
              "performed", 54, 2, stdout);
  PrintVADesc(1, "nfistep", 11, "FileStepCount", 14, "Length of each segment "
              "of the trajectory; the number of segments is nstlim / nfistep",
              54, 2, stdout);
  PrintVADesc(1, "nscm", 11, "ZeroMomentum", 14, "Removal of translational "
              "and rotational center-of-mass (COM) motion at regular "
              "intervals", 54, 2, stdout);
  PrintVADesc(1, "t", 11, "StartTime", 14, "Time at the start of the "
              "simulation (picoseconds); this parameter is for user reference "
              "and otherwise only affects the time entered in the diagnostics "
              "files", 54, 2, stdout);
  PrintVADesc(1, "dt", 11, "TimeStep", 14, "Time step (picoseconds); the "
              "recommended MAXIMUM is 0.0025 if bonds to hydrogen are "
              "constrained, or 0.001 otherwise", 54, 2, stdout);
  PrintVADesc(1, "temp0", 11, "Temperature", 14, "Target temperature for a "
              "constant temperature simulation.  The default is 298.0.  This "
              "value can be used to initialize velocities if a specific "
              "initial temperature is not set.", 54, 2, stdout);
  PrintVADesc(1, "tempi", 11, "StartTemp", 14, "Initial temperature for a "
              "simulation.  The default is -100.0, which commands mdgx to use "
              "temp0 as the initial temperature for things such as velocity "
              "initialization.  If a positive value of tempi is specified, "
              "tempi will be used to initialize velocities.", 54, 2, stdout);
  PrintVADesc(1, "ntt", 11, "Thermostat", 14, "Thermostat specification.  "
              "Numerical values of 0 (default, no thermostat), 1 (Berendsen "
              "thermostat), 2 (Andersen thermostat), 3 (Langevin integrator) "
              "and 4 (Nose-Hoover thermostat) are supported.", 54, 2, stdout);
  PrintVADesc(1, "ntp", 11, "CoordRscl", 14, "Coordinate rescaling "
              "specification.  Numerical values of 0 (default, no rescaling, "
              "constant volume), 1 (isotropic rescaling), and 2 "
              "(anisotropic rescaling) are supported.", 54, 2, stdout);
  PrintVADesc(1, "barostat", 11, "Barostat", 14, "Barostat style.  Numerical "
              "values of 1 (Berendsen, default), and 2 (Monte-Carlo) are "
              "supported.", 54, 2, stdout);
  PrintVADesc(1, "pres0", 11, "Pressure", 14, "Target pressure for a constant "
              "pressure simulation.  The default is 1.0 bar.", 54, 2, stdout);
  PrintVADesc(1, "tautp", 11, "BerendsenTC", 14, "Time constant for "
              "Berendsen temperature coupling.  Default value is 0.4 ps.", 54,
              2, stdout);
  PrintVADesc(1, "gamma_ln", 11, "LangevinFreq", 14, "Langevin collision "
              "frequency in events / ps, when ntt=3.  Default is 0.", 54, 2,
              stdout);
  PrintVADesc(1, "taup", 11, "BerendsenPC", 14, "Compressibility constant for "
              "Berendsen pressure coupling.  Default value is 4.4e-5 / bar.",
              54, 2, stdout);
  PrintVADesc(1, "tauthv", 11, "HooverTC", 14, "Time constant for Hoover "
              "temperature coupling.  Default value is 1.0 ps.", 54, 2,
              stdout);
  PrintVADesc(1, "tauphv", 11, "HooverPC", 14, "Compressibility constant for "
              "Hoover pressure coupling.  Default value is 1.0 / bar.", 54, 2,
              stdout);
  PrintVADesc(1, "mccomp", 11, "MCBarostatPC", 14, "Coordinate rescaling "
              "factor for isotropic Monte-Carlo pressure coupling.  Default "
              "is 2.0e-3, to rescale the volume by +/- 1/10th of 1%.", 54, 2,
              stdout);
  PrintVADesc(1, "mccompx", 11, "MCBarostatPCX", 14, "Coordinate rescaling "
              "factor for anisotropic Monte-Carlo pressure coupling in the X "
              "direction.  Default is 2.0e-3, to rescale the volume by +/- "
              "1/10th of 1%.", 54, 2, stdout);
  PrintVADesc(1, "mccompy", 11, "MCBarostatPCY", 14, "Coordinate rescaling "
              "factor for anisotropic Monte-Carlo pressure coupling in the Y "
              "direction.  Default is to perform only isotropic rescaling "
              "based on the factor for the X direction.", 54, 2, stdout);
  PrintVADesc(1, "mccompz", 11, "MCBarostatPCZ", 14, "Coordinate rescaling "
              "factor for anisotropic Monte-Carlo pressure coupling in the Z "
              "direction.  Default is to perform only isotropic rescaling "
              "based on the factor for the X direction.", 54, 2, stdout);
  PrintVADesc(1, "mcbfrq", 11, "MCBarostatFrq", 14, "Step frequency for "
              "applying the Monte-Carlo barostat, if this barostat is "
              "activated.  Default 100.", 54, 2, stdout);
  PrintVADesc(1, "vrand", 11, "RandomReset", 14, "Time constant for Andersen "
              "temperature coupling.  This is specified as an integer "
              "denoting the number of time steps, and as such is related to "
              "the time step size dt.  Default value is 1000.", 54, 2, stdout);
  PrintVADesc(1, "ig", 11, "RandomSeed", 14, "The random seed for velocity "
              "initialization and thermostats which may require it", 54, 2,
              stdout);
  PrintVADesc(1, "es_cutoff", 11, "ElecCut", 14, "The electrostatic direct "
              "space cutoff (Angstroms)", 54, 2, stdout);
  PrintVADesc(1, "vdw_cutoff", 11, "VdwCut", 14, "The van-der Waals (direct "
              "space) cutoff (Angstroms)", 54, 2, stdout);
  PrintVADesc(1, "cut", 11, "DirectCutoff", 14, "The general (van-der Waals "
              "and electrostatic) direct space cutoff (Angstroms).  This "
              "value, if indicated, will override es_cutoff and vdw_cutoff.",
              54, 2, stdout);
  PrintVADesc(1, "rigidbond", 11, "DoRATTLE", 14, "Set to 1 to activate "
              "RATTLE bond length constraints", 54, 2, stdout);
  PrintVADesc(1, "rigidwat", 11, "DoSETTLE", 14, "Set to 1 to activate "
              "SETTLE water geometry constraints", 54, 2, stdout);
  PrintVADesc(1, "tol", 11, "RattleTol", 14, "Tolerance for RATTLE bond "
              "length constraints", 54, 2, stdout);
  PrintVADesc(1, "scee", 11, "Elec14Fac", 14, "The electrostatic 1-4 "
              "interaction scaling factor", 54, 2, stdout);
  PrintVADesc(1, "scnb", 11, "Vdw14Fac", 14, "The van-der Waals 1-4 "
              "interaction scaling factor", 54, 2, stdout);
  printf("\n");
  PrintParagraph("Thermodynamic integration control variables:", 79, NULL,
                 stdout);
  printf("  Name        Alias        Description\n"
         " ------ ----------------- -----------------------------------------"
         "-------------\n");
  PrintVADesc(1, "icfe", 11, "RunTI", 14, "Flag to activate thermodynamic "
              "integration.  Default 0 (no TI); set to 1 for active.", 54, 2,
              stdout);
  PrintVADesc(1, "clambda", 11, "MixFactor", 14, "The mixing parameter, "
              "(1-L)^k of state 1 and 1-(1-L)^k of state 2", 54, 2, stdout);
  PrintVADesc(1, "klambda", 11, "MixOrder", 14, "The order of the mixing "
              "parameter, (1-L)^k of state 1 and 1-(1-L)^k of state 2", 54, 2,
              stdout);
  PrintVADesc(1, "nsynch", 11, "SynchTI", 14, "Explicit synchronization of "
              "trajectory coordinates will occur every nsynch steps (default "
              "1000).  Set nsynch to 0 to disable this feature.", 54, 2,
              stdout);
  printf("\n");
  PrintParagraph("Output control variables:", 79, NULL, stdout);
  printf("  Name        Alias        Description\n"
         " ------ ----------------- -----------------------------------------"
         "-------------\n");
  PrintVADesc(1, "ntpr", 7, "WriteDiagnostics", 18, "Diagnostics and state "
              "information output frequency", 54, 2, stdout);
  PrintVADesc(1, "ntwx", 7, "WriteCoordinates", 18, "Trajectory coordinates "
              "will be written at this frequency", 54, 2, stdout);
  PrintVADesc(1, "ntwv", 7, "WriteCoordinates", 18, "Trajectory velocities "
              "will be written at this frequency", 54, 2, stdout);
  PrintVADesc(1, "ntwf", 7, "WriteCoordinates", 18, "Trajectory forces will "
              "be written at this frequency", 54, 2, stdout);
  PrintVADesc(1, "ntwr", 7, "WriteCoordinates", 18, "Trajectory forces will "
              "be written at this frequency", 54, 2, stdout);
  PrintVADesc(1, "tchk", 7, "TopologyCheck", 18, "Active by default, set to 0 "
              "to turn off topology and conformation checking at the start of "
              "each run segment", 54, 2, stdout);
}

//-----------------------------------------------------------------------------
// PrintEwaldNamelistVariables: this function provides documentation on &ewald
//                              namelist variables, with whatever aliases may
//                              be available.               
//-----------------------------------------------------------------------------
void PrintEwaldNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("All variables in the &ewald namelist have default values, "
                 "so use of this namelist is optional, but optimization of "
                 "the parameters in this section can be very helpful for "
                 "performing the most efficient molecular simulations.  The "
                 "SANDER and PMEMD programs warn users not to modify "
                 "variables in this section without significant experience; "
                 "what is most important is a clear understanding of all the "
                 "variables and how they will affect the accuracy of "
                 "computed forces.  The state information / diagnostics "
                 "output file will print verbose explanations of the "
                 "consequences of any variables that are changed, so "
                 "modification of these variables, with careful reading of "
                 "the output and checks on the accuracy of computed forces, "
                 "should be safe.", 79, NULL, stdout);
  PrintParagraph("Smooth Particle Mesh Ewald (SPME) Variables:", 79, NULL,
                 stdout);
  printf("    Name        Alias      Description\n"
         " ---------- ------------- ------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "dsum_tol", 11, "DSumTol", 14, "The direct sum tolerance; "
              "at the direct space cutoff es_cutoff (see the &cntrl namelist)"
              ", the ratio between the interaction energy of two point "
              "charges and two Gaussian smeared charges of the same magnitude "
              "will differ from 1 by dsum_tol; this variable controls the "
              "accuracy of the electrostatic direct space sum", 54, 2, stdout);
  PrintVADesc(1, "sigma", 11, "Sigma", 14, "The width (root mean squared "
              "deviation) of spherical Gaussians used to smooth out the "
              "charge density; sigma = ew_coeff / 2; if this value is "
              "supplied by the user it will supercede any entry for ew_coeff "
              "and dsum_tol will be calculated from sigma; otherwise, sigma "
              "and ew_coeff will be calculated from dsum_tol", 54, 2, stdout);
  PrintVADesc(1, "ew_coeff", 11, "EwaldCof", 14, "The Ewald coefficient; this "
              "quantity has a name only because it appears frequently in "
              "the Smooth Particle Mesh Ewald formulas; physically it makes "
              "more sense to consider the Gaussian charge width sigma", 54, 2,
              stdout);
  PrintVADesc(1, "eetbdns", 11, "SplnSpc", 14, "The discretization of the "
              "erfc(beta*r)/r force and energy spline computation tables", 54,
              2, stdout);
  PrintVADesc(1, "rho", 11, "MaxDensity", 14, "The maximum expected density "
              "of the system, g/mL.  Default 2.0, increase to raise the "
              "maximum storage in the direct-space decomposition cell grid.",
              54, 2, stdout);
  PrintVADesc(1, "nfft[123]", 11, "MeshDim[XYZ]", 14, "The number of mesh "
              "points in the X, Y, or Z dimensions, respectively (if the unit "
              "cell / simulaton box is orthorhombic), or otherwise the number "
              "of mesh points along the 1st, 2nd, and 3rd unit cell vectors",
              54, 2, stdout);
  PrintVADesc(1, "ordr[123]", 11, "Order[XYZ]", 14, "The order of particle "
              "<--> mesh interpolation along the 1st, 2nd, and 3rd unit cell "
              "vectors", 54, 2, stdout);
  PrintVADesc(1, "order", 11, "Order", 14, "Sets the interpolation order "
              "isotropically along all unit cell vectors to the specified "
              "value (supercedes ordr[123])", 54, 2, stdout);
  printf("\n");
  PrintParagraph("Long-ranged van-der Waals control parameters:", 79, NULL,
                 stdout);
  printf("    Name        Alias      Description\n"
         " ---------- ------------- ------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "vdwmeth", 11, "vdwMethod", 14, "The method for computing "
              "long-ranged van der Waals interactions.  The default of 1 "
              "implies the inclusion of a homogeneity assumption in the "
              "long-ranged component of the van-der Waals interactions; this "
              "correction changes the computed energy and pressure, but not "
              "forces, and would therefore not affect dynamics in a "
              "simulation at constant volume. A value of zero removes any "
              "such correction and make the van-der Waals energy depend "
              "solely on the pairwise interactions.", 54, 2, stdout);
  printf("\n");
  PrintParagraph("Multi-Level Ewald (MLE) Variables:", 79, NULL, stdout);
  printf("    Name        Alias      Description\n"
         " ---------- ------------- ------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "nlev", 11, "EwaldLevels", 14, "The number of levels in the "
              "multi-level mesh hierarchy (standard Smooth Particle Mesh "
              "Ewald has one level); setting this variable to any value "
              "greater than one will activate MLE (maximum 4)", 54, 2, stdout);
  PrintVADesc(1, "lpad[123]", 11, "Padding[123]", 14, "The number of layers "
              "of \"padding\" for meshes at levels 1, 2, and 3; the "
              "reciprocal space pair potential will be represented exactly on "
              "mesh level 1 up to (and including) lpad1 grid points from the "
              "source, and will be represented to different degrees of "
              "resolution (see cfac[234]) on higher mesh levels, up to "
              "(and including) lpad1 + lpad2 or lpad1 + lpad2 + lpad3 layers "
              "from the source (note that the highest mesh level is not "
              "padded as it involves only one convolution over the entire "
              "simulation box); higher values of lpad[123] will produce more "
              "accurate results (see also ggordr)", 54, 2, stdout);
  PrintVADesc(1, "cfac[234]", 11, "Spread[234]", 14, "The coarsening factor "
              "for higher mesh levels; by definition, cfac1 is 1; generally, "
              "it is advisable to set cfac2 to 1.5-2.0, and to set cfac3 or "
              "cfac4 (if even higher mesh levels are in use) to increasingly "
              "large numbers; although cfac is a real number, but it must be "
              "specified such that the mesh size is an integer multiple of "
              "cfac in every dimension", 54, 2, stdout);
  PrintVADesc(1, "ggordr", 11, "GridOrder", 14, "The order of grid <--> grid "
              "B-Spline interpolation; higher orders of interpolation will "
              "produce more accurate results for given values of lpad[123]",
              54, 2, stdout);
}

//-----------------------------------------------------------------------------
// PrintForceNamelistVariables: the force namelist was added to support 
//                              customization of detailed force report files.
//-----------------------------------------------------------------------------
void PrintForceNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("By default, all available information is printed to the "
                 "force report file; all of the variables listed below are "
                 "set to 1 by default.  Most of this information is not "
                 "needed, however, so much of the output can be suppressed by "
                 "setting these variables to 0.", 79, NULL, stdout);
  printf("   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "var", 11, "VarName", 14, "Information is dumped into a "
              "Matlab-readable force report file; this specifies the name of "
              "the variable that Matlab will use to store the force "
              "information if it reads the report; it may be useful to read "
              "multiple reports at once, and compare the results using Matlab",
              54, 2, stdout);
  PrintVADesc(1, "dumpcrd", 11, "DumpCoord", 14, "Flag to dump the "
              "coordinates of the system", 54, 2, stdout);
  PrintVADesc(1, "dumpbond", 11, "DumpBond", 14, "Flag to dump the forces due "
              "to bonded (1-2) interactions", 54, 2, stdout);
  PrintVADesc(1, "dumpangl", 11, "DumpAngl", 14, "Flag to dump the forces due "
              "to angle interactions", 54, 2, stdout);
  PrintVADesc(1, "dumpdihe", 11, "DumpDihe", 14, "Flag to dump the dihedral "
              "forces", 54, 2, stdout);
  PrintVADesc(1, "dumpdelec", 11, "DumpDElec", 14, "Flag to dump direct sum "
              "electrostatic forces", 54, 2, stdout);
  PrintVADesc(1, "dumprelec", 11, "DumpRElec", 14, "Flag to dump reciprocal "
              "sum electrostatic forces", 54, 2, stdout);
  PrintVADesc(1, "dumpvdw", 11, "DumpVdw", 14, "Flag to dump van-der Waals "
              "forces", 54, 2, stdout);
  PrintVADesc(1, "dumpall", 11, "DumpAll", 14, "Flag to dump total (summed) "
              "forces", 54, 2, stdout);
}

//-----------------------------------------------------------------------------
// PrintFitqNamelistVariables: this function provides documentation on &fitq
//                             namelist varialbes, with whatever aliases may
//                             be available.                
//-----------------------------------------------------------------------------
void PrintFitqNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("The basic concept of fitting charges to reproduce the "
                 "electrostatic potential of a molecule, by finding the "
                 "solution with least squared error, in the presence of "
                 "restraints, is carried over from the original Kollmann "
                 "RESP.  This namelist provides tools for fitting charges in "
                 "a comprehensible manner, with exceptional user control over "
                 "the range of the fitting data.  Because mdgx is unique "
                 "among the current AMBER molecular dynamics engines for its "
                 "ability to use certain types of extra points (virtual "
                 "sites), this namelist is also useful for fitting new charge "
                 "models to accelerate parameter development in mdgx.", 79,
                 NULL, stdout);
  PrintParagraph("Many of these variables can be specified more than once, "
                 "and all instances will accumulate in the result.", 79, NULL,
                 stdout);
  printf("   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "resp", 11, "RespPhi", 14, "File names and numerical weight "
              "of an electrostatic potential to use in fitting.  The format "
              "is <string1> <string2> <real1>, where string1 is a Gaussian "
              "cubegen format file specifying the electrostatic potential "
              "and molecular coordinates and Z-numbers appropriate to the "
              "topology specified by string2 and real1 is the numerical "
              "weight to be given to this conformation in the fit.  This "
              "variable may be specified for all molecular conformations "
              "and systems to be used in the fit.\n", 54, 2, stdout);
  PrintVADesc(1, "ipolq", 11, "IPolQPhi", 14, "File names and numerical "
              "weight of a pair of electrostatic potentials to use in IPolQ "
              "fitting.  The format is <string1> <string2> <string3> <real1>, "
              "where string1 and string2 are Gaussian cubegen format files "
              "relating to the system in vacuum and in a condensed-phase "
              "environment (see the &ipolq namelist).  Note that the "
              "molecular coordinates in both cubegen files must match.  As in "
              "the resp variable format, string3 is the topology, and real1 "
              "is the numerical weight of this conformation.  This extended "
              "REsP procedure supports development of fixed-charge force "
              "fields if one posits that the correct charges of a "
              "non-polarizable model would sit halfway between the charges of "
              "a fully polarized molecule in some solvent reaction field and "
              "the charges of an unpolarized molecule in the gas phase.  This "
              "variable may be specified for all molecular conformations "
              "and systems to be used in the fit.\n", 54, 2, stdout);
  PrintVADesc(1, "conf", 11, "ConfFile", 14, "If specified, mdgx, will output "
              "the first molecular conformation, complete with any added "
              "virtual sites, in PDB format for inspection.  This is useful "
              "for understanding exactly what model is being fitted.", 54, 2,
              stdout);
  PrintVADesc(1, "eprules", 11, "EPRules", 14, "If specified, mdgx, will "
              "output all fitted charges in the form of a Virtual Sites rule "
              "file, which can be given as input to subsequent simulations to "
              "modify the original prmtp and apply the fitted charge model.",
              54, 2, stdout);
  PrintVADesc(1, "hist", 11, "HistFile", 14, "If specified, mdgx will print "
              "a histogram of the distance of all points from the nearest "
              "atom of the molecule.", 54, 2, stdout);
  PrintVADesc(1, "ptrecord", 11, "PtRecFile", 14, "If specified, mdgx will "
              "print a Matlab file of all pointsused in fitting, structure "
              "by structure, with locations, distances to the nearest atom, "
              "and intensities.", 54, 2, stdout);
  PrintVADesc(1, "amblib", 11, "LibOutput", 14, "Name of an Amber-format "
	      "library to write, containing the fitted charges.  This "
	      "requires a template library file which also contains the "
	      "residues (units) that were fitted.", 54, 2, stdout);
  PrintVADesc(1, "ambtmp", 11, "LibTemplate", 14, "Name of the template "
	      "Amber-format library file.  Because mdgx does not understand "
	      "the Amber library format in depth, it must be able to find "
	      "units describing the residues for partial charges that it has "
	      "just fitted.  Only units that match the residue and all atom "
	      "names will be edited.", 54, 2, stdout);
  PrintVADesc(1, "unitmatch", 11, "UnitMatching", 14, "Level of strictness "
	      "by which residue and unit names will be matched in the "
	      "Amber-format output.  Setting this to 'permissive' or "
	      "'PERMISSIVE' will match unit names that contain the full "
	      "residue name (i.e. CARG contains ARG).  Any other setting "
	      "will require complete correspondence.", 54, 2, stdout);
  PrintVADesc(1, "minq", 11, "MinimizeQ", 14, "Restrain the charges of a "
              "group of atoms to zero by the weight given in minqwt.  This "
              "variable may be specified many times for different groups.  "
              "The groups are specified in ambmask format.  In order for "
              "minq directives to have an effect, minqwt must be set to some "
              "non-zero value.", 54, 2, stdout);
  PrintVADesc(1, "equalq", 11, "EqualizeQ", 14, "Restrain the charges of a "
              "group of atoms to have the same values.  Groups are specified "
              "in ambmask format.  This variable may be repeatedly specified.",
              54, 2, stdout);
  PrintVADesc(1, "sumq", 11, "SumQ", 14, "Restrain the total charge of a "
              "group of atoms to have a particular value.  The group is "
              "specified in ambmask format, followed by a real number "
              "indicating the target total charge.  A stiff harmonic "
              "restraint penalty is applied.  This variable may be repeatedly "
              "specified.", 54, 2, stdout);
  PrintVADesc(1, "fixq", 11, "FixQ", 14, "Fix the charge of one or more atoms "
	      "to have a particular value.  The atoms are specified in "
	      "ambmask format, followed by a real number indicating the "
	      "target solvated charge, and optionally a second real number "
	      "indicating the target vacuum charge.  A stiff harmonic "
	      "restraint penalty is applied.  This variable may be repeatedly "
	      "specified.", 54, 2, stdout);
  PrintVADesc(1, "snap", 11, "MaxSnap", 14, "Charges will be adjusted to "
              "eliminate miniscule deviations in the total charge of "
              "specified groups.  This variable sets the maximum adjustment "
              "that can be made, in increments of 1.0e-5 proton charges.", 54,
              2, stdout);
  PrintVADesc(1, "tether", 11, "Tether", 14, "Restrain fitted charges not "
              "zero (small values), but to their values in the input "
              "topologies, to use the original force field as a guide.", 54, 2,
              stdout);
  PrintVADesc(1, "tetherwt", 11, "TetherWeight", 14, "Weight used for "
              "harmonic restraint of fitted charges to values given in the "
              "original force field (in standard REsP), or for restraining "
              "IPolQ charges to their vacuum counterparts (in IPolQ "
              "fitting, when tethering between the charge sets is expected).  "
              "Default zero in both contexts.", 54, 2, stdout);
  PrintVADesc(1, "minqwt", 11, "MinQWeight", 14, "Weight used for restraining "
              "values of charges to zero; as more and more fitting data is "
              "included (either through a higher sampling density of the "
              "electrostatic potential due to each molecular conformation or "
              "additional molecular conformations) higher values of minqwt "
              "may be needed to keep the fitted charges small.  However, with "
              "more data the need to restrain charges may diminish as well.",
              54, 2, stdout);
  PrintVADesc(1, "nfpt", 11, "FitPoints", 14, "The number of fitting points "
              "to select from each electrostatic potential grid.  The points "
              "nearest the molecule, which satisfy the limits set by the "
              "solvent probe and point-to-point distances as defined below, "
              "will be selected for the fit.  Default 1000.", 54, 2, stdout);
  PrintVADesc(1, "psig", 11, "ProbeSig", 14, "The Lennard-Jones sigma "
              "parameter of the solvent probe.  Default 3.16435 (TIP4P "
              "oxygen).", 54, 2, stdout);
  PrintVADesc(1, "peps", 11, "ProbeEps", 14, "The Lennard-Jones parameter of "
              "the solvent probe.  Default 0.16275 (TIP4P oxygen).", 54, 2,
              stdout);
  PrintVADesc(1, "parm", 11, "ProbeArm", 14, "The probe arm; points on the "
              "electrostatic potential grid that would be inaccessible to the "
              "solvent probe may still be included in the fit if they are "
              "within the probe arm's reach.", 54, 2, stdout);
  PrintVADesc(1, "pnrg", 11, "StericLimit", 14, "The maximum Lennard-Jones "
              "energy of the solvent probe at which a point will qualify for "
              "inclusion in the fit.  Default 3.0 kcal/mol.", 54, 2, stdout);
  PrintVADesc(1, "maxmem", 11, "MaxMemory", 14, "The amount of memory "
              "available to this job.  The test is actually the amount "
              "allocated to the fitting and testing matrices, which is the "
              "number of fitting data points times the number of independent "
              "charges being fitted times 16 bytes.  The actual memory usage "
              "will be slightly higher, but the matrices comprise the bulk of "
              "it.\n", 54, 2, stdout);
  PrintVADesc(1, "flim", 11, "Proximity", 14, "The minimum proximity of any "
              "two points to be included in the fit.  Default 0.4A.", 54, 2,
              stdout);
  PrintVADesc(1, "racc", 11, "AcceptAll", 14, "The maximum proximity of any "
              "points from at least one atom of the molecule that will "
              "guarantee its inclusion in the fit, provided that it is "
              "selected within the first nfpt points.  Default 3.0A.", 54, 2,
              stdout);
  PrintVADesc(1, "rmax", 11, "AcceptMax", 14, "The maximum proximity of any "
              "points from at least one atom of the molecule that will "
              "afford a possibility of its inclusion in the fit, provided "
              "that it is selected within the first nfpt points.  Default "
              "6.0A.", 54, 2, stdout);
  PrintVADesc(1, "hbin", 11, "HistogramBin", 14, "If hist is specified, mdgx "
              "will print a histogram reporting the number of fitting points "
              "falling within any particular distance of some atom of the "
              "molecule.  This parameter controls the discretization of the "
              "histogram.", 54, 2, stdout);
  PrintVADesc(1, "verbose", 11, "Verbose", 14, "Print information relating to "
              "progress on the fitting run.  Default is to print such data.  "
              "Set to zero to suppress output.", 54, 2, stdout);
}

//-----------------------------------------------------------------------------
// PrintParamNamelistVariables: the param namelist supports development of
//                              bonded terms, particularly torsion potentials,
//                              with consideration to each parameter possibly
//                              having multiple roles in many different
//                              systems.              
//-----------------------------------------------------------------------------
void PrintParamNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("Bonded terms can be fitted by a linear least squares "
                 "approach, tuning the stiffnesses of various spring-like "
                 "interactions with quantum data as the guide.  The fitting "
                 "functions are fundamentally the same as those used by the "
                 "charge fitting module, but the &param namelist can "
                 "accommodate many more systems.  The basic requirements of "
                 "this module are an Amber parameter file (i.e. parm99.dat), "
                 "an optional frcmod file (i.e. frcmod.ff99SB), and a list of "
                 "system topologies and coordinates coupled to energies "
                 "obtained by a consistent quantum mechanics treatment.  mdgx "
                 "will correlate all adjustable bonded terms found in any of "
                 "the systems by referencing the parameter files, then "
                 "compute a fitting matrix with one adjustable term per "
                 "column.  The fitted parameters are printed to a formatted "
                 "Amber parameter file, while statistics from the run, "
                 "including accuracy in each system and a comprehensive "
                 "analysis of the fitting data, are printed to mdout.", 79,
                 NULL, stdout);
  printf(" System specification and units:\n"
         "   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "sys", 11, "System", 14, "A fitting data point. This keyword "
              "must be followed by three items: the name of a topology file, "
              "the name of a corresponding coordinate file, and the energy of "
              "this system in the stated conformation.", 54, 2, stdout);
  PrintVADesc(1, "printpts", 11, "PrintFitPoints", 14, "Flag to have mdgx "
              "reprint the ENTIRE input file, which for large parameter sets "
              "will consist mostly of the names and energies of fitting data "
              "files.  By default, this feature is suppressed.", 54, 2,
              stdout);
  PrintVADesc(1, "eunits", 11, "EnergyUnits", 14, "Units of the target energy "
              "values.  Default kcal/mol.  Acceptable values include Hartree/"
              "Atomic, kJ/kilojoules, and j/joules.  Case insensitive.", 54,
              2, stdout);

  printf("\n Input quality control settings:\n"
         "   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");  
  PrintVADesc(1, "elimsig", 11, "ElimOutliers", 14, "Flag to activate removal "
              "of molecular conformations whose energies are far outside the "
              "norm for other conformations of the same system.  Default 0 "
              "(do not remove outliers).", 54, 2, stdout);
  PrintVADesc(1, "esigtol", 11 ,"EOutlier", 14, "Tolerance for deviation from "
              "the mean energy value, specified as a function of the standard "
              "deviation for all conformations of the same system.  "
              "Conformations of a system which exceed this threshold will be "
              "reported if verbose is set to 1, and removed from "
              "consideration if elimsig is set to 1.  Default 5.0 sigmas.\n",
              54, 2, stdout);
  PrintVADesc(1, "fsigtol", 11 ,"FOutlier", 14, "Tolerance for deviation from "
              "the fitted energy, a post-hoc analysis of the fitted model.  "
              "High deviations from the resulting force field may indicate "
              "that a data point does not belong in the training set.\n", 54,
              2, stdout);
  PrintVADesc(1, "fdevfloor", 11 ,"OutlierMinDev", 14, "Minimum deviation, in "
              "kcal/mol, needed for a fitted energy to qualify as an outlier.",
              54, 2, stdout); 
  PrintVADesc(1, "ctol", 11, "ConfTol", 14, "Energetic tolerance for strained "
              "bonds or angles in the system.  If any bonded energy term "
              "exceeds this tolerance, mdgx will assume that the atoms have "
              "been entered in the wrong order according to the topology, and "
              "try to rearrange the atoms in order to relax the system.", 54,
              2, stdout);  

  printf("\n Marking parameters for optimization:\n"
         "   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "bonds", 11, "FitBonds", 14, "Requests a linear "
              "least-squares fit for bond stiffnesses in the system.", 54, 2,
              stdout);
  PrintVADesc(1, "angles", 11, "FitAngles", 14, "Requests a linear "
              "least-squares fit for angle stiffnesses in the system.", 54, 2,
              stdout);
  PrintVADesc(1, "torsions", 11, "FitTorsions", 14, "Requests a linear "
              "least-squares fit for torsion stiffnesses in the system.", 54,
              2, stdout);
  PrintVADesc(1, "cmaps", 11, "FitCmaps", 14, "Requests a linear "
              "least-squares fit for all CMAP surfaces in the fitting data.",
	      54, 2, stdout);
  PrintVADesc(1, "fitb", 11, "FitB", 14, "Request that a specific bond "
              "parameter be included in linear least-squares fitting.", 54, 2,
              stdout);
  PrintVADesc(1, "fita", 11, "FitA", 14, "Request that a specific angle "
              "parameter be included in linear least-squares fitting.", 54, 2,
              stdout);
  PrintVADesc(1, "fith", 11, "FitH", 14, "Request that a specific torsion "
              "parameter be included in linear least-squares fitting.", 54, 2,
              stdout);
  PrintVADesc(1, "bondeq", 11, "FitBondEq", 14, "Requests that bond "
              "equilibrium constants be fitted alongside their spring "
              "constants.", 54, 2, stdout);
  PrintVADesc(1, "angleq", 11, "FitAnglEq", 14, "Requests that angle "
              "equilibrium constants be fitted alongside their spring "
              "constants.", 54, 2, stdout);
  PrintVADesc(1, "cmapdens", 11, "CMAPDensity", 14, "Density of fitting "
	      "points for determining CMAP surfaces.  The surfaces will "
	      "always be printed ot the output parameter file with the same "
	      "density as the input parameter file, but bicubic interpolation "
	      "will determine the final result if cmapdens is not the same "
	      "as the input CMAP dimension.  If cmapdens is set to 8 but the "
	      "input force field CMAPs have dimension 24, then the output "
	      "will be interpolated as appropriate.  In practice, the higher "
	      "density CMAPs that come out of such a case will give the same "
	      "results as the low density maps would have.  However, fitting "
	      "with cmapdens set higher than the input force field will "
	      "result in some loss of information.", 54, 2, stdout);
  PrintVADesc(1, "lpost", 11, "BondBasisSep", 14, "Distance from the original "
              "bond equilibrium length to place the equilibrium values of "
              "either of two basis functions for fitting a new bond term.",
              54, 2, stdout);
  PrintVADesc(1, "thpost", 11, "AnglBasisSep", 14, "Distance from the "
              "original angle equilibrium length to place the equilibrium "
              "values of either of two basis functions for fitting a new "
              "bond term.", 54, 2, stdout);
  PrintVADesc(1, "fitscnb", 11, "FitLJ14", 14, "Requests a linear "
              "least-squares fit for Lennard-Jones 1:4 scaling factors.",
              54, 2, stdout);
  PrintVADesc(1, "fitscee", 11, "FitEE14", 14, "Requests a linear "
              "least-squares fit for electrostatic 1:4 scaling factors.",
              54, 2, stdout);

  printf("\n Guiding the fit with parameter restraints or modifications:\n"
         "   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "scnb", 11, "Vdw14Fac", 14, "Sets a universal 1:4 scaling "
              "factor for van-der Waals interactions.  Use this input to "
              "change the scaling on all systems simultaneously.", 54, 2,
              stdout);
  PrintVADesc(1, "scee", 11, "Elec14Fac", 14, "Sets a universal 1:4 scaling "
              "factor for electrostatic interactions.  Use this input to "
              "change the scaling on all systems simultaneously.", 54, 2,
              stdout);
  PrintVADesc(1, "brst", 11 ,"BondRest", 14, "General value for harmonic "
              "restraints on bond stiffness constants.", 54, 2, stdout);
  PrintVADesc(1, "arst", 11 ,"AngleRest", 14, "General value for harmonic "
              "restraints on angle stiffness constants.", 54, 2, stdout);
  PrintVADesc(1, "hrst", 11 ,"DihedralRest", 14, "General value for harmonic "
              "restraints on torsion stiffness constants.", 54, 2, stdout);
  PrintVADesc(1, "sbrst", 11 ,"RestrainB", 14, "Applies a specific restraint "
              "stiffness to the value of a fitted bond, equivalent to "
              "changing brst for that bond alone.  This command takes "
              "subdirectives of atom type names, plus 'Keq' and 'Leq' (each "
              "followed by a positive real number to denote the target "
              "stiffness and equilibrium bond length, respectively).  These "
              "subdirectives may be given in free format.", 54, 2, stdout);
  PrintVADesc(1, "sarst", 11 ,"RestrainA", 14, "Applies a specific restraint "
              "stiffness to the value of a fitted angle, equivalent to "
              "changing arst for that angle alone.  This command takes "
              "subdirectives similar to sbrst, although 'Leq' corresponds "
              "to the equilibrium angle rather than length.", 54, 2, stdout);
  PrintVADesc(1, "shrst", 11 ,"RestrainH", 14, "Applies a specific restraint "
              "stiffness to the value of a fitted torsion amplitude, "
              "equivalent to changing hrst for that torsion alone.  This "
              "command takes subdirectives (in free format) of atom types, "
              "'period' or 'per' followed by a real number for the "
              "periodicity, 'weight' or 'rwt' followed by a real number for "
              "the restraint strength (which scales just like hrst), and "
              "'target' or 'trg' followed by a real number if the target "
              "value of the restraint is non-zero.", 54, 2, stdout);
  PrintVADesc(1, "geom", 11, "Geometry", 14, "Applies a geometry restraint to "
              "the sum of two or more angles, i.e. 360 degrees.  This command "
              "takes subdirectives of 'angl', followed immediately by three "
              "angle atom types, and 'target' or 'trg' followed by a real, "
              "positive number.", 54, 2, stdout);
  PrintVADesc(1, "brstcpl", 11 ,"BondCoupling", 14, "General value for "
              "pricing 1A changes in the fitted equilibrium constant "
              "with kcal/mol-A^2 changes in bond spring constants.  The "
              "default is 5000.0, which penalizes 50 kcal/mol changes in the "
              "stiffness at the same rate as 0.01A changes in the "
              "equilibrium value.  Only relevant with FitBondEq = 1.", 54, 2,
              stdout);
  PrintVADesc(1, "arstcpl", 11 ,"AngleCoupling", 14, "General value for "
              "pricing 1-degree changes in the fitted equilibrium constant "
              "with kcal/mol-rad^2 changes in angle spring constants.  The "
              "default is 114.59, which penalizes 2 kcal/mol changes in the "
              "stiffness at the same rate as 1 degree changes in the "
              "equilibrium value.  Only relevant with FitAnglEq = 1.", 54, 2,
              stdout);
  PrintVADesc(1, "spectrum", 11, "Spectrum", 14, "Request that a particular "
              "bond, angle, or dihedral from among the adjustable parameters "
              "be sampled near various values along a specified range, or "
              "otherwise included as a reoptimizable variable while others "
              "are resampled.  This keyword invokes its own sort of namelist: "
              "the sub-directives can be given in any order, but they must "
              "all be given on the same line.  Words that are not explicitly "
              "sub-directives or values following them may be counted as atom "
              "types, so long as they have fewer than four characters, until "
              "four such types are catalogged.  Sub-directives include "
              "retain (parameters matching this request are reoptimized but "
              "not resampled), sample (parameters matching this request are "
              "resampled), order (followed by a value, 2 = bonds, 3 = angles, "
              "4 = torsions), min and max (followed by values, the resampling "
              "range limits), spc (the resampling discretization), and break "
              "(stop adding new sub-directives to this spectrum command).", 
              54, 2, stdout);

  printf("\n Reporting output (parameters, analysis, and runtime messages):\n"
         "   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "verbose", 11, "ShowProgress", 14, "Alert the user as to the "
              "progress of the fitting procedure.  Runs involving thousands "
              "of molecular conformations and hundreds of parameters can "
              "generally be completed in a few minutes.  Default is 1 (ON).  "
              "Set to zero to suppress output.", 54, 2, stdout);
  PrintVADesc(1, "repall", 11, "ReportAll", 14, "(Deprecated--replaced by the "
              "ParmOutput keyword).  Acceptable values include 0 (ParmOutput "
              "is 'frcmod'), 1 (ParmOutput is 'standard'), and 2 (ParmOutput "
              "is 'stdplus').", 54, 2, stdout);
  PrintVADesc(1, "parmout", 11, "ParmOutput", 14, "This keyword will control "
              "the output of parameters created during the fitting procedure, "
              "including those that were not adjusted by the fit but may have "
              "been part of the original force field.  The default is "
              "'standard': write all parameters from the original force "
              "field (including input frcmod files), with edits made by "
              "the fitting procedure, in a single new Amber parameter file.  "
              "If atom type branching has occured, 'standard' output will "
              "only print parameters for branched atom types that were "
              "involved in the data fitting.  However, it may be useful to "
              "retain the parameters from the original atom types for other "
              "instances in which the branchd atom types may be encountered.  "
              "For this purpose, 'stdplus' may be specified as the output "
              "style, directing all parameters applicable to branched atom "
              "types, including those copied from the original types, to be "
              "included in a new Amber parameter file.  For cases in which "
              "only the edited parameters are desired, an Amber frcmod file "
              "can be written by specifying 'frcmod' here.", 54, 2, stdout);
  PrintVADesc(1, "zeromm", 11, "FittedMMOnly", 14, "Flag to ignore any "
              "contributions from unfitted energy terms during the fit.  "
              "Useful for making force field adjustments.", 54, 2, stdout);
  PrintVADesc(1, "accrep", 11, "AccReport", 14, "Accuracy report on the fit.  "
              "Contains extensive analysis on the resulting parameters, "
              "in MatLab format.", 54, 2, stdout);
  PrintVADesc(1, "title", 11, "ParmTitle", 14, "Parameter file title.  This "
              "is not a file name, but rather the title appearing on the "
              "first line of the printed file.", 54, 2, stdout);
}

//-----------------------------------------------------------------------------
// PrintIPolQNamelistVariables: the ipolq namelist expedites the otherwise
//                              laborious process of computing a solvent
//                              reaction field potential for a molecule and
//                              then setting up quantum calculations in the
//                              vacuum and condensed phases.
//-----------------------------------------------------------------------------
void PrintIPolQNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("The Implicitly Polarized Charge method is laborious to "
                 "implement manually, but this module and associated namelist "
                 "automate much of the process, prevent errors, and make it "
                 "easy to determine convergence in the calculations.  A "
                 "typical REsP procedure requires a series of conformations "
                 "of the molecule of interest; the IPolQ procedure requires "
                 "a series of conformations of the molecule in a condensed "
                 "phase environment, such as water.  With the automation "
                 "afforded by this namelist, users can derive IPolQ charges "
                 "nearly as easily as standard REsP charges.\n", 79, NULL,
                 stdout);
  printf("   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "solute", 11, "SoluteMol", 14, "The solute molecule, "
              "specified by an ambmask string.  This is the molecule of "
              "interest for charge fitting, and will be immobilized during "
              "the simulation.  This must be specified by the user.", 54, 2,
              stdout);
  PrintVADesc(1, "ntqs", 11, "FrameRate", 14, "The rate of charge density "
              "sampling; the number of steps between successive snapshots to "
              "determine the solvent reaction field potential (SRFP).  "
              "Default 1000 (the time step set in the &cntrl or &ipolq "
              "namelists should factor into the ntqs setting).", 54, 2,
              stdout);
  PrintVADesc(1, "nqframe", 11, "FrameCount", 14, "The number of frames used "
              "to compose the SRFP.  Default 10 (this is too low for most "
              "solvent environments).", 54, 2, stdout);
  PrintVADesc(1, "nsteqlim", 11, "EqStepCount", 14, "The number of steps used "
              "to equilibrate the system before charge density collection "
              "begins.  Use this part of the simulation to buffer against any "
              "artifacts that might arise from suddenly freezing the solute "
              "in place.  Default 10000.", 54, 2, stdout);
  PrintVADesc(1, "nblock", 11, "Blocks", 14, "The number of blocks into which "
              "the simulation shall be divided for the purpose of estimating "
              "the convergence of the electrostatic potential.  Default 4.",
              54, 2, stdout);
  PrintVADesc(1, "verbose", 11, "Verbose", 14, "Default 0; set to 1 to "
              "activate step-by-step progress updates printed to the terminal "
              "window.  Useful for monitoring short runs to ensure that the "
              "input successfully completes the SRFP calculation.", 54, 2,
              stdout);
  PrintVADesc(1, "econv", 11, "EConverge", 14, "Convergence tolerance for the "
              "SRFP (convergence checking is not yet implemented).", 54, 2,
              stdout);
  PrintVADesc(1, "nqshell", 11, "QShellCount", 14, "The number of additional "
              "shells of charge to place around the system in order to "
              "approximate the SRFP due to infinite electrostatics in the "
              "confines of an isolated system.  Maximum (and default) is 3, "
              "minimum is 1.", 54, 2, stdout);
  PrintVADesc(1, "nvshell", 11, "VShellCount", 14, "The number of shells "
              "around each atom on which the exact SRFP due to infinite "
              "electrostatics shall be calculated.  Maximum (and default) is "
              "3, minimum is 1.", 54, 2, stdout);
  PrintVADesc(1, "nqphpt", 11, "QSpherePts", 14, "In order to generate the "
              "surface charges that will help in approximating the SRFP, "
              "this number of points is placed equidistant on a sphere.  The "
              "sphere is then rotated randomly and expanded to the radii "
              "indicated by qshell[1,2,3,x].  All points that are on the "
              "sphere due to one atom but within the sphere projected by "
              "another atom are deleted, until only points on the proper "
              "surface remain.  Default 100.", 54, 2, stdout);
  PrintVADesc(1, "nvphpt", 11, "VSpherePts", 14, "Similar to nqphpt above, "
              "but for the shells of SRFP evaluation points.  Default 20.", 54,
              2, stdout);
  PrintVADesc(1, "qshell1", 11, "ExpQBoundary", 14, "The distance at which to "
              "locate the first surface charges, and to stop collecting "
              "charges explicitly from the simulation's non-solute (that is, "
              "solvent) atoms.  Default 5.0.", 54, 2, stdout);
  PrintVADesc(1, "qshell[2,3]", 11, "QShell[2,3]", 14, "The distance at which "
              "to locate the second and third shells of boundary charges.  If "
              "engaged, each shell must be located successively further from "
              "the solute than the previous one.", 54, 2, stdout);
  PrintVADesc(1, "qshellx", 11, "QShellX", 14, "The distance at which to "
              "locate an interior shell of charges, and the minimum distance "
              "at which to take charges explicitly from the simulation.  Use "
              "this feature if there is a risk that QM basis functions might "
              "latch onto the simulation's point charges in absence of any "
              "exclusion effects, and thus distort the wavefunction.", 54, 2,
              stdout);
  PrintVADesc(1, "vshell[1-3]", 11, "VShell[1-3]", 14, "The distances at "
              "which to locate additional shells of exact SRFP evaluation "
              "points.  The SRFP is always evaluated, exactly, at the solute "
              "atom sites.", 54, 2, stdout);
  PrintVADesc(1, "dt", 11, "TimeStep", 14, "The simulation time step.  This "
              "is read in the &ipolq namelist just as if it were present in "
              "the &cntrl namelist, but a value specified in &ipolq overrides "
              "the &cntrl setting.  Default is 0.001ps, set in &cntrl.", 54, 2,
              stdout);
  PrintVADesc(1, "minqwt", 11, "MinQWeight", 14, "The stiffness of harmonic "
              "restraint by which to restrain fitted shell charges to zero.  "
              "Default 0.01.", 54, 2, stdout);
  PrintVADesc(1, "modq", 11, "ModifyQ", 14, "When IPolQ is applied, it is "
              "appropriate to hyper-polarize certain molecules in the SRFP "
              "calculation.  This variable may be specified as many times as "
              "necessary, followed by an ambmask string and a real number "
              "indicating the new charges to be assigned to all atoms in the "
              "mask.  For example, fixed-charge water models should have "
              "their dipoles increased by an amount equal to the original "
              "model's dipole less 1.85 (the dipole of water in vacuum).", 54,
              2, stdout);
  PrintVADesc(1, "prepqm", 11, "QuantumPrep", 14, "Preparatory call for QM "
              "calculations.  This variable may be specified as many times as "
              "necessary.  Each of these calls will be issued, in the order "
              "specified, before executing quantum calculations.", 54, 2,
              stdout);
  PrintVADesc(1, "postqm", 11, "QuantumClean", 14, "Post-processing calls for "
              "QM calculations.  Similar to prepqm directives, called after "
              "QM calculations have been completed.", 54, 2, stdout);
  PrintVADesc(1, "qmprog", 11, "QMPackage", 14, "The quantum package to use.  "
              "Supported packages are \"gaussian\" and \"orca\".", 54, 2,
              stdout);
  PrintVADesc(1, "qmpath", 11, "QMPath", 14, "Path to the primary QM "
              "executable.  This path will be tested, taking into account "
              "prepqm calls, to be sure that the executable exists prior to "
              "running the SRFP calculation.", 54, 2, stdout);
  PrintVADesc(1, "qmcomm", 11, "QMInputFile", 14, "The base name of the QM "
              "input file.  Vacuum and condensed-phase versions will be "
              "written with extensions 'vacu' and 'solv', respectively.  "
              "Default 'IPolQinp'.", 54, 2, stdout);
  PrintVADesc(1, "maxcore", 11, "MaxMemory", 14, "The maximum memory that can "
              "be allocated to arrays for quantum calculations with Orca, or "
              "the maximum total memory that can be allocated for "
              "calculations with Gaussian.", 54, 2, stdout);
  PrintVADesc(1, "qmthreads", 11, "QMThreads", 14, "The number of threads to "
              "be run by the QM package.  If unspecified, the number of "
              "threads run by mdgx itself is taken.  Otherwise any number can "
              "force the QM program to run in a particular configuration.", 54,
              2, stdout);
  PrintVADesc(1, "qmresult", 11, "QMOutputFile", 14, "The base name of the QM "
              "output file, which is given similar extensions to the input "
              "file.  Default 'IPolQout'.", 54, 2, stdout);
  PrintVADesc(1, "ptqfi", 11, "PointQFile", 14, "The name of the point "
              "charges file referenced by orca for including the SRFP into "
              "the condensed-phase calculation.", 54, 2, stdout);
  PrintVADesc(1, "singlept", 11, "SinglePoint", 14, "Rather than compute "
              "electrostatic potentials all around the molecule, perform a "
              "pair of single point calculations in the vacuum and condensed "
              "phases.  This is only supported with ORCA SCF-level "
              "calculations.", 54, 2, stdout);
  PrintVADesc(1, "qmflag", 11, "QMSignal", 14, "The name of the file used to "
              "signal slave processes that the QM calculations launched by "
              "the master are complete.  Default '.mdgx.finqm'.", 54, 2,
              stdout);
  PrintVADesc(1, "qmlev", 11, "QMTheory", 14, "The level of QM theory to "
              "use.  Default MP2.", 54, 2, stdout);
  PrintVADesc(1, "basis", 11, "QMBasis", 14, "The QM basis set to use.  "
              "Default cc-pvTZ.", 54, 2, stdout);
  PrintVADesc(1, "excitation", 11, "ExcitedState", 14, "Target the excitation "
              "state indicated by the specified integer in quantum "
              "calculations.  Default zero (ground state).  This feature may "
              "not give consistent results across different quantum programs!",
              54, 2, stdout);
  PrintVADesc(1, "scrdir", 11, "WorkDirectory", 14, "The scratch directory "
              "to use during QM calculations.  Useful to reduce NFS load.  If "
              "the directory exists, it will be used but not destroyed "
              "following each QM calculation.  If the directory does not "
              "exist at the start of the run, it will be created and later "
              "destroyed.", 54, 2, stdout);
  PrintVADesc(1, "rqminp", 11, "KeepQMInput", 14, "Directive to retain QM "
              "input files after the run.  Default 0 (OFF).", 54, 2, stdout);
  PrintVADesc(1, "rqmchk", 11, "KeepQMCheckPt", 14, "Directive to retain QM "
              "checkpoint file(s) after the run.  Default 0 (OFF).", 54, 2,
              stdout);
  PrintVADesc(1, "rqmout", 11, "KeepQMOutput", 14, "Directive to retain QM "
              "output files after the run.  Default 0 (OFF).", 54, 2, stdout);
  PrintVADesc(1, "rcloud", 11, "KeepQCloud", 14, "Directive to retain the "
              "solvent charge density cloud file after the run.  Default 0 "
              "(OFF).", 54, 2, stdout);
  PrintVADesc(1, "checkex", 11, "CheckExist", 14, "Activates safety checks "
              "for the existence of QM executables (including electrostatic "
              "potential calculators) called at the start of the run.  These "
              "checks attempt to take into account user-specified preparatory "
              "directives (see prepqm above).  Default 1 (ON).  Set to zero "
              "to disable this safeguard, for instance if the checks cannot "
              "find the executables but the preparatory directives, when "
              "fully implemented, are known to result in success.", 54, 2,
              stdout);
  PrintVADesc(1, "unx", 11, "UElecXBin", 14, "The number of grid points on "
              "which to evaluate the electrostatic potential, in the X "
              "direction.  Grid dimensions in Y and Z are set by similar "
              "variables.", 54, 2, stdout);
  PrintVADesc(1, "uhx", 11, "UElecXSpc", 14, "The grid spacing of the "
              "electrostatic potential grid in the X direction.  The grid is "
              "always rectilinear.  Spacings in Y and Z are set by similar "
              "variables.", 54, 2, stdout);
  PrintVADesc(1, "cengrid", 11, "CenterGrid", 14, "Directive to center the "
              "electrostatic potential grid on the location of the molecule "
              "stored in mdgx.  The default behavior varies with each quantum "
              "package: 'orca' activates centering on the molecule whereas "
              "'gaussian' calls for centering on the origin, as Orca does not "
              "reposition the molecule in its output but Gaussian will place "
              "the molecule in a 'Standard Orientation' and leave it there "
              "in the output and checkpoint files used for electrostatic "
              "potential calculations.", 54, 2, stdout);
  PrintVADesc(1, "fmpath", 11, "FormChkPath", 14, "Path to the program called "
              "for converting binary checkpoint files into formatted "
              "checkpoint files.  Needed only if the QM program is 'gaussian'"
              ".", 54, 2, stdout);
  PrintVADesc(1, "uvpath", 11, "UEvalPath", 14, "Path to the program called "
              "for evaluating the electrostatic potential grid.", 54, 2,
              stdout);
  PrintVADesc(1, "grid", 11, "GridFile", 14, "Base name of the electrostatic "
              "potential grid to be written.  As with QM input and output, "
              "this base name is appended 'vacu' or 'solv' for vacuum and "
              "condensed-phase calculations.", 54, 2, stdout);
}

//-----------------------------------------------------------------------------
// PrintConfigsNamelistVariables: the configs module is a tool for hands-on
//                                manipulation of small molecules, giving users
//                                the ability to easily generate many
//                                conformations of a single system sampling
//                                particular degrees of freedom.  This utility
//                                is especially useful for creating force field
//                                data sets and then checking the parameters
//                                afterwards.
//-----------------------------------------------------------------------------
void PrintConfigsNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("The mdgx configuration sampling module offers users rapid, "
                 "hands-on molecular manipulation capabilities.  In "
                 "principle, each operation is no different than a sander "
                 "energy minimization using NMR restraints.  The new feature "
                 "is that mdgx will perform optimizations on any number of "
                 "copies of the molecule at once, subject to NMR restraints "
                 "between the same atoms but with different target settings, "
                 "leading to a diversity of configurations (or, conformations)"
                 "at the end of a single run.  This is effective for cases in "
                 "which all interactions, without regard to a cutoff, must be "
                 "computed.  The results of these operations make good inputs "
                 "to quantum calculations for applications such as force "
                 "field development, and so the output structures can be "
                 "formatted in a number of ways, including input files for "
                 "several quantum packages.", 79, NULL, stdout);
  PrintParagraph("The memory layout of the coordinates is designed to "
                 "expedite this calculation and reduce the time spent looking "
                 "up atom properties.  As some configurations will reach "
                 "energy minima in fewer steps than others, they will be "
                 "removed from circulation to enhance performance.  Adding "
                 "these improvements to the obvious advantage of only booting "
                 "up the program once to do many configurations, this method "
                 "is many times faster than a homemade script running sander "
                 "to make each configuration one by one.  With all of the "
                 "configurations in memory at the end, mdgx is then able to "
                 "make sanity checks, draw conclusions about what aspects of "
                 "the restraint set may be problematic, and even restart "
                 "energy minimizations with different starting configurations "
                 "to try (often successfully) to improve the results.", 79,
                 NULL, stdout);
  PrintParagraph("The keywords for this module are arranged in several "
                 "sections, describing the minimization protocol, output "
                 "format, sanity checking, and finally details of the QM "
                 "protocol.", 79, NULL, stdout);
  printf(" Minimization Protocol:\n   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "verbose", 11, "Verbose", 14, "Sets the verbosity level (0 "
              "is silent, 1 will give frequent updates on the command line.",
              54, 2, stdout);
  PrintVADesc(1, "count", 11, "Replicas", 14, "The number of configurations "
              "to generate, if starting from a single configuration in inpcrd "
              "or restart format.", 54, 2, stdout);
  PrintVADesc(1, "maxcyc", 11, "MaxCycles", 14, "The maximum number of line "
              "minimization steps to attempt in any one round of energy "
              "minimization.", 54, 2, stdout);
  PrintVADesc(1, "ncyc", 11, "SDSteps", 14, "The number of steepest descent "
              "line minimization steps to perform before switching to a "
              "conjugate gradient method.  As with the eponymous keywords in "
              "sander, ncyc must be less than or equal to maxcyc.", 54, 2,
              stdout);
  PrintVADesc(1, "exclmax", 11, "ExclTableSize", 13, "The maximum number of "
              "atoms for which a table of non-bonded scaling factors will be "
              "kept.  For small systems, it is faster to pre-calculate "
              "whether non-bonded interactions will be excluded or attenuated "
              "and store these values in a matrix.  However, this is memory-"
              "intensive and will trash the cache for larger systems.  In "
              "those cases it is better to store a different sort of data "
              "structure that will quickly determine whether two atoms "
              "constitute an exclusion.", 54, 2, stdout);
  PrintVADesc(1, "frctol", 11, "ForceConverge", 14, "Convergence criterion "
              "for the optimization.  This is a quantity of force--if forces "
              "on all particles have lower magnitude than this value, the "
              "energy optimization for that configuration will be deemed "
              "converged.", 54, 2, stdout);
  PrintVADesc(1, "steptol", 11, "StepConverge", 14, "Convergence criterion "
              "for the optimization.  This is a quantity of distance--if the "
              "movement of all particles along the current force vector is "
              "driven lower than this value, the energy optimization for that "
              "configuration will be deemed converged.", 54, 2, stdout);
  PrintVADesc(1, "step0", 11, "InitialStep", 14, "Initial step size for the "
              "energy optimization.  This is a quantity of distance: the "
              "total magnitude of the initial step along the first computed "
              "force vector, that is the square root of the sums of squares "
              "of the displacements of all particles from their original "
              "positions, will be equal to this number (default 0.01A).  The "
              "step size will be iteratively changed throughout optimization "
              "and will be tailored to each configuration.", 54, 2, stdout);
  PrintVADesc(1, "freezeh", 11, "FreezePolarH", 14, "Pertinent only to "
	      "\"FREEZE\" output files, this will casue mdgx to print not "
	      "just constraints for heavy atom rotatable bonds but also "
	      "constraints for rotatable polar hydrogens, useful when "
	      "preventing hydrogen bonding in a QM optimization is desirable.",
              54, 2, stdout);
  PrintVADesc(1, "nshuffle", 11, "ShuffleCount", 14, "The number of times to "
              "restart energy minimizations towards the specified restraint "
              "targets using different initial states.", 54, 2, stdout);
  PrintVADesc(1, "shuffle", 11, "ShuffleStyle", 14, "The type of shuffling to "
              "perform if nshuffle > 0.  Available methods include "
              "\"bootstrap\" (new initial states will be assigned randomly "
              "from existing solutions, with replacement--one solution can "
              "serve as the initial state for more than one configuration), "
              "\"jackknife\" (the default--each existing solution will be "
              "assigned as the initial state for energy reoptimization to one "
              "and only one configuration), and \"proximity\" (every solved "
              "configuration will be evaluated in terms of the restraint "
              "targets of every other, and new initial states will be "
              "randomly chosen from among solutions whose restraint energies "
              "are within a certain threshold of the MINIMUM energy found for "
              "any existing solution with respect to the particular restraint "
              "targets of a given configuration.", 54, 2, stdout);
  PrintVADesc(1, "eprox", 11, "ProximateNrg", 14, "Threshold energy for "
              "taking existing solutions as the intial states for new "
              "attempts at energy minimization if using \"proximity\" "
              "reshuffling (default 5.0 kcal/mol).", 54, 2, stdout);
  PrintVADesc(1, "erep", 11, "ReplacementTol", 14, "Threshold for accepting a "
              "new solution based on a different initial state.  The new "
              "solution must supplant the energy of the existing one by at "
              "least this amount.  Default 1.0e-4 kcal/mol.\n", 54, 2, stdout);
  PrintVADesc(1, "shfdir", 11, "Direction", 14, "The direction to replace "
              "energies when reshuffling energy optimizations.  Choices are "
              "\"up\" and \"down\".  Default is \"down,\" but replacement can "
              "be made to move the energies upwards, finding new local minima "
              "with higher overall energies.", 54, 2, stdout);
  printf("\n Sampling Strategies:\n   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "random", 11, "RandomSample", 14, "Perform random sampling "
              "within a range.  This keyword must be followed by a series of "
              "commands, all on the same line of the input file, but the "
              "order of the sub-commands is flexible.  After seeing this "
              "keyword, mdgx will search the remainder of the line until it "
              "hits another keyword from the &configs namelist; until then it "
              "will associate any input it finds with the previous \"random\" "
              "or \"RandomSample\" keyword.  The range of sampling in this "
              "case is absolute: a flat-bottom harmonic potential will be "
              "constructed, centered on the spot randomly chosen between the "
              "limits \"min\" and \"max\", or given between two { } braces.  "
              "To specify that all configurations be restrained towards a "
              "single target value, the keyword \"center\" may be used in "
              "place of \"min\", \"max\", or { }.  The potential shall be "
              "flat up to a distance \"fbhw\" (flat bottom half width) from "
              "the center, and thereafter rise quadratically with a "
              "coefficient \"Krst\" (stiffness constant K of the restraint) "
              "over a length specified by the \"quadratic\" keyword, or up to "
              "a point at which the restraint force would reach a limit given "
              "by the \"Ftop\" keyword.  Beyond this limit, the force will "
              "be clamped and the restraint potential will be effectively "
              "linear, which helps to ensure that restraints to positions far "
              "from the initial configuration do not break things like "
              "chirality.  Because it is more intuitive to specify a maximum "
              "restraint force than a quadratic window, \"Ftop\" will take "
              "priority over \"quadratic\" if both keywords are given.  The "
              "defaults are to have 64 kcal/mol restraints applied after a "
              "0.5A flat bottom half width, topping out at 32 kcal/mol-A "
              "applied force.", 54, 2, stdout);
  PrintVADesc(1, "uniform", 11, "GridSample", 14, "Perform sampling on "
              "regular intervals within a range.  All of the keywords from "
              "RandomSample apply here as well.", 54, 2, stdout);
  PrintVADesc(1, "rpert", 11, "RandomPerturb", 14, "Perform sampling on "
              "regular intervals within a range based on the arrangement of "
              "atoms in each initial structure.  All of the keywords from "
              "RandomSample apply here as well, except that the range now "
              "specifies minimum and maximum values relative to the initial "
              "arrangement of atoms.  If multiple initial structures are read "
              "in, this will perturb each of them by similar random amounts.",
              54, 2, stdout);
  PrintVADesc(1, "gpert", 11, "GridPerturb", 14, "Perform sampling on "
              "regular intervals within a range based on the initial "
              "arrangement of atoms.  This is to RandomPerturb as GridSample "
              "is to RandomSample.", 54, 2, stdout);
  PrintVADesc(1, "combine", 11, "LinkOperations", 14, "Combine two operations "
              "involving grid-based, interval sampling.  Without any such "
              "combinations, the interval sampling restraints in each "
              "configuration will march from one end of their respective "
              "ranges to the other, in unison--this will generate a line of "
              "configurations in the multi-dimensional space defined by each "
              "restrained coordinate.  To sample two or three dimensions of "
              "the space simultaneously at regular intervals, combine the "
              "operations.  Up to three operations may be combined.  For N "
              "combined operations, mdgx will take the Nth root of the total "
              "number of configurations and take this many samples along each "
              "of the combined restraint dimensions.", 54, 2, stdout);
  PrintVADesc(1, "belly", 11, "MovingAtoms", 14, "Make only the atoms in the "
              "given ambmask string movable during geometry optimization.",
              54, 2, stdout);
  PrintVADesc(1, "nohbond", 11, "PreventHBonds", 14, "Orders mdgx to evaluate "
	      "potential hydrogen bond donors and acceptors in the system, "
	      "and devise a series of restraints to ensure that no such pairs "
	      "come close enough to form hydrogen bonds.  Atoms that are 1-4 "
	      "bonded, or more closely connected, are omitted from the set.",
	      54, 2, stdout);
  PrintVADesc(1, "hblimit", 11, "HBondMargin", 14, "The minimum distance that "
	      "a polar hydrogen will be allowed to approach a heavy atom "
	      "acceptor.  Default 2.5A.", 54, 2, stdout);
  PrintVADesc(1, "fbhwSpinOH", 11, "FBHWSpinOH", 14, "Flat bottom half width "
	      "for restraints randomly rotating hydroxyl groups.  Default "
	      "60 degrees.", 54, 2, stdout);
  PrintVADesc(1, "KSpinOH", 11, "KrstSpinOH", 14, "Restraint stiffness for "
	      "spinning hydroxyl groups.", 54, 2, stdout);
  printf("\n Output format:\n   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "showorig", 11, "ShowOrigins", 14, "Flag to have mdgx show "
              "the original files for each configuration that it solves.  If "
              "all configurations start from a single state in a single file, "
              "the default behavior is to withhold this reporting.  However, "
              "if there are many files, the origin of each configuration may "
              "not be so obvious, and while mdgx does attempt to alphabetize "
              "and organize long lists of files arising from directory "
              "searches or regular expressions the evolution of molecular "
              "configurations may be of interest.  In these cases the default "
              "behavior is to report the origins of each configuration.", 54,
              2, stdout);
  PrintVADesc(1, "outbase", 11, "OutputBase", 14, "The bases of the output "
              "file names for configurations.  The format will be "
              "<base><number>.<suffix>.  Multiple strings may follow this "
              "keyword, so long as they are all on the same line.  Each "
              "string provided will be matched with a suffix and a style "
              "provided, in the order each is given.", 54, 2, stdout);
  PrintVADesc(1, "outsuff", 11, "OutputSuffix", 14, "Suffixes of the output "
              "file names for printed configurations.", 54, 2, stdout);
  PrintVADesc(1, "write", 11, "OutputType", 14, "The type of output to write, "
              "options being \"CRD\" (old Amber .crd format trajectory), "
              "\"CDF\" (Amber netCDF trajectory), \"INPCRD\" (Amber ascii "
              "7-decimal place inpcrd file for individual configurations), "
              "\"PDB\" (PDB format, with descriptions of the way the "
              "configuration was generated in the REMARK section), \"FREEZE\" "
              "(print a file of GeomeTRIC hard constraints for QM "
	      "optimization runs on the OpenFF QC Fractal infrastructure) and "
              "\"ORCA\", \"GAUSSIAN\", \"MOLPRO\", and \"GAMESS\" for input "
              "files to various quantum packages.  If trajectories are being "
              "written, all configurations that pass sanity checks will be "
              "printed to the file.  For the other formats, individual "
              "configurations will be printed to separate files.  More than "
              "one type of output may be written after creating a set of "
              "configurations.", 54, 2, stdout);
  printf("\n Sanity checking:\n   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "bstrain", 11, "BondStrain", 14, "The maximum bond strain "
              "(according to the input force field, as given in the topology "
              "file) that will be tolerated in any configuration that is to "
              "be printed.", 54, 2, stdout);
  PrintVADesc(1, "astrain", 11, "AngleStrain", 14, "The maximum bond angle "
              "strain (according to the input force field, as given in the "
              "topology file) that will be tolerated in any configuration "
              "that is to be printed.", 54, 2, stdout);
  PrintVADesc(1, "strainlim", 11, "StrainLimit", 14, "The maximum restraint "
              "penalty that will be tolerated in any configuration.  Note "
              "that, for any of these sanity checks, convergence of the "
              "energy minimization is NOT an automatic fail--it will simply "
              "be noted in the report file summarizing the process.  "
              "Rather, the sanity checks pertain to features of the "
              "structures that appear well outside the applicable range of "
              "molecular mechanics functions.", 54, 2, stdout);
  printf("\n Quantum controls:\n   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "maxcore", 11, "MaxMemory", 14, "When ordering mdgx to print "
              "configurations as input files to quantum packages, this states "
              "how much memory should be available for QM calculations.", 54,
              2, stdout);
  PrintVADesc(1, "ncpu", 11, "CPUCount", 14, "The number of CPUs to apply in "
              "each QM calculation.", 54, 2, stdout);
  PrintVADesc(1, "spin", 11, "Multiplicity", 14, "The multiplicity to assign "
              "to this system (in all its configurations) for quantum "
              "calculations.  mdgx is not able to calculate this on its own.",
              54, 2, stdout);
  PrintVADesc(1, "qmlev", 11, "QMTheory", 14, "The level of theory to apply "
              "in quantum calculations.", 54, 2, stdout);
  PrintVADesc(1, "basis", 11, "QMBasis", 14, "The basis set to apply in "
              "quantum calculations.", 54, 2, stdout);
  PrintVADesc(1, "chk", 11, "Checkpoint", 14, "The checkpoint file to write "
              "if using Gaussian for QM calculations.", 54, 2, stdout);
  PrintVADesc(1, "optgeom", 11, "OptimizeGeom", 14, "QM optimization will be "
	      "turned on in subsequent QM input files written for the "
	      "resulting configurations.\n", 54, 2, stdout);
}

//-----------------------------------------------------------------------------
// PrintPptdNamelistVariables: print variables associated with the peptide
//                             multi-simulator.
//-----------------------------------------------------------------------------
void PrintPptdNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("The mdgx peptide multi-simulator is the program's first "
		 "CUDA extension for implicit solvent GB and gas-phase "
		 "molecular dynamics.  It treats a GPU as miniature Beowulf "
		 "cluster of streaming multiprocessors (SMPs) and uses the "
		 "device's block scheduler as a queueing system of sorts.  "
		 "This paradigm shift in GPU utilization can backfill idle "
		 "SMPs to reap enormous gains in total throughput on small "
		 "systems (928 atoms maximum) and may even exceed the "
		 "simulation rate of pmemd.cuda for very small systems (less "
		 "than 225 atoms).  While it offers RATTLE, the equivalent of "
		 "SHAKE for mdgx's velocity-Verlet integrator, the module "
		 "also offers a velocity-Verlet I/r-RESPA multi-time stepping "
		 "scheme which performs at least as well as SHAKE in most "
		 "cases, and often considerably better in terms of speed and "
		 "energy conservation.", 79, NULL, stdout);
  PrintParagraph("The sizes of systems served by this module cover a range "
	         "ideal for Generalized Born calculations.  Systems with more "
		 "than 928 atoms will engage the pmemd.cuda GB engine with "
		 "reasonable efficiency.  The mdgx engine is instead designed "
		 "for maximum throughput on one or more systems by simulating "
		 "independent copies on all of the GPU's SMPs.", 79, NULL,
		 stdout);
  PrintParagraph("A reliable approach for getting the best throughput on a "
		 "single system is to call for a number copies of the system "
		 "equal to number of SMPs on the GPU SMP, or twice that "
		 "number if the system has 512 or fewer atoms, or four times "
		 "that number if the system has 256 or fewer atoms.  The "
		 "number of SMPs will be displayed in mdgx output (Section 5,"
		 "'GPU Utilization'), as will the thread block and block grid "
		 "sizes.", 79, NULL, stdout);
  PrintParagraph("For the best throughput on an array of systems with varying "
		 "sizes, the first thing to understand is that simulation "
		 "time will scale as the square of the system size and cause "
		 "each system to finish at a different rate.  If the spread "
		 "of sizes is great, this will create a lot of idle SMPs as "
		 "the GPU works to finish the largest simulation.  However, "
		 "if mdgx has additional systems to run, it can backfill the "
		 "idle SMPs with more work.  The program will automatically "
		 "arrange the systems internally in decreasing order of size, "
		 "to run the largest first and the smallest last.  It is "
		 "therefore advantageous, if trying to simulate many systems "
		 "of disparate sizes, to queue up many more systems than the "
		 "size of the block grid (which will be determined by the "
		 "number of SMPs and the size of the largest sytem).  By "
		 "queueing three to four times as many systems as the size of "
		 "the mdgx block grid, the entire GPU can keep busy.", 79,
		 NULL, stdout);
  PrintParagraph("The main input for this section is the Peptide / oligomer "
		 "keyword, followed by a list of subdirectives reminiscent of "
		 "sander command line input.  Many directives will be carried "
		 "down from the &files and &cntrl namelists, such as DoRATTLE "
		 "/ rigidbonds, thermostat controls, and the time step.  The "
		 "&pptd namelist can override some of these directives for "
		 "specific oligomers.  Other parameters that can influence "
		 "the dynamics, such as the GB style, are native to the &pptd "
		 "namelist as this is the only context in which they can be "
		 "used.", 79, NULL, stdout);
  printf("   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "oligomer", 11, "Peptide", 14, "A system to simulate in "
	      "non-periodic conditions (implicit solvent or vacuum).  After "
	      "seeing this keyword, mdgx will search the remainder of the "
	      "line until it hits another keyword from the &pptd namelist; "
	      "until then it will associate any input it finds with the "
	      "previous \"oligomer\" or \"Peptide\" keyword.  Each oligomer "
	      "requires its own topology and input coordinates, specified by "
	      "the -p and -c flags to mirror sander command line input.  "
	      "Files for mdout, mdcrd, and mdrst can be supplied with "
	      "flags -o, -x, and -r, respectively, again like sander command "
	      "line input.  Multiple copies of the system can be specified by "
	      "including the N-rep flag followed by the number.  It is also "
	      "possible to simulate replicas at a range of temperatures by "
	      "providing the T-rep flag followed by an integer as well as a "
	      "temperature range with the flags Tmin and Tmax (each followed "
	      "by a real number).  Replicas will be simulated at evenly "
	      "spaced intervals of the temperature, inclusive of the two "
	      "end points (i.e. Tmin 100.0 Tmax 200.0 T-rep 11 would create "
	      "replicas at 100.0, 110.0, 120.0, ..., 200.0K).  To simulate "
	      "all replicas at one particular temperature which differs from "
	      "temp0 in &cntrl, temp0 may also be supplied as a flag for a "
	      "specific oligomer.  Also, the -p flag may be replaced by -pi "
	      "and -pf, each followed by a topology file, to create replicas "
	      "based on interpolated topologies.  The two topologies must "
	      "have similar atom counts, names, and bonding patterns, but "
	      "otherwise are just two endpoints.  With two topologies, the "
	      "P-rep flag followed by an integer will specify the number of "
	      "copies to make at regular intervals along a linear "
	      "interpolation between the topologies, again inclusive of the "
	      "end points.", 54, 2, stdout);
  PrintVADesc(1, "igb", 11, "GBStyle", 14, "Type of Generalized Born solvent "
	      "to use.  All standard sander settings, including 7 and 8 (neck "
	      "GB) and 6 (vacuum conditions) are available.", 54, 2, stdout);
  PrintVADesc(1, "offset", 11, "GBOffset", 14, "The offset for GB radii "
	      "calculations.  For igb=8 (Neck GB II), this is 0.195141.  For "
	      "all other models it is 0.09.", 54, 2, stdout);
  PrintVADesc(1, "bondstep", 11, "MinorSteps", 14, "The number of minor steps "
	      "to use in a velocity Verlet I/r-RESPA multiple time-stepping "
	      "scheme.  To say \"bond steps\" is a bit of a misnomer: bond, "
	      "angle, and 1-4 non-bonded interactions are all recalculated "
	      "on each minor step in between major steps where general "
	      "non-bonded and dihedral interactions are calculated.", 54, 2,
	      stdout);
  PrintVADesc(1, "diel", 11, "Dielectric", 14, "Dielectric constant for "
	      "the solvent, whether GB or some continuum homogeneous "
	      "environment.", 54, 2, stdout);
}

//-----------------------------------------------------------------------------
// PrintRuleNamelistVariables: print descriptions of &rule namelist keywords
//                             used to create custom virtual sites (extra
//                             points).
//-----------------------------------------------------------------------------
void PrintRuleNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("Virtual sites, also known as \"extra points\", are a "
		 "powerful means of expanding the monopole distributions of "
		 "molecular models, as well as their surface shapes if the "
		 "virtual sites are given van-der Waals properties.  Through "
		 "Amber20, the application of virtual sites in our premier "
		 "engine, pmemd, has been limited to a handful of famous "
		 "cases: there was little generality beyond the means to "
		 "implement classic TIP-$p and TIP-5P water models, although "
		 "some notable lone pair arrangements were possible to "
		 "produce by hand-editing the prepin files.", 79, NULL,
		 stdout);
  PrintParagraph("The new pmemd engine, in both its CPU and CUDA variants, "
		 "now supports a much more diverse palette of virtual sites "
		 "with two, three, and four atom frames, holding the virtual "
		 "sites at fixed distances or angles from their parent atoms "
		 "or allowing the site positions to stretch along with the "
		 "frame atoms.  What is lacking is a rapid means for adding "
		 "new virtual sites to an existing molecular topology with a "
		 "standard set of simple descriptions.  The CanDo program, "
		 "Amber's upcoming authoritative topology generator, will "
		 "take over this work in time, but until then mdgx provides "
		 "a stopgap means for adding new virtual sites onto an "
		 "existing Amber topology via the &rule namelist.", 79, NULL,
		 stdout);
  PrintParagraph("Virtual sites added by this method will be interwoven into "
		 "the residues of their parent atoms, and pre-existing "
		 "virtual sites, from the custom palette or the original "
		 "frames inferred from bonding patterns, will be incorporated "
		 "as part of the growing system.  The &rule namelist also "
		 "provides a method, frame type NONE, for changing the "
		 "non-bonded properties of a specific atom in the system.  "
		 "This is useful for modifying the charge distribution to "
		 "integrate a new, massless site with a monopole on it.", 79,
		 NULL, stdout);
  PrintParagraph("By far the most important keyword in the &rulenamelist is "
		 "eptype (alias FrameType).  This indicates the basic shape "
		 "of the virtual site's connection to the rest of the "
		 "structure.  mdgx knows seven custom frames in addition to "
		 "the original Amber frames (only the custom frames can be "
		 "specified by this keyword, as the original Amber frames are "
		 "inferred from bonding patterns when no other specifier is "
		 "available).  All of these seven custom frames can be "
		 "specified by two aliases, much like the keywords "
		 "themselves:", 79, NULL, stdout);
  printf("   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "Flex-2", 11, "FlexDis2", 14, "Two atom frame that stretches "
	      "along with the distance between the two parent atoms (which "
	      "are likely bonded to begin with).  The virtual site is placed "
	      "along the line between the parent atom (P) and frame2 atom "
	      "(f2), at a displacement of v12 times the vector (P)->(f2) from "
	      "(P).  Negative values of v12 can form sigma holes from halogen "
	      "atoms (P), while v12=0.5 would place the virtual site at the "
	      "bond midpoint.", 54, 2, stdout);
  PrintVADesc(1, "FD-2", 11, "FixedDis2", 14, "Two atom frame with absolute "
	      "distance v1e maintained between the virtual site and its "
	      "parent atom (P).  The virtual site is placed along the "
	      "direction (P)->(f2).", 54, 2, stdout);
  PrintVADesc(1, "Flex-3", 11, "FlexDis3", 14, "Three atom frame that "
	      "stretches along with the displacements betwee the parent atom "
	      "(P) and frame atoms (f2) and (f3).  It is the two-dimensional "
	      "analog of FD-2 / Flex-2, using proportions v12 and v13 of the "
	      "displacement vectors (P)->(f2) and (P)->(f3) starting again at "
	      "(P).", 54, 2, stdout);
  PrintVADesc(1, "FD-3", 11, "FixedDis3", 14, "Three atom frame maintaining "
	      "an absolute distance v1e between the parent atom (P) and the "
	      "virtual site.  This is the three-dimensional analog of FD-2 /"
	      "FixedDis2, where instead of (f2) we have an imaginary point "
	      "at some proportion v23 of the displacement (f2)->(f3).  In a "
	      "sense, the frame can 'flex' with relative motion of frame "
	      "atoms 2 and 3, (f2) and (f3), but this merely changes the "
	      "direction that the virtual site points coming off of (P).", 54,
	      2, stdout);
  PrintVADesc(1, "FAD-3", 11, "FixAnglDis", 14, "Three atom frame maintaining "
	      "an absolute distane v1e between the parent atom (P) and the "
	      "virtual site (EP), as well as a fixed angle in (EP)-(P)-(f2) "
	      "given by the keyword theta.  TIP4P water can be created by "
	      "setting v1e to 0.125 Angstroms and theta to 0.91185 radians.",
	      54, 2, stdout);
  PrintVADesc(1, "Out-3", 11, "OutOfPlane", 14, "Three atom frame placing the "
	      "virtual site as some linear combination of the parent atom (P) "
	      "to frame atoms (f2) and (f3) displacements, plus a proportion "
	      "of the cross product between these two vectors (normal to the "
	      "plane of the frame atoms).  The virtual site is placed off of "
	      "(P), at v12 times the displacement (P)->(f2) and v13 times the "
	      "displacement (P)->(f3), then out of the plane by v12x13 times "
	      "cross product (P)->(f2) x (P)->(f3).", 54, 2, stdout);
  PrintVADesc(1, "FD-4", 11, "FourPoint", 14, "Four atom frame placing the "
	      "virtual site at a distacnce v1e from the parent atom (P) based "
	      "on a vector computed from cross products of displacements from "
	      "frame atoms (f2), (f3), and (f4).  Relevant inputs include "
	      "v1e, v13, and v23.", 54, 2, stdout);
  PrintVADesc(1, "None", 11, "NONE", 14, "This \"frame\" does not add any "
	      "particles but modifies existing ones.  The charge and "
	      "Lennard-Jones properties of any atom in the topology can be "
	      "changed with a &rule specifying this frame type.n", 54, 2,
	      stdout);
  PrintParagraph("\nStandard keywords for the &rule namelist are as follows:",
		 79, NULL, stdout);
  printf("   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "epname", 11, "ExtraPoint", 14, "The name of the new virtual "
	      "site, akin to atom name CA of residue ARG.", 54, 2, stdout);
  PrintVADesc(1, "atom", 11, "AtomName", 14, "The name of an atom to modify.  "
	      "It is presumed that this atom is original to the topology, but "
	      "technically one can specify the name of a new virtual site "
	      "with this keyword as well.  The intent is to use these "
	      "keywords when modifying an existing atom, and epname / "
	      "ExtraPoint when adding a new particle.", 54, 2, stdout);
  PrintVADesc(1, "residue", 11, "ResidueName", 14, "Name of the residue that "
	      "the parent atom resides in.  Every residue will have unique "
	      "names for all atoms it contains, so specifying the parent atom "
	      "name and the residue is unambiguous identification for a "
	      "virtual site &rule.", 54, 2, stdout);
  PrintVADesc(1, "parent", 11, "ParentAtom", 14, "The parent atom of the "
	      "new virtual site, the primary connection to the rest of the "
	      "molecule.  See also the keywords frame1 and FrameAtom1, also "
	      "pointing to this atom.", 54, 2, stdout);
  PrintVADesc(1, "frame1", 11, "FrameAtom1", 14, "The first frame atom for "
	      "the new virtual site, identical to specifying the parent / "
	      "ParentAtom keywords.", 54, 2, stdout);
  PrintVADesc(1, "frame2", 11, "FrameAtom2", 14, "The second frame atom for "
	      "the new virtual site, needed for all frames.", 54, 2, stdout);
  PrintVADesc(1, "frame3", 11, "FrameAtom3", 14, "The third frame atom for "
	      "the new virtual site, needed for in-plane and out-of-plane "
	      "frames.", 54, 2, stdout);  
  PrintVADesc(1, "frame4", 11, "FrameAtom4", 14, "The fourth frame atom for "
	      "the new virtual site, needed only for fixed-distance pyramidal "
	      "frames.", 54, 2, stdout);
  PrintVADesc(1, "v12", 11, "Vector12", 14, "Displacement multiplier between "
	      "the parent atom (P) and frame atom (f2).  The parent atom is "
	      "also known as frame atom 1, so (f1)->(f2).", 54, 2, stdout);
  PrintVADesc(1, "v13", 11, "Vector13", 14, "Displacement multiplier between "
	      "the parent atom (P) and frame atom (f3).", 54, 2, stdout);
  PrintVADesc(1, "v1e", 11, "Vector1E", 14, "Absolute distance between the "
	      "parent atom (P) and the virtual site / extra point.  This "
	      "specifier is used in fixed distance frame types 2-FD, 3-FD, "
	      "3-FAD, and 4-FD.  Units are Angstroms.", 54, 2, stdout);
  PrintVADesc(1, "theta", 11, "Theta", 14, "Angle between the virtual site, "
	      "parent atom (P), and frame atom 2 (f2).  Units are radians.",
	      54, 2, stdout);
  PrintVADesc(1, "v23", 11, "Vector23", 14, "Proportion of the displacement "
	      "between frame atoms 2 and 3, (f2)->(f3).  This is used in "
	      "frame types 3-FD and 4-FD.", 54, 2, stdout);
  PrintVADesc(1, "v12x13", 11, "Vector12x13", 14, "Proportion of the cross "
	      "product between displacements in the frame atoms.  This is "
	      "needed by the Out-3 frame type.", 54, 2, stdout);
  PrintVADesc(1, "q", 11, "Charge", 14, "Virtual site partial charge (or "
	      "atomic partial charge, if the &rule has frame type NONE and "
	      "modifies an existing atom instead).", 54, 2, stdout);
  PrintVADesc(1, "sig", 11, "Sigma", 14, "Lennard Jones sigma radius for the "
	      "virtual site (or modified atom) when interacting with itself.  "
	      "Lorentz-Berthelot combining rules will determine how these "
	      "Lennard-Jones parameters combine with other atoms in the "
	      "system.", 54, 2, stdout);
  PrintVADesc(1, "eps", 11, "Epsilon", 14, "Lennard Jones Epsilon well depth "
	      "for the virtual site (or modified atom) when interacting with "
	      "itself.  Lorentz-Berthelot combining rules will determine how "
              "these Lennard-Jones parameters combine with other atoms in the "
              "system.", 54, 2, stdout);
  PrintParagraph("\nIn addition to the above keywords, one can find forbidden "
		 "magic hidden in the &rule namelist.  The keywords below "
		 "are for very advanced, niche purposes.", 79, NULL, stdout);
  printf("   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "style", 11, "FrameStyle", 14, "This is a keyword for expert "
	      "use, specifying the internal number of the desired frame.  As "
	      "mdgx and pmemd use different internal number designations for "
	      "the various frames, this can get quite confusing for the "
	      "uninitiated!  The topologies that mdgx prints can be properly "
	      "interpreted by pmemd or other mdgx modules, so the codes can "
	      "talk to each other.", 54, 2, stdout);
  PrintVADesc(1, "excl2", 11, "Exclude2", 14, "Flag to have the virtual site "
	      "inherit exclusions from the second frame atom, not just the "
	      "parent.  The stanadard approach is for a virtual site to be "
	      "counted as 1:1 to its parent and therefore 1:2 to any atoms "
	      "bonded to the parent, and so forth.  For some frames, notably "
	      "a virtual site at a bond midpoint, it is not clear whether "
	      "the exclusions of just one atom are appropriate.  The effects "
	      "of this directive CANNOT be properly interpreted by pmemd, but "
	      "my be of interest for force field development which would then "
	      "spur innovation in the workhorse engine.", 54, 2, stdout);
  PrintVADesc(1, "excl3", 11, "Exclude3", 14, "Flag to have the virtual site "
	      "also inherit exclusions from the third frame atom.", 54, 2,
	      stdout);
  PrintVADesc(1, "excl4", 11, "Exclude4", 14, "Flag to have the virtual site "
	      "also inherit exclusions from the fourth frame atom.", 54, 2,
	      stdout);
}

//-----------------------------------------------------------------------------
// PrintParmEditNamelistVariables: print descriptions of keywords used to edit
//                                 topologies with mdgx, primarily for adding
//                                 custom virtual sites to a system in a way
//                                 that pmemd is prepared ot understand.
//-----------------------------------------------------------------------------
void PrintParmEditNamelistVariables()
{
  PrintSplash(stdout);
  PrintParagraph("As a temporary measure to support the new pmemd capability "
		 "for simulating an expadned palette of virtual sites, mdgx "
		 "has the ability to add virtual sites to a topology via the "
		 "&rule namelist, and then print them with the &parmedit "
		 "module.  The handful of directives in this namelist help "
		 "to perform quality checks and modify standard topologies "
		 "produced by the tLEaP program, until the day when CanDo "
		 "becomes the standard Amber topology builder with these "
		 "virtual site features.", 79, NULL, stdout);
  printf("   Name       Alias      Description\n"
         " -------- ------------- --------------------------------------------"
         "-----------\n");
  PrintVADesc(1, "newprmtop", 11, "NewPrmtop", 14, "The name of the new "
	      "topology to write.  Standard rules for overwriting apply.  It "
	      "is acceptable, if the user assumes the risk to overwrite the "
	      "original topology as all of the needed information is read "
	      "into memory and processed prior to writing a new prmtop.", 54,
	      2, stdout);
  PrintVADesc(1, "newinpcrd", 11, "NewInpcrd", 14, "The name of the new input "
	      "coordinate file to write.  If virtual sites (extra points) are "
	      "inserted, the coordinates will be adjusted as necessary.  New "
	      "virtual sites will be placed in their frame positions prior to "
	      "printing the coordinates, so examination of the new inpcrd "
	      "serves as a convenient sanity check on the list of &rules "
	      "used to build an expanded topology.", 54, 2, stdout);
  PrintVADesc(1, "title", 11, "PrmtopTitle", 14, "Title to give the new "
	      "topology.  A maximum length of 64 characters is recommended.\n",
	      54, 2, stdout);
  PrintVADesc(1, "tchk", 11, "TopologyCheck", 14, "Default 1 (ON), this will "
	      "call for some basic checks on the sanity of the new topology, "
	      "including chirality (for standard amino acids), disulfide "
	      "linkages, and integral net charges on all residues.", 54, 2,
	      stdout);
  PrintVADesc(1, "qtol", 11, "ChargeSanity", 14, "Tolerance for residue "
	      "charge sanity check.  Deviations greater than this amount "
	      "(proton units) for any one residue from an integer net charge, "
	      "or five times this amount from neutrality for the overall "
	      "system, will trigger a warning.", 54, 2, stdout);
}

//-----------------------------------------------------------------------------
// PrintAttributions: users are provided with a straightforward means of seeing
//                    which aspects of this code were taken from other sources.
//                    Attributions are also provided in comments to the source
//                    code, where appropriate.   
//-----------------------------------------------------------------------------
void PrintAttributions()
{
  PrintSplash(stdout);
  PrintParagraph("Attributions:", 79, NULL, stdout);
  PrintParagraph("- Implementation of the SETTLE algorithm (J. Comput. Chem. "
                 "13:952-966, 1992) was adapted from the NAMD program source, "
                 "v2.6, developed by the Theoretical and Computational "
                 "Biophysics Group in the Beckman Institute for Advanced "
                 "Science and Technology at the University of Illinois at "
                 "Urbana-Champaign.", 79, NULL, stdout);
  PrintParagraph("- Implementation of the Smooth Particle Mesh Ewald "
                 "algorithm (J. Chem. Phys. 103, 8577-8593, 1995) includes "
                 "code developed by Thomas A. Darden for optimization of the "
                 "convolution kernel.", 79, NULL, stdout);
  PrintParagraph("- Implementation of the Langevin dynamics integration "
                 "algorithm (Biopolymers 32, 523-535, 1992) is adapted from "
                 "the sff program also distributed with AmberTools.", 79, NULL,
                 stdout);
  PrintParagraph("- Angle fitting, with restraints, is implemented as part "
                 "of the same linear least squares problem that performs "
                 "torsion fitting using concepts put forward by Dr. Kenno "
                 "Vanommeslaeghe (J. Comp. Chem. 36, 1083–1101, 2015).", 79,
                 NULL, stdout);
  PrintParagraph("- Dr. Robert E. Duke is acknowledged for outstanding advice "
                 "and insights into the problem of efficient and scalable "
                 "molecular dynamics.  The mdgx program would not exist "
                 "without his support.", 79, NULL, stdout);
  PrintParagraph("- Rubryc Therapeutics, Inc. provided funds for development "
		 "of the mdgx.cuda engine covering implicit solvent dynamics "
		 "of peptides and other small systems.", 78, NULL, stdout);
}
