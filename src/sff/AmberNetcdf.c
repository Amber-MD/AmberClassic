/*
 * Daniel R. Roe
 * 2010-12-07
 * A C implementation of routines for reading and writing the Amber Netcdf 
 * trajectory format. 
 * Original implementation of netcdf in Amber by Jon Mongan.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h> // For title length
#include "AmberNetcdf.h"

// If MPI define a function that will return current task ID in order to avoid
// printf on non-master threads.
// Otherwise define a function that just returns 0;
int netcdfTaskid();
#ifdef MPI
#include "mpi.h"
int netcdfTaskid() {
  int worldrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldrank);
  return worldrank;
}
#else
int netcdfTaskid() { return 0; }
#endif
    
#ifdef BINTRAJ
#include "netcdf.h"
// ---------- Definitions of Netcdf variable labels ----------------------------
#define NCFRAME "frame"
#define NCSPATIAL "spatial"
#define NCATOM "atom"
#define NCCELL_SPATIAL "cell_spatial"
#define NCCELL_LENGTHS "cell_lengths"
#define NCCELL_ANGULAR "cell_angular"
#define NCCELL_ANGLES "cell_angles"
#define NCCOORDS "coordinates"
#define NCVELO "velocities"
#define NCTEMPERATURE "temp0"
#define NCTIME "time"
#define NCLABEL "label"
#define NCLABELLEN 5
// For debugging make > 0:
#define NCDEBUGVAL 0

/// Macro to initialize the AmberNetcdf structure
#define INIT_AmberNetcdf( _p_ ) \
  _p_->temp0 = 0.0; \
  _p_->restartTime = 0.0; \
  _p_->isNCrestart = 0; \
  _p_->ncid=-1; \
  _p_->frameDID=-1; \
  _p_->ncframe=-1; \
  _p_->currentFrame=0; \
  _p_->atomDID=-1; \
  _p_->ncatom=-1; \
  _p_->ncatom3=-1; \
  _p_->coordVID=-1; \
  _p_->velocityVID=-1; \
  _p_->cellAngleVID=-1; \
  _p_->cellLengthVID=-1; \
  _p_->spatialDID=-1; \
  _p_->labelDID=-1; \
  _p_->cell_spatialDID=-1; \
  _p_->cell_angularDID=-1; \
  _p_->spatialVID=-1; \
  _p_->timeVID=-1; \
  _p_->cell_spatialVID=-1; \
  _p_->cell_angularVID=-1; \
  _p_->TempVID=-1 

// ---------- Private Routines -------------------------------------------------
// checkNCerr()
/** Check error status of err - print error message if not NC_NOERR
  */
static int checkNCerr(int err, const char *message, ...) {
  va_list args;

  if (err!=NC_NOERR) {
    va_start(args,message);
    fprintf(stderr,"NETCDF Error (%s): ",nc_strerror(err));
    vfprintf(stderr,message,args);
    fprintf(stderr,"\n");
    va_end(args);
    return 1;
  }
  return 0;
}

// GetDimInfo()
/** Return the dimension ID of a given attribute in netcdf file ncid.
  * Also set dimension length.
  */ 
static int GetDimInfo(int ncid, const char *attribute, int *length) {
  int dimID;
  size_t slength;
    
  *length = 0;
  slength = 0;
    
  // Get dimid 
  if ( checkNCerr(nc_inq_dimid(ncid, attribute, &dimID),
       "ncGetDimInfo: Getting dimID for attribute %s",attribute)!=0 ) return -1;

  // get Dim length 
  if ( checkNCerr(nc_inq_dimlen(ncid, dimID, &slength),
       "ncGetDimInfo: Getting length for attribute %s",attribute)!=0) return -1;

  *length = (int) slength;
  return dimID;
}

// GetAttrText()
/** Get the information about a netcdf attribute with given vid and 
  * attribute text.
  * Since there is no guarantee that NULL char at the end of retrieved string
  * append one.
  */
static char *GetAttrText(int ncid, int vid, const char *attribute) {
  size_t attlen;
  char *attrText;
  // Get attr length
  if ( checkNCerr(nc_inq_attlen(ncid, vid, attribute, &attlen),
       "ncGetAttrText: Getting length for attribute %s\n",attribute)) return NULL;
  // Allocate space for attr text, plus one for NULL char
  attrText = (char*) malloc( (attlen + 1) * sizeof(char));
  // Get attr text
  if ( checkNCerr(nc_get_att_text(ncid, vid, attribute, attrText),
       "ncGetAttrText: Getting attribute text for %s\n",attribute)) {
    free(attrText);
    return NULL;
  }
  // Append NULL char
  attrText[attlen]='\0';

  return attrText;
}
// -----------------------------------------------------------------------------
#endif

// netcdfDebug()
/** For use in printing various attributes of a previously opened netcdf file.
  */
int netcdfDebug(struct AmberNetcdf *A) {
#ifdef BINTRAJ
  int ndimsp, nvarsp, ngattsp,unlimdimidp;
  int err,i,mytaskid;
  char *varname;

  mytaskid = netcdfTaskid();
  /* ncid:    NetCDF ID, from a previous call to nc open or nc create.
   * ndimsp:  Pointer to location for returned number of dimensions defined for 
   *          this netCDF dataset.
   * nvarsp:  Pointer to location for returned number of variables defined for 
   *         this netCDF dataset.
   * ngattsp: Pointer to location for returned number of global attributes 
   *         defined for this netCDF dataset.
   * unlimdimidp: 
   *  Pointer to location for returned ID of the unlimited dimension, if 
   *  there is one for this netCDF dataset. If no unlimited length 
   *  dimension has been defined, -1 is returned.
   */
  varname=(char*) malloc(1024*sizeof(char));
  err=nc_inq(A->ncid,&ndimsp,&nvarsp,&ngattsp,&unlimdimidp);
  if (mytaskid==0) {
    fprintf(stdout,"========== BEG. NETCDF DEBUG ==========\n");
    fprintf(stdout,"nc_inq returned %i\n",err);
    if (err==NC_NOERR)
      fprintf(stdout,"ndimsp=%i  nvarsp=%i  ngattsp=%i  unlimdimidp=%i\n",
              ndimsp,nvarsp,ngattsp,unlimdimidp);
    else
      fprintf(stdout,"NETCDF Error occurred.\n");
    fprintf(stdout,"NC VARIABLES:\n");
  }
  // Print name of each variable defined in netcdf file
  for (i=0; i<nvarsp; i++) {
    err=nc_inq_varname(A->ncid,i,varname);
    if (mytaskid==0) fprintf(stdout,"  Var %i - ",i);
    if (err==NC_NOERR) {
      if (mytaskid==0) fprintf(stdout,"%s\n",varname);
    } else {
      if (mytaskid==0) fprintf(stdout,"NETCDF Error occured.\n");
    }
  }

  // Struct variables
  if (mytaskid==0) {
    fprintf(stdout,"  TempVID = %i\n",A->TempVID);

    fprintf(stdout,"==========  END NETCDF DEBUG ==========\n");
  }

  free(varname);
  return 0;
#else
  fprintf(stdout,"Error: NAB compiled without Netcdf support.\n");
  fprintf(stdout,"       Recompile with -DBINTRAJ\n");
  return 1;
#endif
}

// netcdfLoad()
/** Load the netcdf trajectory/restart file pointed to by filename and 
  * initialize the AmberNetcdf data structure. Allocate memory for coords 
  * (and velo if restart). Set currentFrame to 0.
  * \return 0 on success, 1 on error.
  */
int netcdfLoad(struct AmberNetcdf *A, char *filename) {
#ifdef BINTRAJ
  int debug;
  char *attrText;

  debug=NCDEBUGVAL;

  // Check netcdf structure
  if (A==NULL) {
    fprintf(stderr,"Error: netcdfLoad: AmberNetcdf is NULL\n");
    return 1;
  }
  // Close if already open
  //netcdfClose(A);

  INIT_AmberNetcdf(A);

  // Do not print error message on error, just return 1
  // This allows one to gracefully check if this is a netcdf file
  if ( nc_open(filename,NC_NOWRITE,&(A->ncid))!=NC_NOERR ) return 1;

  if (debug>0) fprintf(stdout,"Successfully opened %s, ncid=%i\n",filename,A->ncid);
  if (debug>1) netcdfDebug(A);

  // Determine whether this is a netcdf trajectory or restart.
  attrText = GetAttrText(A->ncid,NC_GLOBAL, "Conventions");
  if (attrText==NULL) {
    fprintf(stderr,"Error: netcdfLoad: Could not get netcdf Conventions for %s\n",filename);
    return 1;
  }
  if (strstr(attrText,"AMBERRESTART")!=NULL)
    A->isNCrestart=1;
  else if (strstr(attrText,"AMBER")==NULL) {
    printf("WARNING: Netcdf file %s conventions do not include \"AMBER\" (%s)\n",
            filename, attrText);
  }
  if (debug>0) {
    if (A->isNCrestart)
      fprintf(stdout,"File %s is a netcdf restart file.\n",filename);
    else
      fprintf(stdout,"File %s is a netcdf trajectory file.\n",filename);
  }
  free(attrText);

  // If trajectory, get frame dimension. Not needed for restart
  if (A->isNCrestart) {
    A->ncframe = 1;
  } else {
    A->frameDID=GetDimInfo(A->ncid,NCFRAME,&(A->ncframe));
    if (A->frameDID==-1) return 1;
  } 
  // Get number of atoms
  A->atomDID=GetDimInfo(A->ncid,NCATOM,&(A->ncatom));
  if (A->atomDID==-1) return 1;
  A->ncatom3 = A->ncatom * 3;
  // Get coordinate info
  if (checkNCerr(nc_inq_varid(A->ncid,NCCOORDS,&(A->coordVID)),
      "Getting coordinate VID")!=0) return 1;
  // Get spatial info

  // For restart, get velocity and time info
  if (A->isNCrestart) {
    if ( nc_inq_varid(A->ncid,NCVELO,&(A->velocityVID))==NC_NOERR ){
      if( debug>0 )
         fprintf(stdout,"Netcdf restart %s has velocity info.\n",filename);
    } else {
      if( debug>0 )
          fprintf(stdout,"Netcdf restart %s does not have velocity info.\n",filename);
      A->velocityVID=-1;
    }
    if ( nc_inq_varid(A->ncid, NCTIME, &(A->timeVID)) == NC_NOERR) {
      attrText = GetAttrText(A->ncid,A->timeVID, "units");
      if (attrText==NULL || strcmp(attrText,"picosecond")!=0)
        printf("WARNING: Netcdf restart file %s has time units of %s - expected picosecond.\n",
               filename, attrText);
      if (attrText!=NULL) free(attrText);
      if ( checkNCerr(nc_get_var_double(A->ncid, A->timeVID, &(A->restartTime)),
                      "Getting netcdf restart time.")) return 1;
      if (debug>0) printf("Netcdf restart time= %lf\n",A->restartTime);
    } else {
      fprintf(stderr,"Error: Could not get time from Netcdf restart %s\n",filename);
      return 1;
    }
  }

  // Box info 
  if ( nc_inq_varid(A->ncid,NCCELL_LENGTHS,&(A->cellLengthVID))==NC_NOERR ) {
    if (checkNCerr(nc_inq_varid(A->ncid,NCCELL_ANGLES,&(A->cellAngleVID)),
      "Getting cell angles.")!=0) return 1;
    if (debug>0) fprintf(stdout,"  Netcdf Box information found.\n"); 
  } else {
    A->cellLengthVID=-1;
    A->cellAngleVID=-1;
  }

  // Replica Temperatures
  if ( nc_inq_varid(A->ncid,NCTEMPERATURE,&(A->TempVID)) == NC_NOERR ) {
    if (debug>0) fprintf(stdout,"Netcdf file %s has replica temperatures.\n",filename);
  } else 
    A->TempVID=-1;

  // NOTE: TO BE ADDED
  // spatialDID, labelDID;
  // cell_spatialDID, cell_angularDID;
  // spatialVID, timeVID (for traj), cell_spatialVID, cell_angularVID;

/*  // For restarts, allocate memory for velocities
  if (A->isNCrestart) {
    if (A->velocityVID!=-1) {
      A->Velo=(double*) malloc(A->ncatom3 * sizeof(double));
      if (A->Velo==NULL) {
        fprintf(stderr,"Error: Could not allocate memory for netcdf velocities.\n");
        return 1;
      } 
    }
  // For trajectories, allocate memory for reading single-precision coordinates
  } else {
    A->Coord=(float*) malloc(A->ncatom3 * sizeof(float));
    if (A->Coord==NULL) {
      fprintf(stderr,"Error: Could not allocate memory for netcdf coords.\n");
      return 1;
    }
  }*/

  return 0;
#else
  fprintf(stdout,"Error: NAB compiled without Netcdf support.\n");
  fprintf(stdout,"       Recompile with -DBINTRAJ\n");
  return 1;
#endif
} 

// netcdfClose()
/** Close netcdf trajectory/restart file and free memory */
int netcdfClose(struct AmberNetcdf *A) {
#ifdef BINTRAJ
  if (A==NULL) return 0;
  if (A->ncid!=-1)
    checkNCerr(nc_close(A->ncid),"Closing netcdf file.");
  //if (A->Coord!=NULL) free(A->Coord);
  //if (A->Velo!=NULL) free(A->Velo);
  INIT_AmberNetcdf(A);
#else
  fprintf(stdout,"Error: NAB compiled without Netcdf support.\n");
  fprintf(stdout,"       Recompile with -DBINTRAJ\n");
  return 1;
#endif
  return 0;
}

// netcdfWriteRestart()
/** Create and write to an Amber netcdf restart file with name
  * specified by filename.
  * \param filename File name of netcdf restart
  * \param natom Number of atoms
  * \param X Coordinates
  * \param Velo If not NULL write velocity info
  * \param box If not NULL write box info
  * \param restart_time Restart time
  * \param restart_temp Restart temperature
  * \return 0 on successful write, 1 on error.
  */
int netcdfWriteRestart(char *filename, int natom, double *X, double *Velo, 
                       double *box, double restart_time, double restart_temp) 
{
#ifdef BINTRAJ
  struct AmberNetcdf A;
  int dimensionID[NC_MAX_VAR_DIMS];
  size_t start[2], count[2];
  char xyz[3];
  double restartTime;
  char abc[15] = { 'a', 'l', 'p', 'h', 'a',
                   'b', 'e', 't', 'a', ' ',
                   'g', 'a', 'm', 'm', 'a' };
  char *title;
  const double velocityScale=20.455; // Velocity scaling factor
  int debug=NCDEBUGVAL;

  if (filename==NULL) {
    fprintf(stderr,"Error: netcdfWriteRestart: filename is NULL.\n");
    return 1;
  }
  if (natom<=0) {
    fprintf(stderr,"Error: netcdfWriteRestart: (%s) # atoms <= 0 (%i)\n",filename,natom);
    return 1;
  }
  if (X==NULL) {
    fprintf(stderr,"Error: netcdfWriteRestart: (%s) Input coords are NULL\n",filename);
    return 1;
  }
  title = NULL;
  // Allocate memory for AmberNetcdf structure. Technically not necessary, but
  // makes bookkeeping much easier, as well as consistency with the rest of
  // the code in this file.
  // NOTE: Use static struct for now, otherwise need a free with every return
  //       statement.
  //A = (struct AmberNetcdf*) malloc( sizeof(struct AmberNetcdf) );
  INIT_AmberNetcdf( (&A) );

  // Create file
  if (checkNCerr(nc_create(filename,NC_64BIT_OFFSET,&(A.ncid)),
    "Creating Netcdf restart file %s",filename)) return 1;
  if (debug>0) 
    printf("|   Successfully created Netcdf restart file %s, ncid %i\n",filename,A.ncid);

  // Time variable
  if (checkNCerr(nc_def_var(A.ncid,NCTIME,NC_DOUBLE,0,dimensionID,&(A.timeVID)),
    "Defining time variable.")) return 1;
  if (checkNCerr(nc_put_att_text(A.ncid,A.timeVID,"units",10,"picosecond"),
    "Writing time VID units.")) return 1;
  
  // Spatial dimension and variable
  if (checkNCerr(nc_def_dim(A.ncid,NCSPATIAL,3,&(A.spatialDID)),
    "Defining spatial dimension.")) return 1;
  dimensionID[0] = A.spatialDID;
  if (checkNCerr(nc_def_var(A.ncid,NCSPATIAL,NC_CHAR,1,dimensionID,&(A.spatialVID)),
    "Defining spatial variable.")) return 1;

  // Atom dimension
  A.ncatom = natom;
  A.ncatom3 = natom * 3;
  if (checkNCerr(nc_def_dim(A.ncid,NCATOM,A.ncatom,&(A.atomDID)),
    "Defining atom dimension.")) return 1;

  // Coord variable
  dimensionID[0] = A.atomDID;
  dimensionID[1] = A.spatialDID;
  if (checkNCerr(nc_def_var(A.ncid,NCCOORDS,NC_DOUBLE,2,dimensionID,&(A.coordVID)),
    "Defining coordinates variable.")) return 1;
  if (checkNCerr(nc_put_att_text(A.ncid,A.coordVID,"units",8,"angstrom"),
    "Writing coordinates variable units.")) return 1;

  // Velocity variable
  if (Velo!=NULL) {
    if (checkNCerr(nc_def_var(A.ncid,NCVELO,NC_DOUBLE,2,dimensionID,&(A.velocityVID)),
      "Defining velocities variable.")) return 1;
    if (checkNCerr(nc_put_att_text(A.ncid,A.velocityVID,"units",19,"angstrom/picosecond"),
      "Writing velocities variable units.")) return 1;
    if (checkNCerr(nc_put_att_double(A.ncid,A.velocityVID,"scale_factor",NC_DOUBLE,1, 
                   &velocityScale),
      "Writing velocities scale factor.")) return 1;
  }

  // Cell Spatial
  if (checkNCerr(nc_def_dim(A.ncid,NCCELL_SPATIAL,3,&(A.cell_spatialDID)),
    "Defining cell spatial dimension.")) return 1;
  dimensionID[0]=A.cell_spatialDID;
  if (checkNCerr(nc_def_var(A.ncid,NCCELL_SPATIAL,NC_CHAR,1,dimensionID,&(A.cell_spatialVID)),
    "Defining cell spatial variable.")) return 1;

  // Cell angular
  if (checkNCerr(nc_def_dim(A.ncid,NCLABEL,NCLABELLEN,&(A.labelDID)),
    "Defining label dimension.")) return 1;
  if (checkNCerr(nc_def_dim(A.ncid,NCCELL_ANGULAR,3,&(A.cell_angularDID)),
    "Defining cell angular dimension.")) return 1;
  dimensionID[0] = A.cell_angularDID;
  dimensionID[1] = A.labelDID;
  if (checkNCerr(nc_def_var(A.ncid,NCCELL_ANGULAR,NC_CHAR,2,dimensionID,&(A.cell_angularVID)),
    "Defining cell angular variable.")) return 1;

  // Box Info
  if (box!=NULL) {
    dimensionID[0]=A.cell_spatialDID;
    if (checkNCerr(nc_def_var(A.ncid,NCCELL_LENGTHS,NC_DOUBLE,1,dimensionID,&(A.cellLengthVID)),
      "Defining cell length variable.")) return 1;
    if (checkNCerr(nc_put_att_text(A.ncid,A.cellLengthVID,"units",8,"angstrom"),
      "Writing cell length variable units.")) return 1;
    dimensionID[0]=A.cell_angularDID;
    if (checkNCerr(nc_def_var(A.ncid,NCCELL_ANGLES,NC_DOUBLE,1,dimensionID,&(A.cellAngleVID)),
      "Defining cell angle variable.")) return 1;
    if (checkNCerr(nc_put_att_text(A.ncid,A.cellAngleVID,"units",6,"degree"),
      "Writing cell angle variable units.")) return 1;
  }

  // Set up title
  if (title==NULL) {
    title=(char*) malloc(20*sizeof(char));
    strcpy(title,"AmberNetcdf restart");
  }

  // Attributes
  if (checkNCerr(nc_put_att_text(A.ncid,NC_GLOBAL,"title",strlen(title),title),
    "Writing title.")) return 1;
  if (title!=NULL) free(title);
  if (checkNCerr(nc_put_att_text(A.ncid,NC_GLOBAL,"application",5,"AMBER"),
    "Writing application.")) return 1;
  if (checkNCerr(nc_put_att_text(A.ncid,NC_GLOBAL,"program",11,"AmberNetcdf"),
    "Writing program.")) return 1;
  if (checkNCerr(nc_put_att_text(A.ncid,NC_GLOBAL,"programVersion",3,"1.0"),
    "Writing program version.")) return 1;
  if (checkNCerr(nc_put_att_text(A.ncid,NC_GLOBAL,"Conventions",12,"AMBERRESTART"),
    "Writing conventions.")) return 1;
  if (checkNCerr(nc_put_att_text(A.ncid,NC_GLOBAL,"ConventionVersion",3,"1.0"),
    "Writing conventions version.")) return 1;

  // Replica temperature 
  if (restart_temp>=0) {
    printf("NETCDF: Defining replica temperature in output restart.\n");
    if ( checkNCerr(nc_def_var(A.ncid,NCTEMPERATURE,NC_DOUBLE,0,dimensionID,&(A.TempVID)),
         "Defining replica temperature in netcdf restart.") ) return 1;
    if ( checkNCerr(nc_put_att_text(A.ncid,A.TempVID,"units",6,"kelvin"),
         "Defining replica temperature units in netcdf restart.") ) return 1;
  }

  // Set fill mode
  if (checkNCerr(nc_set_fill(A.ncid, NC_NOFILL, dimensionID),
    "NetCDF setting fill value.")) return 1;

  // End netcdf definitions
  if (checkNCerr(nc_enddef(A.ncid),"NetCDF error on ending definitions."))
    return 1;

  // Specify spatial dimension labels
  start[0] = 0;
  count[0] = 3;

  xyz[0] = 'x'; xyz[1] = 'y'; xyz[2] = 'z';
  if (checkNCerr(nc_put_vara_text(A.ncid, A.spatialVID, start, count, xyz),
    "Error on NetCDF output of spatial VID 'x', 'y' and 'z'")) return 1;

  xyz[0] = 'a'; xyz[1] = 'b'; xyz[2] = 'c';
  if (checkNCerr(nc_put_vara_text(A.ncid, A.cell_spatialVID, start, count, xyz),
    "Error on NetCDF output of cell spatial VID 'a', 'b' and 'c'")) return 1;

  start[0] = 0; start[1] = 0;
  count[0] = 3; count[1] = 5;
  if (checkNCerr(nc_put_vara_text(A.ncid, A.cell_angularVID, start, count, abc),
    "Error on NetCDF output of cell angular VID 'alpha', 'beta ' and 'gamma'"))
    return 1;

  // write coords
  start[0]=0;
  start[1]=0;
  count[0]=A.ncatom;
  count[1]=3;
  if (checkNCerr(nc_put_vara_double(A.ncid,A.coordVID,start,count,X),
      "Netcdf restart Writing coordinates")) return 1;

  // write velocity
  if (Velo!=NULL) {
    if (checkNCerr(nc_put_vara_double(A.ncid,A.velocityVID,start,count,Velo),
        "Netcdf restart writing velocity")) return 1;
  }

  // write box
  if (box!=NULL) { 
    count[0]=3;
    count[1]=0;
    if (checkNCerr(nc_put_vara_double(A.ncid,A.cellLengthVID,start,count,box),
      "Writing cell lengths.")) return 1;
    if (checkNCerr(nc_put_vara_double(A.ncid,A.cellAngleVID,start,count, box+3),
      "Writing cell angles.")) return 1;
  }

  // write time
  if (restart_time>=0) 
    restartTime = restart_time;
  else
    restartTime = 0;
  if (checkNCerr(nc_put_var_double(A.ncid,A.timeVID,&restartTime),
    "Writing restart time.")) return 1;

  // write temperature
  if (restart_temp>=0) {
    if (checkNCerr(nc_put_var_double(A.ncid,A.TempVID,&restart_temp),
      "Writing restart temperature.")) return 1;
  }

  nc_sync(A.ncid); // Necessary? File about to close anyway... 

  checkNCerr(nc_close(A.ncid),"Closing netcdf file.");

  //free(A);
  return 0;
#else
  fprintf(stdout,"Error: NAB compiled without Netcdf support.\n");
  fprintf(stdout,"       Recompile with -DBINTRAJ\n");
  return 1;
#endif
}

// netcdfCreate()
/** Create an Amber netcdf trajectory for the specified #atoms with name 
  * filename. Also add box information if requested.
  */
int netcdfCreate(struct AmberNetcdf *A, char *filename, int natom, int isBox) {
#ifdef BINTRAJ
  int dimensionID[NC_MAX_VAR_DIMS];
  size_t start[3], count[3];
  char *title;
  char xyz[3];
  char abc[15] = { 'a', 'l', 'p', 'h', 'a', 
                   'b', 'e', 't', 'a', ' ',
                   'g', 'a', 'm', 'm', 'a' };
  int debug;

  debug=NCDEBUGVAL;
  title=NULL;

  // Check netcdf structure
  if (A==NULL) {
    fprintf(stderr,"Error: netcdfCreate: AmberNetcdf is NULL\n");
    return 1;
  }
  // Close if already open
  //netcdfClose(A);
  INIT_AmberNetcdf(A); 

  // Create file
  if (checkNCerr(nc_create(filename,NC_64BIT_OFFSET,&(A->ncid)),
    "Creating Netcdf file %s",filename)) return 1;
  if (debug>0) 
    fprintf(stdout,"|   Successfully created Netcdf file %s, ncid %i\n",filename,A->ncid);

  // Frame, Time
  if (checkNCerr(nc_def_dim(A->ncid,NCFRAME,NC_UNLIMITED,&(A->frameDID)),
    "Defining frame dimension.")) return 1;
  dimensionID[0]=A->frameDID;
  if (checkNCerr(nc_def_var(A->ncid,NCTIME,NC_FLOAT,1,dimensionID,&(A->timeVID)),
    "Defining time variable.")) return 1;
  if (checkNCerr(nc_put_att_text(A->ncid,A->timeVID,"units",10,"picosecond"),
    "Writing time VID units.")) return 1;

  // Spatial
  if (checkNCerr(nc_def_dim(A->ncid,NCSPATIAL,3,&(A->spatialDID)),
    "Defining spatial dimension.")) return 1;
  dimensionID[0] = A->spatialDID;
  if (checkNCerr(nc_def_var(A->ncid,NCSPATIAL,NC_CHAR,1,dimensionID,&(A->spatialVID)),
    "Defining spatial variable.")) return 1;

  // Atoms
  if (checkNCerr(nc_def_dim(A->ncid,NCATOM,natom,&(A->atomDID)),
    "Defining atom dimension.")) return 1;
  A->ncatom=natom;
  A->ncatom3 = natom * 3;

  // Coords
  dimensionID[0] = A->frameDID;
  dimensionID[1] = A->atomDID;
  dimensionID[2] = A->spatialDID;
  if (checkNCerr(nc_def_var(A->ncid,NCCOORDS,NC_FLOAT,3,dimensionID,&(A->coordVID)),
    "Defining coordinates variable.")) return 1;
  if (checkNCerr(nc_put_att_text(A->ncid,A->coordVID,"units",8,"angstrom"),
    "Writing coordinates variable units.")) return 1;

  // Cell Spatial
  if (checkNCerr(nc_def_dim(A->ncid,NCCELL_SPATIAL,3,&(A->cell_spatialDID)),
    "Defining cell spatial dimension.")) return 1;
  dimensionID[0]=A->cell_spatialDID;
  if (checkNCerr(nc_def_var(A->ncid,NCCELL_SPATIAL,NC_CHAR,1,dimensionID,&(A->cell_spatialVID)),
    "Defining cell spatial variable.")) return 1;

  // Cell angular
  if (checkNCerr(nc_def_dim(A->ncid,NCLABEL,NCLABELLEN,&(A->labelDID)),
    "Defining label dimension.")) return 1;
  if (checkNCerr(nc_def_dim(A->ncid,NCCELL_ANGULAR,3,&(A->cell_angularDID)),
    "Defining cell angular dimension.")) return 1;
  dimensionID[0] = A->cell_angularDID;
  dimensionID[1] = A->labelDID;
  if (checkNCerr(nc_def_var(A->ncid,NCCELL_ANGULAR,NC_CHAR,2,dimensionID,&(A->cell_angularVID)),
    "Defining cell angular variable.")) return 1;

  // Box Info
  if (isBox>0) {
    dimensionID[0]=A->frameDID;
    dimensionID[1]=A->cell_spatialDID;
    if (checkNCerr(nc_def_var(A->ncid,NCCELL_LENGTHS,NC_DOUBLE,2,dimensionID,&(A->cellLengthVID)),
      "Defining cell length variable.")) return 1;
    if (checkNCerr(nc_put_att_text(A->ncid,A->cellLengthVID,"units",8,"angstrom"),
      "Writing cell length variable units.")) return 1;
    dimensionID[1]=A->cell_angularDID;
    if (checkNCerr(nc_def_var(A->ncid,NCCELL_ANGLES,NC_DOUBLE,2,dimensionID,&(A->cellAngleVID)),
      "Defining cell angle variable.")) return 1;
    if (checkNCerr(nc_put_att_text(A->ncid,A->cellAngleVID,"units",6,"degree"),
      "Writing cell angle variable units.")) return 1;
  }

  // Set up title
  if (title==NULL) {
    title=(char*) malloc(23*sizeof(char));
    strcpy(title,"AmberNetcdf trajectory");
  }

  // Attributes
  if (checkNCerr(nc_put_att_text(A->ncid,NC_GLOBAL,"title",strlen(title),title),
    "Writing title.")) {if (title!=NULL) free(title); return 1;}
  if (title!=NULL) free(title);
  if (checkNCerr(nc_put_att_text(A->ncid,NC_GLOBAL,"application",5,"AMBER"),
    "Writing application.")) return 1;
  if (checkNCerr(nc_put_att_text(A->ncid,NC_GLOBAL,"program",11,"AmberNetcdf"),
    "Writing program.")) return 1;
  if (checkNCerr(nc_put_att_text(A->ncid,NC_GLOBAL,"programVersion",3,"1.0"),
    "Writing program version.")) return 1;
  if (checkNCerr(nc_put_att_text(A->ncid,NC_GLOBAL,"Conventions",5,"AMBER"),
    "Writing conventions.")) return 1;
  if (checkNCerr(nc_put_att_text(A->ncid,NC_GLOBAL,"ConventionVersion",3,"1.0"),
    "Writing conventions version.")) return 1;

  /* Replica temperature 
  if (trajInfo->isREMDTRAJ) {
    fprintf(stdout,"NETCDF: Defining replica temperature in output trajectory.\n");
    dimensionID[0] = NCInfo->frameDID;
    netcdfDefineVariable(NCInfo->ncid, NC_TEMPERATURE, NC_DOUBLE, 1, dimensionID, 
                         &NCInfo->TempVID);
    netcdfPutAttributeText(NCInfo->ncid, NCInfo->TempVarID,"units","kelvin");
  }*/

  // Set fill mode
  if (checkNCerr(nc_set_fill(A->ncid, NC_NOFILL, dimensionID),
    "NetCDF setting fill value.")) return 1;

  // End netcdf definitions
  if (checkNCerr(nc_enddef(A->ncid),"NetCDF error on ending definitions."))
    return 1;

  // Specify spatial dimension labels
  start[0] = 0;
  count[0] = 3;

  xyz[0] = 'x'; xyz[1] = 'y'; xyz[2] = 'z';
  if (checkNCerr(nc_put_vara_text(A->ncid, A->spatialVID, start, count, xyz),
    "Error on NetCDF output of spatial VID 'x', 'y' and 'z'")) return 1;

  xyz[0] = 'a'; xyz[1] = 'b'; xyz[2] = 'c';
  if (checkNCerr(nc_put_vara_text(A->ncid, A->cell_spatialVID, start, count, xyz),
    "Error on NetCDF output of cell spatial VID 'a', 'b' and 'c'")) return 1;

  start[0] = 0; start[1] = 0;
  count[0] = 3; count[1] = 5;
  if (checkNCerr(nc_put_vara_text(A->ncid, A->cell_angularVID, start, count, abc),
    "Error on NetCDF output of cell angular VID 'alpha', 'beta ' and 'gamma'")) 
    return 1;

  // Allocate memory for reading single-precision coordinates
  //A->Coord=(float*) malloc(A->ncatom3 * sizeof(float));
  //if (A->Coord==NULL) {
  //  fprintf(stderr,"Error: Could not allocate memory for netcdf coords.\n");
  //  return 1;
  //}

  return 0;
#else
  fprintf(stdout,"Error: NAB compiled without Netcdf support.\n");
  fprintf(stdout,"       Recompile with -DBINTRAJ\n");
  return 1;
#endif
}

// netcdfGetVelocity
/** If file is a netcdf restart, get the velocities.
  */
int netcdfGetVelocity(struct AmberNetcdf *A, int set, double *V) {
#ifdef BINTRAJ
  int i;
  size_t start[3], count[3];
  float* Velo;

  if (A==NULL) {
    fprintf(stderr,"Error: netcdfGetVelocity: AmberNetcdf structure not allocated.\n");
    return 1;
  }

  if (A->ncid==-1) {
    fprintf(stderr,"Error: netcdfGetVelocity: AmberNetcdf file not open.\n");
    return 1;
  }

  if (V==NULL) {
    fprintf(stderr,"Error: netcdfGetVelocity: Memory for coords not allocated.\n");
    return 1;
  }
  // Set bounds check
  if (set>=A->ncframe || set<0) return 1;

  // ---------- Restart
  if (A->isNCrestart) {
    // Read Velocity
    if (A->velocityVID!=-1) {
      start[0]=0;
      start[1]=0;
      count[0]=A->ncatom;
      count[1]=3;
      if ( checkNCerr(nc_get_vara_double(A->ncid, A->velocityVID, start, count, V),
                      "Getting velocities")!=0 ) return 1;
    } else {
      printf("Warning: No velocities in restart, setting all to 0.0\n");
      for (i=0; i < A->ncatom3; i++)
        V[i] = 0.0;
    }
  // ---------- Trajectory
  } else {
    if (A->velocityVID!=-1) {
      // Allocate space for single-precision velocities
      Velo=(float*) malloc(A->ncatom3 * sizeof(float));
      if (Velo==NULL) {
        fprintf(stderr,"Error: Could not allocate memory for netcdf velocities.\n");
        return 1;
      }
      // Read Velos 
      start[0]=set;
      start[1]=0;
      start[2]=0;
      count[0]=1;
      count[1]=A->ncatom;
      count[2]=3;
      if ( checkNCerr(nc_get_vara_float(A->ncid, A->velocityVID, start, count, Velo),
                      "Getting velocities at frame %i",set)!=0 ) {free(Velo); return 1;}
      // Convert from single precision to double precision
      for (i=0; i < A->ncatom3; i++)
        V[i]=(double) Velo[i];
      free(Velo);
    } else {
      printf("Warning: No velocities in trajectory, setting all to 0.0\n");
      for (i=0; i < A->ncatom3; i++)
        V[i] = 0.0;
    }
  }
  return 0;
#else
  fprintf(stdout,"Error: NAB compiled without Netcdf support.\n");
  fprintf(stdout,"       Recompile with -DBINTRAJ\n");
  return 1;
#endif
}

// netcdfGetFrame()
/** Get the specified frame from amber netcdf trajectory/restart file.
  * Also get the box coords if present and box is not NULL.
  * Also get temperature T if TempVID defined and T is not NULL.
  * Coords are a 1 dimensional array of format X1,Y1,Z1,X2,Y2,Z2,...
  * \return 0 on success, 1 on error.
  */
int netcdfGetFrame(struct AmberNetcdf *A, int set, double *X, double *box) {
#ifdef BINTRAJ
  size_t start[3], count[3];
  int i;
  float *Coord;

  if (A==NULL) {
    fprintf(stderr,"Error: netcdfGetFrame: AmberNetcdf structure not allocated.\n");
    return 1;
  }

  if (A->ncid==-1) {
    fprintf(stderr,"Error: netcdfGetFrame: AmberNetcdf file not open.\n");
    return 1;
  }

  if (X==NULL) {
    fprintf(stderr,"Error: netcdfGetFrame: Memory for coords not allocated.\n");
    return 1;
  }

  // Frame bounds check
  if (set>=A->ncframe || set<0) return 1;

  // ---------- Restart
  if (A->isNCrestart) {
    // Get temperature
    if (A->TempVID!=-1) {
      if ( checkNCerr(nc_get_var_double(A->ncid, A->TempVID, &(A->temp0)), 
                      "Getting restart temperature.")!=0 ) return 1;
      //printf("DEBUG: Restart Temperature %lf\n",A->temp0);
    }
    // Read Coords 
    start[0]=0;
    start[1]=0;
    count[0]=A->ncatom;
    count[1]=3;
    if ( checkNCerr(nc_get_vara_double(A->ncid, A->coordVID, start, count, X),
                    "Getting Coords")!=0 ) return 1;

    // Read box info 
    if (box!=NULL && A->cellLengthVID!=-1) {
      count[0]=3;
      count[1]=0;
      if ( checkNCerr(nc_get_vara_double(A->ncid, A->cellLengthVID, start, count, box),
                      "Getting cell lengths.")!=0 ) return 1;
      if ( checkNCerr(nc_get_vara_double(A->ncid, A->cellAngleVID, start, count, box+3),
                      "Getting cell angles.")!=0 ) return 1;
    }

  // ---------- Trajectory
  } else {
    // Allocate space for single-precision coords
    Coord=(float*) malloc(A->ncatom3 * sizeof(float));
    if (Coord==NULL) {
      fprintf(stderr,"Error: Could not allocate memory for netcdf coords.\n");
      return 1;
    }
    // Read Coords 
    start[0]=set;
    start[1]=0;
    start[2]=0;
    count[0]=1;
    count[1]=A->ncatom;
    count[2]=3;
    if ( checkNCerr(nc_get_vara_float(A->ncid, A->coordVID, start, count, Coord),
                    "Getting frame %i",set)!=0 ) {free(Coord); return 1;}
    // Convert from single precision to double precision
    for (i=0; i < A->ncatom3; i++)
      X[i]=(double) Coord[i];
    free(Coord);
    // Get temperature
    if (A->TempVID!=-1) {
      count[0]=1;
      count[1]=0;
      count[2]=0;
      if ( checkNCerr(nc_get_vara_double(A->ncid, A->TempVID, start,count,&(A->temp0)),
                      "Getting replica temperature.")!=0 ) return 1;
      //fprintf(stderr,"DEBUG: Replica Temperature %lf\n",A->temp0);
    }
    // Read box info 
    if (box!=NULL && A->cellLengthVID!=-1) {
      count[1]=3;
      count[2]=0;
      if ( checkNCerr(nc_get_vara_double(A->ncid, A->cellLengthVID, start, count, box),
                      "Getting cell lengths.")!=0 ) return 1;
      if ( checkNCerr(nc_get_vara_double(A->ncid, A->cellAngleVID, start, count, box+3),
                      "Getting cell angles." )!=0 ) return 1;
    }
  }

  // Set the currentFrame
  A->currentFrame=set; 

  return 0;
#else
  fprintf(stdout,"Error: NAB compiled without Netcdf support.\n");
  fprintf(stdout,"       Recompile with -DBINTRAJ\n");
  return 1;
#endif
}

// netcdfGetNextFrame()
/** Get the frame at currentFrame and increment. Return 0 if no more frames.
  */
int netcdfGetNextFrame(struct AmberNetcdf *A, double *X, double *box) {
#ifdef BINTRAJ
  if ( netcdfGetFrame(A, A->currentFrame, X, box) ) return 0;
  A->currentFrame++;
  return 1;
#else
  fprintf(stdout,"Error: NAB compiled without Netcdf support.\n");
  fprintf(stdout,"       Recompile with -DBINTRAJ\n");
  return 0;
#endif
}

// netcdfWriteFrame()
/** Write coords (and box if specified) to netcdf trajectory file. 
  * \return 0 on success, 1 on error. 
  */
int netcdfWriteFrame(struct AmberNetcdf *A, int set, double *X, double *box) {
#ifdef BINTRAJ
  size_t start[3], count[3];
  int i;
  float *Coord;

  if (A==NULL) {
    fprintf(stderr,"Error: netcdfWriteFrame: AmberNetcdf structure not allocated.\n");
    return 1;
  } 

  if (A->ncid==-1) {
    fprintf(stderr,"Error: netcdfWriteFrame: AmberNetcdf file not open.\n");
    return 1;
  }
    
  if (X==NULL) {
    fprintf(stderr,"Error: netcdfWriteFrame: Memory for coords not allocated.\n");
    return 1;
  }

  // Bounds check
  if (set<0) return 1;

  // ---------- Restart
  if (A->isNCrestart) {
    fprintf(stderr,
      "Error: Called netcdfWriteFrame for netcdf restart; use netcdfWriteRestart instead.\n");
    return 1;
  } 
  // ---------- Trajectory
  // Allocate space for single-precision coords
  Coord=(float*) malloc(A->ncatom3 * sizeof(float));
  if (Coord==NULL) {
    fprintf(stderr,"Error: Could not allocate memory for netcdf coords.\n");
    return 1;
  }
  // Convert to single precision
  for (i=0; i < A->ncatom3; i++)
    Coord[i]=(float) X[i];
  // write coords
  start[0]=set;
  start[1]=0;
  start[2]=0;
  count[0]=1;
  count[1]=A->ncatom;
  count[2]=3;
  if (checkNCerr(nc_put_vara_float(A->ncid,A->coordVID,start,count,Coord),
      "Netcdf Writing frame %i",set)) return 1;
  free(Coord);
  // write box
  if (box!=NULL && A->cellLengthVID!=-1) {
    count[1]=3;
    count[2]=0;
    if (checkNCerr(nc_put_vara_double(A->ncid,A->cellLengthVID,start,count,box),
        "Writing cell lengths.")) return 1;
    if (checkNCerr(nc_put_vara_double(A->ncid,A->cellAngleVID,start,count,box+3),
        "Writing cell angles.")) return 1;
  }
 
  nc_sync(A->ncid); // Necessary after every write??

  // Set the current frame
  A->currentFrame=set;

  return 0;
#else
  fprintf(stdout,"Error: NAB compiled without Netcdf support.\n");
  fprintf(stdout,"       Recompile with -DBINTRAJ\n");
  return 1;
#endif
}  

// netcdfWriteNextFrame()
/** Write coords to currentFrame. Increment currentFrame.
  * Return 1 on error.
  */
int netcdfWriteNextFrame(struct AmberNetcdf *A, double *X, double *box) {
#ifdef BINTRAJ
  if ( netcdfWriteFrame(A, A->currentFrame, X, box) ) return 1;
  A->currentFrame++;
  return 0;
#else
  fprintf(stdout,"Error: NAB compiled without Netcdf support.\n");
  fprintf(stdout,"       Recompile with -DBINTRAJ\n");
  return 1;
#endif
}
  
// netcdfInfo()
/** Print general information about the netcdf file.
  */
int netcdfInfo(struct AmberNetcdf *A) {
#ifdef BINTRAJ
  if (netcdfTaskid()!=0) return 0;
  if (A==NULL) {
    fprintf(stdout,"Error: netcdfInfo: AmberNetcdf is NULL.\n");
    return 1;
  }
  if (A->isNCrestart)
    fprintf(stdout,"  File is a NetCDF AMBER restart");
  else
    fprintf(stdout,"  File is a NetCDF AMBER trajectory");
  if (A->velocityVID!=-1) printf(", with velocity info");
  if (A->cellLengthVID!=-1) fprintf(stdout,", with box info"); 
  if (A->TempVID!=-1) fprintf(stdout,", with temperature info");
  fprintf(stdout,"\n");
  fprintf(stdout,"    %i atoms.\n",A->ncatom);
  fprintf(stdout,"    %i coordinates.\n",A->ncatom3);
  fprintf(stdout,"    %i frames.\n",A->ncframe);

  /*if (debug > 2) {
      if (title != NULL)
        printfone("    title:        \"%s\"\n", p->title);
      if (application != NULL)  
        printfone("    application:  \"%s\"\n", p->application);
      if (program != NULL) 
        printfone("    program:      \"%s\"\n", p->program);
      if (version != NULL) 
        printfone("    version:      \"%s\"\n", p->version);
  }*/
#else
  fprintf(stdout,"Error: NAB compiled without Netcdf support.\n");
  fprintf(stdout,"       Recompile with -DBINTRAJ\n");
  return 1;
#endif
  return 0;
}
