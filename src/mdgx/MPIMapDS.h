#ifndef MPIMapStructs
#define MPIMapStructs

#ifndef PREP_API
#ifdef MPI
#include <mpi.h>
#endif
#endif

#define DSP_MSG_TYPES    16    // The number of direct space message types
#define DSP_CRD_SHARE     0    // Direct space communications get the various
                               //   offsets defined in these macro constants.
#define DSP_EXP_LOC1      3    //   Each function gets a separate offset, or
#define DSP_EXP_LOC2      4    //   perhaps several if it calls multiple pulses
#define DSP_CONSTRAIN     5    //   of communication.
                               //   (0,1,2) = ShareCoordinates,
#define DSP_UPDATE        8    //   (3,4) = ExtraPointLocations,
#define DSP_FRC_MERGE    10    //   5 = Constraints (three variants, skip 6 and
                               //   7), 8 = UpdateCells (two variants, skip 9),
                               //   10 = GatherForcePulse (four variants, skip
                               //   11 through 13), 14 = PostAllMsg in the
#define DSP_RESCALE      14    //   Barostats library (cell grid rescaling)

//-----------------------------------------------------------------------------
// PROCESSCELLGROUP: information about the size and composition of a group of
//                   cells controlled by one process which must communicate
//                   with another process.                  
//-----------------------------------------------------------------------------
struct ProcessCellGroup {
  int ncell;       // The number of cells in this group
  int partner;     // The origin or destination process with which this group
                   //   communicates as referenced in MPI_COMM_WORLD
  int BaseID;      // The message base, as determined by the rank of the
                   //   sending process in the cell grid communicator and the
                   //   direction of the message.  The message base is further
                   //   offset by increments for each function that relies on
                   //   such messages.
  int* cellpt;     // Pointer to the segment of a list of all cells controlled
                   //   by this process; the segment lists cells in this group
};
typedef struct ProcessCellGroup pcgrp;

//-----------------------------------------------------------------------------
// ATOMSHAREPLAN: details the cell-to-cell communication that must occur in
//                order to complete the ShareCoordinates function in the 
//                CellManip library.                                
//-----------------------------------------------------------------------------
struct AtomSharePlan {
  int nsend;       // The number of processes to which this process will
                   //   send information, including itself (this process will
                   //   send ncomm-1 messages over MPI while executing this
                   //   plan)
  int nrecv;       // The number of receives to expect / post
  int* slist;      // The list of all cells controlled by this process (the
                   //   "send from" list)
  int* rlist;      // The list of all cell not controlled by this process but
                   //   which send information to this process (the "receive
                   //   from" list)
  pcgrp* send;     // Detailed list of sends from this process
  pcgrp* recv;     // List of receives for this process
  pcgrp selfrecv;  // The special "self receive" of cells that this plan
                   //   delivered to the same process
};
typedef struct AtomSharePlan ashr;

//-----------------------------------------------------------------------------
// DIRECTCOMMPLAN: holds all plans needed for direct-space communication by
//                 this process.                                     
//-----------------------------------------------------------------------------
struct DirectCommPlan {
  ashr mvshr[3];   // Atom sharing communication
  ashr frcmg[3];   // Force merger communication
};
typedef struct DirectCommPlan dcplan;

//-----------------------------------------------------------------------------
// GRIDSPLICE: details the information present in a grid struct imported from
//             another process.                                    
//-----------------------------------------------------------------------------
struct GridSplice {
  int npag;        // The number of pages controlled, one X value and all Y and
                   //   Z coordinates at that X value
  int ncol;        // The number of pencils / columns controlled, one X and one
                   //   Y value and all Z values
  int npc;         // The number of partial columns controlled, one X and one Y
                   //   value, Z values from Zo to Zf
  int* pagl;       // The list of pages controlled by a process
  int* coll;       // The list of columns controlled by a process, format
                   //   [X1][Y1][X2][Y2] ...
  int* pcl;        // The list of column pieces controlled by a process, format
                   //   [X1][Y1][Z1o][Z1f][X2][Y2][Z2o][Z2f] ...
};
typedef struct GridSplice gsplc;

#endif
