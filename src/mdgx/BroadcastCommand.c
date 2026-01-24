#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef MPI
#include <mpi.h>
#endif
#include "VirtualSites.h"
#include "BroadcastCommand.h"
#include "Matrix.h"

#include "CrdManipDS.h"
#include "pmeDirectDS.h"
#include "pmeRecipDS.h"
#include "TrajectoryDS.h"
#include "ChargeFitDS.h"
#include "CellManipDS.h"
#include "IPolQDS.h"

#ifdef MPI
//-----------------------------------------------------------------------------
// EncodeInteger: encode an integer into the double-precision array for 
//                broadcast.  Integers are encoded with slightly exaggerated
//                values to ensure that, upon conversion of the
//                double-precision reals back to integers the original values
//                will remain intact.                   
//                                                                      
// Arguments:                                                           
//   dscratch:   the array of double-precision reals used to serialize  
//               a struct or other collection of data                   
//   ipos:       position in scratch receiving the next element         
//   val:        the value to write into scratch                        
//-----------------------------------------------------------------------------
static void EncodeInteger(double* dscratch, int ipos, int val)
{
  dscratch[ipos] = (val >= 0) ? val + 1.0e-1 : val - 1.0e-1;
}

//-----------------------------------------------------------------------------
// EncodeLongInt: encode an integer into the double-precision array for 
//                broadcast.  Integers are encoded with slightly exaggerated
//                values to ensure that, upon conversion of the
//                double-precision reals back to integers the original values
//                will remain intact.                   
//                                                                      
// Arguments:                                                           
//   dscratch    the array of double-precision reals used to serialize a struct
//               or other collection of data                   
//   ipos:       position in scratch receiving the next element         
//   val:        the value to write into scratch                        
//-----------------------------------------------------------------------------
static void EncodeLongInt(double* dscratch, int ipos, long int val)
{
  dscratch[ipos] = val + 1.0e-1;
}

//-----------------------------------------------------------------------------
// EncodeLongLongInt: encode a long long int (signed) into two numbers of the 
//                    double-precision array for broadcast.      
//                                                                      
// Arguments:                                                           
//   dscratch:   the array of double-precision reals used to serialize a struct
//               or other collection of data                   
//   ipos:       position in scratch receiving the next element (modified and
//               returned +2 by this function)            
//   val:        the value to write into scratch                        
//-----------------------------------------------------------------------------
static void EncodeLongLongInt(double* dscratch, int *ipos, long long int val)
{
  int val1, val2;
  long long int valtmp2, valtmp1, uval;

  uval = (val < 0) ? -val : val;
  val1 = uval >> 32;
  valtmp1 = val1;
  valtmp2 = uval - (valtmp1 << 32);
  val2 = valtmp2;
  if (val < 0) {
    dscratch[*ipos] = -val1 - 1.0e-1;
    dscratch[*ipos+1] = -val2 - 1.0e-1;
  }
  else {
    dscratch[*ipos] =  val1 + 1.0e-1;
    dscratch[*ipos+1] = val2 + 1.0e-1;
  }
  *ipos += 2;
}

//-----------------------------------------------------------------------------
// DecodeLongLongInt: decode a long long int from two numbers of a      
//                    double-precision array.                           
//                                                                      
// Arguments:                                                           
//   dscratch:   the array of double-precision reals used to serialize  
//               a struct or other collection of data                   
//   ipos:       position in scratch receiving the next element         
//               (modified and returned +2 by this function)            
//-----------------------------------------------------------------------------
static long long int DecodeLongLongInt(double* dscratch, int ipos)
{
  long long int tval;

  tval = dscratch[ipos];
  tval = tval * 4294967296 + dscratch[ipos+1];

  return tval;
}

//-----------------------------------------------------------------------------
// EncodeChar: encode a character string into a much larger one.        
//                                                                      
// Arguments:                                                           
//   cscratch:   the character string to be broadcast                   
//   ipos:       the starting position for encoding cword into cscratch 
//   cword:      the character string to encode                         
//   slen:       the amount of cscratch to be allocated to cword        
//-----------------------------------------------------------------------------
static void EncodeString(char* cscratch, int *ipos, char* cword, int slen)
{
  int i, ii;

  const int wlen = strlen(cword);
  ii = 0;
  const int llim = *ipos;
  const int hlim = llim + slen;
  for (i = llim; i < hlim; i++) {
    cscratch[i] = (ii < wlen) ? cword[ii] : '\0';
    ii++;
  }
  *ipos += slen;
}

//-----------------------------------------------------------------------------
// DecodeString: extract a string from a much longer character array.   
//                                                                      
// Arguments:                                                           
//   dstr:       the string to extract information into                 
//   cscratch:   the array of characters used to serialize a struct or  
//               other collection of data                               
//   ipos:       position in scratch receiving the next element         
//               (modified and returned by this function)               
//   ninc:       the degree to which ipos shall be incremented          
//-----------------------------------------------------------------------------
static void DecodeString(char* dstr, char* cscratch, int *ipos, int ninc)
{
  strncpy(dstr, &cscratch[*ipos], ninc);
  *ipos += ninc;
}

//-----------------------------------------------------------------------------
// EncodeReccon: encode fields of a reccon struct, to the extent that they have
//               been assigned, for broadcast.                
//                                                                      
// Arguments:                                                           
//   rcinp:     the reciprocal space control information struct         
//   dscratch:  array of double precision reals to receive elemental types from
//              rcinp                                        
//   edir:      direction of the encoding (0: encode rcinp -> scratch,  
//              1: encode scratch -> rcinp)                             
//-----------------------------------------------------------------------------
static int EncodeReccon(reccon *rcinp, double* dscratch, int edir)
{
  int i, ipos, nlev;

  if (edir == 0) {
    nlev = rcinp->nlev;
    ipos = 0;
    EncodeInteger(dscratch, ipos++, nlev);
    for (i = 0; i < 3; i++) EncodeInteger(dscratch, ipos++, rcinp->ordr[i]);
    dscratch[ipos++] = rcinp->S;
    if (nlev > 1) {
      EncodeInteger(dscratch, ipos++, rcinp->nslab);
      EncodeInteger(dscratch, ipos++, rcinp->ggordr);
      for (i = 0; i < 4; i++) EncodeInteger(dscratch, ipos++, rcinp->PadYZ[i]);
      for (i = 0; i < 4; i++) dscratch[ipos++] = rcinp->cfac[i];
    }
    for (i = 0; i < 3; i++) EncodeInteger(dscratch, ipos++, rcinp->ng[i]);
    return ipos;
  }
  else {
    ipos = 0;
    nlev = dscratch[ipos++];
    rcinp->nlev = nlev;
    for (i = 0; i < 3; i++) rcinp->ordr[i] = dscratch[ipos++];
    rcinp->S = dscratch[ipos++];
    if (nlev > 1) {
      rcinp->nslab = dscratch[ipos++];
      rcinp->ggordr = dscratch[ipos++];
      for (i = 0; i < 4; i++) rcinp->PadYZ[i] = dscratch[ipos++];
      for (i = 0; i < 4; i++) rcinp->cfac[i] = dscratch[ipos++];
    }
    for (i = 0; i < 3; i++) rcinp->ng[i] = dscratch[ipos++];
    return 0;
  }
}

//-----------------------------------------------------------------------------
// EncodeTrajcon: encode fields of a trajcon struct, to the extent that they
//                have been assigned, for broadcast.               
//                                                                      
// Arguments:                                                           
//   tj:        the trajectory control information struct               
//   dscratch:  array of double precision reals to receive elemental types
//              from the trajcon struct tj                        
//   cscratch:  character string to receive text and file names from tj 
//   edir:      direction of the encoding (0: encode tj -> scratch,     
//              1: encode scratch -> tj)                                
//   tlen:      array containing the number of elements commited to the 
//              double-precision real and character buffers             
//-----------------------------------------------------------------------------
static void EncodeTrajcon(trajcon *tj, double* dscratch, char* cscratch,
			  int edir, int* tlen)
{
  int i, ipos;
  int nipcname, nrstbase, ntrjbase, nvelbase, nfrcbase, ninprow, ninpcol;

  if (edir == 0) {

    // Serialize the numerical information
    ipos = 0;
    EncodeInteger(dscratch, ipos++, tj->mode);
    EncodeInteger(dscratch, ipos++, tj->nsys);
    EncodeInteger(dscratch, ipos++, tj->ntop);
    EncodeLongLongInt(dscratch, &ipos, tj->nstep);
    EncodeLongLongInt(dscratch, &ipos, tj->nfistep);
    EncodeInteger(dscratch, ipos++, tj->RemoveMomentum);
    EncodeInteger(dscratch, ipos++, tj->irest);
    EncodeInteger(dscratch, ipos++, tj->ntwr);
    EncodeInteger(dscratch, ipos++, tj->ntwx);
    EncodeInteger(dscratch, ipos++, tj->ntwv);
    EncodeInteger(dscratch, ipos++, tj->ntwf);
    EncodeInteger(dscratch, ipos++, tj->ntpr);
    EncodeInteger(dscratch, ipos++, tj->ioutfm);
    EncodeInteger(dscratch, ipos++, tj->OverwriteOutput);
    EncodeInteger(dscratch, ipos++, tj->Reckless);
    EncodeInteger(dscratch, ipos++, tj->igseed);
    EncodeInteger(dscratch, ipos++, tj->MaxRattleIter);
    EncodeInteger(dscratch, ipos++, tj->topchk);
    EncodeInteger(dscratch, ipos++, tj->ntt);
    EncodeInteger(dscratch, ipos++, tj->ntp);
    EncodeInteger(dscratch, ipos++, tj->barostat);
    EncodeInteger(dscratch, ipos++, tj->vrand);
    EncodeInteger(dscratch, ipos++, tj->MCBarostatFreq);
    EncodeInteger(dscratch, ipos++, tj->TI);
    EncodeInteger(dscratch, ipos++, tj->mxorder);
    EncodeInteger(dscratch, ipos++, tj->nsynch);
    EncodeLongInt(dscratch, ipos++, tj->rndcon);
    dscratch[ipos++] = tj->mxA;
    dscratch[ipos++] = tj->mxB;
    dscratch[ipos++] = tj->dmxA;
    dscratch[ipos++] = tj->dmxB;
    dscratch[ipos++] = tj->starttime;
    dscratch[ipos++] = tj->currtime;
    dscratch[ipos++] = tj->dt;
    dscratch[ipos++] = tj->rattletol;
    dscratch[ipos++] = tj->Ttarget;
    dscratch[ipos++] = tj->Tinit;
    dscratch[ipos++] = tj->BerendsenTCoupl;
    dscratch[ipos++] = tj->BerendsenPTime;
    dscratch[ipos++] = tj->BerendsenPCoupl;
    dscratch[ipos++] = tj->lnth.gamma_ln;
    for (i = 0; i < 3; i++) dscratch[ipos++] = tj->MCBarostatFac[i];
    dscratch[ipos++] = tj->mcdVmax;
    dscratch[ipos++] = tj->Ptarget;
    dscratch[ipos++] = tj->lambda;
    dscratch[ipos++] = tj->npth.pmass;
    dscratch[ipos++] = tj->npth.qmass;
    dscratch[ipos++] = tj->npth.chi;
    dscratch[ipos++] = tj->npth.eta;
    dscratch[ipos++] = tj->npth.sigma;
    dscratch[ipos++] = tj->npth.TauT;
    dscratch[ipos++] = tj->npth.TauP;
    EncodeInteger(dscratch, ipos++, tj->DMPcrd);
    EncodeInteger(dscratch, ipos++, tj->DMPbond);
    EncodeInteger(dscratch, ipos++, tj->DMPangl);
    EncodeInteger(dscratch, ipos++, tj->DMPdihe);
    EncodeInteger(dscratch, ipos++, tj->DMPdelec);
    EncodeInteger(dscratch, ipos++, tj->DMPrelec);
    EncodeInteger(dscratch, ipos++, tj->DMPvdw);
    EncodeInteger(dscratch, ipos++, tj->DMPall);
    EncodeInteger(dscratch, ipos++, tj->Leash.active);
    EncodeInteger(dscratch, ipos++, tj->Leash.usegrid);
    EncodeInteger(dscratch, ipos++, tj->Leash.usebelly);
    dscratch[ipos++] = tj->Leash.GridScale;
    EncodeInteger(dscratch, ipos++, tj->ipcname.row);
    EncodeInteger(dscratch, ipos++, tj->rstbase.row);
    EncodeInteger(dscratch, ipos++, tj->trjbase.row);
    EncodeInteger(dscratch, ipos++, tj->velbase.row);
    EncodeInteger(dscratch, ipos++, tj->frcbase.row);
    EncodeInteger(dscratch, ipos++, tj->inptext.row);
    EncodeInteger(dscratch, ipos++, tj->inptext.col);
    tlen[0] = ipos;

    // Serialize the character information
    ipos = 0;
    EncodeString(cscratch, &ipos, tj->DMPvar, 32);
    EncodeString(cscratch, &ipos, tj->inpname, MAXNAME);
    EncodeString(cscratch, &ipos, tj->dumpname, MAXNAME);
    EncodeString(cscratch, &ipos, tj->outbase, MAXNAME);
    EncodeString(cscratch, &ipos, tj->outsuff, 32);
    EncodeString(cscratch, &ipos, tj->Leash.GridFile, MAXNAME);
    EncodeString(cscratch, &ipos, tj->Leash.GridDefsFile, MAXNAME);
    EncodeString(cscratch, &ipos, tj->Leash.BellyMask, MAXLINE);
    EncodeString(cscratch, &ipos, tj->Leash.FrozenMask, MAXLINE);
    for (i = 0; i < tj->ipcname.row; i++) {
      EncodeString(cscratch, &ipos, tj->ipcname.map[i], MAXNAME);
    }
    for (i = 0; i < tj->rstbase.row; i++) {
      EncodeString(cscratch, &ipos, tj->rstbase.map[i], MAXNAME);
      EncodeString(cscratch, &ipos, tj->rstsuff.map[i], 32);
    }
    for (i = 0; i < tj->trjbase.row; i++) {
      EncodeString(cscratch, &ipos, tj->trjbase.map[i], MAXNAME);
      EncodeString(cscratch, &ipos, tj->trjsuff.map[i], 32);
    }
    for (i = 0; i < tj->velbase.row; i++) {
      EncodeString(cscratch, &ipos, tj->velbase.map[i], MAXNAME);
      EncodeString(cscratch, &ipos, tj->velsuff.map[i], 32);
    }
    for (i = 0; i < tj->frcbase.row; i++) {
      EncodeString(cscratch, &ipos, tj->frcbase.map[i], MAXNAME);
      EncodeString(cscratch, &ipos, tj->frcsuff.map[i], 32);
    }
    for (i = 0; i < tj->inptext.row*tj->inptext.col; i++) {
      cscratch[ipos] = tj->inptext.data[i];
      ipos++;
    }
    tlen[1] = ipos;
  }
  else {

    // Deconstruct the numerical information
    ipos = 0;
    tj->mode = dscratch[ipos++];
    tj->nsys = dscratch[ipos++];
    tj->ntop = dscratch[ipos++];
    tj->nstep = DecodeLongLongInt(dscratch, ipos);
    tj->nfistep = DecodeLongLongInt(dscratch, ipos+2);
    ipos += 4;
    tj->RemoveMomentum = dscratch[ipos++];
    tj->irest = dscratch[ipos++];
    tj->ntwr = dscratch[ipos++];
    tj->ntwx = dscratch[ipos++];
    tj->ntwv = dscratch[ipos++];
    tj->ntwf = dscratch[ipos++];
    tj->ntpr = dscratch[ipos++];
    tj->ioutfm = dscratch[ipos++];
    tj->OverwriteOutput = dscratch[ipos++];
    tj->Reckless = dscratch[ipos++];
    tj->igseed = dscratch[ipos++];
    tj->MaxRattleIter = dscratch[ipos++];
    tj->topchk = dscratch[ipos++];
    tj->ntt = dscratch[ipos++];
    tj->ntp = dscratch[ipos++];
    tj->barostat = dscratch[ipos++];
    tj->vrand = dscratch[ipos++];
    tj->MCBarostatFreq = dscratch[ipos++];
    tj->TI = dscratch[ipos++];
    tj->mxorder = dscratch[ipos++];
    tj->nsynch = dscratch[ipos++];
    tj->rndcon = dscratch[ipos++];
    tj->mxA = dscratch[ipos++];
    tj->mxB = dscratch[ipos++];
    tj->dmxA = dscratch[ipos++];
    tj->dmxB = dscratch[ipos++];
    tj->starttime = dscratch[ipos++];
    tj->currtime = dscratch[ipos++];
    tj->dt = dscratch[ipos++];
    tj->rattletol = dscratch[ipos++];
    tj->Ttarget = dscratch[ipos++];
    tj->Tinit = dscratch[ipos++];
    tj->BerendsenTCoupl = dscratch[ipos++];
    tj->BerendsenPTime = dscratch[ipos++];
    tj->BerendsenPCoupl = dscratch[ipos++];
    tj->lnth.gamma_ln = dscratch[ipos++];
    for (i = 0; i < 3; i++) tj->MCBarostatFac[i] = dscratch[ipos++];
    tj->mcdVmax = dscratch[ipos++];
    tj->Ptarget = dscratch[ipos++];
    tj->lambda = dscratch[ipos++];
    tj->npth.pmass = dscratch[ipos++];
    tj->npth.qmass = dscratch[ipos++];
    tj->npth.chi = dscratch[ipos++];
    tj->npth.eta = dscratch[ipos++];
    tj->npth.sigma = dscratch[ipos++];
    tj->npth.TauT = dscratch[ipos++];
    tj->npth.TauP = dscratch[ipos++];
    tj->DMPcrd = dscratch[ipos++];
    tj->DMPbond = dscratch[ipos++];
    tj->DMPangl = dscratch[ipos++];
    tj->DMPdihe = dscratch[ipos++];
    tj->DMPdelec = dscratch[ipos++];
    tj->DMPrelec = dscratch[ipos++];
    tj->DMPvdw = dscratch[ipos++];
    tj->DMPall = dscratch[ipos++];
    tj->Leash.active = dscratch[ipos++];
    tj->Leash.usegrid = dscratch[ipos++];
    tj->Leash.usebelly = dscratch[ipos++];
    tj->Leash.GridScale = dscratch[ipos++];

    // Reallocate character matrices
    nipcname = dscratch[ipos++];
    nrstbase = dscratch[ipos++];
    ntrjbase = dscratch[ipos++];
    nvelbase = dscratch[ipos++];
    nfrcbase = dscratch[ipos++];
    ninprow = dscratch[ipos++];
    ninpcol = dscratch[ipos++];
    tj->ipcname = ReallocCmat(&tj->ipcname, nipcname, MAXNAME);
    tj->rstbase = ReallocCmat(&tj->rstbase, nrstbase, MAXNAME);
    tj->trjbase = ReallocCmat(&tj->trjbase, ntrjbase, MAXNAME);
    tj->velbase = ReallocCmat(&tj->velbase, nvelbase, MAXNAME);
    tj->frcbase = ReallocCmat(&tj->frcbase, nfrcbase, MAXNAME);
    tj->rstsuff = ReallocCmat(&tj->rstsuff, nrstbase, 32);
    tj->trjsuff = ReallocCmat(&tj->trjsuff, ntrjbase, 32);
    tj->velsuff = ReallocCmat(&tj->velsuff, nvelbase, 32);
    tj->frcsuff = ReallocCmat(&tj->frcsuff, nfrcbase, 32);
    tj->inptext = CreateCmat(ninprow, ninpcol);
    tj->Leash.GridFile = (char*)malloc(MAXNAME*sizeof(char));
    tj->Leash.GridDefsFile = (char*)malloc(MAXNAME*sizeof(char));
    tj->Leash.BellyMask = (char*)malloc(MAXLINE*sizeof(char));
    tj->Leash.FrozenMask = (char*)malloc(MAXLINE*sizeof(char));

    // Deconstruct the character information
    ipos = 0;
    DecodeString(tj->DMPvar, cscratch, &ipos, 32);
    DecodeString(tj->inpname, cscratch, &ipos, MAXNAME);
    DecodeString(tj->dumpname, cscratch, &ipos, MAXNAME);
    DecodeString(tj->outbase, cscratch, &ipos, MAXNAME);
    DecodeString(tj->outsuff, cscratch, &ipos, 32);
    DecodeString(tj->Leash.GridFile, cscratch, &ipos, MAXNAME);
    DecodeString(tj->Leash.GridDefsFile, cscratch, &ipos, MAXNAME);
    DecodeString(tj->Leash.BellyMask, cscratch, &ipos, MAXLINE);
    DecodeString(tj->Leash.FrozenMask, cscratch, &ipos, MAXLINE);
    for (i = 0; i < tj->ipcname.row; i++) {
      DecodeString(tj->ipcname.map[i], cscratch, &ipos, MAXNAME);
    }
    for (i = 0; i < tj->rstbase.row; i++) {
      DecodeString(tj->rstbase.map[i], cscratch, &ipos, MAXNAME);
      DecodeString(tj->rstsuff.map[i], cscratch, &ipos, 32);
    }
    for (i = 0; i < tj->trjbase.row; i++) {
      DecodeString(tj->trjbase.map[i], cscratch, &ipos, MAXNAME);
      DecodeString(tj->trjsuff.map[i], cscratch, &ipos, 32);
    }
    for (i = 0; i < tj->velbase.row; i++) {
      DecodeString(tj->velbase.map[i], cscratch, &ipos, MAXNAME);
      DecodeString(tj->velsuff.map[i], cscratch, &ipos, 32);
    }
    for (i = 0; i < tj->frcbase.row; i++) {
      DecodeString(tj->frcbase.map[i], cscratch, &ipos, MAXNAME);
      DecodeString(tj->frcsuff.map[i], cscratch, &ipos, 32);
    }
    for (i = 0; i < tj->inptext.row*tj->inptext.col; i++) {
      tj->inptext.data[i] = cscratch[ipos];
      ipos++;
    }
  }
}

//-----------------------------------------------------------------------------
// EncodePrmtop: encode the fields of a prmtop struct NOT read directly from
//               the topology file (but, rather, as part of mdin or command
//               line directives) for broadcast.                
//                                                                      
// Arguments:                                                           
//   tp:        the topology                                            
//   dscratch:  array of double precision reals to receive elemental types from
//              the trajcon struct tj                        
//   cscratch:  character string to receive text and file names from tj 
//   edir:      direction of the encoding (0: encode tp -> scratch,     
//              1: encode scratch -> tp)                                
//   tlen:      array containing the number of elements commited to the 
//              double-precision real and character buffers             
//-----------------------------------------------------------------------------
static void EncodePrmtop(prmtop *tp, double* dscratch, char* cscratch,
			 int edir, int* tlen)
{
  int ipos;

  if (edir == 0) {

    // Commit numerical data
    ipos = 0;
    EncodeInteger(dscratch, ipos++, tp->ljbuck);
    EncodeInteger(dscratch, ipos++, tp->qshape);
    EncodeInteger(dscratch, ipos++, tp->rattle);
    EncodeInteger(dscratch, ipos++, tp->settle);
    dscratch[ipos++] = tp->lj14fac;
    dscratch[ipos++] = tp->elec14fac;
    tlen[0] = ipos;

    // Commit character data
    ipos = 0;
    EncodeString(cscratch, &ipos, tp->norattlemask, MAXNAME);
    EncodeString(cscratch, &ipos, tp->rattlemask, MAXNAME);
    EncodeString(cscratch, &ipos, tp->WaterName, 8);
    EncodeString(cscratch, &ipos, tp->source, MAXNAME);
    EncodeString(cscratch, &ipos, tp->eprulesource, MAXNAME);
    tlen[1] = ipos;
  }
  else {

    // Deconstruct numerical information
    ipos = 0;
    tp->ljbuck = dscratch[ipos++];
    tp->qshape = dscratch[ipos++];
    tp->rattle = dscratch[ipos++];
    tp->settle = dscratch[ipos++];
    tp->lj14fac = dscratch[ipos++];
    tp->elec14fac = dscratch[ipos++];

    // Deconstruct character information
    ipos = 0;
    tp->norattlemask = (char*)malloc(MAXNAME*sizeof(char));
    DecodeString(tp->norattlemask, cscratch, &ipos, MAXNAME);
    tp->rattlemask = (char*)malloc(MAXNAME*sizeof(char));
    DecodeString(tp->rattlemask, cscratch, &ipos, MAXNAME);
    DecodeString(tp->WaterName, cscratch, &ipos, 8);
    DecodeString(tp->source, cscratch, &ipos, MAXNAME);
    DecodeString(tp->eprulesource, cscratch, &ipos, MAXNAME);
  }
}

//-----------------------------------------------------------------------------
// EncodeIpqcon: encode the information for IPolQ calculations.         
//                                                                      
// Arguments:                                                           
//   ipqinp:  IPolQ control data                                        
//   dscratch:  array of double precision reals to receive elemental types from
//              the trajcon struct tj                        
//   cscratch:  character string to receive text and file names from tj 
//   edir:      direction of the encoding (0: encode tp -> scratch,     
//              1: encode scratch -> tp)                                
//   tlen:      array containing the number of elements commited to the 
//              double-precision real and character buffers             
//-----------------------------------------------------------------------------
static void EncodeIpqcon(ipqcon *ipqinp, double* dscratch, char* cscratch,
                         int edir, int* tlen)
{
  int i, ipos;

  if (edir == 0) {

    // Commit numerical data
    ipos = 0;
    EncodeInteger(dscratch, ipos++, ipqinp->ntqs);
    EncodeInteger(dscratch, ipos++, ipqinp->nqframe);
    EncodeInteger(dscratch, ipos++, ipqinp->nQshell);
    EncodeInteger(dscratch, ipos++, ipqinp->nVshell);
    EncodeInteger(dscratch, ipos++, ipqinp->neqstep);
    EncodeInteger(dscratch, ipos++, ipqinp->nQphpt);
    EncodeInteger(dscratch, ipos++, ipqinp->nVphpt);
    EncodeInteger(dscratch, ipos++, ipqinp->nblock);
    EncodeInteger(dscratch, ipos++, ipqinp->verbose);
    EncodeInteger(dscratch, ipos++, ipqinp->retqmout);
    EncodeInteger(dscratch, ipos++, ipqinp->prepcalls.row);
    EncodeInteger(dscratch, ipos++, ipqinp->postcalls.row);
    EncodeInteger(dscratch, ipos++, ipqinp->nqmod);
    EncodeInteger(dscratch, ipos++, ipqinp->excitation);
    dscratch[ipos++] = ipqinp->econv;
    dscratch[ipos++] = ipqinp->minqfac;
    for (i = 0; i < 4; i++) {
      dscratch[ipos++] = ipqinp->Qshell[i];
      dscratch[ipos++] = ipqinp->Vshell[i];
    }
    for (i = 0; i < ipqinp->nqmod; i++) {
      dscratch[ipos++] = ipqinp->QModVal[i];
    }
    tlen[0] = ipos;

    // Commit character data
    ipos = 0;
    EncodeString(cscratch, &ipos, ipqinp->qmprog, MAXNAME);
    EncodeString(cscratch, &ipos, ipqinp->qmpath, MAXNAME);
    EncodeString(cscratch, &ipos, ipqinp->inpfile, MAXNAME);
    EncodeString(cscratch, &ipos, ipqinp->outfile, MAXNAME);
    EncodeString(cscratch, &ipos, ipqinp->ptqfile, MAXNAME);
    EncodeString(cscratch, &ipos, ipqinp->finfile, MAXNAME);
    EncodeString(cscratch, &ipos, ipqinp->grdfile, MAXNAME);
    EncodeString(cscratch, &ipos, ipqinp->qmmeth, MAXNAME);
    EncodeString(cscratch, &ipos, ipqinp->basis, MAXNAME);
    EncodeString(cscratch, &ipos, ipqinp->scrdir, MAXNAME);
    for (i = 0; i < ipqinp->prepcalls.row; i++) {
      EncodeString(cscratch, &ipos, ipqinp->prepcalls.map[i], MAXLINE);
    }
    for (i = 0; i < ipqinp->postcalls.row; i++) {
      EncodeString(cscratch, &ipos, ipqinp->postcalls.map[i], MAXLINE);
    }
    for (i = 0; i < ipqinp->nqmod; i++) {
      EncodeString(cscratch, &ipos, ipqinp->QModMask.map[i], MAXNAME);
    }
    tlen[1] = ipos;
  }
  else {

    // Deconstruct numerical data
    ipos = 0;
    ipqinp->ntqs = dscratch[ipos++];
    ipqinp->nqframe = dscratch[ipos++];
    ipqinp->nQshell = dscratch[ipos++];
    ipqinp->nVshell = dscratch[ipos++];
    ipqinp->neqstep = dscratch[ipos++];
    ipqinp->nQphpt = dscratch[ipos++];
    ipqinp->nVphpt = dscratch[ipos++];
    ipqinp->nblock = dscratch[ipos++];
    ipqinp->verbose = dscratch[ipos++];
    ipqinp->retqmout = dscratch[ipos++];
    ipqinp->prepcalls.row = dscratch[ipos++];
    ipqinp->postcalls.row = dscratch[ipos++];
    ipqinp->nqmod = dscratch[ipos++];
    ipqinp->excitation = dscratch[ipos++];
    ipqinp->econv = dscratch[ipos++];
    ipqinp->minqfac = dscratch[ipos++];
    for (i = 0; i < 4; i++) {
      ipqinp->Qshell[i] = dscratch[ipos++];
      ipqinp->Vshell[i] = dscratch[ipos++];
    }
    ipqinp->QModVal = (double*)malloc(ipqinp->nqmod*sizeof(double));
    ipqinp->QModMask = CreateCmat(ipqinp->nqmod, MAXNAME);
    for (i = 0; i < ipqinp->nqmod; i++) {
      ipqinp->QModVal[i] = dscratch[ipos++];
    }

    // Deconstruct character information
    ipos = 0;
    ipqinp->qmprog = (char*)malloc(MAXNAME*sizeof(char));
    ipqinp->qmpath = (char*)malloc(MAXNAME*sizeof(char));
    ipqinp->inpfile = (char*)malloc(MAXNAME*sizeof(char));
    ipqinp->outfile = (char*)malloc(MAXNAME*sizeof(char));
    ipqinp->ptqfile = (char*)malloc(MAXNAME*sizeof(char));
    ipqinp->finfile = (char*)malloc(MAXNAME*sizeof(char));
    ipqinp->grdfile = (char*)malloc(MAXNAME*sizeof(char));
    ipqinp->qmmeth = (char*)malloc(MAXNAME*sizeof(char));
    ipqinp->basis = (char*)malloc(MAXNAME*sizeof(char));
    ipqinp->scrdir = (char*)malloc(MAXNAME*sizeof(char));
    ipqinp->prepcalls = CreateCmat(ipqinp->prepcalls.row, MAXLINE);
    ipqinp->postcalls = CreateCmat(ipqinp->postcalls.row, MAXLINE);
    DecodeString(ipqinp->qmprog, cscratch, &ipos, MAXNAME);
    DecodeString(ipqinp->qmpath, cscratch, &ipos, MAXNAME);
    DecodeString(ipqinp->inpfile, cscratch, &ipos, MAXNAME);
    DecodeString(ipqinp->outfile, cscratch, &ipos, MAXNAME);
    DecodeString(ipqinp->ptqfile, cscratch, &ipos, MAXNAME);
    DecodeString(ipqinp->finfile, cscratch, &ipos, MAXNAME);
    DecodeString(ipqinp->grdfile, cscratch, &ipos, MAXNAME);
    DecodeString(ipqinp->qmmeth, cscratch, &ipos, MAXNAME);
    DecodeString(ipqinp->basis, cscratch, &ipos, MAXNAME);
    DecodeString(ipqinp->scrdir, cscratch, &ipos, MAXNAME);
    for (i = 0; i < ipqinp->prepcalls.row; i++) {
      DecodeString(ipqinp->prepcalls.map[i], cscratch, &ipos, MAXLINE);
    }
    for (i = 0; i < ipqinp->postcalls.row; i++) {
      DecodeString(ipqinp->postcalls.map[i], cscratch, &ipos, MAXLINE);
    }
    for (i = 0; i < ipqinp->nqmod; i++) {
      DecodeString(ipqinp->QModMask.map[i], cscratch, &ipos, MAXNAME);
    }
  }
}

//-----------------------------------------------------------------------------
// BroadcastInputData: the input command file is read by the master thread.
//                     The information it contains, now encoded in a handful
//                     of structs which are inputs to this function, is
//                     broadcast to all other processes.   
//                                                                      
// Arguments:                                                           
//   dcinp:   direct space control data                                 
//   rcinp:   reciprocal space control data                             
//   tj:      trajectory control data                                   
//   tp:      the topology (not yet filled out, but some fields are assigned
//   myfit:   fitting data control structure                            
//   ipqinp:  IPolQ data control information                            
//-----------------------------------------------------------------------------
void BroadcastInputData(dircon *dcinp, reccon *rcinp, trajcon *tj, prmtop* tp,
			ipqcon *ipqinp)
{
  int i, nelem;
  int tlen[2];
  double* dscratch;
  char* cscratch;

  // Allocate memory for transferring contiguous data
  dscratch = (double*)malloc(256*sizeof(double));
  cscratch = (char*)malloc(1048576*sizeof(char));
  MPI_Bcast(dcinp, 1, tj->MPI_DIRCON, 0, MPI_COMM_WORLD);

  // Broadcast reccon information
  if (tj->tid == 0) {
    nelem = EncodeReccon(rcinp, dscratch, 0);
  }
  MPI_Bcast(&nelem, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(dscratch, nelem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (tj->tid != 0) {
    rcinp->ng = (int*)malloc(3*sizeof(int));
    EncodeReccon(rcinp, dscratch, 1);
  }

  // Broadcast trajcon information
  if (tj->tid == 0) {
    EncodeTrajcon(tj, dscratch, cscratch, 0, tlen);
  }
  MPI_Bcast(tlen, 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(dscratch, tlen[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(cscratch, tlen[1], MPI_CHAR, 0, MPI_COMM_WORLD);
  if (tj->tid != 0) {
    EncodeTrajcon(tj, dscratch, cscratch, 1, tlen);
  }

  // Broadcast topology information
  for (i = 0; i < tj->ntop; i++) {
    if (tj->tid == 0) {
      EncodePrmtop(&tp[i], dscratch, cscratch, 0, tlen);
    }
    MPI_Bcast(tlen, 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(dscratch, tlen[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(cscratch, tlen[1], MPI_CHAR, 0, MPI_COMM_WORLD);
    if (tj->tid != 0) {
      EncodePrmtop(&tp[i], dscratch, cscratch, 1, tlen);
    }
  }

  // Broadcast IPolQ control data
  if (tj->mode == 5) {
    if (tj->tid == 0) {
      EncodeIpqcon(ipqinp, dscratch, cscratch, 0, tlen);
    }
    MPI_Bcast(tlen, 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(dscratch, tlen[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(cscratch, tlen[1], MPI_CHAR, 0, MPI_COMM_WORLD);
    if (tj->tid != 0) {
      EncodeIpqcon(ipqinp, dscratch, cscratch, 1, tlen);
    }
  }

  // FIX ME!  Broadcast fitting control data!

  // Free allocated memory
  free(dscratch);
  free(cscratch);
}

//-----------------------------------------------------------------------------
// BroadcastCoordinates: broadcast an entire coordinates structure to all
//                       processes.                                 
//                                                                      
// Arguments:                                                           
//   crd:     the coordinates to broadcast                              
//   tj:      trajectory control data, which also contains the process ranks in
//            MPI_COMM_WORLD                                   
//   sID:     the ID of the system according to the master trajectory control 
//            struct
//-----------------------------------------------------------------------------
void BroadcastCoordinates(coord *crd, trajcon *tj, int sID)
{
  int i;
  double* buffer;

  buffer = (double*)malloc(32*sizeof(double));

  // The number of atoms in the coordinate set is assumed to be known;
  // this should be true as it is also assumed that space has been    
  // allocated for the coordinates.                                   
  i = 3*crd->natom;
  MPI_Valgrind_bcast(crd->atmid, crd->natom, MPI_INT, 0, tj->SysComm[sID],
		     IDEBUG);
  MPI_Valgrind_bcast(crd->loc, i, MPI_DOUBLE, 0, tj->SysComm[sID], IDEBUG);
  MPI_Valgrind_bcast(crd->prvloc, i, MPI_DOUBLE, 0, tj->SysComm[sID], IDEBUG);
  MPI_Valgrind_bcast(crd->vel, i, MPI_DOUBLE, 0, tj->SysComm[sID], IDEBUG);
  MPI_Valgrind_bcast(crd->prvvel, i, MPI_DOUBLE, 0, tj->SysComm[sID], IDEBUG);
  MPI_Valgrind_bcast(crd->frc, i, MPI_DOUBLE, 0, tj->SysComm[sID], IDEBUG);
  MPI_Valgrind_bcast(crd->prvfrc, i, MPI_DOUBLE, 0, tj->SysComm[sID], IDEBUG);
  if (tj->tid == 0) {
    for (i = 0; i < 3; i++) {
      buffer[i] = crd->gdim[i];
      buffer[i+3] = crd->gdim[i+3];
      buffer[i+6] = crd->hgdim[i];
    }
    for (i = 0; i < 9; i++) {
      buffer[i+9] = crd->U.data[i];
      buffer[i+18] = crd->invU.data[i];
    }
    EncodeInteger(buffer, 27, crd->isortho);
    i = 28;
    EncodeLongLongInt(buffer, &i, tj->currstep);
    EncodeInteger(buffer, 30, tj->currfi);
    EncodeInteger(buffer, 31, tj->irest);
  }
  MPI_Bcast(buffer, 32, MPI_DOUBLE, 0, tj->SysComm[sID]);
  if (tj->tid != 0) {
    for (i = 0; i < 3; i++) {
      crd->gdim[i] = buffer[i];
      crd->gdim[i+3] = buffer[i+3];
      crd->hgdim[i] = buffer[i+6];
    }
    for (i = 0; i < 9; i++) {
      crd->U.data[i] = buffer[i+9];
      crd->invU.data[i] = buffer[i+18];
    }
    crd->isortho = buffer[27];
    tj->currstep = DecodeLongLongInt(buffer, 28);
    tj->currfi = buffer[30];
    tj->irest = buffer[31];
  }
  free(buffer);
}

//-----------------------------------------------------------------------------
// BroadcastEPInfo: broadcast information relating to extra point rules and
//                  extra point definitions within the topology.
//
// Arguments:                                                           
//   tp:    the topology                                                
//-----------------------------------------------------------------------------
void BroadcastEPInfo(prmtop *tp, trajcon *tj)
{
  MPI_Bcast(&tp->EPInserted, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tp->neprule, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (tj->tid != 0) {
    tp->norigatom = tp->natom;
  }
  if (tp->neprule == 0) {
    return;
  }

  // If we're still here, there is at least one
  // EP rule to broadcast and adjustments to   
  // make to other information in the topology.
  int i;
  eprule *tmr;

  if (tj->tid != 0) {
    tp->eprules = (eprule*)malloc(tp->neprule*sizeof(eprule));
  }
  MPI_Bcast(tp->eprules, tp->neprule, tj->MPI_EPRULE, 0, MPI_COMM_WORLD);
  if (tj->tid != 0) {
    for (i = 0; i < tp->neprule; i++) {
      tmr = &tp->eprules[i];
      if (tmr->sig >= 1.0e-8 && tmr->eps >= 1.0e-8) {
	ExpandLJTables(tp, tmr->sig, tmr->eps);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// MPI_Valgrind_bcast: this function wraps MPI_Bcast and ensures that no
//                     broadcast message exceeds a size limit defined at
//                     compile time to help avoid large messages which seem to
//                     overflow some sort of memory limit when using valgrind 
//                     to debug.  This function works for MPI elemental data
//                     types, not derived types.     
//                                                                      
// Arguments:                                                           
//   buffer:   the data buffer (this argument, and also count, datatype, root,
//             and comm, are equivalent to the MPI_bcast arguments in the MPI
//             standard)
//   count:    the number of elements in buffer                         
//   datatype: the type of elements in buffer                           
//   root:     the root process of the broadcast, the one from which all data
//             originates                                          
//   comm:     the communicator                                         
//   idebug:   flag to indicate that debugging is requested             
//-----------------------------------------------------------------------------
int MPI_Valgrind_bcast(void *buffer, int count, MPI_Datatype datatype,
		       int root, MPI_Comm comm, int idebug)
{
  int i, stt, typesize, thismsg, maxmsg, errlen, commsize;
  int *ibuffer;
  float *fbuffer;
  double *dbuffer;
  char *cbuffer;

  // Trivial case: no broadcast necessary
  MPI_Comm_size(comm, &commsize);
  if (commsize == 1) {
    return MPI_SUCCESS;
  }

  // Straightforward case: broadcast the data in one shot
  if (idebug == 0) {
    stt = MPI_Bcast(buffer, count, datatype, root, comm);
    return stt;
  }

  // More complex case: debugging is turned on, so broadcasts
  // of large data should be broken up for valgrind to digest
  if (datatype == MPI_INT) {
    ibuffer = (int*)buffer;
    typesize = sizeof(int);
  }
  else if (datatype == MPI_FLOAT) {
    fbuffer = (float*)buffer;
    typesize = sizeof(float);
  }
  else if (datatype == MPI_DOUBLE) {
    dbuffer = (double*)buffer;
    typesize = sizeof(double);
  }
  else if (datatype == MPI_CHAR) {
    cbuffer = (char*)buffer;
    typesize = sizeof(char);
  }
  else {
    printf("MPI_Valgrind_bcast >> Error.  Data type not recognized.\n");
    exit(1);
  }

  // Send the broadcast all at once
  if (typesize*count < VALGR_MAX_BUFF) {
    stt = MPI_Bcast(buffer, count, datatype, root, comm);
    return stt;
  }

  // Send the broadcast in stages
  i = 0;
  maxmsg = VALGR_MAX_BUFF / typesize;
  while (i < count) {
    thismsg = (i + maxmsg <= count) ? maxmsg : count - i;
    if (datatype == MPI_INT) {
      stt = MPI_Bcast(&ibuffer[i], thismsg, datatype, root, comm);
    }
    else if (datatype == MPI_FLOAT) {
      stt = MPI_Bcast(&fbuffer[i], thismsg, datatype, root, comm);
    }
    else if (datatype == MPI_DOUBLE) {
      stt = MPI_Bcast(&dbuffer[i], thismsg, datatype, root, comm);
    }
    else if (datatype == MPI_CHAR) {
      stt = MPI_Bcast(&cbuffer[i], thismsg, datatype, root, comm);
    }
    i += thismsg;

    if (stt != MPI_SUCCESS) {
      int myrank;
      char errmsg[128];

      MPI_Comm_rank(comm, &myrank);
      printf("MPI_Valgrind_bcast >> Warning, MPI error broadcasting %d "
	     "elements at index\nMPI_Valgrind_bcast >> %d of %d total on "
	     "thread %d.\n", thismsg, i, count, myrank);
      sprintf(errmsg, "MPI_Valgrind_bcast >> Error on thread %d: MPI_ERR_",
	      myrank);
      errlen = strlen(errmsg);
      if (stt == MPI_ERR_COMM) {
	sprintf(&errmsg[errlen], "COMM");
      }
      else if (stt == MPI_ERR_COUNT) {
	sprintf(&errmsg[errlen], "COUNT");
      }
      else if (stt == MPI_ERR_TYPE) {
	sprintf(&errmsg[errlen], "TYPE");
      }
      else if (stt == MPI_ERR_BUFFER) {
	sprintf(&errmsg[errlen], "BUFFER");
      }
      else if (stt == MPI_ERR_ROOT) {
	sprintf(&errmsg[errlen], "ROOT");
      }
      printf("%s\n", errmsg);
    }
  }

  return stt;
}
#endif
