//-----------------------------------------------------------------------------
// This code is included several times by ArraySimulator.cu to make versions
// of the dynamics kernel, with different block sizes and the option of energy
// computations.  These kernels will be called with a minimum of 256 threads
// and a maximum of 1024.
//-----------------------------------------------------------------------------
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 370 || __CUDA_ARCH__ >= 700 || \
			       __CUDA_ARCH__ == 600)
#  define PREFER_AUTO_L1
#endif
{
  // Convenient scalars to import system-specific data onto chip memory
  __shared__ volatile int sysID, atomStart, atomEnd, natom, ntypes, stepidx;
  __shared__ volatile int bondStart, bondPos, bondEnd, randStart, randEnd;
  __shared__ volatile int anglStart, anglPos, anglEnd;
  __shared__ volatile int diheStart, dihePos, diheEnd;
  __shared__ volatile int cnstStart, cnstPos, cnstEnd;
  __shared__ volatile int tileStart, tilePos, tileEnd, gbReffPos, gbDRadPos;
  __shared__ volatile int attnStart, attnPos, attnEnd;
  __shared__ volatile float sysSdfac;
#ifdef GO_RATTLE
  __shared__ volatile float *xprvcrd, *yprvcrd, *zprvcrd;
#endif
  // These arrays will be too big for the compiler to stuff into registers.
  __shared__ volatile int ixfrc[ATOM_LIMIT], iyfrc[ATOM_LIMIT];
  __shared__ volatile int izfrc[ATOM_LIMIT];
#ifdef COMPUTE_ENERGY
  __shared__ volatile int sBondNrg[THREAD_COUNT / GRID];
  __shared__ volatile int sAnglNrg[THREAD_COUNT / GRID];
  __shared__ volatile int sDiheNrg[THREAD_COUNT / GRID];
  __shared__ volatile int sSolvNrg[THREAD_COUNT / GRID];
  __shared__ volatile int sElecNrg[THREAD_COUNT / GRID];
  __shared__ volatile int sVdwNrg[THREAD_COUNT / GRID];
  __shared__ volatile float sKineNrg[THREAD_COUNT / GRID];
#endif
  __shared__ volatile float xcrd[ATOM_LIMIT], xvel[ATOM_LIMIT];
  __shared__ volatile float ycrd[ATOM_LIMIT], yvel[ATOM_LIMIT];
  __shared__ volatile float zcrd[ATOM_LIMIT], zvel[ATOM_LIMIT];
  __shared__ volatile int psi[ATOM_LIMIT], sumdeijda[ATOM_LIMIT];
  __shared__ volatile float reff[ATOM_LIMIT];
  
  // Each block is going to take one system and run with it
  if (threadIdx.x == 0) {
    sysID = blockIdx.x;
#ifdef GO_RATTLE
    xprvcrd = (float*)psi;
    yprvcrd = (float*)sumdeijda;
    zprvcrd = reff;
#endif
  }
  __syncthreads();
  while (sysID < cGms.nsys) {

    // Read the coordinate and parameter bounds of this system from global
    if (threadIdx.x == 0) {
      int2 limits = cGms.DVCatomReadLimits[sysID];
      natom = limits.y - limits.x;
      atomStart = limits.x;
      atomEnd   = limits.y;
      if (cGms.rattle == 1) {
        limits = cGms.DVCcnstReadLimits[sysID];
        cnstStart = limits.x;
        cnstEnd   = limits.y;
      }
    }
    else if (threadIdx.x == GRID) {
      int2 limits = cGms.DVCbondReadLimits[sysID];
      bondStart = limits.x;
      bondEnd   = limits.y;
      ntypes    = cGms.DVCtypeCounts[sysID];
      if (cGms.Tstat.active == 1) {
        int2 limits = cGms.DVCprngReadLimits[sysID];
        randStart = limits.x;
        randEnd   = limits.y;
        sysSdfac = cGms.DVCsdfac[sysID];
      }
    }
    else if (threadIdx.x == GRIDx2) {
      int2 limits = cGms.DVCanglReadLimits[sysID];
      anglStart = limits.x;
      anglEnd   = limits.y;
      limits = cGms.DVCdiheReadLimits[sysID];
      diheStart = limits.x;
      diheEnd   = limits.y;
    }
    else if (threadIdx.x == GRIDx3) {
      int4 limits = cGms.DVCnbReadLimits[sysID];
      tileStart = limits.x;
      tileEnd   = limits.y;
      attnStart = limits.z;
      attnEnd   = limits.w;
    }
    __syncthreads();

    // Read atom coordinates, velocities, and forces from global.
    // Re-center coordinates on energy computation steps.
#ifdef COMPUTE_ENERGY
    int icenx = 0;
    int iceny = 0;
    int icenz = 0;
#endif
    int tgx = (threadIdx.x & GRID_BITS_MASK);
    if (threadIdx.x < natom) {
      int readidx = atomStart + threadIdx.x;
      xcrd[threadIdx.x] = cGms.DVCatomX[readidx];
      ycrd[threadIdx.x] = cGms.DVCatomY[readidx];
      zcrd[threadIdx.x] = cGms.DVCatomZ[readidx];
      xvel[threadIdx.x] = cGms.DVCvelcX[readidx];
      yvel[threadIdx.x] = cGms.DVCvelcY[readidx];
      zvel[threadIdx.x] = cGms.DVCvelcZ[readidx];
      ixfrc[threadIdx.x] = cGms.DVCfrcX[readidx];
      iyfrc[threadIdx.x] = cGms.DVCfrcY[readidx];
      izfrc[threadIdx.x] = cGms.DVCfrcZ[readidx];
#ifdef COMPUTE_ENERGY
      icenx = (int)(xcrd[threadIdx.x] * RCSCALE);
      iceny = (int)(ycrd[threadIdx.x] * RCSCALE);
      icenz = (int)(zcrd[threadIdx.x] * RCSCALE);
#endif
    }
    else if (threadIdx.x < ATOM_LIMIT) {
      xcrd[threadIdx.x] = (float)1.0e5 + (float)100.0*(threadIdx.x - natom);
      ycrd[threadIdx.x] = (float)1.0e5 + (float)100.0*(threadIdx.x - natom);
      zcrd[threadIdx.x] = (float)1.0e5 + (float)100.0*(threadIdx.x - natom);
    }

#ifdef COMPUTE_ENERGY
    // Hijack some energy arrays to do the block-wide reduction
    WarpREDUCE(icenx);
    WarpREDUCE(iceny);
    WarpREDUCE(icenz);
    if (tgx == 0) {
      sBondNrg[threadIdx.x / GRID] = icenx;
      sAnglNrg[threadIdx.x / GRID] = iceny;
      sDiheNrg[threadIdx.x / GRID] = icenz;
    }
    __syncthreads();
    if (threadIdx.x < GRID) {
      if (tgx < THREAD_COUNT / GRID) {
        icenx = sBondNrg[tgx];
      }
      else {
        icenx = 0;
      }
      WarpREDUCE(icenx);
      if (tgx == 0) {
        sBondNrg[0] = __float_as_int((float)icenx / (RCSCALE * (float)natom));
      }
    }
    else if (threadIdx.x < GRIDx2) {
      if (tgx < THREAD_COUNT / GRID) {
        iceny = sAnglNrg[tgx];
      }
      else {
        iceny = 0;
      }
      WarpREDUCE(iceny);
      if (tgx == 0) {
        sAnglNrg[0] = __float_as_int((float)iceny / (RCSCALE * (float)natom));
      }
    }
    else if (threadIdx.x < GRIDx3) {
      if (tgx < THREAD_COUNT / GRID) {
        icenz = sDiheNrg[tgx];
      }
      else {
        icenz = 0;
      }
      WarpREDUCE(icenz);
      if (tgx == 0) {
        sDiheNrg[0] = __float_as_int((float)icenz / (RCSCALE * (float)natom));
      }
    }
    __syncthreads();
    float cenx = __int_as_float(sBondNrg[0]);
    float ceny = __int_as_float(sAnglNrg[0]);
    float cenz = __int_as_float(sDiheNrg[0]);
    if (threadIdx.x < natom) {
      xcrd[threadIdx.x] -= cenx;
      ycrd[threadIdx.x] -= ceny;
      zcrd[threadIdx.x] -= cenz;
    }
    __syncthreads();
#endif
    
    // Dynamics loop
#ifdef COMPUTE_ENERGY
    if (threadIdx.x < THREAD_COUNT / GRID) {
      sBondNrg[threadIdx.x] = 0;
      sAnglNrg[threadIdx.x] = 0;
      sDiheNrg[threadIdx.x] = 0;
      sSolvNrg[threadIdx.x] = 0;
      sElecNrg[threadIdx.x] = 0;
      sVdwNrg[threadIdx.x] = 0;
      sKineNrg[threadIdx.x] = (float)0.0;
    }
    if (threadIdx.x == 0) {
      stepidx = 0;
    }
    __syncthreads();
#else
    if (threadIdx.x == 0) {
      stepidx = 1;
    }
    __syncthreads();
    while (stepidx < cGms.ntpr) {
#endif
      // If a thermostat requiring random numbers is in effect (the
      // Langevin thermostat is the only one currently supported),
      // make a heap of random numbers to use.
      if (cGms.Tstat.active == 1 &&
          (stepidx * cGms.slowFrcMult) % cGms.Tstat.refresh == 0) {
        __syncthreads();
        int prnidx = threadIdx.x + (blockIdx.x * blockDim.x);
        curandState_t *sttptr = (curandState_t*)cGms.prngStates;
        curandState_t locstt = sttptr[prnidx];        

        // First, take the final set of random numbers for this
        // system and put them where the Velocity Verlet coordinate
        // update will find them.
        int rndidx = randStart + threadIdx.x;
        while (rndidx < randStart + 3*natom) {
          cGms.DVCprngHeap[rndidx] =
            cGms.DVCprngHeap[rndidx + (3 * natom * cGms.Tstat.refresh)];
          rndidx += blockDim.x;
        }

        // Generate new random numbers
        rndidx = randStart + 3*natom + threadIdx.x;
        while (rndidx < randEnd) {
          cGms.DVCprngHeap[rndidx] = curand_normal(&locstt);
          rndidx += blockDim.x;
        }
        sttptr[prnidx] = locstt;
        __syncthreads();
      }
      
      // Loop over minor steps, including the slow
      // force computation on the first cycle.
      int i;
      for (i = 0; i < cGms.slowFrcMult; i++) {

        // Prior to moving atoms, store the current locations for
        // constraint applications later.  Use the GB-devoted arrays
        // as these are unneeded at the moment.
#ifdef GO_RATTLE
        if (cGms.rattle == 1 && threadIdx.x < natom) {
          xprvcrd[threadIdx.x] = xcrd[threadIdx.x];
          yprvcrd[threadIdx.x] = ycrd[threadIdx.x];
          zprvcrd[threadIdx.x] = zcrd[threadIdx.x];
        }
#endif        
        // Velocity Verlet half kick with coordinates update.
        // Zero forces when done.
        if (threadIdx.x < natom) {
          float hmdt = cGms.DVCatomHDTM[atomStart + threadIdx.x];          
          if (cGms.Tstat.active == 1) {
            float atmmass = cGms.hdtm2mass / hmdt;
            float rsd = FPSCALEfrc * sysSdfac * sqrt(atmmass);
            int rndidx = randStart +
                         3*(((stepidx * cGms.slowFrcMult) + i) %
                            cGms.Tstat.refresh)*natom +
                         threadIdx.x;
            xvel[threadIdx.x] = (cGms.Tstat.c_explic * xvel[threadIdx.x]) +
                                hmdt * ((float)ixfrc[threadIdx.x] +
                                        rsd * cGms.DVCprngHeap[rndidx]);
            rndidx += natom;
            yvel[threadIdx.x] = (cGms.Tstat.c_explic * yvel[threadIdx.x]) +
                                hmdt * ((float)iyfrc[threadIdx.x] +
                                        rsd * cGms.DVCprngHeap[rndidx]);
            rndidx += natom;
            zvel[threadIdx.x] = (cGms.Tstat.c_explic * zvel[threadIdx.x]) +
                                hmdt * ((float)izfrc[threadIdx.x] +
                                        rsd * cGms.DVCprngHeap[rndidx]);
          }
          else {
            xvel[threadIdx.x] += hmdt * (float)ixfrc[threadIdx.x];
            yvel[threadIdx.x] += hmdt * (float)iyfrc[threadIdx.x];
            zvel[threadIdx.x] += hmdt * (float)izfrc[threadIdx.x];
          }
          xcrd[threadIdx.x] += cGms.dtVF * xvel[threadIdx.x];
          ycrd[threadIdx.x] += cGms.dtVF * yvel[threadIdx.x];
          zcrd[threadIdx.x] += cGms.dtVF * zvel[threadIdx.x];
	}
	
        // Apply geometry constraints.  Use force accumulators
        // to store perturbations.
#ifdef GO_RATTLE
        if (cGms.rattle == 1) {
	  if (threadIdx.x < natom) {
            ixfrc[threadIdx.x] = 0;
            iyfrc[threadIdx.x] = 0;
            izfrc[threadIdx.x] = 0;
	  }
          if (threadIdx.x == 0) {
            cnstPos = cnstStart;
          }
          __syncthreads();
          int warpCnstPos;
          if (tgx == 0) {
            warpCnstPos = atomicAdd((int*)&cnstPos, 32);
          }
          warpCnstPos = __shfl_sync(0xffffffff, warpCnstPos, 0, 32);
          while (warpCnstPos < cnstEnd) {
#ifdef PREFER_AUTO_L1
            uint2 params = cGms.DVCcnstInfo[warpCnstPos + tgx];
#else
            uint2 params = __ldg(&cGms.DVCcnstInfo[warpCnstPos + tgx]);
#endif
            int atma = params.x;
            bool fixlen = (atma != 0xffffffff);
            unsigned int engaged = __ballot_sync(0xffffffff, fixlen);
            if (fixlen) {
              int atmb = (atma & 0xffff);
              atma >>= 16;
#ifdef PREFER_AUTO_L1
              float rma = cGms.DVCatomHDTM[atomStart + atma] * cGms.hdtm2invm;
              float rmb = cGms.DVCatomHDTM[atomStart + atmb] * cGms.hdtm2invm;
#else
              float rma = __ldg(&cGms.DVCatomHDTM[atomStart + atma]) *
                          cGms.hdtm2invm;
              float rmb = __ldg(&cGms.DVCatomHDTM[atomStart + atmb]) *
                          cGms.hdtm2invm;
#endif
              double l0 = (double)params.y / 1.0e8;
              l0 *= l0;
              int niter = 0;
              unsigned int atwork = engaged;
              while (niter < cGms.maxRattleIter && atwork > 0) {
                double rx =
                  ((double)xcrd[atmb] + ((double)ixfrc[atmb] * FPSCALEicn)) -
                  ((double)xcrd[atma] + ((double)ixfrc[atma] * FPSCALEicn));
                double ry =
                  ((double)ycrd[atmb] + ((double)iyfrc[atmb] * FPSCALEicn)) -
                  ((double)ycrd[atma] + ((double)iyfrc[atma] * FPSCALEicn));
                double rz =
                  ((double)zcrd[atmb] + ((double)izfrc[atmb] * FPSCALEicn)) -
                  ((double)zcrd[atma] + ((double)izfrc[atma] * FPSCALEicn));
                double r2 = rx*rx + ry*ry + rz*rz;
                double delta = l0 - r2;
                bool violation = (fabs(delta) > cGms.rattleTol);
                atwork = __ballot_sync(engaged, violation);
                if (violation) {
                  double rrefx = (double)xprvcrd[atmb] - (double)xprvcrd[atma];
                  double rrefy = (double)yprvcrd[atmb] - (double)yprvcrd[atma];
                  double rrefz = (double)zprvcrd[atmb] - (double)zprvcrd[atma];
                  double dot = rx*rrefx + ry*rrefy + rz*rrefz;
                  float term = (float)1.2 * (float)delta /
                               ((float)2.0 * (float)dot * (rma + rmb));
                  float xterm = (float)rrefx * term;
                  float yterm = (float)rrefy * term;
                  float zterm = (float)rrefz * term;
                  atomicAdd((int*)&ixfrc[atma], (int)(-xterm*rma*FPSCALEcn));
                  atomicAdd((int*)&iyfrc[atma], (int)(-yterm*rma*FPSCALEcn));
                  atomicAdd((int*)&izfrc[atma], (int)(-zterm*rma*FPSCALEcn));
                  atomicAdd((int*)&ixfrc[atmb], (int)( xterm*rmb*FPSCALEcn));
                  atomicAdd((int*)&iyfrc[atmb], (int)( yterm*rmb*FPSCALEcn));
                  atomicAdd((int*)&izfrc[atmb], (int)( zterm*rmb*FPSCALEcn));
                }
		__syncwarp(engaged);
                niter++;
              }
            }

            // Increment the warp position
            if (tgx == 0) {
              warpCnstPos = atomicAdd((int*)&cnstPos, 32);
            }
            warpCnstPos = __shfl_sync(0xffffffff, warpCnstPos, 0, 32);
          }
          __syncthreads();

          // Contribute the perturbations back to the coordinates
          if (threadIdx.x < natom) {
            double xmove = (double)(ixfrc[threadIdx.x]) * FPSCALEicn;
            double ymove = (double)(iyfrc[threadIdx.x]) * FPSCALEicn;
            double zmove = (double)(izfrc[threadIdx.x]) * FPSCALEicn;
            xcrd[threadIdx.x] += xmove;
            ycrd[threadIdx.x] += ymove;
            zcrd[threadIdx.x] += zmove;
            xvel[threadIdx.x] += xmove * cGms.invdtVF;
            yvel[threadIdx.x] += ymove * cGms.invdtVF;
            zvel[threadIdx.x] += zmove * cGms.invdtVF;
          }
        }
#endif
        // Reset force accumlators.  Snapshot coordinates and
        // velocities for writing the restart file, if appropriate.
        if (threadIdx.x < natom) {
          ixfrc[threadIdx.x] = 0;
          iyfrc[threadIdx.x] = 0;
          izfrc[threadIdx.x] = 0;
          if (i == 0 && stepidx == cGms.ntpr - 1) {
            cGms.DVCcrdrX[atomStart + threadIdx.x] = xcrd[threadIdx.x];
            cGms.DVCcrdrY[atomStart + threadIdx.x] = ycrd[threadIdx.x];
            cGms.DVCcrdrZ[atomStart + threadIdx.x] = zcrd[threadIdx.x];
            cGms.DVCvelrX[atomStart + threadIdx.x] = xvel[threadIdx.x];
            cGms.DVCvelrY[atomStart + threadIdx.x] = yvel[threadIdx.x];
            cGms.DVCvelrZ[atomStart + threadIdx.x] = zvel[threadIdx.x];
          }
        }
        
        // Set counters for all force components
        if (threadIdx.x == 0) {
          tilePos = tileStart;
          gbReffPos = tileStart;
          gbDRadPos = tileStart;
        }
        else if (threadIdx.x == GRID) {
          bondPos = bondStart;
          anglPos = anglStart;
        }
        else if (threadIdx.x == GRIDx2) {
          dihePos = diheStart;
        }
        else if (threadIdx.x == GRIDx3) {
          attnPos = attnStart;
        }
        __syncthreads();

        // Computation of slowly varying force components.  This happens
        // on the first cycle and the results are multiplied by the cycle
        // count to ensure the correct impulse.
        int warpWorkPos;
        if (i == 0) {

          // Accumulate Generalized Born effective radii
#ifdef GBSOLVENT
#  ifdef COMPUTE_ENERGY
          int usolv = 0;
#  endif
          if (threadIdx.x < natom) {
            psi[threadIdx.x] = 0;
          }
          __syncthreads();
          if (tgx == 0) {
            warpWorkPos = atomicAdd((int*)&gbReffPos, 32);
          }
          warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
          while (warpWorkPos < tileEnd) {

            // Load the pair mask.  The high 16 bits of the
            // mask indicate the position of the tile, and are the same
            // for all 16 of the unsigned int masks in each tile.
#ifdef PREFER_AUTO_L1
            unsigned int pairmask = cGms.DVCpairBitMasks[warpWorkPos + tgx];
#else
            unsigned int pairmask = __ldg(&cGms.DVCpairBitMasks[warpWorkPos +
                                                                tgx]);
#endif	    
            // Initialize atoms and forces
            int katm = (pairmask >> 24)*16 + (tgx & 15);
            int matm = ((pairmask >> 16) & 0xff)*16 + (tgx & 15);
            int loadatm    = (tgx < 16) ? katm : matm;
            float txcrd    = xcrd[loadatm];
            float tycrd    = ycrd[loadatm];
            float tzcrd    = zcrd[loadatm];
#ifdef PREFER_AUTO_L1
            float tbrad    = cGms.DVCrborn[atomStart + loadatm] -
                             cGms.GBOffset;
            float tfs      = cGms.DVCFscreen[atomStart + loadatm];
            int   tneck    = cGms.DVCNeckID[atomStart + loadatm];
#else
            float tbrad    = __ldg(&cGms.DVCrborn[atomStart + loadatm]) -
                             cGms.GBOffset;
            float tfs      = __ldg(&cGms.DVCFscreen[atomStart + loadatm]);
            int   tneck    = __ldg(&cGms.DVCNeckID[atomStart + loadatm]);
#endif
            float tinvrad  = (float)1.0 / tbrad;

            // Detect whether this tile is on the diagonal,
            // so that self-interactions can be eliminated
            int ondiagonal = ((katm >> 4) == (matm >> 4));
              
            // Accumulate psi
            float psiacc = (float)0.0;
            int j;
            for (j = 0; j < 8; j++) {
              int crdSrcLane = (tgx <  16)*(tgx + 16 - j) +
                               (tgx >= 16)*(tgx - 24 + j);
              crdSrcLane += 16 * ((crdSrcLane < 16 && tgx <  16) +
                                  (crdSrcLane <  0 && tgx >= 16));
              float oxcrd = __shfl_sync(0xffffffff, txcrd, crdSrcLane, 32);
              float oycrd = __shfl_sync(0xffffffff, tycrd, crdSrcLane, 32);
              float ozcrd = __shfl_sync(0xffffffff, tzcrd, crdSrcLane, 32);
              float obrad = __shfl_sync(0xffffffff, tbrad, crdSrcLane, 32);
              float ofs   = __shfl_sync(0xffffffff, tfs,   crdSrcLane, 32);
              int   oneck = __shfl_sync(0xffffffff, tneck, crdSrcLane, 32);
              float dx = oxcrd - txcrd;
              float dy = oycrd - tycrd;
              float dz = ozcrd - tzcrd;
              float r2 = dx*dx + dy*dy + dz*dz;
              float invr  = rsqrt(r2);
              float r  = r2 * invr;
                
              // First computation
              float sj  = ofs * obrad;
              float sj2 = sj * sj;
              float psitmp = (float)0.0;
              if (r > (float)4.0 * sj) {
                float invr2  = invr * invr;
                float tmpsd  = sj2 * invr2;
                float dumbo  = TA + tmpsd*(TB +
                                           tmpsd*(TC +
                                                  tmpsd*(TD + tmpsd*TDD)));
                psitmp -= sj * tmpsd * invr2 * dumbo;
              }
              else if (r > tbrad + sj) {
                psitmp -= (float)0.5 * ((sj / (r2 - sj2)) +
                                        ((float)0.5 * invr * log((r - sj) /
                                                                 (r + sj))));
              }
              else if (r > fabs(tbrad - sj)) {
                float theta = (float)0.5 * tinvrad * invr *
                              (r2 + tbrad*tbrad - sj2);
                float uij = (float)1.0 / (r + sj);
                psitmp -= (float)0.25 * (tinvrad*((float)2.0 - theta) - uij +
                                         invr*log(tbrad * uij));
              }
              else if (tbrad < sj) {
                psitmp -= (float)0.5 * ((sj / (r2 - sj2)) +
                                        ((float)2.0 * tinvrad) +
                                        ((float)0.5 * invr * log((sj - r) /
                                                                 (sj + r))));
              }

              // Second computation
              float si  = tfs * tbrad;
              float si2 = si * si;
              float psiret = (float)0.0;
              float oinvrad = (float)1.0 / obrad;
              if (r > (float)4.0 * si) {
                float invr2  = invr * invr;
                float tmpsd  = si2 * invr2;
                float dumbo  = TA + tmpsd*(TB +
                                           tmpsd*(TC +
                                                  tmpsd*(TD + tmpsd*TDD)));
                psiret -= si * tmpsd * invr2 * dumbo;
              }
              else if (r > obrad + si) {
                psiret -= (float)0.5 * ((si / (r2 - si2)) +
                                        ((float)0.5 * invr * log((r - si) /
                                                                 (r + si))));
              }
              else if (r > fabs(obrad - si)) {
                float theta = (float)0.5 * oinvrad * invr *
                              (r2 + obrad*obrad - si2);
                float uij = (float)1.0 / (r + si);
                psiret -= (float)0.25 * (oinvrad*((float)2.0 - theta) - uij +
                                         invr*log(obrad * uij));
              }
              else if (obrad < si) {
                psiret -= (float)0.5 * ((si / (r2 - si2)) +
                                        ((float)2.0 * oinvrad) +
                                        ((float)0.5 * invr * log((si - r) /
                                                                 (si + r))));
              }
                
              // Neck GB contributions
              if ((cGms.igb == 7 || cGms.igb == 8) &&
                  r < tbrad + obrad + cGms.GBNeckCut) {

                // First computation
#ifdef PREFER_AUTO_L1
                float2 neckData = cGms.DVCneckFactors[21*tneck + oneck];
#else
                float2 neckData = __ldg(&cGms.DVCneckFactors[21*tneck +
                                                             oneck]);
#endif
                float mdist  = r - neckData.x;
                float mdist2 = mdist * mdist;
                float mdist6 = mdist2 * mdist2 * mdist2;
                float neck   = neckData.y /
                               ((float)1.0 + mdist2 + (float)0.3*mdist6);
                psitmp -= cGms.GBNeckScale * neck;
          
                // Second computation
#ifdef PREFER_AUTO_L1
                neckData = cGms.DVCneckFactors[21*oneck + tneck];
#else
                neckData = __ldg(&cGms.DVCneckFactors[21*oneck + tneck]);
#endif
                mdist  = r - neckData.x;
                mdist2 = mdist * mdist;
                mdist6 = mdist2 * mdist2 * mdist2;
                neck   = neckData.y /
                         ((float)1.0 + mdist2 + (float)0.3*mdist6);
                psiret -= cGms.GBNeckScale * neck;
              }
                
              // Account for diagonal tiles
              if (ondiagonal) {
                psitmp *= (float)0.5;
                psiret *= (float)0.5;
                if (j == 0 && tgx < 16) {
                  psitmp = (float)0.0;
                  psiret = (float)0.0;
                }
              }

              // Return the accumulated psi value
              int frcRetLane = (tgx <  16)*(tgx + 24 - j) +
                               (tgx >= 16)*(tgx - 16 + j);
              frcRetLane += 16 * ((frcRetLane >= 16 && tgx >= 16) -
                                  (frcRetLane >= 32 && tgx <  16));
              psitmp += __shfl_sync(0xffffffff, psiret, frcRetLane, 32);
              psiacc += psitmp;
            }

            // Contribute to __shared__ psi accumulators
            atomicAdd((int*)&psi[loadatm], (int)(psiacc * FPSCALErad));

            // Increment the tile position
            if (tgx == 0) {
              warpWorkPos = atomicAdd((int*)&gbReffPos, 32);
            }
            warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
          }

          // Final pass to get effective Born radii, initialize GB forces
          __syncthreads();
          if (threadIdx.x < natom) {
            float fipsi = (float)(psi[threadIdx.x]) / FPSCALErad;
            float atmrad;
            int readidx = atomStart + threadIdx.x;
            if (cGms.igb == 1) {

              // Original Hawkins-Craemer-Truhlar effective radii
#ifdef PREFER_AUTO_L1
              float invirad = (float)1.0 / (cGms.DVCrborn[readidx] -
                                            cGms.GBOffset);
#else
              float invirad = (float)1.0 / (__ldg(&cGms.DVCrborn[readidx]) -
                                            cGms.GBOffset);
#endif
              atmrad = (float)1.0 / (invirad + fipsi);
              if (atmrad < (float)0.0) {
                atmrad = (float)30.0;
              }
            }
            else {

              // "GBAO" formulas
#ifdef PREFER_AUTO_L1
              float atmirad = cGms.DVCrborn[readidx] - cGms.GBOffset;
#else
              float atmirad = __ldg(&cGms.DVCrborn[readidx]) - cGms.GBOffset;
#endif
              float invirad = (float)1.0 / atmirad;
              fipsi *= -atmirad;
#ifdef PREFER_AUTO_L1
              atmrad = (float)1.0 /
                       (invirad - tanh((cGms.DVCGBalpha[readidx] -
                                        (cGms.DVCGBbeta[readidx] * fipsi) +
                                        (cGms.DVCGBgamma[readidx] *
                                         fipsi * fipsi)) * fipsi) /
                        cGms.DVCrborn[readidx]);
#else
              atmrad = (float)1.0 /
                       (invirad - tanh((__ldg(&cGms.DVCGBalpha[readidx]) -
                                       (__ldg(&cGms.DVCGBbeta[readidx]) *
                                        fipsi) +
                                       (__ldg(&cGms.DVCGBgamma[readidx]) *
                                        fipsi * fipsi)) * fipsi) /
                        __ldg(&cGms.DVCrborn[readidx]));
#endif
            }
            reff[threadIdx.x] = atmrad;
#ifdef PREFER_AUTO_L1
            float atmq   = cGms.DVCatomQ[readidx];
#else
            float atmq   = __ldg(&cGms.DVCatomQ[readidx]);
#endif
            float expmkf = exp(-cGms.kscale * (cGms.kappa) * atmrad) /
                           cGms.dielectric;
            float dielfac = (float)1.0 - expmkf;
            float atmq2h = (float)0.5 * atmq * atmq;
            float atmqd2h = atmq2h * dielfac;
#  ifdef COMPUTE_ENERGY
            usolv += (int)((-atmqd2h / atmrad) * FPSCALEnrgMP);
#  endif
            sumdeijda[threadIdx.x] = (int)(FPSCALEfrc *
                                           (atmqd2h -
                                            (cGms.kscale * cGms.kappa *
                                             atmq2h * expmkf * atmrad)));
          }
          __syncthreads();
#endif // GBSOLVENT

          // Non-bonded computations.  This uses a 16 x 16 tile
          // framework to solve the typical all-to-all problem.
          // Half efficiency along the diagonal.
          if (tgx == 0) {
            warpWorkPos = atomicAdd((int*)&tilePos, 32);
          }
          warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
#ifdef COMPUTE_ENERGY
          int uelec = 0;
          int uvdw = 0;
#endif
          while (warpWorkPos < tileEnd) {

            // Load pair exclusion bit masks.  The high 16 bits of the
            // mask indicate the position of the tile, and are the same
            // for all 16 of the unsigned int masks in each tile.
#ifdef PREFER_AUTO_L1
            unsigned int pairmask = cGms.DVCpairBitMasks[warpWorkPos + tgx];
#else
            unsigned int pairmask = __ldg(&cGms.DVCpairBitMasks[warpWorkPos +
								tgx]);
#endif

            // Initialize atoms and forces
            int katm = (pairmask >> 24)*16 + (tgx & 15);
            int matm = ((pairmask >> 16) & 0xff)*16 + (tgx & 15);
            int loadatm = (tgx < 16) ? katm : matm; 
            float txcrd = xcrd[loadatm];
            float tycrd = ycrd[loadatm];
            float tzcrd = zcrd[loadatm];
#ifdef PREFER_AUTO_L1
            float tqval = cGms.DVCatomQ[atomStart + loadatm];
            int tljid   = cGms.DVCatomLJID[atomStart + loadatm];
#else
            float tqval = __ldg(&cGms.DVCatomQ[atomStart + loadatm]);
            int tljid   = __ldg(&cGms.DVCatomLJID[atomStart + loadatm]);
#endif
#ifdef GBSOLVENT
            float treff    = reff[loadatm];
            float tdeijda = (float)0.0;
            float diagfac = ((katm >> 4) == (matm >> 4)) ?
                            (float)0.5 : (float)1.0;
#endif            
            float txfrc = (float)0.0;
            float tyfrc = (float)0.0;
            float tzfrc = (float)0.0;
            int j;
            for (j = 0; j < 8; j++) {
              int crdSrcLane = (tgx <  16)*(tgx + 16 - j) +
                               (tgx >= 16)*(tgx - 24 + j);
              crdSrcLane += 16 * ((crdSrcLane < 16 && tgx <  16) +
                                  (crdSrcLane <  0 && tgx >= 16));
              float oxcrd = __shfl_sync(0xffffffff, txcrd, crdSrcLane, 32);
              float oycrd = __shfl_sync(0xffffffff, tycrd, crdSrcLane, 32);
              float ozcrd = __shfl_sync(0xffffffff, tzcrd, crdSrcLane, 32);
              float oqval = __shfl_sync(0xffffffff, tqval, crdSrcLane, 32);
              int   oljid = __shfl_sync(0xffffffff, tljid, crdSrcLane, 32);
              
              // Generalized Born calculation
              float fmag = (float)0.0;
              float qq = tqval * oqval;
              float dx = oxcrd - txcrd;
              float dy = oycrd - tycrd;
              float dz = ozcrd - tzcrd;
              float r2 = dx*dx + dy*dy + dz*dz;
#ifdef GBSOLVENT
              float oreff = __shfl_sync(0xffffffff, treff, crdSrcLane, 32);
              int   oatom = __shfl_sync(0xffffffff, loadatm, crdSrcLane, 32);
              float rb2 = treff * oreff;
              float efac = exp(-r2 / ((float)4.0 * rb2));
              float fgbi = (float)1.0 / sqrt(r2 + rb2 * efac);
              float fgbk = -cGms.kappa * cGms.kscale / fgbi;
              float expmkf = exp(fgbk) / cGms.dielectric;
              float dielfac = (float)1.0 - expmkf;
#  ifdef COMPUTE_ENERGY
              if (oatom != loadatm && oatom < natom && loadatm < natom) {
                usolv += (int)(-qq * dielfac * fgbi * FPSCALEnrgMP * diagfac);
              }
#  endif
              float temp4 = fgbi * fgbi * fgbi;
              float temp6 = qq * temp4 * (dielfac + fgbk * expmkf);
              float temp5 = (float)0.5 * efac * temp6 *
                            (rb2 + (float)0.25*r2) * diagfac;
              float odeijda = (float)0.0;
              if (oatom != loadatm && oatom < natom && loadatm < natom) {
                fmag = temp6 * ((float)1.0 - (float)0.25 * efac) * diagfac;
                tdeijda += treff * temp5;
                odeijda = oreff * temp5;
              }
#endif              
              // Shuffling is done.  Now determine whether
              // this interaction is valid.
              if ((pairmask >> (crdSrcLane - 16*(tgx < 16))) & 0x1) {
                int    pljid  = tljid*ntypes + oljid;
#ifdef PREFER_AUTO_L1
                float2 ljterm = cGms.DVCljFtab[(sysID * cGms.ljABoffset) +
                                               pljid];
#else
                float2 ljterm = __ldg(&cGms.DVCljFtab[(sysID *
						       cGms.ljABoffset) +
						      pljid]);
#endif
                float invr  = rsqrt(r2);
                float invr2 = invr  * invr;
                float invr4 = invr2 * invr2;
                float invr8 = invr4 * invr4;
                fmag += -(qq * invr2 * invr) +
                        invr8*((ljterm.x * invr4 * invr2) + ljterm.y);
#ifdef COMPUTE_ENERGY
#  ifdef PREFER_AUTO_L1
                ljterm = cGms.DVCljUtab[(sysID * cGms.ljABoffset) + pljid];
#  else
                ljterm = __ldg(&cGms.DVCljUtab[(sysID * cGms.ljABoffset) +
					       pljid]);
#  endif
                uelec += (int)(qq * invr * FPSCALEnrgMP);
                uvdw  += (int)(((ljterm.x * invr4 * invr2) + ljterm.y) *
                               invr4 * invr2 * FPSCALEnrgMP);
#endif
              }
              float xPairFrc = fmag * dx;
              float yPairFrc = fmag * dy;
              float zPairFrc = fmag * dz;
              txfrc += xPairFrc;
              tyfrc += yPairFrc;
              tzfrc += zPairFrc;
              int frcRetLane = (tgx <  16)*(tgx + 24 - j) +
                               (tgx >= 16)*(tgx - 16 + j);
              frcRetLane += 16 * ((frcRetLane >= 16 && tgx >= 16) -
                                  (frcRetLane >= 32 && tgx <  16));
              txfrc -= __shfl_sync(0xffffffff, xPairFrc, frcRetLane, 32);
              tyfrc -= __shfl_sync(0xffffffff, yPairFrc, frcRetLane, 32);
              tzfrc -= __shfl_sync(0xffffffff, zPairFrc, frcRetLane, 32);
#ifdef GBSOLVENT
              tdeijda += __shfl_sync(0xffffffff, odeijda, frcRetLane, 32);
#endif
            }
	    
            // Commit forces (and energies) back to __shared__
            atomicAdd((int*)&ixfrc[loadatm], (int)(FPSCALEfrc * txfrc));
            atomicAdd((int*)&iyfrc[loadatm], (int)(FPSCALEfrc * tyfrc));
            atomicAdd((int*)&izfrc[loadatm], (int)(FPSCALEfrc * tzfrc));
#ifdef GBSOLVENT
            atomicAdd((int*)&sumdeijda[loadatm], (int)(FPSCALEfrc * tdeijda));
#endif            
            // Increment the tile position
            if (tgx == 0) {
              warpWorkPos = atomicAdd((int*)&tilePos, 32);
            }
            warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
          }
#ifdef COMPUTE_ENERGY
          WarpREDUCE(uelec);
          WarpREDUCE(uvdw);
#  ifdef GBSOLVENT
          WarpREDUCE(usolv);
#  endif
          if (tgx == 0) {
            sElecNrg[threadIdx.x / GRID] = uelec;
            sVdwNrg[threadIdx.x / GRID] = uvdw;
#  ifdef GBSOLVENT
            sSolvNrg[threadIdx.x / GRID] = usolv;
#  endif
          }
#endif
#ifdef GBSOLVENT
          __syncthreads();
          if (threadIdx.x < natom) {
            int readidx = atomStart + threadIdx.x;
#ifdef PREFER_AUTO_L1
            float atmrad = cGms.DVCrborn[readidx] - cGms.GBOffset;
#else
            float atmrad = __ldg(&cGms.DVCrborn[readidx]) - cGms.GBOffset;
#endif
            float fipsi = (float)(psi[threadIdx.x]) / FPSCALErad;
            if (cGms.igb != 1) {
              fipsi *= -atmrad;
            }
#ifdef PREFER_AUTO_L1
            float thi = tanh((cGms.DVCGBalpha[readidx] -
                             (cGms.DVCGBbeta[readidx] * fipsi) +
                             (cGms.DVCGBgamma[readidx] * fipsi * fipsi)) *
                             fipsi);
            sumdeijda[threadIdx.x] *= (cGms.DVCGBalpha[readidx] -
                                       ((float)2.0 * cGms.DVCGBbeta[readidx] *
                                        fipsi) +
                                       ((float)3.0 * cGms.DVCGBgamma[readidx] *
                                        fipsi * fipsi)) *
                                      ((float)1.0 - thi*thi) * atmrad /
                                      cGms.DVCrborn[readidx];
#else
            float thi = tanh((__ldg(&cGms.DVCGBalpha[readidx]) -
                             (__ldg(&cGms.DVCGBbeta[readidx]) * fipsi) +
                             (__ldg(&cGms.DVCGBgamma[readidx]) * fipsi *
			      fipsi)) * fipsi);
            sumdeijda[threadIdx.x] *= (__ldg(&cGms.DVCGBalpha[readidx]) -
                                       ((float)2.0 *
					__ldg(&cGms.DVCGBbeta[readidx]) *
                                        fipsi) +
                                       ((float)3.0 *
					__ldg(&cGms.DVCGBgamma[readidx]) *
                                        fipsi * fipsi)) *
                                      ((float)1.0 - thi*thi) * atmrad /
                                      __ldg(&cGms.DVCrborn[readidx]);
#endif
          }
          __syncthreads();

          // Compute forces due to Born radii derivatives
          if (tgx == 0) {
            warpWorkPos = atomicAdd((int*)&gbDRadPos, 32);
          }
          warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
          while (warpWorkPos < tileEnd) {

            // Load pair exclusion bit masks, as before
#ifdef PREFER_AUTO_L1
            unsigned int pairmask = cGms.DVCpairBitMasks[warpWorkPos + tgx];
#else
            unsigned int pairmask = __ldg(&cGms.DVCpairBitMasks[warpWorkPos +
								tgx]);
#endif
            int katm = (pairmask >> 24)*16 + (tgx & 15);
            int matm = ((pairmask >> 16) & 0xff)*16 + (tgx & 15);
            int loadatm = (tgx < 16) ? katm : matm;
            float txcrd   = xcrd[loadatm];
            float tycrd   = ycrd[loadatm];
            float tzcrd   = zcrd[loadatm];
#ifdef PREFER_AUTO_L1
            float tfs     = cGms.DVCFscreen[atomStart + loadatm];
            int   tneck   = cGms.DVCNeckID[atomStart + loadatm];
            float tbrad   = cGms.DVCrborn[atomStart + loadatm] -
                            cGms.GBOffset;
#else
            float tfs     = __ldg(&cGms.DVCFscreen[atomStart + loadatm]);
            int   tneck   = __ldg(&cGms.DVCNeckID[atomStart + loadatm]);
            float tbrad   = __ldg(&cGms.DVCrborn[atomStart + loadatm]) -
                            cGms.GBOffset;
#endif
            float tsumda  = (float)(sumdeijda[loadatm]) / FPSCALEfrc;
            float tinvrad = (float)1.0 / tbrad;
            float txfrc = (float)0.0;
            float tyfrc = (float)0.0;
            float tzfrc = (float)0.0;
            int ondiagonal = ((katm >> 4) == (matm >> 4));
            int j;
            for (j = 0; j < 8; j++) {
              int crdSrcLane = (tgx <  16)*(tgx + 16 - j) +
                               (tgx >= 16)*(tgx - 24 + j);
              crdSrcLane += 16 * ((crdSrcLane < 16 && tgx <  16) +
                                  (crdSrcLane <  0 && tgx >= 16));
              float oxcrd   = __shfl_sync(0xffffffff, txcrd,  crdSrcLane, 32);
              float oycrd   = __shfl_sync(0xffffffff, tycrd,  crdSrcLane, 32);
              float ozcrd   = __shfl_sync(0xffffffff, tzcrd,  crdSrcLane, 32);
              float ofs     = __shfl_sync(0xffffffff, tfs,    crdSrcLane, 32);
              int   oneck   = __shfl_sync(0xffffffff, tneck,  crdSrcLane, 32);
              float obrad   = __shfl_sync(0xffffffff, tbrad,  crdSrcLane, 32);
              float osumda  = __shfl_sync(0xffffffff, tsumda, crdSrcLane, 32);
              float oinvrad = (float)1.0 / obrad;
              float dx = oxcrd - txcrd;
              float dy = oycrd - tycrd;
              float dz = ozcrd - tzcrd;
              float r2 = dx*dx + dy*dy + dz*dz;
              float invr  = rsqrt(r2);
              float invr2 = invr * invr;
              float r  = r2 * invr;

              // First computation
              float sj  = ofs * obrad;
              float sj2 = sj * sj;
              float tdatmp;
              if (r > (float)4.0 * sj) {
                float tmpsd  = sj2 * invr2;
                float dumbo  = TE + tmpsd*(TF + tmpsd*(TG +
                                                       tmpsd*(TH +
                                                              tmpsd*THH)));
                tdatmp = tmpsd * sj * invr2 * invr2 * dumbo;
              }
              else if (r > tbrad + sj) {
                float temp1  = (float)1.0 / (r2 - sj2);
                tdatmp = temp1 * sj * ((float)-0.5 * invr2 + temp1) +
                         (float)0.25 * invr * invr2 * log((r - sj) /
                                                          (r + sj));
              }
              else if (r > fabs(tbrad - sj)) {
                float temp1  = (float)1.0 / (r + sj);
                float invr3  = invr2 * invr;
                tdatmp = (float)-0.25 * (((float)-0.5 *
                                          (r2 - tbrad*tbrad + sj2) *
                                          invr3 * tinvrad * tinvrad) +
                                         (invr * temp1 * (temp1 - invr)) -
                                         (invr3 * log(tbrad * temp1)));
              }
              else if (tbrad < sj) {
                float temp1  = (float)1.0 / (r2 - sj2);
                tdatmp = (float)-0.5 * ((sj * invr2 * temp1) -
                                        ((float)2.0 * sj * temp1 * temp1) -
                                        ((float)0.5 * invr2 * invr *
                                         log((sj - r) / (sj + r))));
              }
              else {
                tdatmp = (float)0.0;
              }

              // Second computation
              float si  = tfs * tbrad;
              float si2 = si * si;
              float odatmp;
              if (r > (float)4.0 * si) {
                float tmpsd  = si2 * invr2;
                float dumbo  = TE + tmpsd*(TF + tmpsd*(TG +
                                                       tmpsd*(TH +
                                                              tmpsd*THH)));
                odatmp = tmpsd * si * invr2 * invr2 * dumbo;
              }
              else if (r > obrad + si) {
                float temp1  = (float)1.0 / (r2 - si2);
                odatmp = temp1 * si * ((float)-0.5 * invr2 + temp1) +
                         (float)0.25 * invr * invr2 * log((r - si) /
                                                          (r + si));
              }
              else if (r > fabs(obrad - si)) {
                float temp1  = (float)1.0 / (r + si);
                float invr3  = invr2 * invr;
                odatmp = (float)-0.25 * (((float)-0.5 *
                                          (r2 - obrad*obrad + si2) *
                                          invr3 * oinvrad * oinvrad) +
                                         (invr * temp1 * (temp1 - invr)) -
                                         (invr3 * log(obrad * temp1)));
              }
              else if (obrad < si) {
                float temp1  = (float)1.0 / (r2 - si2);
                odatmp = (float)-0.5 * ((si * invr2 * temp1) -
                                        ((float)2.0 * si * temp1 * temp1) -
                                        ((float)0.5 * invr2 * invr *
                                         log((si - r) / (si + r))));
              }
              else {
                odatmp = (float)0.0;
              }

              // Neck GB contributions
              if ((cGms.igb == 7 || cGms.igb == 8) &&
                  r < tbrad + obrad + cGms.GBNeckCut) {

                // First computation: i -> j
#ifdef PREFER_AUTO_L1
                float2 neckData = cGms.DVCneckFactors[21*tneck + oneck];
#else
                float2 neckData = __ldg(&cGms.DVCneckFactors[21*tneck +
							     oneck]);
#endif
                float mdist = r - neckData.x;
                float mdist2 = mdist * mdist;
                float mdist6 = mdist2 * mdist2 * mdist2;
                float temp1 = (float)1.0 + mdist2 + (float)0.3*mdist6;
                temp1 = temp1 * temp1 * r;
                tdatmp += (((float)2.0*mdist +
                            (float)1.8*mdist2*mdist2*mdist) *
                           neckData.y * cGms.GBNeckScale) / temp1;

                // Second computation: j -> i
#ifdef PREFER_AUTO_L1
                neckData = cGms.DVCneckFactors[21*oneck + tneck];
#else
                neckData = __ldg(&cGms.DVCneckFactors[21*oneck + tneck]);
#endif
                mdist = r - neckData.x;
                mdist2 = mdist * mdist;
                mdist6 = mdist2 * mdist2 * mdist2;
                temp1 = (float)1.0 + mdist2 + (float)0.3*mdist6;
                temp1 = temp1 * temp1 * r;
                odatmp += (((float)2.0*mdist +
                            (float)1.8*mdist2*mdist2*mdist) *
                           neckData.y * cGms.GBNeckScale) / temp1;
              }
              float fmag = (tdatmp * tsumda) + (odatmp * osumda);

              // Account for diagonal tiles
              if (ondiagonal) {
                fmag *= (float)0.5;
                if (j == 0 && tgx < 16) {
                  fmag = (float)0.0;
                }
              }
              float xPairFrc = fmag * dx;
              float yPairFrc = fmag * dy;
              float zPairFrc = fmag * dz;

              // Contribute the derivatives to the forces in registers
              txfrc -= xPairFrc;
              tyfrc -= yPairFrc;
              tzfrc -= zPairFrc;
              int frcRetLane = (tgx <  16)*(tgx + 24 - j) +
                (tgx >= 16)*(tgx - 16 + j);
              frcRetLane += 16 * ((frcRetLane >= 16 && tgx >= 16) -
                                  (frcRetLane >= 32 && tgx <  16));
              txfrc += __shfl_sync(0xffffffff, xPairFrc, frcRetLane, 32);
              tyfrc += __shfl_sync(0xffffffff, yPairFrc, frcRetLane, 32);
              tzfrc += __shfl_sync(0xffffffff, zPairFrc, frcRetLane, 32);
            }
            
            // Commit forces (no energies this time) back to __shared__
            atomicAdd((int*)&ixfrc[loadatm], (int)(FPSCALEfrc * txfrc));
            atomicAdd((int*)&iyfrc[loadatm], (int)(FPSCALEfrc * tyfrc));
            atomicAdd((int*)&izfrc[loadatm], (int)(FPSCALEfrc * tzfrc));

            // Increment the tile position
            if (tgx == 0) {
              warpWorkPos = atomicAdd((int*)&gbDRadPos, 32);
            }
            warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
          }
#endif          
          // Dihedral computation.  This introduces some funky pre-processor
          // directives, but it's critical to keep the kernel from locking.
          // In the simple force computation case, each dihedral can get
          // computed without regard to any of the others, but if energies
          // are being computed then the entire warp must pass into the
          // while loop together for the shuffle in WarpREDUCE().
          if (tgx == 0) {
            warpWorkPos = atomicAdd((int*)&dihePos, 32);
          }
          warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
#ifdef COMPUTE_ENERGY
          int udihe = 0;
#endif
          while (warpWorkPos < diheEnd) {
            if (warpWorkPos + tgx < diheEnd) {
#ifdef PREFER_AUTO_L1
              uint2 myDiheID = cGms.DVCdiheIDs[warpWorkPos + tgx];
#else
              uint2 myDiheID = __ldg(&cGms.DVCdiheIDs[warpWorkPos + tgx]);
#endif
              unsigned int jatm = myDiheID.x;
              unsigned int katm = ((jatm >> 10) & 0x3ff);
              unsigned int matm = (jatm & 0x3ff);
              jatm = (jatm >> 20);
              unsigned int patm = myDiheID.y;
              float Npdc = (float)(patm & 0xffff);
              patm = (patm >> 16);

              // Compute displacements
              float ab[3], bc[3], cd[3], crabbc[3], crbccd[3];
              ab[0] = xcrd[katm] - xcrd[jatm];
              ab[1] = ycrd[katm] - ycrd[jatm];
              ab[2] = zcrd[katm] - zcrd[jatm];
              bc[0] = xcrd[matm] - xcrd[katm];
              bc[1] = ycrd[matm] - ycrd[katm];
              bc[2] = zcrd[matm] - zcrd[katm];
              cd[0] = xcrd[patm] - xcrd[matm];
              cd[1] = ycrd[patm] - ycrd[matm];
              cd[2] = zcrd[patm] - zcrd[matm];

              // A device function would be nice here, but it would only be
              // needed three times and the register count absolutely must stay
              // below 64 per thread.  Macro to reach minimal register usage.
              DvcCROSSPf(ab, bc, crabbc);
              DvcCROSSPf(bc, cd, crbccd);
              float fmag;
              if ((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] +
                   crabbc[2]*crabbc[2] >= (float)1.0e-4) &&
                  (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] +
                   crbccd[2]*crbccd[2] >= (float)1.0e-4)) {
                float costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] +
                                 crabbc[2]*crbccd[2];
                costheta /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] +
                                  crabbc[2]*crabbc[2]) *
                                 (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] +
                                  crbccd[2]*crbccd[2]));
                float scr[3];
                DvcCROSSPf(crabbc, crbccd, scr);
                costheta = (costheta < (float)-1.0) ?
                           (float)-1.0 : (costheta > (float)1.0) ?
                                         (float)1.0 : costheta;
                float theta;
                if (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > (float)0.0) {
                  theta = acosf(costheta);
                }
                else {
                  theta = -acosf(costheta);
                }
#ifdef PREFER_AUTO_L1
                float2 dhinfo = cGms.DVCdiheInfo[warpWorkPos + tgx];
#else
                float2 dhinfo = __ldg(&cGms.DVCdiheInfo[warpWorkPos + tgx]);
#endif
                float sangle = Npdc*theta - dhinfo.y;
                fmag = FPSCALEfrc * dhinfo.x * Npdc * sinf(sangle);
                float mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
                float invab = (float)1.0 / mgab;
                float mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
                float invbc = (float)1.0 / mgbc;
                float mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
                float invcd = (float)1.0 / mgcd;
                float cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]) *
                             invab*invbc;
                float isinb2 = (cosb*cosb < (float)0.9999) ?
                               (float)1.0 / ((float)1.0 - cosb*cosb) :
                               (float)0.0;
                float cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) *
                             invbc * invcd;
                float isinc2 = (cosc*cosc < (float)0.9999) ?
                               (float)1.0 / ((float)1.0 - cosc*cosc) :
                               (float)0.0;
                isinb2 *= fmag;
                isinc2 *= fmag;
                float invabc = invab * invbc;
                float invbcd = invbc * invcd;
                crabbc[0] *= invabc;
                crabbc[1] *= invabc;
                crabbc[2] *= invabc;
                crbccd[0] *= invbcd;
                crbccd[1] *= invbcd;
                crbccd[2] *= invbcd;

                // Transform the dihedral forces to cartesian coordinates
                float fa = -invab * isinb2;
                float fb1 = (mgbc - mgab*cosb) * invabc * isinb2;
                float fb2 = cosc * invbc * isinc2;
                float fc1 = (mgbc - mgcd*cosc) * invbcd * isinc2;
                float fc2 = cosb * invbc * isinb2;
                float fd = -invcd * isinc2;

                // Apply the dihedral forces.  The integer conversion here is
                // going to cause rounding that violates Newton's third law,
                // but it's not clear how to do any better.
                atomicAdd((int*)&ixfrc[jatm], (int)(crabbc[0] * fa));
                atomicAdd((int*)&ixfrc[katm],
                          (int)(fb1*crabbc[0] - fb2*crbccd[0]));
                atomicAdd((int*)&ixfrc[matm],
                          (int)(-fc1*crbccd[0] + fc2*crabbc[0]));
                atomicAdd((int*)&ixfrc[patm], (int)(-fd * crbccd[0]));
                atomicAdd((int*)&iyfrc[jatm], (int)(crabbc[1] * fa));
                atomicAdd((int*)&iyfrc[katm],
                          (int)(fb1*crabbc[1] - fb2*crbccd[1]));
                atomicAdd((int*)&iyfrc[matm],
                          (int)(-fc1*crbccd[1] + fc2*crabbc[1]));
                atomicAdd((int*)&iyfrc[patm], (int)(-fd * crbccd[1]));
                atomicAdd((int*)&izfrc[jatm], (int)(crabbc[2] * fa));
                atomicAdd((int*)&izfrc[katm],
                          (int)(fb1*crabbc[2] - fb2*crbccd[2]));
                atomicAdd((int*)&izfrc[matm],
                          (int)(-fc1*crbccd[2] + fc2*crabbc[2]));
                atomicAdd((int*)&izfrc[patm], (int)(-fd * crbccd[2]));
#ifdef COMPUTE_ENERGY
                udihe += FPSCALEnrgHP * dhinfo.x * ((float)1.0 + cos(sangle));
#endif
              }
            }

            // Increment the work counter
            if (tgx == 0) {
              warpWorkPos = atomicAdd((int*)&dihePos, 32);
            }
            warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
          }
#ifdef COMPUTE_ENERGY
          WarpREDUCE(udihe);
          if (tgx == 0) {
            sDiheNrg[threadIdx.x / GRID] = udihe;
          }
#endif
          // The most stable thing to do is to compute the forces at
          // nominal strength, accumulate them completely, then scale
          // them here.  Fewer mults in the inner loops, but synchronization.
          if (cGms.slowFrcMult > 1) {
            __syncthreads();
            if (threadIdx.x < natom) {
              ixfrc[threadIdx.x] *= cGms.slowFrcMult;
              iyfrc[threadIdx.x] *= cGms.slowFrcMult;
              izfrc[threadIdx.x] *= cGms.slowFrcMult;
            }
            __syncthreads();          
          }
        }
        // ^^ End slow force computations inside the cycle zero conditional ^^

        // Compute the fast force components, beginning with 1:4 attenuations.
        // Compute the energy only once, on cycle zero of the MTS for this
        // step.  In this manner, all energies will correspond to the
        // coordinates used by the dihedral and non-bonded routines.
#ifdef COMPUTE_ENERGY
        int uelecAttn = 0;
        int uvdwAttn = 0;
#endif
        if (tgx == 0) {
          warpWorkPos = atomicAdd((int*)&attnPos, 32);
        }
        warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
        while (warpWorkPos < attnEnd) {
          if (warpWorkPos + tgx < attnEnd) {
#ifdef PREFER_AUTO_L1
            int jatm = cGms.DVCattnIDs[warpWorkPos + tgx];
#else
            int jatm = __ldg(&cGms.DVCattnIDs[warpWorkPos + tgx]);
#endif
            int katm = (jatm & 0xffff);
            jatm >>= 16;
#ifdef PREFER_AUTO_L1
            float2 scaling = cGms.DVCattnFactors[warpWorkPos + tgx];
#else
            float2 scaling = __ldg(&cGms.DVCattnFactors[warpWorkPos + tgx]);
#endif
            float dx = xcrd[katm] - xcrd[jatm];
            float dy = ycrd[katm] - ycrd[jatm];
            float dz = zcrd[katm] - zcrd[jatm];
            float r2 = dx*dx + dy*dy + dz*dz;
            float invr = rsqrt(r2);
            float invr2 = invr * invr;
            float invr4 = invr2 * invr2;
            float invr6 = invr4 * invr2;
            float invr8 = invr4 * invr4;
#ifdef PREFER_AUTO_L1
            float qq = cGms.DVCatomQ[atomStart + katm] *
                       cGms.DVCatomQ[atomStart + jatm];
            int pljid = (cGms.DVCatomLJID[atomStart + jatm] * ntypes) +
                        cGms.DVCatomLJID[atomStart + katm];
            float2 ljterm = cGms.DVCljFtab[(sysID * cGms.ljABoffset) + pljid];
#else
            float qq = __ldg(&cGms.DVCatomQ[atomStart + katm]) *
                       __ldg(&cGms.DVCatomQ[atomStart + jatm]);
            int pljid = (__ldg(&cGms.DVCatomLJID[atomStart + jatm]) * ntypes) +
                        __ldg(&cGms.DVCatomLJID[atomStart + katm]);
            float2 ljterm = __ldg(&cGms.DVCljFtab[(sysID * cGms.ljABoffset) +
						  pljid]);
#endif
            float fmag = -(scaling.x * qq * invr2 * invr) +
                         (scaling.y * invr8 * ((ljterm.x * invr6) + ljterm.y));
            fmag *= FPSCALEfrc;
#ifdef COMPUTE_ENERGY
#ifdef PREFER_AUTO_L1
            ljterm = cGms.DVCljUtab[(sysID * cGms.ljABoffset) + pljid];
#else
            ljterm = __ldg(&cGms.DVCljUtab[(sysID * cGms.ljABoffset) + pljid]);
#endif
            uelecAttn += (int)(FPSCALEnrgMP * scaling.x * qq * invr);
            uvdwAttn  += (int)(FPSCALEnrgMP * scaling.y *
                               ((ljterm.x * invr6) + ljterm.y) * invr6);
#endif            
            int ifx = (int)(fmag * dx);
            int ify = (int)(fmag * dy);
            int ifz = (int)(fmag * dz);
            atomicAdd((int*)&ixfrc[jatm],  ifx);
            atomicAdd((int*)&iyfrc[jatm],  ify);
            atomicAdd((int*)&izfrc[jatm],  ifz);
            atomicAdd((int*)&ixfrc[katm], -ifx);
            atomicAdd((int*)&iyfrc[katm], -ify);
            atomicAdd((int*)&izfrc[katm], -ifz);          
          }
        
          // Increment the tile position
          if (tgx == 0) {
            warpWorkPos = atomicAdd((int*)&attnPos, 32);
          }
          warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
        }
#ifdef COMPUTE_ENERGY
        if (i == 0) {
          WarpREDUCE(uelecAttn);
          WarpREDUCE(uvdwAttn);
          if (tgx == 0) {
            sElecNrg[threadIdx.x / GRID] += uelecAttn;
            sVdwNrg[threadIdx.x / GRID] += uvdwAttn;
          }
        }
#endif
        // Bond stretching interactions
#ifdef COMPUTE_ENERGY
        int ubond = 0;
#endif
        if (tgx == 0) {
          warpWorkPos = atomicAdd((int*)&bondPos, 32);
        }
        warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
        while (warpWorkPos < bondEnd) {
          if (warpWorkPos + tgx < bondEnd) {
#ifdef PREFER_AUTO_L1
            int jatm      = cGms.DVCbondIDs[warpWorkPos + tgx];
            float2 basics = cGms.DVCbondBasics[warpWorkPos + tgx];
            float4 augs   = cGms.DVCbondAugs[warpWorkPos + tgx];
#else
            int jatm      = __ldg(&cGms.DVCbondIDs[warpWorkPos + tgx]);
            float2 basics = __ldg(&cGms.DVCbondBasics[warpWorkPos + tgx]);
            float4 augs   = __ldg(&cGms.DVCbondAugs[warpWorkPos + tgx]);
#endif
            int katm = (jatm >> 16);
            jatm = (jatm & 0xffff);
            float dx = xcrd[katm] - xcrd[jatm];
            float dy = ycrd[katm] - ycrd[jatm];
            float dz = zcrd[katm] - zcrd[jatm];
            float r2 = dx*dx + dy*dy + dz*dz;
            float r = sqrt(r2);
            float dl = (basics.y >= (float)0.0) ? basics.y - r : r + basics.y;
            float dlpu = (r > augs.y) ? augs.y - r : (float)0.0;
            float dlpr = (r < augs.w) ? augs.w - r : (float)0.0;
            float fmag = -2.0 * (basics.x*dl + augs.x*dlpu + augs.z*dlpr) / r;
            fmag *= FPSCALEfrc;
            int ifx = (int)(fmag * dx);
            int ify = (int)(fmag * dy);
            int ifz = (int)(fmag * dz);
            atomicAdd((int*)&ixfrc[jatm],  ifx);
            atomicAdd((int*)&iyfrc[jatm],  ify);
            atomicAdd((int*)&izfrc[jatm],  ifz);
            atomicAdd((int*)&ixfrc[katm], -ifx);
            atomicAdd((int*)&iyfrc[katm], -ify);
            atomicAdd((int*)&izfrc[katm], -ifz);          
#ifdef COMPUTE_ENERGY
            ubond += FPSCALEnrgHP * ((basics.x * dl * dl) +
                                     (augs.x * dlpu * dlpu) +
                                     (augs.z * dlpr * dlpr));
#endif
          }
          
          // Increment the work counter
          if (tgx == 0) {
            warpWorkPos = atomicAdd((int*)&bondPos, 32);
          }
          warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
        }
#ifdef COMPUTE_ENERGY
        if (i == 0) {
          WarpREDUCE(ubond);
          if (tgx == 0) {
            sBondNrg[threadIdx.x / GRID] = ubond;
          }
        }
#endif
        // Angle bending interactions
#ifdef COMPUTE_ENERGY
        int uangl = 0;
#endif
        if (tgx == 0) {
          warpWorkPos = atomicAdd((int*)&anglPos, 32);
        }
        warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
        while (warpWorkPos < anglEnd) {
          if (warpWorkPos + tgx < anglEnd) {
#ifdef PREFER_AUTO_L1
            unsigned int jatm = cGms.DVCanglIDs[warpWorkPos + tgx];
            float2 ainfo      = cGms.DVCanglInfo[warpWorkPos + tgx];
#else
            unsigned int jatm = __ldg(&cGms.DVCanglIDs[warpWorkPos + tgx]);
            float2 ainfo      = __ldg(&cGms.DVCanglInfo[warpWorkPos + tgx]);
#endif
            unsigned int katm = ((jatm >> 10) & 0x3ff);
            unsigned int matm = (jatm & 0x3ff);
            jatm = (jatm >> 20);
            float bax = xcrd[jatm] - xcrd[katm];
            float bay = ycrd[jatm] - ycrd[katm];
            float baz = zcrd[jatm] - zcrd[katm];
            float bcx = xcrd[matm] - xcrd[katm];
            float bcy = ycrd[matm] - ycrd[katm];
            float bcz = zcrd[matm] - zcrd[katm];
            float mgba = bax*bax + bay*bay + baz*baz;
            float mgbc = bcx*bcx + bcy*bcy + bcz*bcz;
            float invbabc = 1.0/sqrt(mgba*mgbc);
            float costheta = (bax*bcx + bay*bcy + baz*bcz) * invbabc;
            costheta = (costheta < (float)-1.0) ? (float)-1.0 :
                       (costheta > (float)1.0)  ? (float)1.0  : costheta;
            float theta = acos(costheta);
            float dtheta = theta - ainfo.y;
            float dA = -2.0 * ainfo.x * dtheta /
                       sqrt((float)1.0 - costheta*costheta);
            dA *= FPSCALEfrc;
            float sqba = dA / mgba;
            float sqbc = dA / mgbc;
            float mbabc = dA * invbabc;
            int iadfx = (bcx*mbabc - costheta*bax*sqba);
            int icdfx = (bax*mbabc - costheta*bcx*sqbc);
            int iadfy = (bcy*mbabc - costheta*bay*sqba);
            int icdfy = (bay*mbabc - costheta*bcy*sqbc);
            int iadfz = (bcz*mbabc - costheta*baz*sqba);
            int icdfz = (baz*mbabc - costheta*bcz*sqbc);
            atomicAdd((int*)&ixfrc[jatm], -iadfx);
            atomicAdd((int*)&iyfrc[jatm], -iadfy);
            atomicAdd((int*)&izfrc[jatm], -iadfz);
            atomicAdd((int*)&ixfrc[katm], iadfx + icdfx);
            atomicAdd((int*)&iyfrc[katm], iadfy + icdfy);
            atomicAdd((int*)&izfrc[katm], iadfz + icdfz);
            atomicAdd((int*)&ixfrc[matm], -icdfx);
            atomicAdd((int*)&iyfrc[matm], -icdfy);
            atomicAdd((int*)&izfrc[matm], -icdfz);
#ifdef COMPUTE_ENERGY
            uangl += FPSCALEnrgHP * ainfo.x * dtheta * dtheta;
#endif
          }
          
          // Increment the work counter
          if (tgx == 0) {
            warpWorkPos = atomicAdd((int*)&anglPos, 32);
          }
          warpWorkPos = __shfl_sync(0xffffffff, warpWorkPos, 0, 32);
        }
#ifdef COMPUTE_ENERGY
        if (i == 0) {
          WarpREDUCE(uangl);
          if (tgx == 0) {
            sAnglNrg[threadIdx.x / GRID] = uangl;
          }
        }
#endif
        __syncthreads();
        // ^^ End of fast force computations, and all force   ^^
        //    computations for this MTS minor step.
        
        // Velocity Verlet half kick
#ifdef COMPUTE_ENERGY
        float ukine;
#endif
        if (threadIdx.x < natom) {
#ifdef PREFER_AUTO_L1
          float hmdt = cGms.DVCatomHDTM[atomStart + threadIdx.x];
#else
          float hmdt = __ldg(&cGms.DVCatomHDTM[atomStart + threadIdx.x]);
#endif
          if (cGms.Tstat.active == 1) {
            float atmmass = cGms.hdtm2mass / hmdt;
            float rsd = FPSCALEfrc * sysSdfac * sqrt(atmmass);
            int rndidx = randStart +
                         3*((((stepidx * cGms.slowFrcMult) + i) %
                             cGms.Tstat.refresh) + 1)*natom +
                         threadIdx.x;
            xvel[threadIdx.x] = (hmdt*((float)ixfrc[threadIdx.x] +
                                       (rsd * cGms.DVCprngHeap[rndidx])) +
                                 xvel[threadIdx.x]) * cGms.Tstat.c_implic;
            rndidx += natom;
            yvel[threadIdx.x] = (hmdt*((float)iyfrc[threadIdx.x] +
                                       (rsd * cGms.DVCprngHeap[rndidx])) +
                                 yvel[threadIdx.x]) * cGms.Tstat.c_implic;
            rndidx += natom;
            zvel[threadIdx.x] = (hmdt*((float)izfrc[threadIdx.x] +
                                       (rsd * cGms.DVCprngHeap[rndidx])) +
                                 zvel[threadIdx.x]) * cGms.Tstat.c_implic;
          }
          else {
            xvel[threadIdx.x] += hmdt * (float)ixfrc[threadIdx.x];
            yvel[threadIdx.x] += hmdt * (float)iyfrc[threadIdx.x];
            zvel[threadIdx.x] += hmdt * (float)izfrc[threadIdx.x];
          }
        }

        // Apply velocity constraints
#ifdef GO_RATTLE
        if (cGms.rattle == 1) {

          // Stash forces in registers
          int TLixf, TLiyf, TLizf;
          if (threadIdx.x < natom) {
            TLixf = ixfrc[threadIdx.x];
            TLiyf = iyfrc[threadIdx.x];
            TLizf = izfrc[threadIdx.x];
            ixfrc[threadIdx.x] = 0;
            iyfrc[threadIdx.x] = 0;
            izfrc[threadIdx.x] = 0;
          }
          if (threadIdx.x == 0) {
            cnstPos = cnstStart;
          }
          __syncthreads();

          // Perform constraints, again using the integer force arrays
          // to accumulate perturbations
          int warpCnstPos;
          if (tgx == 0) {
            warpCnstPos = atomicAdd((int*)&cnstPos, 32);
          }
          warpCnstPos = __shfl_sync(0xffffffff, warpCnstPos, 0, 32);
          while (warpCnstPos < cnstEnd) {
#ifdef PREFER_AUTO_L1
            uint2 params = cGms.DVCcnstInfo[warpCnstPos + tgx];
#else
            uint2 params = __ldg(&cGms.DVCcnstInfo[warpCnstPos + tgx]);
#endif
            int atma = params.x;
            bool fixlen = (atma != 0xffffffff);
            unsigned int engaged = __ballot_sync(0xffffffff, fixlen);
            if (fixlen) {
              int atmb = (atma & 0xffff);
              atma >>= 16;
#ifdef PREFER_AUTO_L1
              float rma = cGms.DVCatomHDTM[atomStart + atma] * cGms.hdtm2invm;
              float rmb = cGms.DVCatomHDTM[atomStart + atmb] * cGms.hdtm2invm;
#else
              float rma = __ldg(&cGms.DVCatomHDTM[atomStart + atma]) *
                          cGms.hdtm2invm;
              float rmb = __ldg(&cGms.DVCatomHDTM[atomStart + atmb]) *
                          cGms.hdtm2invm;
#endif
              double l0 = (double)params.y / 1.0e8;
              l0 *= l0;
              float l0f = (float)l0;
              int niter = 0;
              unsigned int atwork = engaged;
              while (niter < cGms.maxRattleIter && atwork > 0) {
                double rx = xcrd[atmb] - xcrd[atma];
                double ry = ycrd[atmb] - ycrd[atma];
                double rz = zcrd[atmb] - zcrd[atma];
                double vx = (xvel[atmb] +
                             ((double)ixfrc[atmb] * FPSCALEicnv)) -
                            (xvel[atma] + ((double)ixfrc[atma] * FPSCALEicnv));
                double vy = (yvel[atmb] +
                             ((double)iyfrc[atmb] * FPSCALEicnv)) -
                            (yvel[atma] + ((double)iyfrc[atma] * FPSCALEicnv));
                double vz = (zvel[atmb] +
                             ((double)izfrc[atmb] * FPSCALEicnv)) -
                            (zvel[atma] + ((double)izfrc[atma] * FPSCALEicnv));
                double dot = rx*vx + ry*vy + rz*vz;
                float term = (float)-1.2 * (float)dot / ((rma + rmb) * l0f);
                bool violation = (fabs(term) > cGms.rtoldt);
                atwork = __ballot_sync(engaged, violation);
                if (violation) {
                  float xterm = (float)rx * term;
                  float yterm = (float)ry * term;
                  float zterm = (float)rz * term;
                  atomicAdd((int*)&ixfrc[atma], (-xterm * rma * FPSCALEcnv));
                  atomicAdd((int*)&iyfrc[atma], (-yterm * rma * FPSCALEcnv));
                  atomicAdd((int*)&izfrc[atma], (-zterm * rma * FPSCALEcnv));
                  atomicAdd((int*)&ixfrc[atmb], ( xterm * rmb * FPSCALEcnv));
                  atomicAdd((int*)&iyfrc[atmb], ( yterm * rmb * FPSCALEcnv));
                  atomicAdd((int*)&izfrc[atmb], ( zterm * rmb * FPSCALEcnv));
                }
                __syncwarp(engaged);
                niter++;
              }
            }

            // Increment the counter
            if (tgx == 0) {
              warpCnstPos = atomicAdd((int*)&cnstPos, 32);
            }
            warpCnstPos = __shfl_sync(0xffffffff, warpCnstPos, 0, 32);
          }
          __syncthreads();

          // Add the perturbations back to the velocity arrays,
          // then un-stash the forces.
          if (threadIdx.x < natom) {
            xvel[threadIdx.x] += (double)ixfrc[threadIdx.x] * FPSCALEicnv;
            yvel[threadIdx.x] += (double)iyfrc[threadIdx.x] * FPSCALEicnv;
            zvel[threadIdx.x] += (double)izfrc[threadIdx.x] * FPSCALEicnv;
            ixfrc[threadIdx.x] = TLixf;
            iyfrc[threadIdx.x] = TLiyf;
            izfrc[threadIdx.x] = TLizf;
          }
        }
#endif
        // Compute kinetic energy and temperature
#ifdef COMPUTE_ENERGY
        if (i == 0) {
          if (threadIdx.x < natom) {
#ifdef PREFER_AUTO_L1
            float hmdt = cGms.DVCatomHDTM[atomStart + threadIdx.x];
#else
            float hmdt = __ldg(&cGms.DVCatomHDTM[atomStart + threadIdx.x]);
#endif
            ukine = (xvel[threadIdx.x] * xvel[threadIdx.x]) +
                    (yvel[threadIdx.x] * yvel[threadIdx.x]) +
                    (zvel[threadIdx.x] * zvel[threadIdx.x]);
            if (cGms.Tstat.active == 1) {
              ukine *= (float)1.0 + ((float)0.5 * cGms.Tstat.gamma_ln *
                                     cGms.dt);
            }
            ukine *= ((float)0.25 * sqrt((float)418.4) * cGms.dt /
                      FPSCALEfrc) / hmdt;
          }
          else {
            ukine = (float)0.0;
          }
          WarpREDUCE(ukine);
          if (tgx == 0) {
            sKineNrg[threadIdx.x / GRID] = ukine;
          }
        }
#endif
        __syncthreads();
      }

      // Increment the step counter
#ifndef COMPUTE_ENERGY
      __syncthreads();
      if (threadIdx.x == 0) {
        stepidx++;
      }
      __syncthreads();
    }
#else
    if (threadIdx.x < GRID) {
      int utmp = (tgx < THREAD_COUNT / GRID) ? sBondNrg[tgx] : 0;
      WarpREDUCE(utmp);
      if (tgx == 0) {
        cGms.Ubond.DevcData[(sysID * cGms.nsgmdout) + sgc] =
          (double)utmp / FPSCALEnrgHP;
      }
      float ukine = (tgx < THREAD_COUNT / GRID) ? sKineNrg[tgx] : (float)0.0;
      WarpREDUCE(ukine);
      if (tgx == 0) {
        cGms.Ukine.DevcData[(sysID * cGms.nsgmdout) + sgc] = (double)ukine;
      }
    }
    else if (threadIdx.x < GRIDx2) {
      int utmp = (tgx < THREAD_COUNT / GRID) ? sAnglNrg[tgx] : 0;
      WarpREDUCE(utmp);
      if (tgx == 0) {
        cGms.Uangl.DevcData[(sysID * cGms.nsgmdout) + sgc] =
          (double)utmp / FPSCALEnrgHP;
      }
    }
    else if (threadIdx.x < GRIDx3) {
      int utmp = (tgx < THREAD_COUNT / GRID) ? sDiheNrg[tgx] : 0;
      WarpREDUCE(utmp);
      if (tgx == 0) {
        cGms.Udihe.DevcData[(sysID * cGms.nsgmdout) + sgc] =
          (double)utmp / FPSCALEnrgHP;
      }
    }
    else if (threadIdx.x < GRIDx4) {
      int utmp = (tgx < THREAD_COUNT / GRID) ? sElecNrg[tgx] : 0;
      utmp >>= MP2LP_BITS;
      WarpREDUCE(utmp);
      if (tgx == 0) {
        cGms.Uelec.DevcData[(sysID * cGms.nsgmdout) + sgc] =
          (double)utmp / FPSCALEnrgLP;
      }
    }
    else if (threadIdx.x < GRIDx5) {
      int utmp = (tgx < THREAD_COUNT / GRID) ? sVdwNrg[tgx] : 0;
      utmp >>= MP2LP_BITS;
      WarpREDUCE(utmp);
      if (tgx == 0) {
        cGms.Uvdw.DevcData[(sysID * cGms.nsgmdout) + sgc] =
          (double)utmp / FPSCALEnrgLP;
      }
    }
    else if (threadIdx.x < GRIDx6) {
      int utmp = (tgx < THREAD_COUNT / GRID) ? sSolvNrg[tgx] : 0;
      utmp >>= MP2LP_BITS;
      WarpREDUCE(utmp);
      if (tgx == 0) {
        cGms.Usolv.DevcData[(sysID * cGms.nsgmdout) + sgc] =
          (double)utmp / FPSCALEnrgLP;
      }
    }
#endif

    // Store coordinates, velocities, and forces in global memory
    if (threadIdx.x < natom) {
      int storeidx = atomStart + threadIdx.x;
      cGms.DVCatomX[storeidx] = xcrd[threadIdx.x];
      cGms.DVCatomY[storeidx] = ycrd[threadIdx.x];
      cGms.DVCatomZ[storeidx] = zcrd[threadIdx.x];      
      cGms.DVCvelcX[storeidx] = xvel[threadIdx.x];
      cGms.DVCvelcY[storeidx] = yvel[threadIdx.x];
      cGms.DVCvelcZ[storeidx] = zvel[threadIdx.x];
      cGms.DVCfrcX[storeidx]  = ixfrc[threadIdx.x];
      cGms.DVCfrcY[storeidx]  = iyfrc[threadIdx.x];
      cGms.DVCfrcZ[storeidx]  = izfrc[threadIdx.x];
    }
    
    // Increment the system counter
    __syncthreads();
    if (threadIdx.x == 0) {
#ifdef COMPUTE_ENERGY
      sysID = atomicAdd(&cGms.DVCsystemPos[sgc], 1);
#else
      sysID = atomicAdd(&cGms.DVCsystemPos[sgc + cGms.nsgmdout], 1);
#endif
    }
    __syncthreads();
  }
}
#ifdef PREFER_AUTO_L1
#  undef PREFER_AUTO_L1
#endif
