/* SGCT-Based-Fault-Tolerant-Taxila-LBM Code source file.
   Copyright (c) 2015, Md Mohsin Ali. All rights reserved.
   Licensed under the terms of the BSD License as described in the LICENSE_FT_CODE file.
   This comment must be retained in any redistributions of this source file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <mpi.h>
#include <assert.h>
#ifdef TAU_PROF
#include <Profile/Profiler.h>
#endif

#include "FTGridCombine3D.h"
#include "GcpGeneral.cpp"    // solves general coefficient problem

#define GSTAG 5 /*message tag for gather & scatter comms*/
#define NOT_SET -1

#define FAILURE_ON_EXTRA_LAYER   3
#define MAX_GRID_FAILED          nGrids()
#define THREE_D                  3 
#define TWO_D                    2

//#define PRINT_gxG_getGridId

// return which grid id global rank is in
int FTGridCombine3D::getGid(int globalRank) {
   int rank = localCombRank(globalRank);   
   assert(rank != -1);
   if (!(0 <= rank && rank < nprocs)) {
      printf("%d: getGid(bad rank=%d), nprocs=%d\n", myrank, rank, nprocs);
   }   
   assert (0 <= rank && rank < nprocs); 
   int v = gIds[rank];
   return (v);
}

int FTGridCombine3D::nGrids() {
  if (haveExtraGrids) {
     if (is2D())
        return (4*level - 6);                 // 2 extra layers
     else
        return (2*level*level - 4*level + 4); // 1 extra plane
  }
  else
     return (GridCombine3D::nGrids());        // call base routine
  printf("Error in FTGridCombine3D::nGrids().\n");
  exit(1);
}

int FTGridCombine3D::diagRank(int gid) {
  if (is2D()) {
     int nGridsTwoDiags = (2*level - 1);
     int nGridsThreeDiags = (3*level - 3);
     if (gid >= 0 && gid < nGridsTwoDiags) 
        return (gid >= level); 
     else if (haveExtraGrids) { // 2 extra layers
        if (gid >= nGridsTwoDiags && gid < nGridsThreeDiags)
           return 2;      
        else if (gid >= nGridsThreeDiags && gid < nGrids())
           return 3;      
        else // gid >= nGrids() which should happen if separate group of processes
             // used for combination
           return 4;      
     }        
  }
  else {
     if (gid < (level+1)*level/2)
        return 0;
     else if (gid < level*level)
        return 1;
     else if (gid < (level*level + (level-1)*(level-2)/2))  
        return 2;
     else if (haveExtraGrids && gid >= (level*level + (level-1)*(level-2)/2)) // 1 extra plane
        return 3;      
  }
  printf("Error in FTGridCombine3D::diagRank().\n");
  exit(1);
}

double FTGridCombine3D::getCombCoeff(int g) {
  if (haveExtraGrids) {
     double c = listCoeffs[g];
     if (g != gidSingle  &&  gidSingle != GID_ALL)
       c = 0.0;
     else if (g == gidSingle  &&  gidSingle != GID_ALL)
       c = 1.0;
     if (isOverrideCoeff && c != 0.0)
       c = overrideCoeff;
     //printf("getCombCoeff(%d)=%f diagRank=%d is2d=%d\n", g, c, diagRank(g), is2D());
     return (c); 
  }
  else {
     double c = (diagRank(g)!=1)? +1.0: (is2D()? -1.0: -2.0); // default
     if (g != gidSingle  &&  gidSingle != GID_ALL)
       c = 0.0;
     else if (g == gidSingle  &&  gidSingle != GID_ALL)
       c = 1.0;
     if (isOverrideCoeff && c != 0.0)
       c = overrideCoeff;
     // printf("getCombCoeff(%d)=%f diagRank=%d is2d=%d\n", g, c, diagRank(g), is2D());
     return (c); 
  }
}

int FTGridCombine3D::getGrank(int rank) {
  return (GridCombine3D::getGrank(rank));
}

// enumerate grids starting from x intercept, then towards y, then z
Vec3D<int> FTGridCombine3D::gridIx(int gid) {
  assert (0 <= gid  &&  gid <= nGrids());
  Vec3D<int> g;
  int l = level, d = diagRank(gid);
  int dIx = 0, dIO = gid;
  for (int i=0; i < d; ++i) // subtract grids from upper hyperplanes 
    dIO -= (is2D()? l: (l+1)*l/2) , l--; 
  while (dIO >= l) { // l is length of the dIx'th diagonal in the hyperplane
    assert (l > 0); assert (!is2D());
    dIO -= l; l--; ++dIx;
  }
  // grid gid is the dIO'th grid in the dIx'th diagonal in the d'th hyperplane
  g.x = gxS.x - d - dIx - dIO;
  g.y = gxS.y - level + 1 + dIO;
  g.z = is2D()? 0: gxS.z - level + 1 + dIx;
  if (is2D()) { // check against previously derived 2D formulas
     if (!haveExtraGrids) {
        assert (g.x == gxS.x - gid + d*(level - 1));
        assert (g.y == gxS.y + gid - (d+1)*level + 1); 
     }
     else {
        assert (g.x == gxS.x - d - dIO);
        assert (g.y == gxS.y - level + 1 + dIO);         
     }
  }
  return g; 
}


int FTGridCombine3D::localCombRank(int rank) {
   // Converting global non-failed and failed ranks into the ranks of combination
   // communicators
   for (int dp = 0; dp < degOfParall; ++dp) {
      // nprocs is the processes count for a single combination communicator
      if (rank >= (dp*nprocs) && rank <= ((dp+1)*nprocs - 1)) {
        return (rank - (dp*nprocs));
      }
   }
   printf("Error in FTGridCombine3D::localCombRank().\n");
   exit(1);
}//localCombRank()

FTGridCombine3D::FTGridCombine3D(int lv, int pd[], bool fixedP, int sgProcs, 
			     bool is2dim, Vec3D<int> gridS, Timer *t, int verb, bool dbg,
                             int degP, MPI_Comm comm, bool haveExtraGrids_) 
                             : GridCombine3D() {
  listCoeffs = NULL;
  haveExtraGrids = haveExtraGrids_;
  pgSComm = comm;
  degOfParall = degP;  
  level = lv; pD[0] = pd[0]; pD[1] = pd[1]; pD[2] = pd[2]; pD[3] = pd[3];
  fixedProcs = fixedP;
  is2d = is2dim, gxS = gridS; // needed here to set up is2D()
  assert (!is2d || gxS.z == 0);
  nprocs = nProcs(is2D(), level, pD, fixedProcs, haveExtraGrids);
  gRanks = new int[nprocs]; gIds = new int[nprocs];
  timer = t;

  verbosity = verb; debug = dbg;
  MPI_Comm_rank(pgSComm, &myrank);

  isOverrideCoeff = false;
  gidSingle = GID_ALL;

  pgs = new ProcGrid3D* [nGrids()];
  int pCount = 0; // count of global processes covered so far
  for (int i = 0; i < nGrids(); ++i) { // define process grid for each grid
    int pInc = pD[diagRank(i)]; // assumes fixedProcs
    for (int j = pCount; j < pCount + pInc; ++j) { // setup process maps
      gRanks[j] = j - pCount;
      gIds[j] = i;
    }
    pgs[i] = new ProcGrid3D((pCount<=myrank && myrank<pCount+pInc)? myrank: -1,
			    getProcConfig(pInc, is2d, gridIx(i)),
			    pgSComm, pCount);
    pCount += pInc;
  }
  pgU = pgs[getGid(myrank)]; // grid for this process       
  pgS = new ProcGrid3D(myrank < (1<<sgProcs)? myrank: -1, 
		       getProcConfig(1 << sgProcs, is2d, gridS), pgSComm);

  gidU = getGid(myrank);
  gxU = gridIx(gidU);

  // check that s.g. process (0,0,0) receives some points from this grid;
  // otherwise sgProcs is too large for simple interpolation scheme used.
  // The following is equivalent to gxU - level - 1 >= log2(gS->P)
  assert((pgS->G2L(gridSz1(gxU), Vec3D<int>(0))).prod() > 0); 

  selfBuf = 0;
#ifdef USE_BSEND
  sendBuf = 0; sendBufSz = 0;
#else
  msgSendRequest = 0; msgSendBuff = 0;
  nSendMsgs = 0; maxSendMsgs = 0;
#endif
} //FTGridCombine3D()

void FTGridCombine3D::gatherScatter(HaloArray3D* u, FTSparseGrid3D* usg) {
#ifdef TAU_PROF
  TAU_PROFILE("void FTGridCombine3D::gatherScatter(HaloArray3D*, FTSparseGrid3D*)", " ", TAU_USER);
#endif
  gather(u, usg);
  scatter(u, usg);
} //gatherScatter()

void FTGridCombine3D::gather(HaloArray3D* u, FTSparseGrid3D* usg) {
  sendInit(true, u->B);
  gatherSend(u);
  gatherRecv(usg);
  sendWait(true);
} //gather()

void FTGridCombine3D::scatter(HaloArray3D* u, FTSparseGrid3D* usg) {
  sendInit(false, usg->B);
  scatterSend(usg);
  scatterRecv(u);
  sendWait(false);
} //scatter()

void FTGridCombine3D::sendInit(bool gather, int B) {
   if (gather) {

#ifdef USE_BSEND
     sendBufSz = 0;
#else
     maxSendMsgs = 0;
#endif
     if (getCombCoeff(gidU) != 0.0 && pgU->myrank != -1) { // participating
       Vec3D<int> nMsgs = pgS->P / pgU->P + pgS->P.min(Vec3D<int> (2));
#ifdef USE_BSEND
       Vec3D<int> us = pgS->G2L(gridSz1(gxS), pgS->P - 1) + 1;
       sendBufSz = nMsgs.prod() * (MPI_BSEND_OVERHEAD +
				  B * us.prod() * sizeof(double)); 
#else
       maxSendMsgs = nMsgs.prod();
#endif   
     }
   } else {
#ifdef USE_BSEND
     sendBufSz = 0; sendBuf = 0; 
#else
     maxSendMsgs = 0;   nSendMsgs = 0;
     msgSendRequest = 0; msgSendBuff = 0; 
#endif
     if (pgS->myrank < 0) 
       return;
     for (int g=0; g < nGrids(); ++g) {
       Vec3D<int> nMsgs = pgs[g]->P/pgS->P + pgs[g]->P.min(Vec3D<int>(2));
#ifdef USE_BSEND
       Vec3D<int> us = pgs[g]->G2L(gridSz1(gridIx(g)), pgs[g]->P - 1);
       sendBufSz += nMsgs.prod() * (MPI_BSEND_OVERHEAD +
				    B * us.prod() * sizeof(double));
#else
       maxSendMsgs += nMsgs.prod();
#endif
     }
  }
#ifdef USE_BSEND
  if (verbosity > 1)
    printf("%d: BSend %s buffer size %ld\n", pgU->myrank, 
	   gather? "gather": "scatter", sendBufSz); 
  if (sendBufSz > 0) {
    sendBuf = new char [sendBufSz];
    assert (sendBuf != 0);
    MPI_Buffer_attach((void *) sendBuf, sendBufSz);
  }
#else
  if (maxSendMsgs > 0) {
    msgSendRequest = new MPI_Request [maxSendMsgs];
    msgSendBuff = new double* [maxSendMsgs];
  }
#endif
} //sendInit()


void FTGridCombine3D::sendWait(bool gather) {
  timer->start("sendWait", 0, 1);
  if (verbosity > 1)
    printf("%d: %s send wait\n", pgU->myrank, gather? "gather": "scatter");
#ifdef USE_BSEND
  MPI_Barrier(pgS->comm);   // ensure all recvs have completed
#endif
  if (gather || pgS->myrank >= 0) {
#ifdef USE_BSEND
    if (sendBuf != 0) { // this process has participated in the sendInit()
      int bSz; void *bPtr;
      MPI_Buffer_detach(&bPtr, &bSz);
      assert (bPtr == (void*) sendBuf);
      assert (bSz == sendBufSz);  
      delete[] sendBuf;
    }
    sendBuf = 0;
    sendBufSz = 0;
#else
    if (nSendMsgs > 0) {
      MPI_Status *sendReqStatus = new MPI_Status [nSendMsgs];
      MPI_Waitall(nSendMsgs, msgSendRequest, sendReqStatus);      
      delete[] sendReqStatus; 
    }
    for (int i = 0; i < nSendMsgs; ++i)
      delete[] msgSendBuff[i];
    delete[] msgSendRequest; delete[] msgSendBuff; 
    msgSendRequest = 0; msgSendBuff = 0;
    maxSendMsgs = nSendMsgs = 0; 
#endif
  }
  timer->stop("sendWait");
} //sendWait()

// current process sends its grid to the respective sparse grid process.
// There may be > 1 of these when more processes are used for
// the sparse grid than the current grid (the expected case). 
// nb. grid dimensions are always <= those of the sparse grid (=> / rSU>={1,1})

void FTGridCombine3D::gatherSend(HaloArray3D* u) {
  gatherSendScatterRecv(true, u);
}
void FTGridCombine3D::scatterRecv(HaloArray3D* u) {
  gatherSendScatterRecv(false, u);
}

#define CHECK_RECV_OK(s, len) \
       do { \
	 int count, rv  = MPI_Get_count(&(s), MPI_DOUBLE, &count); \
	 assert(rv != MPI_ERR_TRUNCATE); assert(count == (len)); \
       } while (0) 
   
void FTGridCombine3D::gatherSendScatterRecv(bool send, HaloArray3D* u) {
  // check that s.g. process (0,0,0) receives some points from this grid;
  // otherwise sgProcs is too large for simple interpolation scheme used.
  // The following is equivalent to gxU - level - 1 >= log2(gS->P)
  assert((pgS->G2L(gridSz(gxU), Vec3D<int>(0))).prod() > 0);

  if (send && getCombCoeff(gidU) == 0.0)
    return;
  Vec3D<int> rSU = gridSz1(gxS) / gridSz1(gxU); 
  assert (!is2D() || rSU.z == 1);
  Vec3D<int> NU = gridSz1(gxU), NS = gridSz1(gxS);
  // now find SG co-ordinates and offsets corresponding to our first point
  Vec3D<int> srcG0 = pgU->L2G0(NU), destG0 = srcG0 * rSU; 
  Vec3D<int> p0S = pgS->getP0(destG0, NS), pS;  
  Vec3D<int> offs0S = pgS->getOffs0(destG0, NS); 
  Vec3D<int> nU = pgU->G2L(NU);      // we will send away all of our points
  int i=0; 
  int B = u->B;                      // element block size
#ifndef USE_BSEND // define recv message buffers and associated data
  int maxMsgs = (pgS->P / pgU->P + 2).prod();
  MPI_Request *msgRequest = new MPI_Request[maxMsgs];
  double** msgBuff = new double* [maxMsgs];
  int nMsgs = 0;
  Vec3D<int> *msgI0 = new Vec3D<int> [maxMsgs], // recv msg ijk offset
             *msgDn = new Vec3D<int> [maxMsgs]; // recv msg size
#endif
  pS.x = p0S.x;
  while (i < nU.x) {
    int j=0, dnx = nU.x;    // force termination if j or k loop is empty   
    pS.y = p0S.y;
    while (j < nU.y) {
      int k=0, dny = nU.y;  // force termination if j or k loop is empty      
      pS.z = p0S.z;
      while (k < nU.z) { //send next msg to SG process pS, offset offsS in SG
	Vec3D<int> ijk = Vec3D<int>(i, j, k);
	Vec3D<int> offsS = Vec3D<int>(i==0? offs0S.x: 0, j==0? offs0S.y: 0, 
				      k==0? offs0S.z: 0);
	assert (pS <= pgS->P); //check if not stopped when should, bad getP0()
	//following case not yet implemented; will require sending prev. points
	// Should not occur if pgS process grid sizes are powers of 2
	assert (pgS->L2G0(NS, pS) % rSU == Vec3D<int>(0));

      	Vec3D<int> nS = pgS->G2L(NS, pS) - offsS; // num. points on dest
	Vec3D<int> dn = nS / rSU;         // corresp. number of points here
	assert (pgS->ownsData(pS, destG0 + rSU*ijk, rSU*dn, NS));//sanity check
	//now add last row/col (grid length=2^k+1), if truncated by / rSU above
	//o.w. the case nS % rSU != (0,0) is covered as we pass next row/col
	//dn = dn + (pgS->lastProc(pS) % rSU); // add if dest has them // commented for 2^k grids version
	dn = dn.min(nU - ijk);             // truncate to the end of our points
	int rankS = pgS->getRank(pS);
	if (verbosity > 1)
	  printf("%d: %s %dx%dx%d from pt (%d,%d,%d) @proc %d=(%d,%d,%d)\n",
		 pgU->myrank, send? "send": "recv", dn.x*B, dn.y, dn.z,
	         srcG0.x*B+i, srcG0.y+j, srcG0.y+k ,rankS, pS.x, pS.y, pS.z);
	assert (dn.prod() > 0);  // as we still have points to send

        // Process with non-zero grid coeff sends (gather-send) their value 
	if (send) {
	  // we need to include next row/col for interpolation where rSU>1
	  // note: interpolate must avoid accessing them for rSU=1 if we don't
	  //       in any dimesnion, u must have a halo accomodate inclusion
	  timer->start("pack", u->l.prod(), 2);
	  Vec3D<int> dnS = dn + Vec3D<int>(1, rSU.y > 1, rSU.z > 1); 
	  double *buf = u->pack(ijk, dnS);
	  timer->stop("pack");
	  if (rankS != pgU->myrank) {
#ifdef USE_BSEND
	    MPI_Bsend(buf, dnS.prod()*B, MPI_DOUBLE, rankS, GSTAG, pgS->comm);
	    delete[] buf;
#else
	    MPI_Isend(buf, dnS.prod()*B, MPI_DOUBLE, rankS, GSTAG, 
		      pgS->comm, &msgSendRequest[nSendMsgs]);
	    msgSendBuff[nSendMsgs] = buf; 
	    ++nSendMsgs; assert (nSendMsgs <= maxSendMsgs);
#endif
	  } else 
	    selfBuf = buf;
	} else { // recv
	  double *uR = new double [dn.prod()*B];
	  if (rankS != pgS->myrank) {
#ifdef USE_BSEND
	    MPI_Status s; 
	    MPI_Recv(uR, dn.prod()*B, MPI_DOUBLE, rankS, GSTAG, pgS->comm, &s);
	    CHECK_RECV_OK(s, dn.prod()*B);
#else
	    MPI_Irecv(uR, dn.prod()*B, MPI_DOUBLE, rankS, GSTAG, pgS->comm, 
		      &msgRequest[nMsgs]);
	    msgBuff[nMsgs] = uR; msgI0[nMsgs] = ijk; msgDn[nMsgs] = dn;
	    ++nMsgs; assert(nMsgs <= maxMsgs);	    
#endif
	  } else {
	    delete[] uR;
	    assert (selfBuf!=0); // set by scatterSend()
	    uR = selfBuf;
	    selfBuf = 0;
	  }
#ifndef USE_BSEND
	  if (rankS == pgS->myrank) {
#endif
	    timer->start("unpack", dn.prod(), 2);
	    u->unpack(uR, ijk, dn);
	    timer->stop("unpack");
	    delete[] uR;
#ifndef USE_BSEND
	  }
#endif	    
	} 
	k += dn.z; dnx = dn.x; dny = dn.y; ++pS.z;
      } //for (k...)
      j += dny; ++pS.y;
    } //for (j...)
    i += dnx; ++pS.x;
  } //for (i...)
  
#ifndef USE_BSEND 
  if (!send) {
    for (int i=0; i < nMsgs; ++i) {
      int msgIx; MPI_Status s; 
      MPI_Waitany(nMsgs, msgRequest, &msgIx, &s); 
      CHECK_RECV_OK(s, msgDn[msgIx].prod()*B);
      timer->start("unpack", msgDn[msgIx].prod(), 2);
      u->unpack(msgBuff[msgIx], msgI0[msgIx], msgDn[msgIx]);
      timer->stop("unpack");
    }
  }
  for (int i = 0; i < nMsgs; ++i)
    delete[] msgBuff[i];
  delete[] msgRequest; delete[] msgBuff; delete[] msgDn; delete[] msgI0;
#endif
} //gatherSendScatterRecv()

// for each component grid, gather contributions from respective processes. 
// There may be > 1 of these when fewer processes are used for
// the sparse grid than the component grid. 
void FTGridCombine3D::gatherRecv(FTSparseGrid3D* uS) {
  gatherRecvScatterSend(true, uS);
}
void FTGridCombine3D::scatterSend(FTSparseGrid3D* uS) {
  gatherRecvScatterSend(false, uS);
}

void FTGridCombine3D::gatherRecvScatterSend(bool recv, FTSparseGrid3D* uS) { 
  if (pgS->myrank < 0) // this process does not hold part of sparse grid
    return;
#ifndef USE_BSEND
  int maxMsgs = 0;
  for (int g=0; g < nGrids(); ++g) {
    if (recv && getCombCoeff(g) != 0.0)
      maxMsgs += (pgs[g]->P / pgS->P + 2).prod();
  }
  MPI_Request* msgRequest = new MPI_Request [maxMsgs];
  int nMsgs = 0;
  Vec3D<int> *msgI0 = new Vec3D<int> [maxMsgs], // recv msg s..g. ijk offset 
    *msgDn = new Vec3D<int> [maxMsgs], //  s.g. size corresp. recv msg  
    *msgRSU = new Vec3D<int> [maxMsgs]; // recv msg  interpolation ratio
  HaloArray3D** msgUR = new HaloArray3D* [maxMsgs]; // recv msg data
  double *msgCoeff = new double [maxMsgs]; 
#endif

  for (int g=0; g < nGrids(); ++g) {
    double coeff = getCombCoeff(g);
    if (recv && coeff == 0.0) // this grid has been excluded from the combination
      continue;
    Vec3D<int> gxU = gridIx(g);
    ProcGrid3D * pgU = pgs[g];
    Vec3D<int> rSU = gridSz1(gxS) / gridSz1(gxU); 
    assert (!is2D() || rSU.z == 1);
    Vec3D<int> NU = gridSz1(gxU), NS = gridSz1(gxS);
    // grid U process co-ordinates and offsets corresponding to our first point
    Vec3D<int> destG0 = pgS->L2G0(NS), srcG0 = destG0 / rSU; 
    Vec3D<int> p0U = pgU->getP0(srcG0, NU), pU;
    Vec3D<int> offs0U = pgU->getOffs0(srcG0, NU);  
    Vec3D<int> nS = pgS->G2L(NS);
    int B = uS->B;                      // element block size
    int i = 0;
    pU.x = p0U.x;
    while (i < nS.x) {
      int j=0, dnx = nS.x;
      pU.y = p0U.y;
      while (j < nS.y) {
	int k=0, dny = nS.y;
	pU.z = p0U.z;
	while (k < nS.z) {
	  Vec3D<int> ijk = Vec3D<int>(i, j, k);
	  Vec3D<int> offsU = Vec3D<int>(i==0? offs0U.x: 0, j==0? offs0U.y: 0, 
					k==0? offs0U.z: 0);	  
	  Vec3D<int> nU = pgU->G2L(NU, pU) - offsU; // no. of points on sender
	  Vec3D<int> dn = nU * rSU;  // corresponding no. of points here
	  dn = dn.min(nS - ijk);     // cut off at our boundary 
	  nU = nU.min(dn / rSU);     // do not expect extra points to be sent
	  Vec3D<int> lastPts = (pgU->lastProc(pU)).min(pgS->lastProc()) % rSU;
	  //nU = nU + lastPts;         // unless they're the last ones // must commented for 2^k grids version  
	  int rankU = pgU->getRank(pU);
	  if (verbosity > 1)
	    printf("%d: %s %dx%dx%d pts grid %d proc %d=%d,%d,%d corresp. pt (%d,%d,%d) %dx%dx%d\n", 
		   pgS->myrank, recv? "recv": "send", nU.x*B, nU.y, nU.z, g, 
		   rankU, pU.x, pU.y, pU.z, destG0.x*B+i, destG0.y+j, destG0.z+k,
		   dn.x*B, dn.y, dn.z);
	  assert (pU <= pgU->P); //check if not stopped when should / bad getP0
	  assert (pgU->ownsData(pU, srcG0 + ijk/rSU, nU, NU)); // sanity check
	  assert (dn.prod() > 0); // as we have more points to receive

	  if (getCombCoeff(getGid(rankU)) != 0 && recv) { // include next points for interpolation
	    nU = nU + Vec3D<int>(1, rSU.y > 1, rSU.z > 1); 
	    HaloArray3D *uR = new HaloArray3D(nU, Vec3D<int>(0), B);
	    if (rankU != pgS->myrank) {
#ifdef USE_BSEND
	      MPI_Status s; 
	      MPI_Recv(uR->u, nU.prod()*B, MPI_DOUBLE, rankU, GSTAG, 
		       pgS->comm, &s);
	      CHECK_RECV_OK(s, nU.prod()*B);
#else
	      MPI_Irecv(uR->u, nU.prod()*B, MPI_DOUBLE, rankU, GSTAG,
			pgS->comm, &msgRequest[nMsgs]);
	      msgI0[nMsgs] = ijk; msgDn[nMsgs] = dn; msgUR[nMsgs] = uR;
	      msgRSU[nMsgs] = rSU; msgCoeff[nMsgs] = coeff;
	      ++nMsgs; assert(nMsgs <= maxMsgs);
#endif
	    } else {
	      delete[] uR->u;
	      assert (selfBuf!=0); // set by gatherSend()
	      uR->u = selfBuf;
	      selfBuf = 0;
	    }
#ifndef USE_BSEND
	    if (rankU == pgS->myrank) {
#endif
	      if (verbosity >= 4)
		uR->print(pgS->myrank, "uR");
	      timer->start("interpolate", dn.prod(), 2);
	      uS->interpolate(coeff, uR, rSU, ijk, dn);
	      timer->stop("interpolate");
	      delete uR;
#ifndef USE_BSEND
	    }
#endif	      
	  } else { // send
	    timer->start("sample", uS->numElts(), 2);
	    double *buf = uS->sample(ijk, rSU, nU);
	    timer->stop("sample");
	    if (rankU != pgS->myrank) {
#ifdef USE_BSEND
	      MPI_Bsend(buf, nU.prod()*B, MPI_DOUBLE, rankU, GSTAG, pgS->comm);
	      delete[] buf;
#else
	      MPI_Isend(buf, nU.prod()*B, MPI_DOUBLE, rankU, GSTAG, pgS->comm,
			&msgSendRequest[nSendMsgs]);
	      msgSendBuff[nSendMsgs] = buf;
	      ++nSendMsgs; assert (nSendMsgs <= maxSendMsgs);
#endif
	    } else 
	      selfBuf = buf;
	  }
	  k += dn.z; dny = dn.y; dnx = dn.x; ++pU.z;
	} //while (k...)
        j += dny; ++pU.y;
      } //while (j...)
      i += dnx; ++pU.x;
    } //while (i...)
  } //while (g...)

#ifndef USE_BSEND
  if (recv) {
    for (int i=0; i < nMsgs; ++i) {
      int msgIx; MPI_Status s;
      MPI_Waitany(nMsgs, msgRequest, &msgIx, &s);
      CHECK_RECV_OK(s, msgUR[msgIx]->l.prod());
      timer->start("interpolate", msgDn[msgIx].prod(), 2);
      uS->interpolate(msgCoeff[msgIx], msgUR[msgIx], msgRSU[msgIx], 
	              msgI0[msgIx], msgDn[msgIx]);
      timer->stop("interpolate");
      delete msgUR[msgIx];
    }
  }
  delete[] msgRequest; delete[] msgDn; delete[] msgI0;
  delete[] msgUR; delete[] msgRSU; delete[] msgCoeff; 
#endif
} //gatherRecvScatterSend()

void FTGridCombine3D::setGridId(int x, int y, int z, int gId) {
   int ny, nz;
   ny = (level+1);
   nz = is2D()? 2: (level+1);
   gxGtoGridId[z + nz * (y + ny * x)] = gId;
} //setGridId()

int FTGridCombine3D::getGridId(int x, int y, int z) {
   int ny, nz;
   ny = (level+1);
   nz = is2D()? 2: (level+1);
   return (gxGtoGridId[z + nz * (y + ny * x)]);
} //getGridId()

void FTGridCombine3D::setGridCoordinatesGridIds() {
   // Allocate memory for listCoeffs
   listCoeffs = new double [nGrids()];
      
   // Determining 3D coordinates (indices started from 1) of grid of grids
   int nx, ny, nz;
   nx = ny = (level+1);
   nz = is2D()? 2: (level+1);
   gxGtoGridId = new int [nx*ny*nz]; // 1 more (in each dimension) since indices started from 1

   for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
         for (int k = 0; k < nz; ++k) {
            setGridId(i, j, k, -1); // init gxGtoGridId with -1
         }
      }
   }

   gxG = new Vec3D<int> [nGrids()];
   int gridId = 0, zSize = 0;
   for (int planeId = 0; planeId <= diagRank(nGrids()-1); ++planeId) { // nGrids()-1 is the max gridId
      zSize = is2D()? 1: (level - planeId);
      for (int z = 1; z <= zSize; ++z) {
         for (int y = 1; y <= level - (z-1) - planeId; ++y) {
            gxG[gridId] = Vec3D<int>(level - planeId - (z-1) - (y-1), y, z);
            setGridId(level - planeId - (z-1) - (y-1), y, z, gridId); // set gxGtoGridId with gridId
            ++gridId;
         }
      }
   }
#ifdef PRINT_gxG_getGridId
   // Testing coordinates of grid of grids, and grid id of coordinates of grid of grids
   if (myrank == 0) {
      for (int t = 0; t < nGrids(); ++t) {
         printf("[%d] ===== gxG[%d] = (%d,%d,%d), getGridId(%d,%d,%d) = %d\n", myrank, t, 
                gxG[t].x, gxG[t].y, gxG[t].z, gxG[t].x, gxG[t].y, gxG[t].z, 
                getGridId(gxG[t].x, gxG[t].y, gxG[t].z));
      }
   }
#endif   
} //setGridCoordinatesGridIds()

FTGridCombine3D::~FTGridCombine3D() {
   // nGrids() of GridCombine3D (base) and FTGridCombine3D (derived) classes are different.
   // In derived class, it may be larger.
   // Base destructor destructs only pgs[0] to pgs[GridCombine3D::nGrids()-1].
   // Derived destructor destructs pgs[GridCombine3D::nGrids()] to pgs[nGrids()-1]  
   for (int i=GridCombine3D::nGrids(); i < nGrids(); ++i) {
      //delete[] pgs[i];
      delete pgs[i];
   }
   delete[] gxG;
   delete[] gxGtoGridId;
   if (listCoeffs != NULL)
      delete[] listCoeffs;
} //~FTGridCombine3D()


void FTGridCombine3D::alternateCombCoeffs(int* failedList, int numFailed) {
   double sTime = 0.0;
   int globalRank, globalSize;
   MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
   MPI_Comm_rank(MPI_COMM_WORLD, &globalSize);
   if(globalRank == 0){     
      sTime = MPI_Wtime();
   }  
   if (is2D())
      coeffSelectUpdate2D(failedList, numFailed);
   else
      coeffSelectUpdate3D(THREE_D, failedList, numFailed);     
   if (globalRank == 0 && verbosity > 0) {
      printf("\n");  
      for(int i = 0; i < nGrids(); ++i){                                             
	 printf("[%d] ===== Coefficient of grid %d = %0.1f =====\n", 
                globalRank, i, listCoeffs[i]);
      }
      printf("\n"); 
   }
   if (globalRank == 0)
      printf("[%d]----- Creating coefficient list takes %0.6f Sec (MPI_Wtime) -----\n", 
             globalRank, MPI_Wtime() - sTime);    
} //alternateCombCoeffs()

void FTGridCombine3D:: coeffSelectUpdate2D(int * failedList, int numFailed) {
   double * coeffList = new double[nGrids()];
   double * coeffListAlt = new double[nGrids()];
   int * isGridFailed = new int[nGrids()];
   int * maxGridList = new int [level]; 
   int fid, gid, levId, levId2, globalRank, 
       failedGridCounter = 0, 
       dimSizeY, dimSizeX, 
       totMax = 0, maxFound;   
   
   int ** maxGridList2D = new int *[level+2];
   for(levId = 0; levId < (level+2); ++levId)
      maxGridList2D[levId] = new int [level+2-levId];

   MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);

   // Initialization
   #pragma omp parallel for default(shared)
   for(gid = 0; gid < nGrids(); ++gid){
      coeffList[gid] = 0.0;
      coeffListAlt[gid] = 0.0;
      isGridFailed[gid] = 0;
   }

   #pragma omp parallel for default(shared)
   for(levId = 0 ; levId < level; ++levId)
      maxGridList[levId] = NOT_SET; 

   #pragma omp parallel for default(shared)
   for(levId = 0; levId < (level+2); ++levId){
      for(levId2 = 0; levId2 < (level+2-levId); ++levId2)       
         maxGridList2D[levId][levId2] = NOT_SET;
   }    

   // Determine which grids are failed
   #pragma omp parallel for default(shared)
   for(fid = 0; fid < numFailed; ++fid)
      isGridFailed[getGid(failedList[fid])] = 1;

   // Determine number of failed grids
   for(gid = 0; gid < nGrids(); ++gid){
      if(isGridFailed[gid] == 1) {
         if(globalRank == 0)
            printf("====================== grid %d is failed\n", gid);
         ++failedGridCounter; // I have to use it for testing the condition of how many grid is possible to fail
      }
   }
   
   // Check if all grids are failed. exit if so
   if(failedGridCounter == MAX_GRID_FAILED){
       if(myrank == 0)
          printf("\n[%d] ***** All grids are failed. Exiting the application. *****\n\n", myrank);
       exit(1);
   }

   // Determine a list which contains max
   for(gid = 0; gid < nGrids(); ++gid){
      // For each grid, check all the possible grids on upper set for max.
      // If there is no such grid found, flag that non-failed grid as max.
      maxFound = 0;//reset 
      dimSizeY = gxG[gid].y;         
      dimSizeX = gxG[gid].x;               
      for(levId2 = dimSizeX; levId2 < (level + 2 - dimSizeY); ++levId2){
          for(levId = dimSizeY; levId < (level + 2 - levId2); ++levId){
              if(maxGridList2D[levId2][levId] != NOT_SET){
                 maxFound = 1;
                 break;
              }
          }
          if(maxFound == 1)
             break;          
      }
      if(!isGridFailed[gid] && !maxFound){          
         maxGridList[dimSizeY - 1] = gid;    
         maxGridList2D[dimSizeX][dimSizeY] = gid;       
         ++totMax; 
      }
   } 

   // Printing max grid list
   if(myrank == 0 && verbosity > 0){
      for(levId = 0 ; levId < level; ++levId)
         printf("[%d] ***** maxGridList[%d] = %d *****\n", myrank, levId, maxGridList[levId]);
   }
   
   // Calculating coefficients
   // Scanning the max values from bottom-to-top and top-to-bottom order
   // Generate two sets of coefficient values
   // Choose the one which contains both the first and last grid of upper layer (if possible)
   int onlyContribGrid = 0, minusCoeffCounter = 0;
   if(totMax > 0){
      int nextLevId, failedTest, intersectId, xSize, ySize;  
      // Scanning from bottom-to-top
      for(levId = 0; levId < level; ++levId){
         while(maxGridList[levId] == NOT_SET && levId < level)//this does not hold max value, try for next grid
            ++levId;
         if(levId < level)
            onlyContribGrid = maxGridList[levId];

         for(nextLevId = (levId+1); nextLevId < level; ++nextLevId){
            while(maxGridList[nextLevId] == NOT_SET && nextLevId < level){//this does not hold max value, try for next grid
               ++nextLevId;
            }           
      
            if(nextLevId < level){ //at least two max values found
               xSize = gxG[maxGridList[nextLevId]].x;
               ySize = gxG[maxGridList[levId]].y;

               intersectId = getGridId(xSize, ySize, 1); // z index 1 indicates 2D
               failedTest = (intersectId == NOT_SET)? 1: isGridFailed[intersectId];
 
               if(failedTest == 1){ //intersection is failed or invalid
                  continue;
               }
               else if(failedTest == 0){ //intersection is not failed or not invalid
                  coeffList[maxGridList[levId]] = +1.0;
                  coeffList[maxGridList[nextLevId]] = +1.0;
                  coeffList[intersectId] = -1.0;
	
                  ++minusCoeffCounter;	          

                  levId = nextLevId-1; // replace levId with nextLevId in outer loop (after increment)
                  nextLevId = level-1; //break; //breaking inner loop (condition will fail after the increment in loop)
               }
            }//end of if(nextLevId < level...
         }//end of for(nextLevId = ...
      }//end of for(levId = 0...

      // Scanning from top-to-bottom     
      for(levId = level-1; levId >= 0; levId--){
         while(maxGridList[levId] == NOT_SET && levId >= 0)//this does not hold max value, try for next grid
            levId--;
         if(levId >= 0)
            onlyContribGrid = maxGridList[levId];
         for(nextLevId = (levId-1); nextLevId >=  0; nextLevId--){
            while(maxGridList[nextLevId] == NOT_SET && nextLevId >= 0){//this does not hold max value, try for next grid
               nextLevId--;
            }           

            if(nextLevId >= 0){ //at least two max values found
               xSize = gxG[maxGridList[levId]].x;
               ySize = gxG[maxGridList[nextLevId]].y;

               intersectId = getGridId(xSize, ySize, 1); // z index 1 indicates 2D
               failedTest = (intersectId == NOT_SET)? 1: isGridFailed[intersectId];

               if(failedTest == 1) //intersection is failed or invalid
                  continue;
               else if(failedTest == 0){ //intersection is not failed or not invalid
                  coeffListAlt[maxGridList[levId]] = +1.0;
                  coeffListAlt[maxGridList[nextLevId]] = +1.0;
                  coeffListAlt[intersectId] = -1.0;
	
                  ++minusCoeffCounter;	          

                  levId = nextLevId+1; // replace levId with nextLevId in outer loop (after decrement)
                  nextLevId = 0; //break; //breaking inner loop (condition will fail after the decrement in loop)
               }
            }//end of if(nextLevId < level...
         }//end of for(nextLevId = ...
      }//end of for(levId = level-1...
   }//end of if(totMax > 0 ...        

   // coeffListAlt is the first priority
   if(totMax == 1 || minusCoeffCounter == 0)
      coeffListAlt[onlyContribGrid] = +1.0;

   // Choose the one which contains both the first and last grid of upper layer (if possible)
   // Copy the values
   if(coeffListAlt[0] == 1 && coeffListAlt[level-1] == 1){
      #pragma omp parallel for default(shared)
      for(gid = 0; gid < nGrids(); ++gid)
         listCoeffs[gid] = coeffListAlt[gid];
   }
   else if(coeffList[0] == 1 && coeffList[level-1] == 1){
      #pragma omp parallel for default(shared)
      for(gid = 0; gid < nGrids(); ++gid)
         listCoeffs[gid] = coeffList[gid];
   }
   else{
      #pragma omp parallel for default(shared)
      for(gid = 0; gid < nGrids(); ++gid)
         listCoeffs[gid] = coeffListAlt[gid];
   }

   // Memory release 
   for(levId = 0; levId < level; ++levId)
      delete[]  maxGridList2D[levId];
   delete[] maxGridList2D;
   
   delete[] coeffList;
   delete[] coeffListAlt;
   delete[] isGridFailed;
   delete[] maxGridList;
} // coeffSelectUpdate2D()

void FTGridCombine3D::coeffSelectUpdate3D(int dimension, int * failedList, int numFailed) {
   int * isGridFailed = new int [nGrids()];    // an array with value 1 at index i (i>=0) means
                                               // grid with id i fails
   int globalRank,                             // global MPI rank
       gid,                                    // a counter
       fid,                                    // a counter
       numGridFailureExtraLayer = 0,           // grid failure count on extra layer
       failedGridCounter = 0,                  // count number of failed grids           
       * failedGridList = NULL;                // store an array of grids that fails
   
   MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);  

   // Initialization
   #pragma omp parallel for default(shared)
   for(gid = 0; gid < nGrids(); ++gid){
      // Default coefficients
      if (diagRank(gid) == 0)
         listCoeffs[gid] = +1.0;      
      else if (diagRank(gid) == 1)
         listCoeffs[gid] = -2.0;      
      else if (diagRank(gid) == 2)
         listCoeffs[gid] = +1.0;      
      else
         listCoeffs[gid] = 0.0;      
      // Grid is not failed
      isGridFailed[gid] = 0;     
   }

   if (numFailed > 0) {
      // Determine which grids fail
      #pragma omp parallel for default(shared)
      for(fid = 0; fid < numFailed; ++fid){
         isGridFailed[getGid(failedList[fid])] = 1;
      }

      // Determine number of failed grids
      for (gid = 0; gid < nGrids(); ++gid) {
         if (isGridFailed[gid] == 1) {
            if(globalRank == 0)
               printf("===== Grid %d fails =====\n", gid);
            ++failedGridCounter;
         }
      }

      // Application will exit if MAX_GRID_FAILED fails
      if (failedGridCounter >= MAX_GRID_FAILED) {
         if (globalRank == 0) 
            printf("\n[%d] ***** failedGridCounter fails its constraint. Exiting "
                   "the application. *****\n\n", globalRank);
         exit(1);
      }    

      // Determine a list of grids that fails
      failedGridList = new int[failedGridCounter];
      failedGridCounter = 0; //reuse failedGridCounter variable
      for (gid = 0; gid < nGrids(); ++gid) {
         if (isGridFailed[gid] == 1) {
            failedGridList[failedGridCounter++] = gid;
            if (diagRank(gid) == FAILURE_ON_EXTRA_LAYER)  // failure on forth layer (extra layer)
               ++numGridFailureExtraLayer;
         }   
      }   

      // Updating coefficients
      if (numGridFailureExtraLayer == failedGridCounter) { // all failures on forth layer (extra layer)
         if (globalRank == 0) 
            printf("\n[%d] ***** All failures are on the grids on extra plane. Return default "
                   "coefficients. *****\n", globalRank);
      }      
      else // failure on first, and/or second, and/or third layer(s)
         // Selecting alternate combination coefficients for multiple grid failures
         // by definition, level on coeffSelectUpdateAlgorithm(...) is one less than the definition on other functions
         coeffSelectUpdateAlgorithm3D(gxG, dimension, level-1, failedGridList, failedGridCounter, listCoeffs, 0);          
   } // end of if (numFailed > 0)
   
   // Memory release 
   delete[] isGridFailed;
   if (numFailed > 0)
      delete[] failedGridList;
} //coeffSelectUpdate3D()


