/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// Parallel 3D Sparse Grid Combination class
//   handles a parallel 3D combination algorithm
//   For 4+D grids, parallelization of dimensions outside of the 3D being
//   used for the combination is possible; if p is the total ||ization factor,
//   there will be p groups of processes performing independent combination
//   algorithms, requiring p separate SGCT communicators (comm) passed to this
//   class.
// written by Peter Strazdins, Jun 14

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <mpi.h>
#include <assert.h>
#include "GridCombine3D.h"
#include <cmath>        /* log2() */

#define GSTAG 5 /*message tag for gather & scatter comms*/
#define GSWTAG 6 /*message tag for gather & scatter wait acks*/

bool GridCombine3D::is2D() {
  return is2d;
}

// return size corresponding to co-ordinate gix                               
int gridSz(int gix) {
  return (1 + (1 << gix));
}

Vec3D<int> GridCombine3D::gridSz(Vec3D<int> gix) {
  Vec3D<int> g;
  g.x = 1 + (1 << gix.x);
  g.y = 1 + (1 << gix.y);
  g.z = (is2d == 0) + (1 << gix.z); 
  // n.b. gix.z == 0 doesn't neccessarily indicate 2D case
  return (g);
}

// return which grid id global rank is in
int GridCombine3D::getGid(int rank) {
   assert (0 <= rank && rank < nprocs); 
   int v = gIds[rank];
   return (v);
}

int GridCombine3D::nGrids() {
  if (is2D()) 
    return (2*level - 1);        
  else
    return (level*level + (level-1)*(level-2)/2);    
}

int GridCombine3D::diagRank(int gid) {
  assert (0 <= gid  && gid < nGrids());
  if (is2D())
    return (gid >= level);
  else if (gid < (level+1)*level/2)
    return 0;
  else if (gid < level*level)
    return 1;
  else
    return 2;
}


double GridCombine3D::getCombCoeff(int g) {
  double c = (diagRank(g)!=1)? +1.0: (is2D()? -1.0: -2.0); // default
  if (g != gidSingle  &&  gidSingle != GID_ALL)
    c = 0.0;
  if (isOverrideCoeff && c != 0.0)
    c = overrideCoeff;
  return (c); 
}

void GridCombine3D::overrideCombCoeff(double coeff) {
  isOverrideCoeff = true;
  overrideCoeff = coeff;
}
void GridCombine3D::restoreCombCoeff() {
  isOverrideCoeff = false;
}
void GridCombine3D::singleGridCombine(int g) {
  gidSingle = g;
}

int GridCombine3D::getGrank(int rank) {
  assert (0 <= rank && rank < nprocs); 
  int v = gRanks[rank];
  return (v);
}

// enumerate grids starting from x intercept, then towards y, then z
Vec3D<int> GridCombine3D::gridIx(int gid) {
  assert (0 <= gid  &&  gid <= nGrids());
  Vec3D<int> g;
  int l = level, d = diagRank(gid);
  int dIx = 0, dIO = gid;
  for (int i=0; i < d; i++) // subtract grids from upper hyperplanes 
    dIO -= (is2D()? l: (l+1)*l/2) , l--; 
  while (dIO >= l) { // l is length of the dIx'th diagonal in the hyperplane
    assert (l > 0); assert (!is2D());
    dIO -= l; l--; dIx++;
  }
  Vec3D<int> gS = hMode? gxSsave: gxS; // ensure we get size of original SG
  // grid gid is the dIO'th grid in the dIx'th diagonal in the d'th hyperplane
  g.x = gS.x - d - dIx - dIO;
  g.y = gS.y - level + 1 + dIO;
  g.z = is2D()? 0: gS.z - level + 1 + dIx;
  if (is2D()) { // check against previously derived 2D formulas
    assert (g.x == gS.x - gid + d*(level - 1));
    assert (g.y == gS.y + gid - (d+1)*level + 1); 
  }
  return g; 
}

GridCombine3D::GridCombine3D(int lv, int pd[], bool fixedP, int sgProcs, 
			     bool is2dim, Vec3D<int> gridS, Timer *t, 
			     MPI_Comm comm, int verb) {
  level = lv; pD[0] = pd[0]; pD[1] = pd[1]; pD[2] = pd[2]; 
  fixedProcs = fixedP; pgSComm = comm;
  is2d = is2dim ,gxS = gridS; // needed here to set up is2D()  
  assert (!is2d || gxS.z == 0);
  hMode = false; // need here for gridIx() to work!
  nprocs = nProcs(is2D(), level, pD, fixedProcs);
  timer = t;

  int np; MPI_Comm_size(pgSComm, &np);
  assert (np == nprocs); // o.w. pgSComm's size not consistent with level, pD 
  verbosity = verb; 
  MPI_Comm_rank(pgSComm, &myrank);

  isOverrideCoeff = false;
  gidSingle = GID_ALL;

  pgs = new ProcGrid3D* [nGrids()];
  gRanks = new int[nprocs]; gIds = new int[nprocs];
  gRank0 = new int[nGrids()];
  int pCount = 0; // count of global processes covered so far
  for (int i = 0; i < nGrids(); i++) { // define process grid for each grid
    int pInc = pD[diagRank(i)]; // assumes fixedProcs
    assert (pCount + pInc <= nprocs);
    gRank0[i] = pCount;
    for (int j = pCount; j < pCount + pInc; j++) { // setup process maps
      gRanks[j] = j - pCount;
      gIds[j] = i;
    }
    pgs[i] = new ProcGrid3D((pCount<=myrank && myrank<pCount+pInc)? myrank: -1,
			    getProcConfig(pInc, is2d,  gridIx(i)),
			    pgSComm, pCount);
    pCount += pInc;
  }

  pgU = pgs[getGid(myrank)]; // grid for this process       
  pgS = new ProcGrid3D(myrank < (1<<sgProcs)? myrank: -1, 
		       getProcConfig(1 << sgProcs, is2d, gridS), pgSComm);

  gidU = getGid(myrank);
  gxU = gridIx(gidU);

  selfBuf = 0;
#ifdef USE_BSEND
  sendBuf = 0; sendBufSz = 0; 
  sendRanks = 0;
#else
  msgSendRequest = 0; msgSendBuff = 0;
#endif
  nSendMsgs = 0; maxSendMsgs = 0;
} //GridCombine3D()

GridCombine3D::~GridCombine3D() {
  for (int i=0; i < nGrids(); i++)
      delete pgs[i];
  delete[] pgs; delete pgS;
  delete[] gRanks; delete[] gIds;
  //delete[] gRank0;
}

void GridCombine3D::gatherScatter(HaloArray3D* u, SparseGrid3D* usg) {
  gather(u, usg);
  scatter(u, usg);
}

void GridCombine3D::gather(HaloArray3D* u, SparseGrid3D* usg) {
  assert (u->B >= 1  &&  u->B == usg->B);
  timer->start("gather", 0, 0);
  sendInit(true, u->B);
  gatherSend(u); 
  gatherRecv(usg);
  sendWait(true);
  timer->stop("gather");
} //gatherScatter()

void GridCombine3D::scatter(HaloArray3D* u, SparseGrid3D* usg) {
  assert (u->B >= 1  &&  u->B == usg->B);
  timer->start("scatter", 0, 0);
  sendInit(false, usg->B);
  scatterSend(usg);
  scatterRecv(u);
  sendWait(false);
  timer->stop("scatter");
} //scatterSingle()

void GridCombine3D::sendInit(bool gather, int B) {
  maxSendMsgs = 0; nSendMsgs = 0; 
#ifdef USE_BSEND
  sendBufSz = 0; sendBuf = 0; 
  sendRanks = 0;
#else
  msgSendRequest = 0; msgSendBuff = 0; 
#endif
  if (gather) {
    if (getCombCoeff(gidU) != 0.0  &&  pgU->myrank != -1) { // participating
      Vec3D<int> nMsgs = pgS->P / pgU->P + pgS->P.min(Vec3D<int> (2));
      maxSendMsgs = nMsgs.prod();
#ifdef USE_BSEND
      Vec3D<int> us = pgS->G2L(gridSz(gxS), pgS->P - 1) + 1;
      sendBufSz = nMsgs.prod() * (MPI_BSEND_OVERHEAD +
				  B * us.prod() * sizeof(double)); 
#endif   
    }
  } else {
     if (pgS->myrank < 0) 
       return;
     for (int g=0; g < nGrids(); g++) {
       if (hMode && !(hx <= gridIx(g))) // g too small for hier. basis comb 
	 continue;
       Vec3D<int> nMsgs = pgs[g]->P/pgS->P + pgs[g]->P.min(Vec3D<int>(2));
       maxSendMsgs += nMsgs.prod();
#ifdef USE_BSEND
       // for hMode, this will over-estimate as any combination of 
       // hierarchical surpluses on the current grid is no larger than the grid
       Vec3D<int> us = pgs[g]->G2L(gridSz(gridIx(g)), pgs[g]->P - 1);
       sendBufSz += nMsgs.prod() * (MPI_BSEND_OVERHEAD +
				    B * us.prod() * sizeof(double));
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
  sendRanks = new int [maxSendMsgs]; 
#else
  msgSendRequest = new MPI_Request [maxSendMsgs];
  msgSendBuff = new double* [maxSendMsgs];
#endif
} //sendInit()

void GridCombine3D::sendWait(bool gather) {
  timer->start("sendWait", 0, 1);
  if (verbosity > 1)
    printf("%d: %s send wait\n", pgU->myrank, gather? "gather": "scatter");
#ifdef USE_BSEND
  MPI_Barrier(pgS->comm);   // ensure all recvs have completed
#endif
  if (gather || pgS->myrank >= 0) {
#ifdef USE_BSEND
    MPI_Status s; double b;
    for (int i=0; i < nSendMsgs; i++)
      MPI_Recv(&b, 0, MPI_DOUBLE, sendRanks[i], GSWTAG, pgS->comm, &s);
    delete[] sendRanks; sendRanks = 0;
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
    for (int i = 0; i < nSendMsgs; i++)
      delete[] msgSendBuff[i];
    delete[] msgSendRequest; delete[] msgSendBuff; 
    msgSendRequest = 0; msgSendBuff = 0;
#endif
    maxSendMsgs = nSendMsgs = 0; 
  }
  timer->stop("sendWait");
}

// current process sends its grid to the respective sparse grid process.
// There may be > 1 of these when more processes are used for
// the sparse grid than the current grid (the expected case). 
// nb. grid dimensions are always <= those of the sparse grid (=> / rSU>={1,1})

void GridCombine3D::gatherSend(HaloArray3D* u) {
  timer->start("gatherSend", u->l.prod(), 1);
  gatherSendScatterRecv(true, u);
  timer->stop("gatherSend");
}
void GridCombine3D::scatterRecv(HaloArray3D* u) {
  timer->start("scatterRecv", u->l.prod(), 1);
  gatherSendScatterRecv(false, u);
  timer->stop("scatterRecv");
}

#define CHECK_RECV_OK(s, len) \
       do { \
	 int count, rv  = MPI_Get_count(&(s), MPI_DOUBLE, &count); \
	 assert(rv != MPI_ERR_TRUNCATE); assert(count == (len)); \
       } while (0) 
   
void GridCombine3D::gatherSendScatterRecv(bool send, HaloArray3D* u) {
  if (pgU->myrank == -1) {//not part of subgrid from current (small) surplus
    assert (hMode);
    return;
  }
  // check that s.g. process (0,0,0) receives some points from this grid;
  // otherwise sgProcs is too large for simple interpolation scheme used.
  // The following is equivalent to gxU - level - 1 >= log2(gS->P)
  assert((pgS->G2L(gridSz(gxU), Vec3D<int>(0))).prod() > 0); 

  if (send && getCombCoeff(gidU) == 0.0)
    return;
  assert (!hMode || gxU <= gxUsave); 
  Vec3D<int> rSU = gridSz1(gxS) / gridSz1(gxU); 
  assert (!is2D() || rSU.z == 1);
  assert (!hMode  ||  rSU == Vec3D<int>(1));
  Vec3D<int> NU = gridSz(gxU), NS = gridSz(gxS);
  // now find SG co-ordinates and offsets corresponding to our first point
  Vec3D<int> srcG0 = pgU->L2G0(NU), destG0 = srcG0 * rSU; 
  Vec3D<int> p0S = pgS->getP0(destG0, NS), pS;  
  Vec3D<int> offs0S = pgS->getOffs0(destG0, NS); 
  Vec3D<int> nU = pgU->G2L(NU);      // we will send away all of our points
  int i=0; 
  int B = u->B;                      // element block size
#ifndef USE_BSEND /* define recv message buffers and associated data*/
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
#if 0
	if (!pgS->ownsData(pS, destG0 + rSU*ijk, rSU*dn, NS))
	  printf("%d: pgU->id=" V3DFMT " pgU->P=" V3DFMT
		 " pS->P=" V3DFMT " NU= " V3DFMT " Ns= " V3DFMT 
		 " nU=" V3DFMT " ijk=" V3DFMT " offs0S=" V3DFMT  "\n",  
		 pgU->myrank,  V3DLST(pgU->id), V3DLST(pgU->P), 
		 V3DLST(pgS->P),V3DLST(NU),  V3DLST(nS), 
		 V3DLST(nU),  V3DLST(ijk), V3DLST(offs0S)); 
#endif
	assert (pgS->ownsData(pS, destG0 + rSU*ijk, rSU*dn, NS));//sanity check
	//now add last row/col (grid length=2^k+1), if truncated by / rSU above
	//o.w. the case nS % rSU != (0,0) is covered as we pass next row/col
	dn = dn + (pgS->lastProc(pS) % rSU); // add if dest has them
	dn = dn.min(nU - ijk);             // truncate to the end of our points
	int rankS = pgS->getRank(pS);
	if (verbosity > 2)
	  printf("%d: %s %dx%dx%d from pt (%d,%d,%d) @proc %d=(%d,%d,%d)\n",
		 pgU->myrank, send? "send": "recv", dn.x*B, dn.y, dn.z,
	         srcG0.x*B+i, srcG0.y+j, srcG0.y+k ,rankS, pS.x, pS.y, pS.z);
	assert (dn.prod() > 0);  // as we still have points to send

        // Process with non-zero grid coeff sends (gather-send) their value 
	if (send) {
	  // we need to include next row/col for interpolation where rSU>1
	  // note: interpolate must avoid accessing them for rSU=1 if we don't.
	  //       in any dimension, u must have a halo accommodate inclusion
	  timer->start("pack", u->l.prod(), 2);
	  Vec3D<int> dnS = dn;
	  if (!hMode) // we do not interpolate in hMode
	    dnS = dnS + Vec3D<int>(1, rSU.y > 1, rSU.z > 1); 
	  double *buf = u->pack(ijk, dnS);
	  timer->stop("pack");
	  if (rankS != pgU->myrank) {
#ifdef USE_BSEND
	    MPI_Bsend(buf, dnS.prod()*B, MPI_DOUBLE, rankS, GSTAG, pgS->comm);
	    delete[] buf;
	    sendRanks[nSendMsgs] = rankS;
#else
	    MPI_Isend(buf, dnS.prod()*B, MPI_DOUBLE, rankS, GSTAG, 
		      pgS->comm, &msgSendRequest[nSendMsgs]);
	    msgSendBuff[nSendMsgs] = buf; 
#endif
	    nSendMsgs++; assert (nSendMsgs <= maxSendMsgs);
	  } else 
	    selfBuf = buf;
	} else { // recv
	  double *uR = new double [dn.prod()*B];
	  if (rankS != pgS->myrank) {
#ifdef USE_BSEND
	    MPI_Status s; double b;
	    MPI_Recv(uR, dn.prod()*B, MPI_DOUBLE, rankS, GSTAG, pgS->comm, &s);
	    CHECK_RECV_OK(s, dn.prod()*B);
	    MPI_Bsend(&b, 0, MPI_DOUBLE, rankS, GSWTAG, pgS->comm);
#else
	    MPI_Irecv(uR, dn.prod()*B, MPI_DOUBLE, rankS, GSTAG, pgS->comm, 
		      &msgRequest[nMsgs]);
	    msgBuff[nMsgs] = uR; msgI0[nMsgs] = ijk; msgDn[nMsgs] = dn;
	    nMsgs++; assert(nMsgs <= maxMsgs);	    
#endif
	  } else {
	    delete[] uR;
	    assert (selfBuf!=0); // set by scatterSend()
	    uR = selfBuf;
	    selfBuf = 0;
	  }
#ifndef USE_BSEND
	  if (rankS == pgS->myrank) { // receive and unpack later
#endif
	    timer->start("unpack", dn.prod(), 2);
	    u->unpack(uR, ijk, dn);
	    timer->stop("unpack");
	    delete[] uR;
#ifndef USE_BSEND
	  }
#endif	    
	} 
	k += dn.z; dnx = dn.x; dny = dn.y; pS.z++;
      } //for (k...)
      j += dny; pS.y++;
    } //for (j...)
    i += dnx; pS.x++;
  } //for (i...)
  
#ifndef USE_BSEND 
  if (!send) {
    for (int i=0; i < nMsgs; i++) {
      int msgIx; MPI_Status s; 
      MPI_Waitany(nMsgs, msgRequest, &msgIx, &s); 
      CHECK_RECV_OK(s, msgDn[msgIx].prod()*B);
      timer->start("unpack", msgDn[msgIx].prod(), 2);
      u->unpack(msgBuff[msgIx], msgI0[msgIx], msgDn[msgIx]);
      timer->stop("unpack");
    }
  }
  for (int i = 0; i < nMsgs; i++)
    delete[] msgBuff[i];
  delete[] msgRequest; delete[] msgBuff; delete[] msgDn; delete[] msgI0;
#endif

} //gatherSendScatterRecv()

// for each component grid, gather contributions from respective processes. 
// There may be > 1 of these when fewer processes are used for
// the sparse grid than the component grid. 
void GridCombine3D::gatherRecv(SparseGrid3D* uS) {
  timer->start("gatherRecv", uS->numElts(), 1);
  if (gidSingle == GID_ALL) {
    timer->start("zero", uS->numElts(), 2);
    uS->zero();
    timer->stop("zero");
  }
  gatherRecvScatterSend(true, uS);
  timer->stop("gatherRecv");
}
void GridCombine3D::scatterSend(SparseGrid3D* uS) {
  timer->start("scatterSend", uS->numElts(), 1);
  gatherRecvScatterSend(false, uS);
  timer->stop("scatterSend");
}

void GridCombine3D::gatherRecvScatterSend(bool recv, SparseGrid3D* uS) { 
  assert (!hMode || gxU <= gxUsave); 
  if (pgS->myrank < 0) // this process does not hold part of sparse grid
    return;
#ifndef USE_BSEND
  int maxMsgs = 0;
  for (int g=0; g < nGrids(); g++) {
    if (hMode && !(hx <= gridIx(g))) // grid too small for hier. basis comb 
      continue;
    if (recv && getCombCoeff(g) != 0.0)
      maxMsgs += (pgs[g]->P / pgS->P + 2).prod();
  }
  MPI_Request* msgRequest = new MPI_Request [maxMsgs]; // recv msg handle
  int nMsgs = 0;
  Vec3D<int> *msgI0 = new Vec3D<int> [maxMsgs], // recv msg s..g. ijk offset 
    *msgDn = new Vec3D<int> [maxMsgs], //  s.g. size corresp. recv msg  
    *msgRSU = new Vec3D<int> [maxMsgs]; // recv msg  interpolation ratio
  HaloArray3D** msgUR = new HaloArray3D* [maxMsgs]; // recv msg data
  double *msgCoeff = new double [maxMsgs]; 
#endif

  assert (!hMode || uS->useFullGrid);
  for (int g=0; g < nGrids(); g++) {
    double coeff = getCombCoeff(g);
    if (recv && coeff == 0.0) //this grid is excluded from the combination
      continue;
    if (hMode && !(hx <= gridIx(g))) // grid too small for hier. basis comb.
      continue;
    Vec3D<int> gxU = hMode? gxS: gridIx(g);
    ProcGrid3D *pgU = pgs[g];
    Vec3D<int> rSU = gridSz1(gxS) / gridSz1(gxU); 
    assert (!is2D() || rSU.z == 1);
    assert (!hMode  ||  rSU == Vec3D<int>(1));
    Vec3D<int> NU = gridSz(gxU), NS = gridSz(gxS);
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
	  nU = nU + lastPts;         // unless they're the last ones  
	  int rankU = pgU->getRank(pU);
	  if (verbosity > 2)
	    printf("%d: %s %dx%dx%d pts grid %d proc %d=%d,%d,%d corresp. pt (%d,%d,%d) %dx%dx%d\n", 
		   pgS->myrank, recv? "recv": "send", nU.x*B, nU.y, nU.z, g, 
		   rankU, pU.x, pU.y, pU.z, destG0.x*B+i, destG0.y+j, destG0.z+k,
		   dn.x*B, dn.y, dn.z);
	  assert (pU <= pgU->P); //check if not stopped when should / bad getP0
	  assert (pgU->ownsData(pU, srcG0 + ijk/rSU, nU, NU)); // sanity check
	  assert (dn.prod() > 0); // as we have more points to receive
	  assert (g == getGid(rankU)); // if passes, replace with g below

	  if (getCombCoeff(getGid(rankU)) != 0 && recv) { 
	    if (!hMode) // include next points for interpolation
	      nU = nU + Vec3D<int>(1, rSU.y > 1, rSU.z > 1); 
	    HaloArray3D *uR = new HaloArray3D(nU, Vec3D<int>(0), B);
	    if (rankU != pgS->myrank) {
#ifdef USE_BSEND
	      MPI_Status s; double b;
	      MPI_Recv(uR->u, nU.prod()*B, MPI_DOUBLE, rankU, GSTAG, 
		       pgS->comm, &s);
	      CHECK_RECV_OK(s, nU.prod()*B);
	      MPI_Bsend(&b, 0, MPI_DOUBLE, rankU, GSWTAG, pgS->comm);
#else
	      MPI_Irecv(uR->u, nU.prod()*B, MPI_DOUBLE, rankU, GSTAG,
			pgS->comm, &msgRequest[nMsgs]);
	      msgI0[nMsgs] = ijk; msgDn[nMsgs] = dn; msgUR[nMsgs] = uR;
	      msgRSU[nMsgs] = rSU; msgCoeff[nMsgs] = coeff;
	      nMsgs++; assert(nMsgs <= maxMsgs);
#endif
	    } else {
	      delete[] uR->u;
	      assert (selfBuf!=0); // set by gatherSend()
	      uR->u = selfBuf;
	      selfBuf = 0;
	    }
#ifndef USE_BSEND /* unless uR comes from here we recv & interp. it later */ 
	    if (rankU == pgS->myrank) {
#endif
	      if (verbosity > 4)
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
	      sendRanks[nSendMsgs] = rankU;
#else
	      MPI_Isend(buf, nU.prod()*B, MPI_DOUBLE, rankU, GSTAG, pgS->comm,
			&msgSendRequest[nSendMsgs]);
	      msgSendBuff[nSendMsgs] = buf;
#endif
	      nSendMsgs++; assert (nSendMsgs <= maxSendMsgs);
	    } else 
	      selfBuf = buf;
	  }
	  k += dn.z; dny = dn.y; dnx = dn.x; pU.z++;
	} //while (k...)
        j += dny; pU.y++;
      } //while (j...)
      i += dnx; pU.x++;
    } //while (i...)
  } //while (g...)

#ifndef USE_BSEND
  if (recv) {
    if (verbosity > 4) {
      uS->print(pgS->myrank, "uS");
    }
    for (int i=0; i < nMsgs; i++) {
      int msgIx; MPI_Status s;
      MPI_Waitany(nMsgs, msgRequest, &msgIx, &s);
      CHECK_RECV_OK(s, msgUR[msgIx]->l.prod());
      timer->start("interpolate", msgDn[msgIx].prod(), 2);
#if 0
      if (verbosity > 4) {
	printf("%d: %d of %d\n", pgS->myrank, i, nMsgs);
	printf("%d: rS="V3DFMT" i0="V3DFMT" dn="V3DFMT"\n", pgS->myrank,
	       V3DLST(msgRSU[msgIx]), V3DLST(msgI0[msgIx]), V3DLST(msgDn[msgIx]));
      }
#endif
      uS->interpolate(msgCoeff[msgIx], msgUR[msgIx], msgRSU[msgIx], 
			msgI0[msgIx], msgDn[msgIx]);
#if 0
      if (verbosity > 4) {
	uS->print(pgS->myrank, "uS");
      }
#endif
      timer->stop("interpolate");
      delete msgUR[msgIx];
    }
  }
  delete[] msgRequest; delete[] msgDn; delete[] msgI0;
  delete[] msgUR; delete[] msgRSU; delete[] msgCoeff; 
#endif
} //gatherRecvScatterSend()

HaloArray3D* GridCombine3D::setHierComb(Vec3D<int> hx_, Vec3D<int> hgx) {
  gxUsave = gxU; gxSsave = gxS; pgSsave = pgS;
  hx = hx_; gxU = gxS = hgx; 
  hMode = true;
  // create new pgS for procs from grid g such that gridIx(g) <= hg
  int Ptot = 0;
  for (int g=0; g < nGrids(); g++) {
    if (hx <= gridIx(g)) {
      Ptot += pgs[g]->P.prod();
    }
  }
  Ptot = 1 << ((int) log2(Ptot));
  Ptot = std::min(Ptot, gridSz1(hgx).prod());
  int *r0 = new int [nGrids()]; int *Ps = new int [nGrids()]; 
  int nP = 0; int Ptotg = 0; 
  int hrank = -1;
  for (int g=0; g < nGrids() && Ptotg < Ptot; g++) {
    //    printf("g=%d grid " V3DFMT " gU=%d gxS=" V3DFMT "\n", 
    //	   g, V3DLST(gridIx(g)), gidU, V3DLST(gxS));
    if (hx <= gridIx(g)) {
      if (g == gidU)
	hrank = myrank;
      r0[nP] = gRank0[g]; Ps[nP] = pgs[g]->P.prod(); 
      Ptotg += Ps[nP];
      nP++;
    }
  }
  Vec3D<int> PS = getProcConfig(Ptot, is2d, hgx);
#if 0
  if (!(PS <= gridSz1(hgx)))
    printf("%d: hgx=" V3DFMT " PS=" V3DFMT "\n", pgSsave->myrank,
	   V3DLST(hgx), V3DLST(PS)); 
#endif
  assert (PS <= gridSz1(hgx)); //otherwise block distribution breaks down
  pgS = new ProcGrid3D(hrank, PS, pgSComm, nP, r0, Ps);
  delete[] r0; delete[] Ps; 
  return (new HaloArray3D(pgS->G2L(gridSz(hgx)), Vec3D<int>(0)));
} //setHierComb() 

void GridCombine3D::resetHierComb() {
  delete pgS;
  gxU = gxUsave; gxS = gxSsave; pgS = pgSsave;
  hMode = false;
} //resetHierComb 
