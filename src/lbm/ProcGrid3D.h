/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// 3D Process Grid Class for handling 3D data grids distributed in a 
// block fashion.  
// Written by Peter Strazdins, Jun 14

#ifndef PROCGRID3D_INCLUDED
#define PROCGRID3D_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

class ProcGrid3D {
public:
  int myrank;    // this process's MPI rank in comm (-1 if not part of this grid)
  int *rank0;    // rank0[0] = MPI rank of logical process (0,0,0) in this grid
  Vec3D<int> P;  // processor grid size
  Vec3D<int> id; // 3D rank in processor grid
  MPI_Comm comm;
  // extension for hierarchical basis combination, where contiguous sub-ranks
  // in comm, representing a selection of component grids, are used for 
  // a series of applications of the SGCT
  int nP;        // rank0[0:nP-1] are ranks in comm of contiguous sub-ranks 
  int nSz;       // of component process grids of sizes Psz[0:nSz-1],   
  int *Psz, *nPsz; // where there are nPsz[i] grids of size Psz[i].
  // we will have nP = sum (i: 0 <= i < nSz) nPsz[i];

  /* extension to support process sub-grids for the component grids for
     coalesced hierarchical surpluses whose size in any dimensions is smaller
     than P, thus violating the block distribution. Here we
     temporarily set rank0, id, P to the process sub-grid which will
     hold 1 element in any `too small' dimension, still in a block
     distribution. Processes not in this sub-grid have myrank==-1.
  */
  Vec3D<int> offsP, strP; //process offset and stride for sub-grid
  Vec3D<int> zeroDel; //a `delta' on the mapping for the last process for
                      //coalesced supluses (containing surplus 0)
  int myrankOrig/*, rank0Orig*/; // stores original values
  Vec3D<int> POrig, idOrig;

  ProcGrid3D(int rank, Vec3D<int> P_, MPI_Comm c, int r0=0) {
    P = P_;
    myrank = rank;
    comm = c;
    nP = nSz = 1; 
    rank0 = new int[nP]; Psz = new int[nSz]; nPsz = new int[nSz]; 
    rank0[0] = r0;
    Psz[0] = P.prod();  nPsz[0] = 1;
    if (myrank == -1) 
      id.x = id.y = id.z = -1;
    else { 
      id.x = (myrank - rank0[0]) % P.x;
      id.y = ((myrank - rank0[0]) / P.x) % P.y;
      id.z = ((myrank - rank0[0]) / P.x) / P.y;
    }
    strP = Vec3D<int> (1); offsP = zeroDel = Vec3D<int> (0);
    myrankOrig = myrank; POrig = P, idOrig = id; //rank0Orig = rank0[0];
  } //ProcGrid3D()
  
  // pre: r0, Ps have nP elements
  ProcGrid3D(int rank, Vec3D<int> P_, MPI_Comm c, int nP_, int r0[], int Ps[]){
    P = P_;
    myrank = rank;
    comm = c;
    nP = nP_;
    rank0 = new int[nP]; Psz = new int[nP]; nPsz = new int[nP];
    rank0[0] = r0[0]; Psz[0] = Ps[0]; nPsz[0] = 1;
    nSz = 1; 

    strP = Vec3D<int> (1); offsP = zeroDel = Vec3D<int> (0);
    myrankOrig = myrank; POrig = P, idOrig = id; //rank0Orig = rank0[0];

    for (int i=1; i < nP; i++) {
      rank0[i] = r0[i];
      assert (rank0[i] - rank0[i-1] >= Ps[i-1]); //check r0[] correctly formed
      if (Ps[i] == Ps[i-1]) 
	nPsz[nSz-1]++;
      else 
	Psz[nSz] = Ps[i], nPsz[nSz] = 1, nSz++;
      assert (nSz <= nP);      
    }
    int Ptot = 0, ntot = 0;
    for (int i=0; i < nSz; i++) {
      Ptot += Psz[i] * nPsz[i], ntot += nPsz[i];
    }
    assert (Ptot >= P.prod()); // check P,Ps[] are consistent
#if 0
    if (ntot != nP)
      printf("%d: ProcGrid3D(%d, " V3DFMT ") ntot=%d nP=%d\n",
	     myrankOrig, myrank, V3DLST(P), ntot, nP);
#endif
    assert (ntot == nP); // sanity check

    id = Vec3D<int>(-1);
    if (myrank == -1)
      return;
    int myrankIx = -1, pid0 = 0;
    for (int i=0; i < nP  &&  myrankIx == -1; i++) 
      if (rank0[i] <= myrank && myrank < rank0[i] + Ps[i]) 
	myrankIx = i, pid0 += myrank - rank0[i]; 
      else
	pid0 += Ps[i];
    assert (0 <= myrankIx && myrankIx < nP); // o.w. invalid r0[], Ps[]

    if (pid0 >= P.prod()) // already have enough processes
      myrank = -1;
    else { 
      id.x = pid0 % P.x;
      id.y = (pid0 / P.x) % P.y;
      id.z = (pid0 / P.x) / P.y;
    }
    idOrig = id; 
  } //ProcGrid3D()

  ~ProcGrid3D() {
    delete rank0; delete Psz; delete nPsz;
  }

  void setSubgrid(Vec3D<int> sgPoffs, Vec3D<int> sgPstr, Vec3D<int>sgP,
		  Vec3D<int> zDel) {
    assert (nP == 1); // only implemented for grids made up a single rank set
#if 0
    if (!(sgPoffs + sgPstr * (sgP - 1) - zDel <= P - 1)) {
      int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      printf("%d: %d sgP:offs=" V3DFMT " str=" V3DFMT " sgP= " V3DFMT " zD= "
	     V3DFMT " P=" V3DFMT "\n", rank,
	     myrankOrig, V3DLST(sgPoffs), V3DLST(sgPstr),
	     V3DLST(sgP), V3DLST(zDel), V3DLST(P));
    }
#endif
    assert (sgPoffs + sgPstr * (sgP - 1) - zDel <= P - 1); //a valid sub-grid
    assert (Vec3D<int>(0) <= zDel  &&  zDel <= Vec3D<int>(1));

    offsP = sgPoffs; strP = sgPstr; P = sgP; zeroDel = zDel;
    int pid0 = myrank - rank0[0];
    id = Vec3D<int> (-1);
    if (pid0 >= 0) {
      for (int d=0; d < 3; d++) {
	//deltas are only on coalesced (0..) surpluses, which have zero offsets
	assert (zeroDel.v[d] == 0  ||  offsP.v[d] == 0); 
	int idv = (pid0 - offsP.v[d]) % POrig.v[d];
	int dv = (idv == POrig.v[d]-1)? zeroDel.v[d]: 0;
 	if ((idv + dv) % strP.v[d] == 0) // this process is in the sub-grid
	  id.v[d] = (idv + dv) / strP.v[d];
	assert (dv == 0 || id.v[d] == P.v[d]-1);
	pid0 = pid0 / POrig.v[d];      
      }
    }  
    if (id.x == -1 || id.y == -1 || id.z == -1)
      myrank = -1; 
  } //setSubgrid()
  
  void resetSubgrid() {
    strP = Vec3D<int> (1); offsP = zeroDel = Vec3D<int> (0);
    myrank = myrankOrig; id = idOrig;  P = POrig; // rank0[0] = rank0Orig;
  }  

  // get rank of process (pidx, pidy) in this grid
  int getRank(Vec3D<int> pid) {
    assert (Vec3D<int>(0) <= pid  &&  pid <= P); 
    for (int d=0; d < 3; d++) // scale to place in original grid
      pid.v[d] = offsP.v[d] + strP.v[d] * pid.v[d] 
	- ((pid.v[d] == P.v[d]-1)? zeroDel.v[d]: 0);
    assert (Vec3D<int>(0) <= pid  &&  pid <= POrig); 
    int pid0 = pid.x + POrig.x*(pid.y + POrig.y*pid.z);
    int myrankIx = 0, i=0;
    do {
      int d = std::min(pid0 / Psz[i], nPsz[i]);
      myrankIx += d; pid0 -= d*Psz[i];
      if (d == nPsz[i])
	i++;
      else 
	break;
    } while (1);
    assert (i < nSz);
    assert (myrankIx < nP); assert (pid0 < Psz[i]);
    return (rank0[myrankIx] + pid0);
  } //getRank()

  // return global starting index in dim. d for global vector of length N
  // for this process 
  int L2G0(int d, int N) {
    return L2G0(d, N, id.v[d]);
  }
  // ... and for process index in pid
  int L2G0(int d, int N, int pid) {
    assert (0 <= d && d <= 2);  
    if (!(0 <= pid && pid < P.v[d])) {
      printf("%d: L2G0(%d,%d,%d) P=%d error\n", myrank, d, N, pid, P.v[d]);
      fflush(stdout); exit(1);
    }
    assert (0 <= pid && pid < P.v[d]);
    return ((N / P.v[d]) * pid); 
  }
  Vec3D<int> L2G0(Vec3D<int> N, Vec3D<int> pid) {
    return Vec3D<int>(L2G0(0, N.x, pid.x), L2G0(1, N.y, pid.y),
		      L2G0(2, N.z, pid.z));
  }
  Vec3D<int> L2G0(Vec3D<int> N) {
    return L2G0(N, id);    
  } 
  
  // return local length corresponding to N in dim d in for this process 
  int G2L(int d, int N) {
    if (myrank == -1) // this process is not part of grid
      return 0;
    return G2L(d, N, id.v[d]);
  }
  // and for process index pid
  int G2L(int d, int N, int pid) {
    assert (0 <= d && d <= 2);  
    if (!(0 <= pid && pid < P.v[d])) {
      printf("%d: G2L(%d,%d,%d) P=%d ERROR!\n", myrank, d, N, pid, P.v[d]);
      fflush(stdout);
    }
    assert (0 <= pid && pid < P.v[d]);
    int n = N / P.v[d];
    if (pid == P.v[d]-1)
      n += N % P.v[d];
    return n;
  }
  Vec3D<int> G2L(Vec3D<int> N, Vec3D<int> pid) {
    return Vec3D<int>(G2L(0, N.x, pid.x), G2L(1, N.y, pid.y), 
		      G2L(2, N.z, pid.z));
  }
  Vec3D<int> G2L(Vec3D<int> N) {
    if (myrank == -1) // this process is not part of grid
      return Vec3D<int>(0, 0, 0);
    return G2L(N, id);    
  } 

  // returns whether process pid is last in dimension d (has extra elements)
  Vec3D<int> lastProc() {
    return lastProc(id);
  }
  int lastProc(int d, int pid) {
    return (P.v[d]-1 == pid);
  }
  Vec3D<int> lastProc(Vec3D<int> pid) {
    return Vec3D<int>(lastProc(0, pid.x), lastProc(1, pid.y), 
		      lastProc(2, pid.z));
  }

  // return process index for element N0 in dimension d corresp. to N
  int getP0(int d, int N0, int N) {
    assert (0 <= d && d <= 2);
    assert (0 <= N0 && N0 <= N);
    assert (P.v[d] <= N);
    // it is possible, e.g. with P=7 N=33 N0=28, that the last process 
    // has >= 2x the points of others
    return (std::min(N0 / (N / P.v[d]), P.v[d]-1));
  }
  Vec3D<int> getP0(Vec3D<int> N0, Vec3D<int> N) {
    return Vec3D<int>(getP0(0, N0.x, N.x), getP0(1, N0.y, N.y),
		      getP0(2, N0.z, N.z));
  } 
 
  // return offset in process for element N0 in dimension d corresp. to N
  int getOffs0(int d, int N0, int N) {
    assert (0 <= d && d <= 2);
    assert (0 <= N0 && N0 <= N);
    int n = N / P.v[d];
    int o = (N0/n < P.v[d]-1)? N0 % n: N0 - n*(P.v[d]-1);
    if (!(o == N0 % n || (N0/n >= P.v[d]-1 &&  o >= N0 % n)))
      printf("%d: getOffs0(%d,%d,%d): P=%d o=%d, o'=%d\n", myrank, d, N0, N, 
	     P.v[d], o, N0 % n);
    assert (o == N0 % n || (N0/n >= P.v[d]-1 && o >= N0 % n));
    return o;
    // return (N0 % n + ((N0==N-1 && N>1)? n: 0));
  }
  Vec3D<int> getOffs0(Vec3D<int> N0, Vec3D<int> N) {
    return Vec3D<int>(getOffs0(0, N0.x, N.x), getOffs0(1, N0.y, N.y),
		      getOffs0(2, N0.z, N.z));
  } 

  // returns if process pid in grid owns points N0..N0+dN-1. for grid length N 
  bool ownsData(int d, int pid, int N0, int  dN, int N) {
    return (L2G0(d, N, pid) <= N0 && N0+dN <= L2G0(d, N, pid)+G2L(d, N, pid)); 
  }
  bool ownsData(Vec3D<int> pid, Vec3D<int> N0, Vec3D<int> dN, Vec3D<int> N) {
    bool rv = ownsData(0, pid.x, N0.x, dN.x, N.x) &&
      ownsData(1, pid.y, N0.y, dN.y, N.y) &&
      ownsData(2, pid.z, N0.z, dN.z, N.z);
    if (rv == 0)
      printf("%d: (%d) !ownsData(%d %d %d, %d %d %d, %d %d %d, %d,%d,%d)\n", 
	     myrankOrig, myrank, pid.x,pid.y,pid.z, N0.x,N0.y,N0.z, dN.x,dN.y,dN.z,  
	     N.x,N.y,N.z);
    return (rv);
  }
  
  // return rank of process in direction dir and dimension d 
  // left: dir=-1,d=0: right: dir=+1,d=0; up: dir=-1,d=1; down: dir=+1,d=1 
  int neighbour(int dir, int d) {
    assert (0 <= d && d <= 2);  
    assert (dir==-1 || dir==1);  
    Vec3D<int> nid = id;
    nid.v[d] = (P.v[d] + id.v[d] + dir) % P.v[d];
    return (getRank(nid));
  }
};

#endif /*PROCGRID3D_INCLUDED*/
