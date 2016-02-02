/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// Parallel 3D Sparse Grid Combination via Hierarchical Basis class
//   forms each hierarchy, uses GridCombine3D to do the combination on the
//   hierarchy and the restores the combined grid.
//   Credits to Brendan Harding MSI ANU for providing a serial implementation

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>
#include "HierachDecomp3D.h"

// #define NO_COALESC /* combine on surpluses without any coalescing */

// #define NO_SCHED_COMP /* no load balancing on coalesced hierarchies */

// maximize coalescing to increase range of solvable problems 
// (+ reduce startup) at expense of some parallelization
// #define MAX_COALESC 

// 3d hierarchical update/restore on a single a plane (within node) 
void HierachDecomp3D::hierUpdate(bool formHier, HaloArray3D *u, 
				 Vec3D<int> nU, int d, int j0, int dl){
  assert (0 <= d && d < 3); assert (nU.v[d] == 1);
  double sc = formHier? -0.5: +0.5;
  Vec3D<int> ic = Vec3D<int>(0, 0, 0); ic.v[d] = j0;
  Vec3D<int> il = ic, ir = ic; il.v[d] -= dl, ir.v[d] += dl; 
  assert (0 <= il.v[d]);
  for (int k = 0; k < nU.z; k++) 
    for (int j = 0; j < nU.y; j++) 
      for (int i = 0; i < nU.x; i++) {
	Vh(u, i+ic.x, j+ic.y, k+ic.z) += sc * (Vh(u, i+il.x, j+il.y, k+il.z) +
					       Vh(u, i+ir.x, j+ir.y, k+ir.z));
      } 
} //hierUpdate()

// version of hierUpdate() with one exterior plane (left: dl<0, right: dl>0)
void HierachDecomp3D::hierUpdateB(bool formHier, HaloArray3D *u, 
				  Vec3D<int> nU, int d, int j0, int dl, double *b) {
  assert (0 <= d && d < 3); assert (nU.v[d] == 1);
  double sc = formHier? -0.5: +0.5;
  Vec3D<int> ic = Vec3D<int>(0, 0, 0); ic.v[d] = j0;
  Vec3D<int> il = ic; il.v[d] += dl; 
  assert (0 <= il.v[d]);
  for (int k = 0; k < nU.z; k++) 
    for (int j = 0; j < nU.y; j++) 
      for (int i = 0; i < nU.x; i++) {
	Vh(u, i+ic.x, j+ic.y, k+ic.z) += sc * (Vh(u, i+il.x, j+il.y, k+il.z) +
						*b);
	b++;
      } 
} //hierUpdateB()

// version of hierUpdate() with 2 exterior planes 
void HierachDecomp3D::hierUpdateBB(bool formHier, HaloArray3D *u, 
					  Vec3D<int> nU, 
		  int d, int j0, double *bl, double *br) {
  assert (0 <= d && d < 3); assert (nU.v[d] == 1);
  double sc = formHier? -0.5: +0.5;
  Vec3D<int> ic = Vec3D<int>(0, 0, 0); ic.v[d] = j0;
  for (int k = 0; k < nU.z; k++) 
    for (int j = 0; j < nU.y; j++) 
      for (int i = 0; i < nU.x; i++) {
	Vh(u, i+ic.x, j+ic.y, k+ic.z) += sc * (*bl + *br);
	bl++; br++;
      } 
} //hierUpdateBB()

// extract local part of hierarchy from u[hOffs:hStride:hN] to uH 
void HierachDecomp3D::packHier(Vec3D<int> hOffs, Vec3D<int> hStr, 
			       Vec3D<int> hN, HaloArray3D *u, HaloArray3D *uH){
  assert (u->B == 1); //TODO
  for (int k = 0; k < hN.z; k++) 
    for (int j = 0; j < hN.y; j++) 
      for (int i = 0; i < hN.x; i++) 
	Vh(uH, i, j, k) = 
	  Vh(u, i*hStr.x + hOffs.x, j*hStr.y + hOffs.y, k*hStr.z + hOffs.z);
}

// update uH into local part of hierarchy in u[hOffs:hStride:hN] 
void HierachDecomp3D::unpackHier(Vec3D<int> hOffs, Vec3D<int> hStr, 
				 Vec3D<int> hN, HaloArray3D *u, HaloArray3D *uH) {
  assert (u->B == 1); //TODO
  for (int k = 0; k < hN.z; k++) 
    for (int j = 0; j < hN.y; j++) 
      for (int i = 0; i < hN.x; i++) 
	Vh(u, i*hStr.x + hOffs.x, j*hStr.y + hOffs.y, k*hStr.z + hOffs.z) =
	  Vh(uH, i, j, k);
}


#define HTAG 7
#define CHECK_RECV_OK(s, len) \
  do { \
    int count, rv  = MPI_Get_count(&(s), MPI_DOUBLE, &count);	\
    assert(rv != MPI_ERR_TRUNCATE); assert(count == (len));	\
  } while (0)

void HierachDecomp3D::hierarchize3D(bool formHier, HaloArray3D *u) {
  Vec3D<int> NU = gc->gridSz(gxU), nU = pgU->G2L(NU);
  Vec3D<int> NL = pgU->L2G0(NU), NR = NL + nU - 1;
  assert (u->B == 1); //TODO
  for (int d=0; d < 3; d++) {
    Vec3D<int> nUd = nU; nUd.v[d] = 1; 
    int nSz = nUd.prod();
    double *urL = new double [nSz], *urR =new double [nSz]; 
    // in the 2D case, gxU.v[2]==0 and this loop drops out
    for (int l = formHier? 0: gxU.v[d]-1; 
	 formHier? (l < gxU.v[d]): (l >= 0); l += formHier? 1: -1) {
      MPI_Status s; 
      int nMsgs= 0; MPI_Request msgRequest[2]; double *msgBuf[2];
      int dl = (1 << l), dl2 = 2*dl;
      //calc. local indices of left & right-most source points in hierarchy
      int jLm = NL.v[d] % dl2, jRm = NR.v[d] % dl2; 
      int jsL = (jLm > 0)? dl2 - jLm: 0, 
        jsR = nU.v[d]-1 - jRm;   
      //printf("%d: l=%d jsL=%d jsR=%d NL=%d NR=%d n=%d\n", pgU->myrank, l, 
      //     jsL, jsR, NL.v[d], NR.v[d], nU.v[d]);  
      assert ((NL.v[d] + jsL) % dl2 == 0 && (dl2 + NL.v[d] + jsR) % dl2 == 0);

      // perform local updates on dest points
      for (int j = jsL+dl; j <= jsR-dl; j+= dl2)
	hierUpdate(formHier, u, nUd, d, j, dl);  

      // calculate local indices of left & right-most dest points in hierarchy
      int jdL = jsL + ((jsL >= dl)? -dl: +dl), 
	jdR = jsR + ((jsR + dl < nU.v[d])? +dl: -dl);
      //printf("%d: l=%d jdL=%d jdR=%d NL=%d NR=%d n=%d\n", pgU->myrank, l, 
      //     jdL, jdR, NL.v[d], NR.v[d], nU.v[d]);  
      assert ((NL.v[d] + jdL) % dl2 == dl && (dl2 + NL.v[d] + jdR) % dl2 == dl);

      if (jsL < nU.v[d]  &&  jsL < jdL &&  NL.v[d] + jsL >= dl) { 
	//we own jsL and there is an external point to the left that needs it
	Vec3D<int> isL = Vec3D<int>(0); isL.v[d] = jsL;
	double* uS = u->pack(isL, nUd); 
	Vec3D<int> NdL = NL; NdL.v[d] -= dl - jsL;
	assert(NdL.v[d] >= 0); assert(NdL.v[d] % dl2 == dl);
	Vec3D<int> dst = pgU->getP0(NdL, NU); 
	assert (dst.v[d] <  pgU->id.v[d]);
	if (verbosity > 2)
	  printf("%d: %s(%d): d=%d l=%d send L g%d=l%d "V3DFMTSZ" to %d="V3DFMT"\n",
		 pgU->myrank, "hierarchize3D", formHier, d, l, NL.v[d], jsL, 
		 V3DLST(nUd), pgU->getRank(dst), V3DLST(dst));
	MPI_Isend(uS, nSz, MPI_DOUBLE, pgU->getRank(dst), HTAG, pgU->comm,
		  &msgRequest[nMsgs]);
	msgBuf[nMsgs] = uS; nMsgs++;
      }

      if (0 <= jsR  &&  jdR < jsR  &&  NL.v[d] + jsR + dl < NU.v[d]) { 
	//we own jsR and there is an external point to the right that needs it
	Vec3D<int> isR = Vec3D<int> (0); isR.v[d] = jsR;
	double* uS = u->pack(isR, nUd);
	Vec3D<int> NdR = NR;  NdR.v[d] += dl - (NR.v[d] - (NL.v[d] + jsR)); 
	assert(NdR.v[d] < NU.v[d]);  assert(NdR.v[d] % dl2 == dl);
	Vec3D<int> dst = pgU->getP0(NdR, NU);  
	assert (pgU->id.v[d] < dst.v[d]);
	if (verbosity > 2)
	  printf("%d: %s(%d) d=%d l=%d sendR g%d=l%d "V3DFMTSZ" to %d="V3DFMT"\n",
		 pgU->myrank, "hierarchize3D", formHier, d, l, NR.v[d], jsR, 
		 V3DLST(nUd), pgU->getRank(dst),V3DLST(dst));
	MPI_Isend(uS, nSz, MPI_DOUBLE, pgU->getRank(dst), HTAG, pgU->comm,
		  &msgRequest[nMsgs]);
	msgBuf[nMsgs] = uS; nMsgs++;
      }

      if (jdL < nU.v[d] && jdL < jsL) { //jdL is here & needs points from left 
	Vec3D<int> NsL = NL; NsL.v[d] -= dl - jdL;  
	assert (NsL.v[d] >= 0); assert (NsL.v[d] % dl2 == 0);	
	Vec3D<int> src = pgU->getP0(NsL, NU); assert (src.v[d] < pgU->id.v[d]);
	if (verbosity > 2)
	  printf("%d: %s(%d) d=%d l=%d recv L g%d=l%d "V3DFMTSZ" from %d=%d,%d,%d\n",
		 pgU->myrank, "hierarchize3D", formHier, d, l, NL.v[d], jdL, 
		 V3DLST(nUd), pgU->getRank(src), V3DLST(src));
	MPI_Recv(urL, nSz, MPI_DOUBLE, pgU->getRank(src), HTAG, pgU->comm, &s);
	CHECK_RECV_OK(s, nSz);
	if (jdL + dl < nU.v[d])
	  hierUpdateB(formHier, u, nUd, d, jdL, dl, urL);
	else
	  assert (jdL == jdR); // this process holds a single dest point
      }

      if (0 <= jdR && jsR < jdR) { //jdR is here & needs points from right
	Vec3D<int> NsR = NR;  NsR.v[d] += dl - (NR.v[d] - (NL.v[d] + jdR));  
	assert(NsR.v[d] < NU.v[d]); assert(NsR.v[d] % dl2 == 0);
	Vec3D<int> src = pgU->getP0(NsR, NU);  
	assert (pgU->id.v[d] < src.v[d]);
	if (verbosity > 2)
	  printf("%d: %s(%d) d=%d l=%d recv R g%d=l%d "V3DFMTSZ" from %d="V3DFMT"\n",
		 pgU->myrank, "hierarchize3D", formHier, d, l, NR.v[d], jdR, 
		 V3DLST(nUd), pgU->getRank(src), V3DLST(src));

	MPI_Recv(urR, nSz, MPI_DOUBLE, pgU->getRank(src), HTAG, pgU->comm, &s);
	CHECK_RECV_OK(s, nSz);
	if (jdR - dl >= 0)
	  hierUpdateB(formHier, u, nUd, d, jdR, -dl, urR);
	else {
	  assert (jdL == jdR); // this process holds a single dest point
	  hierUpdateBB(formHier, u, nUd, d, jdR, urL, urR);
	}
      }
      for (int m = 0; m < nMsgs; m++) {
	  MPI_Status s;
	  MPI_Wait(&msgRequest[m], &s);
	  delete[] msgBuf[m]; 
      }
	
    } //for(l...)     
  delete[] urL; delete[] urR;
  }  //for(d...)
} //formHierarchy3D

// combine contiguous range of hierarchical subspaces hgMin:hG
void HierachDecomp3D::combHierUpSet3D(Vec3D<int> hgMin,
				      Vec3D<int> hg, HaloArray3D *u) {
  if (verbosity > 1)
    printf("%d: combHierUpSet(" V3DFMT "; " V3DFMT ")\n",
	   pgU->myrank, V3DLST(hgMin), V3DLST(hg));
  hCombCntCheckSet(hgMin, hg);
  if (!(hg <= gxU)) // this process's grid does not hold hg
    return;

  Vec3D<int> hgC = hg; //grid index used for combination
  Vec3D<int> NU = gc->gridSz(gxU), nU = pgU->G2L(NU), NL = pgU->L2G0(NU); 
  //global offset is 2^(gxU.x-hg.x), stride is 2^(gxU.x-hg.x+1) [similarly y,z]
  //exception: if hgMin.x==0, we include all hierarchies from 0: 
  //           offset is 0, stride is 2^(gxU.x-hg.x), have 2^hg.x+1 global elts
  Vec3D<int> hOffs = Vec3D<int>(0), hStride = Vec3D<int>(1);
  for (int d=0; d < 3; d++) {
    assert (hgMin.v[d] == 0 || hgMin.v[d] == hg.v[d]); // only cases handled
    if (d == 2 && gc->is2D()) {  // 2D case
      continue;
    }
    if (hg.v[d] == 0) { // single zero hierarchy (boundary) is a special case
      hStride.v[d] = 1 << gxU.v[d];
      hgC.v[d]++;
    } else if (hgMin.v[d] == hg.v[d]) { // combine on single index
      hOffs.v[d] = 1 << (gxU.v[d] - hg.v[d]);
      hStride.v[d] = hOffs.v[d] * 2;
      hgC.v[d]--; //grid size of hierarchy is one less than index in hierarchy 
    } else { // combine on indices 0..hg
      hStride.v[d] = 1 << (gxU.v[d] - hg.v[d]);
    }
    if (NL.v[d] >= hOffs.v[d]) { // local offset != global offset
      int m = (NL.v[d] - hOffs.v[d]) % hStride.v[d];
      hOffs.v[d] = (m > 0)? hStride.v[d] - m: 0;
    }
  }
  Vec3D<int> hNLoc = (nU - hOffs + hStride - 1) / hStride;
  if (verbosity > 1)
    printf("%d: combHierUpSet: o=" V3DFMT " str= " V3DFMT " n =" V3DFMT "\n",
	   pgU->myrank, V3DLST(hOffs), V3DLST(hStride), V3DLST(hNLoc));  

  // this call sets the combination grid's process grid, which must be 
  // comprised of the participating component grids' process grids before 
  // sub-griding is applied
  HaloArray3D *uHs = gc->setHierComb(hg, hgC);
  uHs->zero(); 

  for (int g=0; g < gc->nGrids(); g++) { 
    Vec3D<int> sgP = gc->pgs[g]->P, hgP = gc->gridSz1(hgC); 
    if (!(sgP <= hgP)) {//process grid is not always smaller than subspace size
      //we need to take a sub-grid of this process grid to maintain block dist.
      Vec3D<int> sgPoffs, sgPstr, zeroDel; 
      for (int d=0; d < 3; d++) {
	sgPstr.v[d] = (sgP.v[d] <= hgP.v[d])? 1: (hg.v[d] == 0)? sgP.v[d]:
	  sgP.v[d]/hgP.v[d];
	sgPoffs.v[d] = (hgMin.v[d] == 0)?  0: sgPstr.v[d]/2;
	zeroDel.v[d] = (hgMin.v[d] == 0) && (sgPstr.v[d] > 1);
	sgP.v[d] = std::min(sgP.v[d], hgP.v[d] + (hg.v[d]>0? zeroDel.v[d]: 0));
      }
      gc->pgs[g]->setSubgrid(sgPoffs, sgPstr, sgP, zeroDel);
    }
  }
  //note that, from the above, *pgU may now point to a sub-grid of the original

  //pad uH for single hierarchy: size is 2^k, GridCombine3D expects 2^k+1
  Vec3D<int> pad = pgU->lastProc() * hg.whereEq(hgMin); 
  if (gc->is2D()) pad.z = 0; // 2D case 

#if 1
  if (!(pgU->myrank==-1 || hNLoc + pad == pgU->G2L(gc->gridSz(hgC))) ||
      !(pgU->myrank!=-1 || (hNLoc + pad).prod() == 0))
    printf("%d: hgC=" V3DFMT " hN=" V3DFMT " pad=" V3DFMT " n= " V3DFMT
	   " N=" V3DFMT " id= " V3DFMT  " P=" V3DFMT "\n",  
	   pgU->myrankOrig, V3DLST(hgC), V3DLST(hNLoc), V3DLST(pad), 
	   V3DLST(gc->gridSz(hgC)),
	   V3DLST(pgU->G2L(gc->gridSz(hgC))), V3DLST(pgU->id), V3DLST(pgU->P));
  fflush(stdout);
#endif
  //check local size matches the expected value, whether or not participating
  assert (pgU->myrank!=-1 || (hNLoc + pad).prod() == 0);
  assert (pgU->myrank==-1 || hNLoc + pad == pgU->G2L(gc->gridSz(hgC))); 
  HaloArray3D *uH = new HaloArray3D(hNLoc + pad, u->halo); 
	  uH->zero(); // ensures padded elements are initialized
  packHier(hOffs, hStride, hNLoc, u, uH);

  SparseGrid3D *uHsg = new SparseGrid3D(uHs, gc->is2D(), hgC, gc->pgS);
  gc->gatherScatter(uH, uHsg); 
  gc->resetHierComb();
  delete uHsg;
  delete uHs;

  for (int g=0; g < gc->nGrids(); g++) {
    gc->pgs[g]->resetSubgrid();
  }

  unpackHier(hOffs, hStride, hNLoc, u, uH);
  if (verbosity > 4)
    uH->print(pgU->myrank, "hierUpSetComb");  
  delete uH;
} //combHierUpSet3D()


#ifndef MAX_COALESC
#define DH 0
#else
#define DH 1 
#endif

void HierachDecomp3D::hierarchComb2D(HaloArray3D *u) {
  hierarchize3D(true, u);
  if (verbosity > 4)
    u->print(pgU->myrank, "hier field");
#ifdef NO_COALESC
  for (int hy = 0; hy < gxS.y; hy++) 
    for (int hx = 0; hx < gxS.x; hx++) {
      if (hy + hx <= gxS.y + gxS.x - level)
	combHierUpSet3D(Vec3D<int>(hx, hy, 0), Vec3D<int>(hx, hy, 0), u); 
    }
#else
#ifndef NO_SCHED_COMP /* try to schedule for parallel computation */
  int schedStride = 2;
#else
  int schedStride = 1; 
#endif
  int hx = gxS.x - level + 1 + DH, hy = gxS.y - 1 - DH;
  while (hy >= gxS.y - level + 1 + DH) {
    for (int s = 0; s < schedStride; s++)
      // these should be able to run in parallel when schedStride > 1
      for (int i = s;  hy-i >= gxS.y - level + 1 + DH;  i+= schedStride)
	combHierUpSet3D(Vec3D<int>(hx+i,hy-i,0), Vec3D<int>(hx+i,hy-i, 0), u); 
#ifndef NO_SCHED_COMP
    schedStride++;
#endif
    hy--;
  } 
  for (int i=0; i < level-1-DH; i++) {
    int hy = gxS.y - 1 - i;
#ifndef NO_SCHED_COMP
    int hx = gxS.x - 1 - i; // enables limited ||ism of the 2 calls, i<level/2
#else
    int hx = gxS.x - level + 1 + i; 
#endif
    combHierUpSet3D(Vec3D<int>(0,hy,0), Vec3D<int>(gxS.x-level+DH,hy,0), u);
    combHierUpSet3D(Vec3D<int>(hx,0,0), Vec3D<int>(hx,gxS.y-level+DH,0), u); 
  }
  combHierUpSet3D(Vec3D<int>(0), Vec3D<int>(gxS.x-level+DH, gxS.y-level+DH, 0),
		  u);
#endif /* else NO_COALESC */
  if (verbosity > 4)
    u->print(pgU->myrank, "hier field comb");
  hierarchize3D(false, u);
} //hierarchComb2D()


// formula to determine whether (i,j,k) is a h.s. index common to > 1 grid
bool HierachDecomp3D::isCommonHierSurplus3D(int i, int j, int k) {
  Vec3D<int> ijk = Vec3D<int>(i, j, k);
  assert (i < gxS.x && j < gxS.y && k < (gc->is2D()? 1: gxS.z));
  int sijk = ijk.sum(), sg = gxS.sum();
  if (gc->is2D()) 
    return (sijk <= sg - level);
  for (int d=0; d<3; d++) // check if any 2 indices fail the 2D criteria
    if (sijk - ijk.v[d] > sg - gxS.v[d] - level)
      return false;
  return (sijk <= sg - 2*level+1);
}

void HierachDecomp3D::hierarchComb3D(HaloArray3D *u) {
  hierarchize3D(true, u);
  if (verbosity > 4)
    u->print(pgU->myrank, "hier field");
#ifdef NO_COALESC
  for (int hz = 0; hz < std::max(gxS.z, 1); hz++) 
    for (int hy = 0; hy < gxS.y; hy++) 
      for (int hx = 0; hx < gxS.x; hx++) {
	if (isCommonHierSurplus3D(hx, hy, hz))
	  combHierUpSet3D(Vec3D<int>(hx, hy, hz), Vec3D<int>(hx, hy, hz), u); 
    }
#else
#ifndef NO_SCHED_COMP /* try to schedule for parallel computation */
  int schedStride = 2;
#else
  int schedStride = 1; 
#endif
  Vec3D<int> ht = gxS - (level - 1 - DH);
  for (int strCol = 0;  strCol < schedStride*schedStride; strCol++) {
    for (int l0 = level-1-DH; l0 > 0; l0--) {
      int i=0, j=0;
      for (int k=l0-1;  k >= 0;  k--,j++) {
	for (int di=0; di <= k; di++) {
	  int ijCol = (i+di)%schedStride + schedStride * (j%schedStride);
	  assert (i+di+j+k-di == l0-1);//traverse along tri plane at dist. l0-1
	  Vec3D<int> h = ht + Vec3D<int> (i+di, j, k-di); 
	  if (ijCol == strCol) //same colors can run in || 
	    combHierUpSet3D(h, h, u);
	}
      }
    }
  }

  Vec3D<int> hl = gxS - (level - DH);
  for (int i = 0; i < level-1-DH; i++)
    for (int j = 0; j < level-1-DH - i; j++) {
      // printf("%d: hierarchComb3D %d,%d\n", pgU->myrank, i, j);
      combHierUpSet3D(Vec3D<int>(0,      ht.y+i,  ht.z+j), 
		      Vec3D<int>(hl.x,   ht.y+i,  ht.z+j), u);
      combHierUpSet3D(Vec3D<int>(ht.x+i, 0,       ht.z+j), 
		      Vec3D<int>(ht.x+i, hl.y,    ht.z+j), u);
      combHierUpSet3D(Vec3D<int>(ht.x+i, ht.y+j,  0), 
		      Vec3D<int>(ht.x+i, ht.y+j,  hl.z), u);      
    }

  for (int i=0; i < level-1-DH; i++) {
    // have limited ||ism of 3 calls, where i<level/2
    Vec3D<int> h = gxS - (i + 1);
    combHierUpSet3D(Vec3D<int>(0,  0,  h.z), Vec3D<int>(hl.x,hl.y,h.z), u);
    combHierUpSet3D(Vec3D<int>(0,  h.y,0),   Vec3D<int>(hl.x,h.y, hl.z), u);
    combHierUpSet3D(Vec3D<int>(h.x,0,  0),   Vec3D<int>(h.x, hl.y,hl.z), u); 
  }
  combHierUpSet3D(Vec3D<int>(0), hl, u);
#endif /* else NO_COALESC */
  if (verbosity > 4)
    u->print(pgU->myrank, "hier field comb");
   hierarchize3D(false, u);
} //hierarchComb3D()





HierachDecomp3D::HierachDecomp3D(GridCombine3D *gc_, int verb) {
  gc = gc_;  
  level = gc->level;
  gxU = gc->gxU;
  pgU = gc->pgU;
  gxS = gc->gxS;
  verbosity = verb;
  hCombCnt = new int** [gxS.z+1];
  for (int k=0; k <= gxS.z; k++) {
    hCombCnt[k] = new int* [gxS.y+1];
    for (int j=0; j <= gxS.y; j++) {
      hCombCnt[k][j] = new int [gxS.x+1];
      for (int i=0; i <= gxS.x; i++)
	hCombCnt[k][j][i] = 0;
    }
  }
} //HierachDecomp3D()

void HierachDecomp3D::hCombCntCheckSet(Vec3D<int> h0, Vec3D<int> h1) {
  for (int k=h0.z; k <= h1.z; k++) {
    assert (0 <= k && k < (gc->is2D()? 1: gxS.z));
    for (int j=h0.y; j <= h1.y; j++) {
      assert (0 <= j && j < gxS.y);
      for(int i=h0.x; i <= h1.x; i++) {
	assert (0 <= i && i < gxS.x);
	if (!isCommonHierSurplus3D(i, j, k))
	  printf("%d: hCombCntCheckSet i=%d,j=%d,k=%d\n",  pgU->myrankOrig,   
		 i, j, k);
	assert (isCommonHierSurplus3D(i, j, k));
	if (hCombCnt[k][j][i] != 0)
	  printf("%d: hCombCntCheckSet done i=%d,j=%d,k=%d\n",  pgU->myrankOrig, 
		 i, j, k);
	assert (hCombCnt[k][j][i] == 0);
	hCombCnt[k][j][i]++;
      }
    }
  }
}  

void HierachDecomp3D::hCombCntCheck() {
  for (int k=0; k < (gc->is2D()? 1: gxS.z); k++) 
    for (int j=0; j < gxS.y; j++) 
      for(int i=0; i < gxS.x; i++) {
	if (isCommonHierSurplus3D(i, j, k)) {
	  if (hCombCnt[k][j][i] != 1)
	    printf("%d: hCombCntCheck:i=%d,j=%d,k=%d v=%d\n",  pgU->myrankOrig, 
		   i, j, k, hCombCnt[k][j][i]);
	  assert (hCombCnt[k][j][i] == 1);
	} else {
	  if (hCombCnt[k][j][i] != 0)
	    printf("%d: hCombCntCheck:i=%d,j=%d,k=%d v=%d\n",  pgU->myrankOrig, 
		   i, j, k, hCombCnt[k][j][i]);
	  assert (hCombCnt[k][j][i] == 0);
	}
	hCombCnt[k][j][i] = 0; // clear for repeated combinations
      }
}  

void HierachDecomp3D::gatherScatter(HaloArray3D* u) {
  if (gc->is2D())
    hierarchComb2D(u);
  else
    hierarchComb3D(u);
  hCombCntCheck();
} //gatherScatter()

HierachDecomp3D::~HierachDecomp3D() {
  for (int k=0; k <= gxS.z; k++) {
    for (int j=0; j <= gxS.y; j++)
      delete[] hCombCnt[k][j];
    delete[] hCombCnt[k];
  }
  delete[] hCombCnt;
}  

