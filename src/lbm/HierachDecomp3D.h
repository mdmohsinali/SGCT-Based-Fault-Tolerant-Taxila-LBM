/* ParSGCT Code source file.
Copyright (c) 2015 Peter Edward Strazdins. All rights reserved.
Licensed under the terms of the BSD License as described in the LICENSE_SGCT file.
This comment must be retained in any redistributions of this source file.
*/

// Parallel 3D Sparse Grid Combination via Hierarchical Basis class
// written by Peter Strazdins, Dec 14

#include "Vec3D.h"
#include "ProcGrid3D.h"
#include "HaloArray3D.h"
#include "GridCombine3D.h"

class HierachDecomp3D {
 protected:
  // helper functions manipulating HaloArray3D data
  static void hierUpdate(bool formHier, HaloArray3D *u, Vec3D<int> nU,
			 int d, int j0, int dl);
  static void hierUpdateB(bool formHier, HaloArray3D *u, Vec3D<int> nU,
			  int d, int j0, int dl, double *b) ;
  static void hierUpdateBB(bool formHier, HaloArray3D *u, Vec3D<int> nU, 
			   int d, int j0, double *bl, double *br);
  static void packHier(Vec3D<int> hOffs, Vec3D<int> hStr, Vec3D<int> hN,
		       HaloArray3D *u, HaloArray3D *uH);
  static void unpackHier(Vec3D<int> hOffs, Vec3D<int> hStr, Vec3D<int> h,
			 HaloArray3D *u, HaloArray3D *uH);

  // higher level helper functions
  void hierarchize3D(bool formHier, HaloArray3D *u);
  void combHierUpSet3D(Vec3D<int> hgMin, Vec3D<int> hg, HaloArray3D *u);
  void hierarchComb2D(HaloArray3D *u); 
  bool isCommonHierSurplus3D(int i, int j, int k);
  void hierarchComb3D(HaloArray3D *u);   

  void hCombCntCheckSet(Vec3D<int> h0, Vec3D<int> h1) ;
  void hCombCntCheck();

 protected:
  int verbosity; 
  Timer *timer;
  GridCombine3D *gc;
  int ***hCombCnt;      //3D array of size gxS testing if eahc suplus coverd 1x

 public:
  // cached versions of gc's fields
  Vec3D<int> gxU; // this process' component grid's co-ordinates
  ProcGrid3D* pgU; // this process' component grid's process grid
  Vec3D<int> gxS; // combination grid co-ordinates
  int level;

  HierachDecomp3D(){};
  // create object using gc's grid parameters
  HierachDecomp3D(GridCombine3D *gc, int verbosity=0);
  ~HierachDecomp3D();

  void gatherScatter(HaloArray3D* u);

}; //HierachDecomp3D

