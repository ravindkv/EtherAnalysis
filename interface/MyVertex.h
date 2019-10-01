#ifndef __MYVERTEX_H__
#define __MYVERTEX_H__

#include "TROOT.h"
#include <map>
#include <vector>
#include <string>

#include "MomentumVec.h"

class MyVertex
{
 public:
  MyVertex();
  ~MyVertex();
  
  void Reset();
  double chi2;
  double normalizedChi2;
  int  ndof;
  int totVtx;
  double rho;  // for maxd0 selection (<2)
  bool isValid;
  Point3D        XYZ;      //  position
  Point3D        ErrXYZ;   //  error position
  //fixedGridRhoAll: 
  //https://github.com/cms-analysis/flashgg/blob/e2fac35487f23fe05b20160d7b51f34bd06b0660/ Taggers/python/globalVariables_cff.py
  double rhoAll;  
  
 private :

};
#endif
