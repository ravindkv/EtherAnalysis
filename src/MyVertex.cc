#include "interface/MyVertex.h"

MyVertex::MyVertex():
 chi2(-9.0),
 normalizedChi2(-9.0),
 ndof(0),
 totVtx(0),
 rho(-9.0),
 isValid(true),
 rhoAll(-9)
{
}

MyVertex::~MyVertex()
{
}

void MyVertex::Reset()
{
 chi2 = -9.0;
 normalizedChi2 = -9.0;
 ndof = 0;
 totVtx =0;
 rho = -9.0;
 ErrXYZ.SetCoordinates(0.0,0.0,0.0);
 isValid = true;
 XYZ.SetCoordinates(0.0,0.0,0.0);
 rhoAll = -9;
}
