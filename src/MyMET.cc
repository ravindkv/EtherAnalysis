#include "interface/MyMET.h"

MyMET::MyMET():
  
  metName(""),
  sumEt(-9.)
{
}

MyMET::~MyMET()
{
}
void MyMET::Reset()
{
  metName = "";
  sumEt = -9.;
  p4.SetCoordinates(0.0,0.0,0.0,0.0);
}
