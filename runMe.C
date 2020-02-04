#include "TString.h"
#include "TROOT.h"
//check
void runMe(char const *arg){
  gROOT->ProcessLine(TString::Format(".L %s.C+",arg));
  gROOT->ProcessLine(TString::Format("%s t",arg));
  gROOT->ProcessLine("t.processEvents()");
}
