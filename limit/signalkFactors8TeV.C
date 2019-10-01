
double MT = 174.3; 
//---------------------------------
//the momentum transfer scale Q2 at which to evaluate alphaEM
//https://github.com/pinkenburg/pythiaeRHIC/blob/048fc8f55d8ef1c8efce0ede8877e317325f5945/pythia/pyalem.f
//https://github.com/lovedeepkaursaini/lstarGI/blob/master/Theory/calc.cc
//---------------------------------
//See R. Kleiss et al, CERN 89-08, vol. 3, pp. 129-131.
double pyalem( double Q2){
  ////use pythia default
  double MSTU101 = 1; 
  double PARU104 = 1; 
  double PARU1 = 3.141592653589793; 
  double PARU101 = 0.00729735; 
  double AEMPI = PARU101/(3*PARU1);
  double PRIGG = 0; 
  double PARU103 = 0.007764;   /// 1/128.8; 
  if( MSTU101 <= 0 || Q2 < 2E-6 ) 
    PRIGG = 0; 
  else if (MSTU101 == 2 && Q2 < PARU104) 
    PRIGG = 0; 
  else if (MSTU101 == 2 ) 
    PRIGG = 1- PARU101/PARU103; 
  else if(Q2 < 0.09) 
    PRIGG = AEMPI*(13.4916+log(Q2))+0.00835*log(1+Q2) ;
  else if(Q2 < 9)
    PRIGG = AEMPI*(16.32 + 2*log(Q2))+ 0.00238*log(1+3.927*Q2); 
  else if(Q2 <1E4)
    PRIGG = AEMPI*(13.4955 + 3*log(Q2))+ 0.00165+0.00299*log(1+Q2); 
  else PRIGG = AEMPI *(13.4955+3*log(Q2))+0.00221+0.00293*log(1+Q2);
  double PYALEM = PARU101/(1-PRIGG); 
  return PYALEM; 
}

//---------------------------------
//decay width
//---------------------------------
//from pywidth.f
double wid(double M, double L ,int cha){
  double SH = M*M; 
  double SHR = M; 
  double AEM = 0.0079772474; /// about 1/137.035
  ///running aem. 
  AEM =  pyalem(SH);
  double FAC = SH/(L*L)*SHR; 
  double QF; 
  double f = 1.0;   //// f
  double ff = 1.0;   /// f'
  double xw = 0.23122;  ///dpg 2007 constant. sin^2(theta_w). pythia use 0.232; 
  double xw1 = 1-xw; 
  double WDTP; 
  double MZ = 91.1876; 
  double MW = 80.403; 
  double RM1; 
  if(cha ==1){  ///mu->mug
    QF = -f/2.0 - ff/2.0; 
    WDTP = FAC*AEM*QF*QF*1/4.0; 
  }else if( cha==2){///mu*->Zmu
    QF = -f*xw1/2.0+ff*xw/2.0; 
    RM1 = MZ*MZ/SH; 
    WDTP = FAC*AEM*QF*QF/(8*xw*xw1)*(1-RM1)*(1-RM1)*(2+RM1);
  }else if( cha==3){ ///mu*->Wnu
    //QF = f/sqrt(2.0);
    RM1 = MW*MW/SH; 
    WDTP = FAC*AEM*f*f/(16*xw)*(1-RM1)*(1-RM1)*(2+RM1);
  }else{
    cout<<"channel: "<<cha<<"  not defined"<<endl;
    exit(1);
  }
  return WDTP; 
}

//---------------------------------
////for e*->e+t+tbar needs to calculate again( MT= 175, not small for small e*)
//---------------------------------
double width_ci_ttbar(double *x,double *par){
  double Lambda = par[0];
  double M = par[1];
  double m23 = x[0]; ///
  //m23 = sqrt(m23);
  double PI = 3.141592653589793; 
  double g = sqrt(4*PI);
  double m2 = 0; ///lepton mass set zero
  double m3 = MT; 
  double m4 = MT;
  double X = (-pow(m2,2) + pow(m23,2) + pow(m3,2))/(2.*m23);
  double Y = (pow(M,2) - pow(m23,2) - pow(m4,2))/(2.*m23);
  double A = pow(X+Y,2)- pow( sqrt(X*X-m3*m3)+sqrt(Y*Y-m4*m4),2);
  double B = pow(X+Y,2)- pow( sqrt(X*X-m3*m3)-sqrt(Y*Y-m4*m4),2);
  double MSquare = (2*pow(g,4)*(-pow(m2,2) + pow(m23,2) - pow(m3,2))*(pow(M,2) - pow(m23,2) + pow(m4,2)))/
    pow(Lambda,4);
  double val = 1/pow(2*PI,3) * 1/(32*pow(M,3)) * MSquare * (B-A);
  val *= 2*m23; 
  return val; 
}

//---------------------------------
// width for contact interaction
//---------------------------------
double widci(double isCorrection, double M,double L){
  //// 1/(96*pi)*NC*S*M^5/L^4
  //double SHR = M; 
  //double NC = 
  double PI = 3.141592653589793; 
  double WDTP = 1/(96*PI)*pow(M,5)/pow(L,4); 
  double sum; 
  double NC = 1; 
  double S = 2; 
  ///mu*->mu+mumu
  sum = WDTP*NC*S; 
  NC = 1; 
  S = 1; 
  ///mu*->mu+ e/tau + E/Tau  or 3 neutrino pairs
  sum += WDTP*NC*S*5; 
  NC = 3; 
  S = 1; 
  ///mu*->mu+ qQ; 
  //sum += WDTP*NC*S*6; 
  sum += WDTP*NC*S*5; // mass set to be 0 
  if(M>2*MT && isCorrection){
    double m23Square_l = pow(MT,2);
    double m23Square_u = pow(M-MT,2);
    ///top
    TF1 *ftmp = new TF1("ftmp",width_ci_ttbar,sqrt(m23Square_l),sqrt(m23Square_u),2);
    ftmp->SetParameter(0,L);
    ftmp->SetParameter(1,M);
    double tmp = ftmp->Integral(sqrt(m23Square_l),sqrt(m23Square_u));
    sum += tmp*S*NC; 
  }
  else if(M>2*MT) sum += WDTP*S*NC;
  return sum; 
}

///Shilpi - May 29, 2011
/////e*->eee
double widci_ee(double M,double L){
  //// 1/(96*pi)*NC*S*M^5/L^4
  double PI = 3.141592653589793; 
  double WDTP = 1/(96*PI)*pow(M,5)/pow(L,4); 
  double sum = 0; 
  double NC = 1; 
  double S = 2; 
  ///mu*->mu+mumu
  sum = WDTP*NC*S; 
  return sum; 
}

//e*->e mu/tau mu/tau or 2 neutrinoes
double widci_mutau(double M,double L){
  //// 1/(96*pi)*NC*S*M^5/L^4
  double PI = 3.141592653589793; 
  double WDTP = 1/(96*PI)*pow(M,5)/pow(L,4); 
  double sum = 0; 
  double NC = 1; 
  double S = 2; 
  ///mu*->mu+mumu
  //sum = WDTP*NC*S; 
  NC = 1; 
  S = 1; 
  ///mu*->mu+ e/tau + E/Tau  or 3 neutrino pairs
  sum += WDTP*NC*S*5; 
  return sum; 
}

//e*->eqq
double widci_qq(double M,double L){
  //// 1/(96*pi)*NC*S*M^5/L^4
  //double SHR = M; 
  //double NC = 
  double PI = 3.141592653589793; 
  double WDTP = 1/(96*PI)*pow(M,5)/pow(L,4); 
  double sum = 0; 
  double NC = 1; 
  double S = 2; 
  NC = 3; 
  S = 1; 
  ///mu*->mu+ qQ; 
  //sum += WDTP*NC*S*6; 
  sum += WDTP*NC*S*5; // mass set to be 0 
  if(M>2*MT){
    double m23Square_l = pow(MT,2);
    double m23Square_u = pow(M-MT,2);
    ///top
    TF1 *ftmp = new TF1("ftmp",width_ci_ttbar,sqrt(m23Square_l),sqrt(m23Square_u),2);
    ftmp->SetParameter(0,L);
    ftmp->SetParameter(1,M);
    double tmp = ftmp->Integral(sqrt(m23Square_l),sqrt(m23Square_u));
    sum += tmp*S*NC; 
  }
  return sum; 
}
//////////////////////////////////////////////////////////////////////////
///simplified branching ratio vs M/L
void BF_factor(double res[]){
  ////from pywidth.f
  double AEM = 0.0079772474; /// about 1/137.035
  //double QF; 
  double f = 1.0;   //// f
  double ff = 1.0;   /// f'
  double xw = 0.23122;  ///sin^2(theta_w). 
  double xw1 = 1-xw; 
  double QF1 = -f/2.0 - ff/2.0; 
  double QF2 = -f*xw1/2.0+ff*xw/2.0; 
  double A = 1/4.0*AEM*QF1*QF1; 
  double B = 1/(8*xw*xw1)*AEM*QF2*QF2*2;
  double C = AEM*f*f/(16*xw)*2; 
  double D = 25/(96*3.1415926535);
  res[0] = A; 
  res[1] = B; 
  res[2] = C; 
  res[3] = D; 
}

double pyalps( double Q2){
  ///order to run alpahS
  int MSTU111 = 1; 
  ///default quark mass.
  ///starting from 1. 
  double PMAS[7]={0,0.33,0.33,0.5,1.5,4.8,175};
  /////coeff
  double STEPDN[7] ={0,0,0,0.10568,0.13398,0.17337,0};
  double STEPUP[7]={0,0,0,0,-0.11413,-0.14563,-0.18988};
  int NF = 4; 
  double PARU112 =0.192;
  double PARU1 = 3.141592653589793; 
  double PARU111 = 0.2; 
  //// 0-order alpahS
  if(MSTU111 <=0 ){
    double val = PARU111; 
    return val; 
  }
  double Q2EFF = Q2; 
  double ALAM2 = PARU112*PARU112; 
  /// C... MSTU(115), default 0. MSTU(112) default 5, the nominal # of flavours assumed in//
  ///C... alphas expression
  ////C... paru(112)  D = 0.25GeV Lambad used in running alphaS, acutally print out  0.192 
  ////C... paru(113) D = 1, the flavour threshold for the effecitve number of flavours Nf to use
  ///C... in the alphaS expression, Q2 = paru113*mq^2. 
  ///C... paru(114) D = 4GeV^2  Q2 values below which alphaS is assumed to be constant 
  ////C... MSTU(113)  default 3: minimum number of flavours that may be assumed in alphas expression.
  ///C... MSTU(114), default 5, maximum....
  ////C... MSTU(112), the nominal number of ...
  int MSTU113 = 3; 
  int PARU113 = 1; 
  int MSTU114 = 5; 
  ///// treatment of alphS singularity for Q2->0. 
  int MSTU115 = 0;  ///allow it to diverge like 1/log(Q2/Lambda2);
  while(NF >  max(3,MSTU113)){
    double Q2THR = PARU113*pow(PMAS[NF],2);
    if(Q2EFF < Q2THR){
      NF = NF-1; 
      double Q2RAT = Q2THR/ALAM2; 
      ALAM2 = ALAM2*pow(Q2RAT,2.0/(33-2*NF));
      if(MSTU111 == 2) ALAM2 = ALAM2*pow(log(Q2RAT),STEPDN[NF]);
      ////continue; 
    }else break; 
  }
  cout<<"after 1st loop. "<<ALAM2<<" "<<Q2EFF<<endl;
  while( NF < min(6,MSTU114)){
    double Q2THR = PARU113*pow(PMAS[NF+1],2);
    if(Q2EFF> Q2THR){
      NF = NF +1; 
      double Q2RAT = Q2THR/ALAM2; 
      ALAM2 = ALAM2*pow(Q2RAT,-2.0/(33-2*NF));
      if(MSTU111 ==2) ALAM2 = ALAM2*pow(log(Q2RAT),STEPUP[NF]);
    }else break; 
  }
  cout<<"after 2nd loop. "<<ALAM2<<" "<<Q2EFF<<endl;
  if(MSTU115 == 1) Q2EFF = Q2EFF + ALAM2; 
  double tmp = Q2EFF/ALAM2;
  cout<<"Q2EFF/AlAM2: "<<Q2EFF<<" "<<ALAM2<<endl;
  cout<<"log: "<<log(tmp)<<endl;
  if(tmp<1.0001) tmp  = 1.0001; 
  double ALGQ = log(tmp);
  double B0 = (33-2*NF)/6.0; 
  double PARU115 = 10; 
  double val = 2*PARU1/(B0*ALGQ);
  if(val>PARU115) val = PARU115; 
  if(MSTU111 ==1)
    return val;
  ///2nd order

  // B1=(153D0-19D0*NF)/6D0
  //      PYALPS=MIN(PARU(115),PARU(2)/(B0*ALGQ)*(1D0-B1*LOG(ALGQ)/
  //   &  (B0**2*ALGQ)))
  double B1 = (153-19*NF)/6.0; 
  val = val*(1-B1*log(ALGQ)/(B0*B0*ALGQ)); 
  return val; 
}

///brachcing ratio of mu*->mu+ga
float getBranchRatioMuStarMG(float M,float Lambda){
  if(Lambda<1000){
    cout<<"warning Lambda<1TEV.."<<endl;
  }
  float br =wid(M,Lambda,1)/(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1)+widci(false, M,Lambda));
  return br; 
}

///////////////////////////////////SHILPI/////////////////////////////////
///brachcing ratio of l*->l+ga
float brlg(float M,float Lambda){
  if(Lambda<1000){
    cout<<"warning Lambda<1TEV.."<<endl;
  }
  float br =wid(M,Lambda,1)/(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1)+widci(false, M,Lambda));
  //float br =wid(M,Lambda,1)/(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1));
  return br; 
}

///brachcing ratio of l*->l+Z
float corrFacbrlz(float M,float Lambda){
  if(Lambda<1000){
    cout<<"warning Lambda<1TEV.."<<endl;
  }
  float br =(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1))/(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1)+widci(false, M,Lambda));
  return br; 
}

float brpho(float M,float Lambda){
  if(Lambda<1000){
    cout<<"warning Lambda<1TEV.."<<endl;
  }
  float br =wid(M,Lambda,1)/(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1)+widci(false, M,Lambda));
  return br; 
}

float brlz(float M,float Lambda){
  if(Lambda<1000){
    cout<<"warning Lambda<1TEV.."<<endl;
  }
  float br =wid(M,Lambda,2)/(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1)+widci(false, M,Lambda));
  return br; 
}

float brCI(float M,float Lambda){
  if(Lambda<1000){
    cout<<"warning Lambda<1TEV.."<<endl;
  }
  float br =widci(false, M,Lambda)/(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1)+widci(false, M,Lambda));
  return br; 
}

///brachcing ratio of l*->W+nu
float brnw(float M,float Lambda){
  if(Lambda<1000){
    cout<<"warning Lambda<1TEV.."<<endl;
  }
  float br =wid(M,Lambda,3)/(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1)+widci(false, M,Lambda));
  return br; 
}

///brachcing ratio of contact interaction
float brci(float M,float Lambda){
  if(Lambda<1000){
    cout<<"warning Lambda<1TEV.."<<endl;
  }
  float br =widci(false, M,Lambda)/(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1)+widci(false, M,Lambda));
  return br; 
}

///brachcing ratio of contact interaction: e*->eee
float brci_ee(float M,float Lambda){
  if(Lambda<1000){
    cout<<"warning Lambda<1TEV.."<<endl;
  }
  float br =widci_ee(M,Lambda)/(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1)+widci(false, M,Lambda));
  return br; 
}

///brachcing ratio of contact interaction: e*->e mu/tau mu/tau neutrinoes
float brci_mutau(float M,float Lambda){
  if(Lambda<1000){
    cout<<"warning Lambda<1TEV.."<<endl;
  }
  
  float br =widci_mutau(M,Lambda)/(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1)+widci(false, M,Lambda));
  return br; 
}

///brachcing ratio of contact interaction: e*->e qq
float brci_qq(float M,float Lambda){
  if(Lambda<1000){
    cout<<"warning Lambda<1TEV.."<<endl;
  }
  float br =widci_qq(M,Lambda)/(wid(M,Lambda,2)+wid(M,Lambda,3)+wid(M,Lambda,1)+widci(false, M,Lambda));
  return br; 
}

////////SHILPI////////////////////////////////////////////////////////
float Lambda0 = 10; 
float Xsec_Sig[10] = {1.610e-2,8.717e-3,4.806e-3,2.657e-3,1.474e-3,3.236e-4,6.589e-5};
float Mas_mustar[10]={200,400,600,800,1000,1500,2000};
//input the mass point j >=0 && j<=6 
///xlim is the cross section limit value 
double searchCS(int j,double xlimt){
  int useCorrDecay = 0; 
  ///first find L1 and L2
  double L1 = 3; 
  double L2 = 3; 
  double M = Mas_mustar[j];
  int found = 0; 
  double br ;
  double Lstart = 0.1; 
  if(j>=5) Lstart = 0.05; 
  for(double L = Lstart; L<=50; L += 1){
    L1 = L ;
    L2 = L+1; 
    double xs = Xsec_Sig[j]*pow(Lambda0,4)/pow(L1,4);
    ////corrected due to top
    if(useCorrDecay==1) br= wid(M,L1*1000,1)/(wid(M,L1*1000,2)+wid(M,L1*1000,3)+wid(M,L1*1000,1)+widci(false, M,L1*1000));
    else br = wid(M,L1*1000,1)/(wid(M,L1*1000,2)+wid(M,L1*1000,3)+wid(M,L1*1000,1)+widci(true, (M,L1*1000));
    double xsbr1 = xs*br; 
    xs = Xsec_Sig[j]*pow(Lambda0,4)/pow(L2,4);
    if(useCorrDecay==1) br= wid(M,L2*1000,1)/(wid(M,L2*1000,2)+wid(M,L2*1000,3)+wid(M,L2*1000,1)+widci(false, M,L2*1000));
    else br = wid(M,L2*1000,1)/(wid(M,L2*1000,2)+wid(M,L2*1000,3)+wid(M,L2*1000,1)+widci(true, (M,L2*1000));
    double xsbr2 = xs*br; 
    if(xsbr1>xlimt && xsbr2<xlimt){
      found =1; 
      break; 
    }
  } 
  if(found ==0) {
    cout<<"eroor not found "<<endl; 
    return -1 ;
  }
  double stepsize = 1; ///inGeV
  double errMin = 1;
  L1 *= 1000;
  L2 *= 1000;
  double Lm = 3;
  for(double L = L1; L<= L2; L += stepsize){
    double xsecmux =  Xsec_Sig[j] * pow(Lambda0,4)/pow(L/1000.,4);
    if(useCorrDecay==1) br= wid(M,L,1)/(wid(M,L,2)+wid(M,L,3)+wid(M,L,1)+widci(false, M,L));
    else br = wid(M,L,1)/(wid(M,L,2)+wid(M,L,3)+wid(M,L,1)+widci(true, (M,L));
    double xmuxmug = xsecmux*br;
    double err = fabs(xmuxmug - xlimt);
    if(err<errMin){
      errMin = err;
      Lm = L;
    }
  }
  if( errMin > 0.002){
    cout<<"wrong errMin "<< errMin <<endl; 
    exit(1);
  }
  errMin = 1;
  double Lm1 = Lm; 
  //low finally do once again search with useCorrDeacy =1 
  useCorrDecay = 1; 
  for(double L = Lm1 - 100; L<= Lm1 + 100; L += stepsize){
    double xsecmux =  Xsec_Sig[j] * pow(Lambda0,4)/pow(L/1000.,4);
    if(useCorrDecay==1) br= wid(M,L,1)/(wid(M,L,2)+wid(M,L,3)+wid(M,L,1)+widci(false, M,L));
    else br = wid(M,L,1)/(wid(M,L,2)+wid(M,L,3)+wid(M,L,1)+widci(true, (M,L));
    double xmuxmug = xsecmux*br;
    double err = fabs(xmuxmug - xlimt);
    if(err<errMin){
      errMin = err;
      Lm = L;
    }
  }
  if( errMin > 0.002){
    cout<<"wrong errMin2 "<< errMin <<endl; 
    exit(1);
  }
  return Lm; 
}
