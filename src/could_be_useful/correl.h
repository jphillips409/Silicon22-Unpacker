#include "solution.h"
#include <iostream>


/**
 * !\brief multi-particle correlations
 *
 * This class histograms the kinetic energy of multi-particle
 * events in their center-of-mass frame and looks for the 
 * correlations between the fragments
 */


class correl
{
 public:
#ifdef rel
  CEinstein Kinematics;
#else
  CNewton Kinematics;
#endif
  solution *proton[7];
  solution *alpha[5];
  solution *H2[5];
  solution *He3[5];
  solution *H3[5];
  solution *Li6[5];
  solution *Be7[5];
  solution *B8[5];
  solution *B10[5];
  solution *C9[5];
  solution *C9Fake[5];
  solution *C9Fake_Channeling[5];
  solution *C9Fake_CsIreaction[5];
  solution *C9Fake_CsIreaction_C11[5];
  solution *C10Fake_CsIreaction[5];
  solution *C10[5];
  solution *C11[5];
  solution *C12[5];
  solution *O13[5];
  solution *frag[7];
  solution *gamma[5];
  int N;
  float energyTot;

  void load(solution*);
  void Reset();
  float findErel();
  void setMask();
  void zeroMask();
  void zeroProton();
  void zeroAlpha();
  void zero3He(); 
  void zeroH2();
  void printP();
  void makeArray(bool flagMask = 0);

  int multProton;
  int multAlpha;
  int mult3He;
  int multH2;
  int mult6Li;
  int mult7Be;
  int multH3;
  int multGamma;
  int mult10B;
  int mult8B;
  int mult9C;
  int mult9CFake;
  int mult9CFake_Channeling;
  int mult9CFake_CsIreaction;
  int mult9CFake_CsIreaction_C11;
  int mult10CFake_CsIreaction;
  int mult10C;
  int mult11C;
  int mult12C;
  int mult13O;

  float Vcmframevector[3];
  float Vcmframe;
  float VcmframePerp;
  float VcmframePara;
  float TKEL;
  float Erel;
  float KEToT;
  float thetaCM;
  float phiCM;
  float thetaInReactCM;

  float check_ke;
  float check_mass;

  int maskProton[5];
  int maskAlpha[5];
  int mask3He[5];
  int mask2H[5];
  int mask3H[5];
  int mask6Li[5];
  int mask7Be[5];
  int mask9C[5];
  int mask9CFake[5];
  int mask9CFake_Channeling[5];
  int mask9CFake_CsIreaction[5];
  int mask9CFake_CsIreaction_C11[5];
  int mask10CFake_CsIreaction[5];
  int mask8B[5];
  int mask10B[5];
  int mask10C[5];
  int mask11C[5];
  int mask12C[5];
  int mask13O[5];

  int ifront;
  int istrip;
  int itele;
  // int Norder;
  float Vlab;

  int multAlp[14];
  int multPro[14];
  int mult3he[14];
  int multh2[14];
  int mult6li[14];
  int mult7be[14];
  int mult9c[14];
  int mult9cFake[14];
  int mult9cFake_Channeling[14];
  int mult9cFake_CsIreaction[14];
  int mult9cFake_CsIreaction_C11[14];
  int mult10cFake_CsIreaction[14];
  int mult8b[14];
  int mult10b[14];
  int mult10c[14];
  int mult11c[14];
  int mult12c[14];
  int mult13o[14];
  int multh3[14];
  int mult_Gamma[14];

  float thetaAlpha;
  float phiAlpha;
  float Thetaopen;
  float ppThetaRel;
  float cosThetapp;
  // float protvel1; 
  // float protvel2; 
  float pTheta0;
  float pTheta2;
  float pTheta3;
  float pTheta1;
  float aTheta0;
  float pThetaMax;
  float aTheta1;
  float thetadip;
  float Epp;
  float CosAngle;

  float CosAngleY1;
  float CosAngleY2;
  float Ey1;
  float Ey2;


  float Ey1_10;
  float Ey2_10;
  float CosAngleAlphaY1_10;
  float CosAngleAlphaY2_10;

  float Epp_10;
  float CosAngle_10;

  float CosAngleAlpha;
  float CosAngleAlphaY1;
  float CosAngleAlphaY2;
  float Eaa;
  float ThetaP;
  float ThetaA1;
  float ThetaA2;
  float PhiP;
  float PhiA1;
  float PhiA2;
  float cosAA;
 
  float aarelbeamangle;
  float prelbeamangle;

  float Eap_min, Eap_max;

  float ThetaA1_P, ThetaA2_P;

  float cosAlpha[2];
  float cosProton0[2];
  float cosProton1[2];
  float energy3;

  float velocityCM;
  float momentumCM;


  //jacobi
  void getJacobi();
  float cosThetaT;
  float cosThetaY[2];
  float cosThetaV;
  float getAlphaMom();
  float x_T;
  float x_Y[2];
  void rotate();


  float momC[3];
  float Mtot[3];
  float cosThetaC;
  float PperpC;
  float PparaC;
  float PtotC;
  float cosAlphaQ;

};
