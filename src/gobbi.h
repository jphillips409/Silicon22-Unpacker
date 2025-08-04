#ifndef Gobbi_H
#define Gobbi_H
#include <iostream>
#include <string>
#include "TRandom.h"
#include "histo_sort.h"
#include "telescope.h"
#include "calibrate.h"
#include "HINP.h"
#include "caen.h"
#include "mtdc.h"
#include "s800_results.h"

//two structures for matching up CsI times and Energies
struct dataE
{
  int ienergy;
  int id;
  int itime;
  float energy;
  float time;
	float qdcmatch; //psd parameter, filled in psd calc
  float qdc; //psd parameter, filled in psd calc
};
struct dataT
{
  int itime;
  int id;
  float time;
  int itdc;
  int ichan;
};
struct dataQ
{
  int ienergy;
  int id;
  int itime;
  float energy;
  float time;
};

class gobbi
{
 public:
  float TargetThickness;

  gobbi(TRandom* ran, histo_sort * Histo1);
  ~gobbi();
  histo_sort * Histo;
  telescope * Telescope[4];
  telescope * S800;
  HINP * SiADC;
  caen * CsIADC;
  mtdc * CsITDC;
  caen * CsIQDC;

  calibrate * FrontEcal;
  calibrate * BackEcal;
  calibrate * CsIEcal;
  calibrate * CsIEcal90;

  calibrate * FrontTimecal;
  calibrate * BackTimecal;
  calibrate * CsITimecal;

  calibrate * FrontLowEcal;
  calibrate * BackLowEcal;

  void SetTarget(double Targetdist, float TargetThickness);
  void reset();

  bool unpack(unsigned short*& point,unsigned short*& tdcpoint,int runno);
  bool unpackSi_HINP4(unsigned short*& point);
  bool unpackCsI_ADC(unsigned short*& point);
  bool unpackCsI_TDC(unsigned short*& point);
  bool unpackCsI_QDC(unsigned short*& point);
  void addFrontEvent(int quad, unsigned short chan, unsigned short high, unsigned short low, unsigned short time);
  void addBackEvent(int quad, unsigned short chan, unsigned short high, unsigned short low, unsigned short time);
  void addCsIEvent(int quad, unsigned short chan, unsigned short high, unsigned short low, unsigned short time);

  int analyze(s800_results S800_results, int runno);
  void SiNeighbours();
  int matchTele();

  int multidEE;
  int multiECsI;
  int NsimpleECsI = 0;
  int NmultiECsI = 0;
  int NSameCsI = 0;

  int Gfvb_cnt[4];

  int NE;
  int NT;
  int NQ;
  dataE DataE[56];
  dataT DataT[56];
	dataQ DataQ[56];
  void MatchCsIEnergyTime(int runno);
	void CsIQDCMatch();

  int GrunNum = 0;

  int prtswtch;

  int counter=0; //counts unpack numbers

  int csiTimeMin[20];
  int csiTimeMax[20];

  int num_protons = 0;

  int num1p = 0;
  bool flag1p = 0;
  int num2p = 0;
  bool flag2p = 0;
  bool flag3p = 0;
  bool flagalpha = 0;

  //Tracking types of S800 coincidence events
  bool S800coinc_Si = false; //Only Si
  bool S800coinc_CsI = false; //Only CsI
  bool S800coinc_both = false; //Si and CsI
  bool S800_complete = false; //Matched Si/CsI

};


#endif
