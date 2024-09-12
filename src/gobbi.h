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
#include "TDC1190.h"

//two structures for matching up CsI times and Energies
struct dataE
{
  int ienergy;
  int id;
  int itime;
  float energy;
  float time;
};
struct dataT
{
  int itime;
  int id;
};

//Struct needed to pass s800 results in for time gates
struct s800_results
{
  bool trig_coin;
  bool trig_singles;
  bool trig_s800_singles;
  int Zbeam;
  int Abeam;
  int Zresidue;
  int Aresidue;
  void Reset()
  {
    trig_coin = false;
    trig_singles = false;
    trig_s800_singles = false;
    Zbeam =-1;
    Abeam=-1;
    Zresidue=-1;
    Aresidue=-1;
  }
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
  TDC1190 * CsITDC;

  calibrate * FrontEcal;
  calibrate * BackEcal;
  calibrate * CsIEcal;

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
  void addFrontEvent(int quad, unsigned short chan, unsigned short high, unsigned short low, unsigned short time);
  void addBackEvent(int quad, unsigned short chan, unsigned short high, unsigned short low, unsigned short time);
  void addCsIEvent(int quad, unsigned short chan, unsigned short high, unsigned short low, unsigned short time);

  int analyze(s800_results S800_results);
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
  dataE DataE[56];
  dataT DataT[56];
  void MatchCsIEnergyTime();

  int GrunNum = 0;

  int prtswtch;

  int counter=0; //counts unpack numbers

  int csiTimeMin[20];
  int csiTimeMax[20];

};


#endif
