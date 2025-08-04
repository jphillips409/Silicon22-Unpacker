#include "histo_sort.h"
#include "histo_read.h"
#include "TRandom.h"
#include "caen.h"
#include "mtdc.h"
#include "calibrate.h"
#include "TVector3.h"
#include <cmath>


struct dataEC
{
  int ienergy;
  int id;
  int itime;
  float time;
  int iRing;
  int iLoc;
  int iadc;
  int ichan;
  int addbackmult;
  float theta;
  float energy;
  float Total;
  float dop_energy;
  float phi;
  TVector3 pos;
};

struct dataTC
{
  int itime;
  int id;
  float time;
  int itdc;
  int ichan;
};

struct mapC
{
  int iRing;
  int iLoc;
  int iTDC;
  int iTDCChan;
};

class ceasar
{
 private:
  TRandom * ran;
  histo_sort * Histo;
  histo_read * Histo_read;
  //doppler *Doppler;
  caen Caen[6];
  mtdc Mtdc[6];

  mapC MapC[6][32];
  calibrate * calCeasar;
  calibrate * calCeasarT;

 public:
  ceasar(TRandom* ran0, histo_sort * Histo0,histo_read * Histo1, float shift0);
  ~ceasar();
//  bool unpack(unsigned short *point);
  bool unpack(unsigned short *point, int runno);

  void init();

  string chipmap;
  string calfile;
  string posmap;
  
  dataEC DataEC[192];
  dataTC DataTC[192];
  
  dataEC select[192];
  dataEC added[192];
  dataEC tempAdd; //For summing and then passing to added
  dataEC NoTadded[192];
  int N_NoTaddback;
  int Nselect;
  int Nadded;

  float mag[10][24]; //10 rings and 24 max detectors
  float angle[10][24]; //10 rings and 24 max detectors
  float angle2[10][24];
  float Txfp[3];
  int Nxfp;
  float TRF;

  int NE;
  int NT;
  int NTSelect; //Number of select gammas


  void Reset();

};
