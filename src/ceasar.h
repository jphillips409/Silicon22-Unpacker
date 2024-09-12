#include "histo_sort.h"
#include "histo_read.h"
#include "TRandom.h"
#include "caen.h"
#include "TDC1190.h"
#include "calibrate.h"
#include <cmath>


struct dataEC
{
  int ienergy;
  int id;
  int itime;
  int iRing;
  int iLoc;
  int iqdc;
  int ichan;
  int addbackmult;
  float theta;
  float energy;
  float Total;
  float dop_energy;
  float phi;
  float time;
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
  caen Caen;
  TDC1190 ** tdc;

  mapC MapC[6][32];
  calibrate * calCeasar;
  calibrate * calCeasarT;

 public:
  ceasar(TRandom* ran0, histo_sort * Histo0,histo_read * Histo1, float shift0);
  ~ceasar();
//  bool unpack(unsigned short *point);
  bool unpack(unsigned short *point, unsigned short *tdcpoint, int runno);

  void init(float shift);

  dataEC DataEC[192];
  dataTC DataTC[192];
  
  dataEC select[192];
  dataEC added[192];
  dataEC NoTadded[192];
  int N_NoTaddback;
  int Nselect;
  int Nadded;

  float angle[10][24]; //10 rings and 24 max detectors
  float angle2[10][24];
  float Txfp[3];
  int Nxfp;
  float TRF;

  int NE;
  int NT;


  void reset();

};
