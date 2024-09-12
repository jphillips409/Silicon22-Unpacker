#ifndef det_
#define det_
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include "TRandom.h"
#include "TFile.h"
//#include "hira.h"
#include "gobbi.h"
#include "histo_sort.h"
#include "histo_read.h"
#include "correl2.h"
#include "doppler.h"
#include "constants.h"
#include "S800.h"
#include "ceasar.h"
#include "losses.h"
#include "solution.h"
#include "janus.h"

using namespace std;

/**
 *!\brief detector readout
 *
 * responsible for unpacking the data stream for physics events
 */
class det
{
 public:
  //  det(histo_sort * Histo_sort, forest * Forest);
  det(histo_sort* Histo0, histo_read* Histo1, int setting);
  ~det();
  void setRunno(int runno);
  bool unpack(ifstream *point,int runno,int sourceID,int fragmentsize);
  void Reset();

  TRandom * ran;
  gobbi *Gobbi;
  S800 *s800;
  janus *Janus;

  ceasar * Ceasar;
  //doppler is private

  CLosses *losses_fiber;
  CLosses *losses_target;
  
  //forest * Forest;
  
  void analyze(int event, int run);

  int Eventt;
  int treelevel;

  float Egamma0;
  float thetarel0;
  float res_Vel0;
  float thetagamma;
  float thetares;
  int detID;
  int Ring;
  int Loc;
  float Edop;
 
  int release0;
  int release3;
  int release4;
  int release5;
  int release6;
  int Ar33_37Cabeam;
  int Ar33_36Kbeam;

  int N_s800_singles=0;
  int N_coin=0;
  int N_singles=0;

  int Nresidue;
  int Nbadresidue;

  int solnZ;
  int solnnoZ;
  int NinvMassK35;
  int NinvMassK36;
  int justK35;
  int justK36;
  int K35twoprot;
  int K35alpha;
  int K36twoprot;
  int K36alpha;
  int K35other;
  int K36other;
  int K35withnozsecondary;
  int NAr32;

  int cntCa37 =0;
  int cntCa38 =0;

  int Nmaxmix;
  int EventMixCounter = 0;
  solution * EventMixerP[100];
  solution * EventMixer33Ar[100];
  int EventSingleMixCounter = 0;
  int NmaxmixSingles;
  solution * EventMixerSingleP[3000];

  correl2 Correl;

  //void treeGrow();
  //void loadTree(Event * event);


  int CsImult;
  int count = 0;

  int S800ID;
  int SiID;
  int JanusID;
  int CAESARID;

  bool LoadS800toSolution();
  
  
 private:
  doppler *Doppler;
  histo_sort * Histo_sort;
  histo_read * Histo_read;
  

};
#endif
