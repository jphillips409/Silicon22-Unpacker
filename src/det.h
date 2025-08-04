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
#include "TTree.h"
#include "wood.h"
#include "corrcomb.h"
#include "wood_gammas.h"

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
  CLosses *losses_fiberAl;

  corrcomb *Corrcomb;
  
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

  //proton coincidence with S800 beam
  int N_coin_Si23 = 0;
  int N_coin_Si23_1p = 0;
  int N_coin_Si23_Si23_1p = 0;
  int N_coin_Si23_2p = 0;
  int N_coin_Si23_3p = 0;
  int N_coin_Si22_1p = 0;


  int N_coin_Si23_Mg20 = 0;
  int N_coin_Si23_Mg20_1p = 0;
  int N_coin_Si23_Mg20_2p = 0;
  bool Si23_Mg20_2p_flag = false;
  bool Si23_Si22_p_flag = false;
  bool Si23_Si23_p_flag = false;

  int N_coin_Ne17_3p = 0;
  bool Ne17_3p_flag = false;
  int N_coin_Ne18_3p = 0;
  bool Ne18_3p_flag = false;

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

  int s800coinc_cnt;

  int Nmaxmix;
  int EventMixCounter = 0;
  solution * EventMixerP[100];
  solution * EventMixer33Ar[100];
  int EventSingleMixCounter = 0;
  int NmaxmixSingles;
  solution * EventMixerSingleP[3000];

  correl2 Correl;
  void corr_1H();
  void corr_12N();
  void corr_13N();
  void corr_14O();
  void corr_15O();
  void corr_15F();
  void corr_16F();
  void corr_17F();
  void corr_16Ne();
  void corr_17Ne();
  void corr_18Ne();
  void corr_19Ne();
  void corr_19Na();
  void corr_20Na();
  void corr_21Na();
  void corr_19Mg();
  void corr_20Mg();
  void corr_21Mg();
  void corr_22Mg();
  void corr_20Al();
  void corr_21Al();
  void corr_22Al();
  void corr_23Al();
  void corr_22Si();
  void corr_23Si();
  void corr_24Si();
  void corr_23P();
  void corr_3p18Ne();
  void corr_24P();

  //alphas
  void corr_a15O();

  //Initialize gamma-ray tree
  wood_gammas * tgammas;

  //Initialize corr trees
  wood * p11C;
  wood * p12C;
  wood * p13N;
  wood * p14N;
  wood * p14O;
  wood * pp14O;
  wood * p15O;
  wood * pp15O;
  wood * p16O;
  wood * pp16O;
  wood * p18F;
  wood * pp18F;
  wood * p17F;
  wood * pp17F;
  wood * p18Ne;
  wood * p19Ne;
  wood * p20Ne;
  wood * pp17Ne;
  wood * pp18Ne;
  wood * p20Na;
  wood * pp19Ne;
  wood * p21Na;
  wood * ppp17Ne;
  wood * p20Mg;
  wood * ppp18Ne;
  wood * p21Mg;
  wood * p22Mg;
  wood * pp20Mg;
  wood * p22Al;
  wood * pp21Mg;
  wood * p23Al;
  wood * ppp20Mg;
  wood * p22Si;
  wood * p23Si;
  wood * pppp18Ne;

  //alphas
  wood * a15O;

  int corr_23P_counter = 0;

  //void treeGrow();
  //void loadTree(Event * event);


  int CsImult;
  int count = 0;

  int S800ID;
  int SiID;
  int JanusID;
  int CAESARID;

  int runnum;

  bool LoadS800toSolution();

  //Heavy frag unit vector
  double xH;
  double yH;
  double zH;
  double EH; //Euclidian vector length
  
  
 private:
  doppler *Doppler;
  histo_sort * Histo_sort;
  histo_read * Histo_read;
  

};
#endif
