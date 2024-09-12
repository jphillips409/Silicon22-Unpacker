#include "TFile.h"
#include "TTree.h"
#include <string>
using namespace std;
class wood
{
 public:
  wood(int N,bool gamma,string* name);
  ~wood();
  TFile * file;
  TTree * T;

  int N;
  bool gamma;

  int id1;
  int id2;
  int id3;
  int id4;
  int id5;
  int id6;

  float M1[3];
  float M2[3];
  float M3[3];
  float M4[3];
  float M5[3];
  float M6[3];

  float et1;
  float et2;
  float et3;
  float et4;
  float et5;
  float et6;

  int ifront1;
  int ifront2;
  int ifront3;
  int ifront4;
  int ifront5;
  int ifront6;

  int iback1;
  int iback2;
  int iback3;
  int iback4;
  int iback5;
  int iback6;

  


  float Ex;
  float Vcm;
  float thetaCM;

  int Ngamma;
  float Egamma[10];


  float energy_p1;
  float energy_p2;
  float energy_p3;
  float energy_p4;
  float energy_p5;
  float energy_p6;


  float denergy1;
  float denergy2;
  float denergy3;
  float denergy4;
  float denergy5;
  float denergy6;


};

