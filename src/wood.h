#include "TFile.h"
#include "TTree.h"
#include <string>
using namespace std;
class wood
{
 public:
  wood(int N,bool gamma,string name);
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

  int itele1;
  int itele2;
  int itele3;
  int itele4;
  int itele5;
  int itele6;

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
  
  float Erel;
  float Ex;
  float Vcm;
  float thetaCM;
  float cos_thetaH;

  int Ngamma;
  int Ngamma_Select;
  float Egamma[15];
  float Egamma_Select[15];
  float Tgamma[15];
  float Tgamma_Select[15];
  int Chgamma[15];
  int Chgamma_Select[15];

  float energy_p1;
  float energy_p2;
  float energy_p3;
  float energy_p4;
  float energy_p5;
  float energy_p6;

  float energy_p1_R;
  float energy_p2_R;
  float energy_p3_R;
  float energy_p4_R;
  float energy_p5_R;
  float energy_p6_R;

  float denergy1_R;
  float denergy2_R;
  float denergy3_R;
  float denergy4_R;
  float denergy5_R;
  float denergy6_R;

  float time1;
  float time2;
  float time3;
  float time4;
  float time5;
  float time6;

  float theta1;
  float theta2;
  float theta3;
  float theta4;
  float theta5;
  float theta6;

  float phi1;
  float phi2;
  float phi3;
  float phi4;
  float phi5;
  float phi6;

  int runnum;
  int beamZ;

  float theta_s800;
  float phi_s800;


};

