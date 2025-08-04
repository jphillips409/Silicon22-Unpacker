#include "TFile.h"
#include "TTree.h"
#include <string>
#include <vector>
using namespace std;
class wood_gammas
{
 public:
  wood_gammas(string name);
  ~wood_gammas();
  void reset();
  TFile * file;
  TTree * T;

  int Ngamma;
  int Ngamma_Select;
  vector<float> Egamma;
  vector<float> Egamma_Select;
  vector<float> Thetagamma;
  vector<float> Thetagamma_Select;
  vector<float> Phigamma;
  vector<float> Phigamma_Select;
  vector<float> Tgamma;
  vector<float> Tgamma_Select;
  vector<int> Chgamma;
  vector<int> Chgamma_Select;

  //S800 only gives one heavy fragment
  int Zres;
  int Ares;
  float theta_res;
  float phi_res;
  float ekin_res;

  int beamZ;

};

