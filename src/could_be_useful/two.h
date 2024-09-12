#include "TFile.h"
#include "TTree.h"

class two
{
 public:
  two();
  ~two();
  TFile * file;
  TTree * T;

  int run;
  int ida;
  int idCore;
  int ifrontCore;
  int ibackCore;
  int ifronta;
  int ibacka;

  float energyR;
  float energy;
  float denergy;
  float theta;
  float phi;

  float MMa;
  float MMCore;
  float Ma[3];
  float MCore[3];
  float eta;
  float etCore;

  


};
