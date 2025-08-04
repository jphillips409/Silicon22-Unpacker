#include "wood_gammas.h"
wood_gammas::wood_gammas(string  name)
{

  //gROOT->ProcessLine("#include <vector>")

  file = new TFile (name.c_str(),"RECREATE");
  T = new TTree("T","T");


  T->Branch("Ngamma",&Ngamma);
  T->Branch("Ngamma_Select",&Ngamma_Select);
  T->Branch("Egamma",&Egamma);
  T->Branch("Egamma_Select",&Egamma_Select);
  T->Branch("Thetagamma",&Thetagamma);
  T->Branch("Thetagamma_Select",&Thetagamma_Select);
  T->Branch("Phigamma",&Phigamma);
  T->Branch("Phigamma_Select",&Phigamma_Select);
  T->Branch("Tgamma",&Tgamma);
  T->Branch("Tgamma_Select",&Tgamma_Select);
  T->Branch("Chgamma",&Chgamma);
  T->Branch("Chgamma_Select",&Chgamma_Select);

  T->Branch("Zres",&Zres,"Zres/I");
  T->Branch("Ares",&Ares,"Ares/I");
  T->Branch("theta_res",&theta_res,"theta_res/F");
  T->Branch("phi_res",&phi_res,"phi_res/F");
  T->Branch("ekin_res",&ekin_res,"ekin_res/F");

  T->Branch("beamZ",&beamZ,"beamZ/I");

}

//Reset everything or else the vectors go crazy
void wood_gammas::reset()
{

  Ngamma = -1;
  Ngamma_Select = -1;
  Egamma.clear();
  Egamma_Select.clear();
  Thetagamma.clear();
  Thetagamma_Select.clear();
  Phigamma.clear();
  Phigamma_Select.clear();
  Tgamma.clear();
  Tgamma_Select.clear();
  Chgamma.clear();
  Chgamma_Select.clear();

  Zres = -1;
  Ares = -1;
  theta_res = -1;
  phi_res = -1;
  ekin_res = -1;

  beamZ = -1;

}

wood_gammas::~wood_gammas()
{
  file->Write();
  delete T;
  delete file;
}
