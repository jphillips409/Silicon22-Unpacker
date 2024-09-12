#include "two.h"
two::two()
{
  file = new TFile ("two.root","RECREATE");
  T = new TTree("T","T");

  T->Branch("run",&run,"run/I");
  T->Branch("idCore",&idCore,"idCore/I");
  T->Branch("ida",&ida,"ida/I");
  T->Branch("ifrontCore",&ifrontCore,"ifrontCore/I");
  T->Branch("ibackCore",&ibackCore,"ibackCore/I");
  T->Branch("ifronta",&ifronta,"ifronta/I");
  T->Branch("ibacka",&ibacka,"ibacka/I");

  T->Branch("MMa",&MMa,"MMa/F");
  T->Branch("MMCore",&MMCore,"MMCore/F");
  T->Branch("Ma",Ma,"Ma[3]/F");
  T->Branch("MCore",MCore,"MCore[3]/F");
  T->Branch("eta",&eta,"eta/F");
  T->Branch("etCore",&etCore,"etCore/F");

//   T->Branch("denergy",&denergy,"denergy/F");  
//   T->Branch("energyR",&energyR,"energyR/F");  
//   T->Branch("energy",&energy,"energy/F");  
//   T->Branch("theta",&theta,"theta/F");  
//   T->Branch("phi",&phi,"phi/F");  

}

two::~two()
{
  file->Write();
  delete T;
  delete file;
}
