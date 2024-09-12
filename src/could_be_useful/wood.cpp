#include "wood.h"
wood::wood(int N0,bool gamma0, string* name)
{
  N = N0;
  gamma = gamma0;
  file = new TFile (name->c_str(),"RECREATE");
  T = new TTree("T","T");

  T->Branch("id1",&id1,"id1/I");
  T->Branch("id2",&id2,"id2/I");
  if (N >= 3)  T->Branch("id3",&id3,"id3/I");
  if (N >= 4)  T->Branch("id4",&id4,"id4/I");
  if (N >= 5)  T->Branch("id5",&id5,"id5/I");
  if (N >= 6)  T->Branch("id6",&id6,"id6/I");

  T->Branch("ifront1",&ifront1,"ifront1/I");
  T->Branch("ifront2",&ifront2,"ifront2/I");
  if (N >= 3)  T->Branch("ifront3",&ifront3,"ifront3/I");
  if (N >= 4)  T->Branch("ifront4",&ifront4,"ifront4/I");
  if (N >= 5)  T->Branch("ifront5",&ifront5,"ifront5/I");
  if (N >= 6)  T->Branch("ifront6",&ifront6,"ifront6/I");

  T->Branch("iback1",&iback1,"iback1/I");
  T->Branch("iback2",&iback2,"iback2/I");
  if (N >= 3)  T->Branch("iback3",&iback3,"iback3/I");
  if (N >= 4)  T->Branch("iback4",&iback4,"iback4/I");
  if (N >= 5)  T->Branch("iback5",&iback5,"iback5/I");
  if (N >= 6)  T->Branch("iback6",&iback6,"iback6/I");

  T->Branch("M1",M1,"M1[3]/F");
  T->Branch("M2",M2,"M2[3]/F");
  if (N >= 3)  T->Branch("M3",M3,"M3[3]/F");
  if (N >= 4)  T->Branch("M4",M4,"M4[3]/F");
  if (N >= 5)  T->Branch("M5",M5,"M5[3]/F");
  if (N >= 6)  T->Branch("M6",M6,"M6[3]/F");

  T->Branch("et1",&et1,"et1/F");
  T->Branch("et2",&et2,"et2/F");
  if (N>=3)  T->Branch("et3",&et3,"et3/F");
  if (N>=4)  T->Branch("et4",&et4,"et4/F");
  if (N>=5)  T->Branch("et5",&et5,"et5/F");
  if (N>=6)  T->Branch("et6",&et6,"et6/F");

  T->Branch("Ex",&Ex,"Ex/F");
  T->Branch("Vcm",&Vcm,"Vcm/F");
  T->Branch("thetaCM",&thetaCM,"thetaCM/F");

  if (gamma)
    {
     T->Branch("Ngamma",&Ngamma,"Ngamma/I");
     T->Branch("Egamma",Egamma,"Egamma/F");
    }

  T->Branch("energy_p1",&energy_p1,"energy_p1/F");
  T->Branch("energy_p2",&energy_p2,"energy_p2/F");
  if (N>=3)   T->Branch("energy_p3",&energy_p3,"energy_p3/F");
  if (N>=4)   T->Branch("energy_p4",&energy_p4,"energy_p4/F");
  if (N>=5)   T->Branch("energy_p5",&energy_p5,"energy_p5/F");
  if (N>=6)   T->Branch("energy_p6",&energy_p6,"energy_p6/F");


  T->Branch("denergy1",&denergy1,"denergy1/F");
  T->Branch("denergy2",&denergy2,"denergy2/F");
  if (N>=3)   T->Branch("denergy3",&denergy3,"denergy3/F");
  if (N>=4)   T->Branch("denergy4",&denergy4,"denergy4/F");
  if (N>=5)   T->Branch("denergy5",&denergy5,"denergy5/F");
  if (N>=6)   T->Branch("denergy6",&denergy6,"denergy6/F");



}

wood::~wood()
{
  file->Write();
  delete T;
  delete file;
}
