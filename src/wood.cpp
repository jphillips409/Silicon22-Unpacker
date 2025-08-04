#include "wood.h"
wood::wood(int N0,bool gamma0, string  name)
{
  N = N0;
  gamma = gamma0;
  file = new TFile (name.c_str(),"RECREATE");
  T = new TTree("T","T");

  T->Branch("id1",&id1,"id1/I");
  if (N >= 2 ) T->Branch("id2",&id2,"id2/I");
  if (N >= 3)  T->Branch("id3",&id3,"id3/I");
  if (N >= 4)  T->Branch("id4",&id4,"id4/I");
  if (N >= 5)  T->Branch("id5",&id5,"id5/I");
  if (N >= 6)  T->Branch("id6",&id6,"id6/I");

  T->Branch("itele1",&itele1,"itele1/I");
  if (N >= 2)  T->Branch("itele2",&itele2,"itele/I");
  if (N >= 3)  T->Branch("itele3",&itele3,"itele/I");
  if (N >= 4)  T->Branch("itele4",&itele4,"itele/I");
  if (N >= 5)  T->Branch("itele5",&itele5,"itele/I");
  if (N >= 6)  T->Branch("itele6",&itele6,"itele/I");

  T->Branch("ifront1",&ifront1,"ifront1/I");
  if (N >= 2)  T->Branch("ifront2",&ifront2,"ifront2/I");
  if (N >= 3)  T->Branch("ifront3",&ifront3,"ifront3/I");
  if (N >= 4)  T->Branch("ifront4",&ifront4,"ifront4/I");
  if (N >= 5)  T->Branch("ifront5",&ifront5,"ifront5/I");
  if (N >= 6)  T->Branch("ifront6",&ifront6,"ifront6/I");

  T->Branch("iback1",&iback1,"iback1/I");
  if (N >= 2)  T->Branch("iback2",&iback2,"iback2/I");
  if (N >= 3)  T->Branch("iback3",&iback3,"iback3/I");
  if (N >= 4)  T->Branch("iback4",&iback4,"iback4/I");
  if (N >= 5)  T->Branch("iback5",&iback5,"iback5/I");
  if (N >= 6)  T->Branch("iback6",&iback6,"iback6/I");

  T->Branch("M1",M1,"M1[3]/F");
  if (N >= 2)T->Branch("M2",M2,"M2[3]/F");
  if (N >= 3)  T->Branch("M3",M3,"M3[3]/F");
  if (N >= 4)  T->Branch("M4",M4,"M4[3]/F");
  if (N >= 5)  T->Branch("M5",M5,"M5[3]/F");
  if (N >= 6)  T->Branch("M6",M6,"M6[3]/F");

  T->Branch("et1",&et1,"et1/F");
  if (N>=2)  T->Branch("et2",&et2,"et2/F");
  if (N>=3)  T->Branch("et3",&et3,"et3/F");
  if (N>=4)  T->Branch("et4",&et4,"et4/F");
  if (N>=5)  T->Branch("et5",&et5,"et5/F");
  if (N>=6)  T->Branch("et6",&et6,"et6/F");


  T->Branch("time1",&time1,"time1/F");
  if (N>=2)  T->Branch("time2",&time2,"time2/F");
  if (N>=3)  T->Branch("time3",&time3,"time3/F");
  if (N>=4)  T->Branch("time4",&time4,"time4/F");
  if (N>=5)  T->Branch("time5",&time5,"time5/F");
  if (N>=6)  T->Branch("time6",&time6,"time6/F");


  T->Branch("Erel",&Erel,"Erel/F");
  T->Branch("Ex",&Ex,"Ex/F");
  T->Branch("Vcm",&Vcm,"Vcm/F");
  T->Branch("thetaCM",&thetaCM,"thetaCM/F");
  T->Branch("cos_thetaH",&cos_thetaH,"cos_thetaH/F");

  if (gamma)
    {
     T->Branch("Ngamma",&Ngamma,"Ngamma/I");
     T->Branch("Ngamma_Select",&Ngamma_Select,"Ngamma_Select/I");
     T->Branch("Egamma",Egamma,"Egamma[15]/F");
     T->Branch("Egamma_Select",Egamma_Select,"Egamma_Select[15]/F");
     T->Branch("Tgamma",Tgamma,"Tgamma[15]/F");
     T->Branch("Tgamma_Select",Tgamma_Select,"Tgamma_Select[15]/F");
     T->Branch("Chgamma",Chgamma,"Chgamma[15]/I");
     T->Branch("Chgamma_Select",Chgamma_Select,"Chgamma_Select[15]/I");
    }

  T->Branch("energy_p1",&energy_p1,"energy_p1/F");
  if (N>=2)  T->Branch("energy_p2",&energy_p2,"energy_p2/F");
  if (N>=3)  T->Branch("energy_p3",&energy_p3,"energy_p3/F");
  if (N>=4)  T->Branch("energy_p4",&energy_p4,"energy_p4/F");
  if (N>=5)  T->Branch("energy_p5",&energy_p5,"energy_p5/F");
  if (N>=6)  T->Branch("energy_p6",&energy_p6,"energy_p6/F");

  T->Branch("energy_p1_R",&energy_p1_R,"energy_p1_R/F");
  if (N>=2)  T->Branch("energy_p2_R",&energy_p2_R,"energy_p2_R/F");
  if (N>=3)  T->Branch("energy_p3_R",&energy_p3_R,"energy_p3_R/F");
  if (N>=4)  T->Branch("energy_p4_R",&energy_p4_R,"energy_p4_R/F");
  if (N>=5)  T->Branch("energy_p5_R",&energy_p5_R,"energy_p5_R/F");
  if (N>=6)  T->Branch("energy_p6_R",&energy_p6_R,"energy_p6_R/F");

  T->Branch("denergy1_R",&denergy1_R,"denergy1_R/F");
  if (N>=2)  T->Branch("denergy2_R",&denergy2_R,"denergy2_R/F");
  if (N>=3)  T->Branch("denergy3_R",&denergy3_R,"denergy3_R/F");
  if (N>=4)  T->Branch("denergy4_R",&denergy4_R,"denergy4_R/F");
  if (N>=5)  T->Branch("denergy5_R",&denergy5_R,"denergy5_R/F");
  if (N>=6)  T->Branch("denergy6_R",&denergy6_R,"denergy6_R/F");

  T->Branch("theta1",&theta1,"theta1/F");
  if (N>=2)  T->Branch("theta2",&theta2,"theta2/F");
  if (N>=3)  T->Branch("theta3",&theta3,"theta3/F");
  if (N>=4)  T->Branch("theta4",&theta4,"theta4/F");
  if (N>=5)  T->Branch("theta5",&theta5,"theta5/F");
  if (N>=6)  T->Branch("theta6",&theta6,"theta6/F");

  T->Branch("phi1",&phi1,"phi1/F");
  if (N>=2)  T->Branch("phi2",&phi2,"phi2/F");
  if (N>=3)  T->Branch("phi3",&phi3,"phi3/F");
  if (N>=4)  T->Branch("phi4",&phi4,"phi4/F");
  if (N>=5)  T->Branch("phi5",&phi5,"phi5/F");
  if (N>=6)  T->Branch("phi6",&phi6,"phi6/F");

  T->Branch("runnum",&runnum,"runnum/I");
  T->Branch("beamZ",&beamZ,"beamZ/I");

  T->Branch("theta_s800",&theta_s800,"theta_s800/F");
  T->Branch("phi_s800",&phi_s800,"phi_s800/F");

}

wood::~wood()
{
  file->Write();
  delete T;
  delete file;
}
