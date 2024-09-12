#include "histo_read.h"
#include "histOn.h"

histo_read::histo_read(string suffix)
{

  ostringstream outstring;
  string name;

  Ntele = 14;
  Nstrip = 32;
  NCsI = 56;

  //create root file
  string filename = "read"+suffix+".root";
  cout << "writing histo_read to: " << filename << endl;
  file_read = new TFile (filename.c_str(),"RECREATE");
  file_read->cd();


  //directories for Correl
  dirCorr = new TDirectoryFile("corr","corr");

  
  dirGamma = new TDirectoryFile("gamma", "gamma");
  dirCRDC_proj = new TDirectoryFile("CRDC_proj", "CRDC_proj");


  //delta E vs E
  dirdEE = new TDirectoryFile("dEE","dEE");

  dirProjections = new TDirectoryFile("projections","projections");
  dirProjectionsRaw = new TDirectoryFile("projections_raw","projections_raw");
 
  //directory for Front back mixed spectra
  dirFB = new TDirectoryFile("fb","fb");

  //directory for CsI Spectra
  dirCsI = new TDirectoryFile("CsI","CsI");
  dirCsIGate = dirCsI->mkdir("CsIGate","CsIGate");
  dirSi = new TDirectoryFile("Si","Si");
  dirELoss = dirSi->mkdir("ELoss","ELoss");


  CorrelationTable = new TH2I("CorrelationTable","",6,-0.5,5.5,24,0,25.5);


  if (On_dEE) dEE = new TH2I*[NCsI];

  if (On_dEE_pro) dEE_Projections = new TH1I**[Ntele];
  if (On_dEE_pro) dEE_ProjectionsRaw = new TH1I**[Ntele];


  Light = new TH1I*[NCsI];
  if (On_Eloss) ELoss = new TH1I*[NCsI];

  for(int icsi =0;icsi <NCsI;icsi++)
  {
    outstring.str("");
    outstring << "dEE_" << icsi;
    name = outstring.str();
    dirdEE->cd();
    if (On_dEE)
	  {
      dEE[icsi] = new TH2I(name.c_str(),"",1024,0,4095,5000,0,500);//4096
      dEE[icsi]->GetXaxis()->SetTitle("CsI Light [ADC Channel number]");
      dEE[icsi]->GetYaxis()->SetTitle("Si Energy [MeV]");
  	}
    outstring.str("");
    outstring << "Light_" << icsi;
    name = outstring.str();
    dirCsIGate->cd();
    //Light[icsi] = new TH1I(name.c_str(),"",4096,0,4096);
    Light[icsi] = new TH1I(name.c_str(),"",1500,0,600);
    Light[icsi]->GetXaxis()->SetTitle("Csi Light [arb]");
    Light[icsi]->GetXaxis()->CenterTitle();
    
    Light[icsi]->GetYaxis()->SetTitle("Counts");
    Light[icsi]->GetYaxis()->CenterTitle();

    outstring.str("");
    outstring << "ELoss_"<< icsi;
    name = outstring.str();
    dirELoss->cd();
    if (On_Eloss) ELoss[icsi] = new TH1I(name.c_str(),"",5000,0,500);
  }
  dirdEE->cd();
  HitMap = new TH2I("HitMap","",1200,-30,30,800,-20,20);

  
  dirProjections->cd();
  for(int itele=0; itele<Ntele; itele++)
    {
      if (On_dEE_pro)dEE_Projections[itele] = new TH1I*[Nstrip];
      for(int istrip=0; istrip<Nstrip; istrip++)
	{
	  outstring.str("");
	  outstring << "tele_" << itele << "_strip_" << istrip;
	  name = outstring.str();
	  if (On_dEE_pro)dEE_Projections[itele][istrip] = new TH1I(name.c_str(),"",2000,0,450);
	}
    }
  
  dirProjectionsRaw->cd();
  for(int itele=0; itele<Ntele; itele++)
    {
      if (On_dEE_pro)dEE_ProjectionsRaw[itele] = new TH1I*[Nstrip];
      for(int istrip=0; istrip<Nstrip; istrip++)
	{
	  outstring.str("");
	  outstring << "tele_" << itele << "_strip_" << istrip;
	  name = outstring.str();
	  if (On_dEE_pro)dEE_ProjectionsRaw[itele][istrip] = new TH1I(name.c_str(),"",4000,0,16000);
	}
    }

  dirFB->cd();

  if (On_FB)
    {
     FBDiff = new TH1I*[Ntele];
     FBDiffLG = new TH1I*[Ntele];
     FB = new TH2I*[Ntele];
     FBMult = new TH2I*[Ntele];



     for (int itele = 0;itele<Ntele;itele++)
       {
         outstring.str("");
         outstring << "FBMult_" << itele;
         name = outstring.str();
         dirFB->cd();
         FBMult[itele] = new TH2I(name.c_str(),"",10,-0.5,9.5,10,-0.5,9.5);
         FBMult[itele]->GetXaxis()->SetTitle("Front Multiplicity");
         FBMult[itele]->GetYaxis()->SetTitle("Back Multiplicity");

         outstring.str("");
         outstring << "FBDiff_" << itele;
         name = outstring.str();
         dirFB->cd();
         FBDiff[itele] = new TH1I(name.c_str(),"",100,0,20);
         FBDiff[itele]->GetXaxis()->SetTitle("Front - Back");

         outstring.str("");
         outstring << "FB_" << itele;
         name = outstring.str();
         dirFB->cd();

         FB[itele] = new TH2I(name.c_str(),"",4096,0,500,4096,0,500);
         FB[itele]->GetXaxis()->SetTitle("Front");
         FB[itele]->GetYaxis()->SetTitle("Back");



         outstring.str("");
         outstring << "FBDiffLG_" << itele;
         name = outstring.str();
         dirFB->cd();
         FBDiffLG[itele] = new TH1I(name.c_str(),"",100,0,20);
         FBDiffLG[itele]->GetXaxis()->SetTitle("Front - Back");
          

       }
  
    }

  int Nbins = 900;


  string correlation;



  Rigidity_protongated2plus = new TH1I("Rigidity_protongated2plus","",100,1.85,2.25);
  Rigidity_gammagated2plus = new TH1I("Rigidity_gammagated2plus","",100,1.85,2.25);
  Rigidity_gammagated2plus_strict = new TH1I("Rigidity_gammagated2plus_strict","",100,1.85,2.25);
  cosbeamCMtoHF_Ex_p35K = new TH2I("cosbeamCMtoHF_Ex_p35K","",100,-2,15,100,-1,1);
  cosbeamCMtoHF_Erel_p35K = new TH2I("cosbeamCMtoHF_Erel_p35K","",120,0,12,100,-1,1);

  Energy_protongated2plus = new TH1I("Energy_protongated2plus","",100,900,1100);
  EnergyKin_protongated2plus = new TH1I("EnergyKin_protongated2plus","",100,50,100);
  costheta_protongated2plus = new TH1I("costheta_protongated2plus","",100,-1.2,1.2);

  Ex_ineachCsI = new TH2I("Ex_ineachCsI","",20,0,20,100,0,8);


  correlation = "_36Ca_2p34Ar";
  name = "Erel" + correlation;
  Erel_36Ca_2p34Ar = new TH1I(name.c_str(),"",Nbins,0,18);
  name = "Ex" + correlation;
  Ex_36Ca_2p34Ar = new TH1I(name.c_str(),"",Nbins,-2,16);
  name = "ThetaCM" + correlation;
  ThetaCM_36Ca_2p34Ar = new TH1I(name.c_str(),"",50,0,2);
  name = "VCM" + correlation;
  VCM_36Ca_2p34Ar = new TH1I(name.c_str(),"",200,7,14);



  
}
//*********************************************
void histo_read::write()
{
  file_read->Write();
  cout << "histo written" << endl;
  file_read->Close();
  /*
    for (int i=0;i<Ntele;i++)
    {
    delete red[i];
    }
    delete [] red;
  */
}
