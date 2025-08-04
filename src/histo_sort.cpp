#include "histo_sort.h"


histo_sort::histo_sort(string suffix)
{
  ostringstream outstring;
  string filename = "sort"+suffix+".root";
  cout << "writing histo_sort to: " << filename << endl;
  //create root file
  file = new TFile (filename.c_str(), "RECREATE");

  //tree = new TTree("scalars", "run scalars");
  
  file->cd();

  //**************************************************************************************
  //Directories
  //**************************************************************************************
    
  //Gobbi directories
  dirSummary = new TDirectoryFile("GobbiSum", "GobbiSum"); //name, title
  //make subdirectories
  //Energies and DeltaEnergies, Raw+Calibrated
  dir1dFrontE_R = dirSummary->mkdir("1dFrontE_R","1dFrontE_R");
  dir1dBackE_R = dirSummary->mkdir("1dBackE_R","1dBackE_R");

  dir1dFrontlowE_R = dirSummary->mkdir("1dFrontlowE_R","1dFrontlowE_R");
  dir1dBacklowE_R = dirSummary->mkdir("1dBacklowE_R","1dBacklowE_R");

  dir1dFrontE_cal = dirSummary->mkdir("1dFrontE_cal","1dFrontE_cal");
  dir1dBackE_cal = dirSummary->mkdir("1dBackE_cal","1dBackE_cal");

  dir1dFrontlowE_cal = dirSummary->mkdir("1dFrontlowE_cal","1dFrontlowE_cal");
  dir1dBacklowE_cal = dirSummary->mkdir("1dBacklowE_cal","1dBacklowE_cal");

  //raw time summaries
  dir1dFrontTime_R = dirSummary->mkdir("1dFrontTime_R","1dFrontTime_R");
  dir1dBackTime_R = dirSummary->mkdir("1dBackTime_R","1dBackTime_R");

  dir1dCsI_Energy = dirSummary->mkdir("1dCsI_Energy","1dCsI_Energy");
  dir1dCsI_Time = dirSummary->mkdir("1dCsI_Time","1dCsI_Time");
  dir1dCsI_QDC = dirSummary->mkdir("1dCsI_QDC","1dCsI_QDC");

  //directory for DeltaE-E plots
  dirDEEplots = new TDirectoryFile("DEEplots","DEEplots");
  dirPSD = dirDEEplots->mkdir("PSDplots","PSDplots");
  dirhitmaps = new TDirectoryFile("hitmaps","hitmaps");

  //fiber
  dirFiber = new TDirectoryFile("Fiber","Fiber");
  dir1dBoard0TOA_R = dirFiber->mkdir("Board0TOA_R","Board0TOA_R");
  dir1dBoard0TOT_R = dirFiber->mkdir("Board0TOT_R","Board0TOT_R");
  dir1dBoard1TOA_R = dirFiber->mkdir("Board1TOA_R","Board1TOA_R");
  dir1dBoard1TOT_R = dirFiber->mkdir("Board1TOT_R","Board0TOT_R");
  dirHitMap = dirFiber->mkdir("FiberHitMap","FiberHitMap");

  //Inv mass directories
  dirInvMass = new TDirectoryFile("InvMass", "InvMass");
  dirCorrCombs = dirInvMass->mkdir("CorrCombs","CorrCombs");
  dir1H = dirInvMass->mkdir("1H","1H");
  dir16F = dirInvMass->mkdir("16F","16F");
  dir17F = dirInvMass->mkdir("17F","17F");
  dir18Ne = dirInvMass->mkdir("18Ne","18Ne");
  dir19Ne = dirInvMass->mkdir("19Ne","19Ne");
  dir19Na = dirInvMass->mkdir("19Na","19Na");
  dir20Na = dirInvMass->mkdir("20Na","20Na");
  dir21Na = dirInvMass->mkdir("21Na","21Na");
  dir19Mg = dirInvMass->mkdir("19Mg","19Mg");
  dir20Mg = dirInvMass->mkdir("20Mg","20Mg");
  dir21Mg = dirInvMass->mkdir("21Mg","21Mg");
  dir2p19Ne = dir21Mg->mkdir("2p19Ne","2p19Ne");
  dir22Mg = dirInvMass->mkdir("22Mg","22Mg");
  dir21Al = dirInvMass->mkdir("21Al","21Al");
  dir22Al = dirInvMass->mkdir("22Al","22Al");
  dir23Al = dirInvMass->mkdir("23Al","23Al");
  dir3p18Ne = dir21Al->mkdir("3p18Ne","3p18Ne");
  dir22Si = dirInvMass->mkdir("22Si","22Si");
  dir23Si = dirInvMass->mkdir("23Si","23Si");
  dirp22Al = dir23Si->mkdir("p22Al","p22Al");
  dir2p21Mg = dir23Si->mkdir("2p21Mg","2p21Mg");
  dir24Si = dirInvMass->mkdir("24Si","24Si");
  dir23P = dirInvMass->mkdir("23P","23P");
  dir3p20Mg = dir23P->mkdir("3p20Mg","3p20Mg");
  dirp22Si = dir23P->mkdir("22Si","22Si");
  dir24P = dirInvMass->mkdir("24P","24P");
  dir3p17Ne = dirInvMass->mkdir("20Al","20Al");
  
  //S800 directory
  dirS800 = new TDirectoryFile("S800","S800");
  dirS800Raw = dirS800->mkdir("Raw","Raw");
  dirS800Cal = dirS800->mkdir("Cal","Cal");
  dirS800PID = dirS800->mkdir("PID","PID");
  dirS800Trig = dirS800->mkdir("Trig","Trig");
  dirS800veldist = dirS800->mkdir("veldist","veldist");


  //for calibrated spectra
  int Nbin = 5000;
  float Ecal_Emax = 50.0;
  float CsI_Emax = 200;

  DeltaT_VME_Janus = new TH1D("DeltaT_VME_Janus","",1000,-10000,10000);
  DeltaT_VME_S800 = new TH1D("DeltaT_VME_S800","",1000,-10000,10000);
  DeltaT_S800_Janus = new TH1D("DeltaT_S800_Janus","",1000,-10000,10000);

  //T diff vs event number
  DeltaT_VME_Janus_vsEvtCnt = new TH2D("DeltaT_VME_Janus_vsEvtCnt","",1000,0,1000000,500,-3000,3000);
  DeltaT_VME_S800_vsEvtCnt = new TH2D("DeltaT_VME_S800_vsEvtCnt","",1000,0,1000000,500,-3000,3000);
  DeltaT_S800_Janus_vsEvtCnt = new TH2D("DeltaT_S800_Janus_vsEvtCnt","",1000,0,1000000,500,-3000,3000);

  //Histogram the number of words in S800 events
  S800_NumWords_R = new TH1D("S800_NumWords_R","",5000,0,30000);
  S800_NumWords_Accepted = new TH1D("S800_NumWords_Accepted","",5000,0,30000);
  
  //2D for e1up time and obj time vs gobbi copies
  DB5T_vs_gCopy = new TH2D("DB5T_vs_gCopy","",501,-500,500,500,-4000,4000);
  objT_vs_gCopy = new TH2D("objT_vs_gCopy","",501,-1000,1000,500,-4000,4000);
  RFT_vs_gCopy = new TH2D("RFT_vs_gCopy","",501,-1000,1000,500,-4000,4000);
  S800_RFTime = new TH1I("S800_RFTime","",1000,0,1000);
  S800_RFTime_coin = new TH1I("S800_RFTime_coin","",1000,0,1000);
  S800_RFTime_wbeam = new TH1I("S800_RFTime_wbeam","",1000,0,1000);
  S800_RFTime_Sibeam = new TH1I("S800_RFTime_Sibeam","",1000,0,1000);
  S800_RFTime_Albeam = new TH1I("S800_RFTime_Albeam","",1000,0,1000);
  S800_RFTime_Mgbeam = new TH1I("S800_RFTime_Mgbeam","",1000,0,1000);
  S800_RFTime_Nabeam = new TH1I("S800_RFTime_Nabeam","",1000,0,1000);
  Gobbi_RFTime = new TH1I("Gobbi_RFTime","",2000,-1000,1000);
  
  NCsI = 16;
  Ntele = 4;
  Nstrip = 32; //strips per silicon
  Ncaesar = 192;
  
  //**************************************************************************************
  //Gobbi Summary spectra

  dirSummary->cd();

  
  ostringstream name;

  //All Gobbi summaries

  int Nchan = 4*Nstrip;
  sumFrontE_R = new TH2I("sumFrontE_R","",Nchan,0,Nchan,1024,0,8192);
  sumFrontE_R->SetOption("colz");
  sumBackE_R = new TH2I("sumBackE_R","",Nchan,0,Nchan,1024,0,8192);
  sumBackE_R->SetOption("colz");
  sumFrontlowE_R = new TH2I("sumFrontlowE_R","",Nchan,0,Nchan,1024,0,8192);
  sumBacklowE_R = new TH2I("sumBacklowE_R","",Nchan,0,Nchan,1024,0,8192);

  b1s2Zline = new TH2I("b1s2Zline","",32, -0.5, 31.5,1024,0,8192);
  b1s2Zline->SetOption("colz");

  b1s2FvsB = new TH2I("b1s2FvsB","",1000, 0, 100,1000,0,100);
  b1s2FvsB->SetOption("colz");

  sumFrontE_cal = new TH2I("sumFrontE_cal","",Nchan,0,Nchan,Nbin,0,250);
  sumFrontE_cal->SetOption("colz");
  sumBackE_cal = new TH2I("sumBackE_cal","",Nchan,0,Nchan,Nbin,0,Ecal_Emax);
  sumBackE_cal->SetOption("colz");
  sumFrontlowE_cal = new TH2I("sumFrontlowE_cal","",Nchan,0,Nchan,Nbin,0,250);
  sumFrontlowE_cal_thetacorr = new TH2I("sumFrontlowE_cal_thetacorr","",Nchan,0,Nchan,Nbin,0,250);
  sumBacklowE_cal = new TH2I("sumBacklowE_cal","",Nchan,0,Nchan,Nbin,0,250);
  sumBacklowE_cal_thetacorr = new TH2I("sumBacklowE_cal_thetacorr","",Nchan,0,Nchan,Nbin,0,250);

  sumCsIE_R = new TH2I("sumCsIE_R","",16,0,16,1024,0,4095);
  sumCsIE_R->SetOption("colz");
  sumCsIE_cal = new TH2I("sumCsIE_cal","",16,0,16,500,0,CsI_Emax);
  sumCsIE_cal->SetOption("colz");

  //Gobbi times
  sumFrontTime_R = new TH2I("sumFrontTime_R","",Nchan,0,Nchan,512,0,16383);
  sumFrontTime_R->SetOption("colz");
  sumFrontTime_cal = new TH2I("sumFrontTime_cal","",Nchan,0,Nchan,512,0,16383);
  sumFrontTime_cal->SetOption("colz");
  sumBackTime_R = new TH2I("sumBackTime_R","",Nchan,0,Nchan,512,0,16383);
  sumBackTime_R->SetOption("colz");
  sumBackTime_cal = new TH2I("sumBackTime_cal","",Nchan,0,Nchan,512,0,16383);
  sumBackTime_cal->SetOption("colz");
  sumCsITime_R = new TH2I("sumCsITime_R","",16,0,16,3000,-15000,15000);
  sumCsITime_R->SetOption("colz");
  sumCsITime_cal = new TH2I("sumCsITime_cal","",16,0,16,2000,-1000,1000);
  sumCsITime_cal->SetOption("colz");

  //**************************************************************************************
  //**************************************************************************************
  //**************************************************************************************
  //**************************************************************************************
  //1D spectra for silicons
  //Proton KE for different CsI arrangments
  /*for (int csi=0;csi<4;csi++)
  {
    for (int i=0;i<4;i++)
    {
      dir1dCsI_Linearity->cd();
      name.str("");
      name << "pKE_AllCsI_" << csi << "_" << i;
      CsI_ProtonKE[csi][i] = new TH1I(name.str().c_str(),"",1024,0,4095);

      name.str("");
      name << "CsI_Linearity_" << csi << "_" << i;
      CsI_Linearity[csi][i] = new TH2I(name.str().c_str(),"",1024,0,4095,375,0,75);

      name.str("");
      name << "pKE_Theta_AllCsI_" << csi << "_" << i;
      CsI_ProtonKE_Angle[csi][i] = new TH2I(name.str().c_str(),"",4096,0,4096,180,0,90);
    }
  }*/


  //create all Gobbi 1d spectra
  for (int board_i=0; board_i<Ntele; board_i++)
  {
    for (int chan_i=0; chan_i<Nstrip; chan_i++)
    {
      //individual Front Energy
      dir1dFrontE_R->cd();
      name.str("");
      name << "FrontE_R" << board_i << "_" << chan_i;
      FrontE_R[board_i][chan_i] = new TH1I(name.str().c_str(),"",512,0,2024);

      dir1dFrontlowE_R->cd();
      name.str("");
      name << "FrontElow_R" << board_i << "_" << chan_i;
      FrontElow_R[board_i][chan_i] = new TH1I(name.str().c_str(),"",2048,0,8192);

      dir1dFrontTime_R->cd();
      name.str("");
      name << "FrontTime_R" << board_i << "_" << chan_i;
      FrontTime_R[board_i][chan_i] = new TH1I(name.str().c_str(),"",1024,0,16383);

      dir1dFrontE_cal->cd();
      name.str("");
      name << "FrontE_cal" << board_i << "_" << chan_i;

      double Emax2 = Ecal_Emax;
      if (board_i == 4) Emax2 = 200;

      FrontE_cal[board_i][chan_i] = new TH1I(name.str().c_str(),"",Nbin,0,Emax2);

      dir1dFrontlowE_cal->cd();
      name.str("");
      name << "FrontlowE_cal" << board_i << "_" << chan_i;
      FrontlowE_cal[board_i][chan_i] = new TH1I(name.str().c_str(),"",Nbin,0,Ecal_Emax);


      //individual Back Energy
      dir1dBackE_R->cd();
      name.str("");
      name << "BackE_R" << board_i << "_" << chan_i;
      BackE_R[board_i][chan_i] = new TH1I(name.str().c_str(),"",512,0,2024);

      dir1dBacklowE_R->cd();
      name.str("");
      name << "BackElow_R" << board_i << "_" << chan_i;
      BackElow_R[board_i][chan_i] = new TH1I(name.str().c_str(),"",1024,0,4095);

      dir1dBackTime_R->cd();
      name.str("");
      name << "BackTime_R" << board_i << "_" << chan_i;
      BackTime_R[board_i][chan_i] = new TH1I(name.str().c_str(),"",1024,0,16383);

      dir1dBackE_cal->cd();
      name.str("");
      name << "BackE_cal" << board_i << "_" << chan_i;
      BackE_cal[board_i][chan_i] = new TH1I(name.str().c_str(),"",Nbin,0,Ecal_Emax);

    }
  }

  //create all the CsI 1d spectra
  for (int ichan=0; ichan<NCsI; ichan++)
  {
    dir1dCsI_Energy->cd();
    name.str("");
    name << "CsI_Energy_" << ichan << "R_unmatched";
    CsI_Energy_R_um[ichan] = new TH1I(name.str().c_str(),"",1024,0,4095);

    name.str("");
    name << "CsI_Energy_" << ichan << "R_matched";
    CsI_Energy_R[ichan] = new TH1I(name.str().c_str(),"",1024,0,4095);

    name.str("");
    name << "CsI_Energy_" << ichan << "cal";
    CsI_Energy_cal[ichan] = new TH1I(name.str().c_str(),"",200,0,200);

    name.str("");
    name << "CsI_Energyp_" << ichan << "cal";
    CsI_Energyp_cal[ichan] = new TH1I(name.str().c_str(),"",200,0,200);

    //Take the center of each CsI
    name.str("");
    name << "CsI_Energy_" << ichan << "R_center";
    CsI_Energy_R_center[ichan] = new TH1I(name.str().c_str(),"",1024,0,4095);

    //Proton Beam Energy
    name.str("");
    name << "CsI_Energy_" << ichan << "protonBE";
    CsI_Energy_protonBE[ichan] = new TH1I(name.str().c_str(),"",200,0,200);

    dir1dCsI_Time->cd();
    name.str("");
    name << "CsI_Time_" << ichan << "R_unmatched";
    CsI_Time_R_um[ichan] = new TH1I(name.str().c_str(),"",3000,-15000,15000);

    name.str("");
    name << "CsI_Time_" << ichan << "R_matched";
    CsI_Time_R[ichan] = new TH1I(name.str().c_str(),"",500,-1500,1000);

    name.str("");
    name << "CsI_Time_" << ichan << "cal";
    CsI_Time_cal[ichan] = new TH1I(name.str().c_str(),"",2000,-1000,1000);

  }

  //Create the S800 TDC channels in our VME for cross-checking
  CsI_Time_R_um[16] = new TH1I("S800 DB5 Time","",5000,-15000,15000);
  CsI_Time_R_um[17] = new TH1I("S800 Object Time","",5000,-15000,15000);
  CsI_Time_R_um[18] = new TH1I("S800 RF Time","",5000,-15000,15000);

  CsI_Time_cal[16] = new TH1I("S800 DB5 Time Cal","",2000,-1000,1000);
  CsI_Time_cal[17] = new TH1I("S800 Object Time Cal","",2000,-1000,1000);
  CsI_Time_cal[18] = new TH1I("S800 RF Time Cal","",2000,-10000,10000);
  
  for (int i =0;i<32;i++) //32 qdc channels
  {
    dir1dCsI_QDC->cd();
    name.str("");
    name << "CsI_QDC_" << i << "R";
    CsI_QDC_R[i] = new TH1I(name.str().c_str(),"",1024,0,4095);

    dir1dCsI_QDC->cd();
    name.str("");
    name << "CsI_QDC_" << i << "matched";
    CsI_QDC_matched[i] = new TH1I(name.str().c_str(),"",1024,0,4095);
  }

  //create all spectra based on telescopes
  dirDEEplots->cd();

  for (int tele=0; tele<4; tele++)
  {
    name.str("");
    name << "DEE" << tele;
    //float Emax = 80; //Defines likely max energy loss in gobbi dE-E and bins
    //int Ebins = 500;
    //float dEmax = 22;
    //float dEbins = 800;
    float Emax = 80*2; //Plots expanded to include heavier particles
    int Ebins = 500*2;
    float dEmax = 22*6;
    float dEbins = 800*3;

    for (int icsi=0; icsi<4; icsi++)
    {
  		dirDEEplots->cd();    
      name.str("");
      name << "DEE_CsI" << tele << "_" << icsi;   
      DEE_CsI[tele][icsi] = new TH2I(name.str().c_str(),"",512,0,4096,500,0,50);

      name.str("");
      name << "DEE_CsIS800Coinc" << tele << "_" << icsi;   
      DEE_CsI_S800Coinc[tele][icsi] = new TH2I(name.str().c_str(),"",512,0,4096,500,0,50);

      name.str("");
      name << "DEE_CsIS800Coinc_Si23" << tele << "_" << icsi;   
      DEE_CsI_S800Coinc_Si23[tele][icsi] = new TH2I(name.str().c_str(),"",512,0,4096,500,0,50);

      name.str("");
      name << "DEE_CsI_0deg" << tele << "_" << icsi;   
      DEE_CsI_0deg[tele][icsi] = new TH2I(name.str().c_str(),"",4096,0,4096,400,0,100);

      name.str("");
      name << "DEE_CsI_lowgain" << tele << "_" << icsi;   
      DEE_CsI_lowgain[tele][icsi] = new TH2I(name.str().c_str(),"",4096,0,4096,400,0,100);

      name.str("");
      name << "timediff" << tele << "_" << icsi;    
      timediff_CsI[tele][icsi] = new TH1I(name.str().c_str(),"",1000,-2000,2000);

			dirPSD->cd();
			name.str("");
			name << "PSD_CsI_unmatched" << tele << "_" << icsi;   
      PSD_CsI_unmatched[tele][icsi] = new TH2I(name.str().c_str(),"",4096,0,4096,4096,0,4096);

			dirPSD->cd();
			name.str("");
			name << "PSD_CsI_matched" << tele << "_" << icsi;   
      PSD_CsI_matched[tele][icsi] = new TH2I(name.str().c_str(),"",4096,0,4096,4096,0,4096);
    } 
  }
    

  dirhitmaps->cd();
  testinghitmap = new TH2I("testinghitmap","", 100,-10,10,100,-10,10);
  xyhitmap = new TH2I("xyhitmap","", 100,-10,10,100,-10,10);
  
  protonhitmap = new TH2I("protonhitmap","", 100,-10,10,100,-10,10);
  tphitmap = new TH2I("tphitmap","",1000,-40,40,1000,-40,-40);

  //Front and back strips vs CsI id
  for (int i=0;i<4;i++)
  {
     name.str("");
     name << "GFrontStrip_CsI_unmatched_" << i;
     GFrontStrip_CsI[i] = new TH2I(name.str().c_str(),"",4,-0.5,3.5,32,-0.5,31.5);

     name.str("");
     name << "GBackStrip_CsI_matched_" << i;
     GBackStrip_CsI[i] = new TH2I(name.str().c_str(),"",4,-0.5,3.5,32,-0.5,31.5);
  }

  //Hit map for different types of CsI
  for (int i=0;i<4;i++)
  {
    name.str("");
    name << "CsIHitMap_" << i;
    CsIHitMap[i] = new TH2I(name.str().c_str(),"",100,-10,10,100,-10,10);
  }


  dirdee = new TDirectoryFile("dee","dee");
  dirdee->cd();
  dee = new TH2I*[NCsI];
  dee_S800 = new TH2I*[NCsI];

  //string name = "";

  for(int icsi =0;icsi <NCsI;icsi++)
  {
    name.str("");
    name << "DEE_" << icsi;
    dirdee->cd();

    dee[icsi] = new TH2I(name.str().c_str(),"",1024,0,4095,1024,0,100);//4096


    name.str("");
    name << "DEE_S800_" << icsi;
    dirdee->cd();

    dee_S800[icsi] = new TH2I(name.str().c_str(),"",1024,0,4095,1024,0,100);//4096

  }


  file->cd();


  //fiber
  dirFiber->cd();
	tot_summary_blue = new TH2I("tot_summary_blue", "Blue Fibers ToT Summary (gain-matched)", 64, 0, 64, 512, 0, 512);
	tot_summary_red = new TH2I("tot_summary_red", "Red Fibers ToT Summary (gain-matched)", 64, 0, 64, 512, 0, 512);
	
	tot_summary_blue_matched = new TH2I("tot_summary_blue_matched", "Blue Fibers ToT Summary Matched (gain-matched)", 64, 0, 64, 512, 0, 512);
	tot_summary_red_matched = new TH2I("tot_summary_red_matched", "Red Fibers ToT Summary Matched (gain-matched)", 64, 0, 64, 512, 0, 512);

	// singles histograms
  int fnum = 64;
  //Make 1d plots
  for (int i=0;i<fnum;i++)
  {
    dir1dBoard0TOA_R->cd();
    name.str("");
    name << "Board0_TOA_" << i << "_R";
    B0TOA_R[i] = new TH1I(name.str().c_str(),"",250,0,500);

    dir1dBoard0TOT_R->cd();
    name.str("");
    name << "Board0_TOT_" << i << "_R";
    B0TOT_R[i] = new TH1I(name.str().c_str(),"",2048,0,4095);

    dir1dBoard1TOA_R->cd();
    name.str("");
    name << "Board1_TOA_" << i << "_R";
    B1TOA_R[i] = new TH1I(name.str().c_str(),"",250,0,500);

    dir1dBoard1TOT_R->cd();
    name.str("");
    name << "Board1_TOT_" << i << "_R";
    B1TOT_R[i] = new TH1I(name.str().c_str(),"",2048,0,4095);
  }


  toa_hist = new TH1I("toa_hist", "Time of Arrival", 4096, 0, 4096);

	// matched events histograms
  dirHitMap->cd();
  Fiber_ixiy = new TH2I("Fiber_ixiy","",64,0,64,64,0,64);
  Fiber_xy = new TH2I("Fiber_xy","",64,-16,16,64,-16,16);

	Fiber_totx = new TH1F("Fiber_totx","",64,-16,16);
  Fiber_toty = new TH1F("Fiber_toty","",64,-16,16);
  Fiber_postotx = new TH1F("Fiber_postotx","",64,0,64);
  Fiber_postoty = new TH1F("Fiber_postoty","",64,0,64);
  Fiber_postoax = new TH1F("Fiber_postoax","",2048,0,2048);
  Fiber_postoay = new TH1F("Fiber_postoay","",2048,0,2048);
  Fiber_toax = new TH1I("Fiber_toax","",2048,0,2048);
  Fiber_toay = new TH1I("Fiber_toay","",2048,0,2048);
  

  Fiber_X = new TH1I("Fiber_X","",200,-10,10);
  Fiber_Y = new TH1I("Fiber_Y","",200,-10,10);
  Fiber_XY = new TH2I("Fiber_XY","",64,-8,8,64,-8,8);
  Fiber_Xid = new TH1I("Fiber_Xid","",64,0,64);
  Fiber_Yid = new TH1I("Fiber_Yid","",64,0,64);
  Fiber_XYid = new TH2I("Fiber_XYid","",64,0,64,64,0,64);
  Fiber_XBeam = new TH1I("Fiber_XBeam","",200,-10,10);
  Fiber_YBeam = new TH1I("Fiber_YBeam","",200,-10,10);
  Fiber_XYBeam =new TH2I("Fiber_XYBeam","",64,-8,8,64,-8,8);


  file->cd();


  dirS800Trig->cd();
  singles_trig_time = new TH1I("singles_trig_time","",1024,0,2000);
  S800_Csi_time = new TH1I("S800_CsI_time","",1024,0,2000);
  S800_Csi_time_with_proton = new TH1I("S800_CsI_time_with_proton","",1024,-1000,2000);

  file->cd();
  //CAESAR
 caeDir = file->mkdir("CAESAR");
  caeDir->cd();
  summaryDir = caeDir->mkdir("Summary Plots");
  summaryDir->cd();
  energyiSum = new TH2F("energyiSum","ADC Summary UnCalibrated",Ncaesar, 0, Ncaesar, 2048, 0, 8192);
  energySum = new TH2F("energySum","ADC Summary Calibrated",Ncaesar, 0, Ncaesar, 2048, 0, 8192);
  energyTSum = new TH2F("energyTSum","ADC Summary Calibrated (TDC Channel Multiplicity > 0)",Ncaesar, 0, Ncaesar, 2048, 0, 8192);
  timeSum = new TH2F("timeSum","TDC Summary",Ncaesar,0,Ncaesar,5000,0,5000);
  caeDir->cd();

  sumDir = caeDir->mkdir("Sum Plots");
  sumDir->cd();
  energyTot = new TH1F("energyTot","Total Energy", 2048, 0, 8192);
  energyTot->GetXaxis()->SetTitle("Caesar Energy (keVish)");
  energyTTot = new TH1F("energyTTot","Total Energy (TDC Channel Multiplicity > 0", 2048, 0, 8192);
  energyTTot->GetXaxis()->SetTitle("Caesar Energy (keVish)");
  energyTTotvsT = new TH2F("energyTTotvsT","Total Energy (TDC Channel Multiplicity > 0 vs iTime", 2048, 0, 8192,5000,0,5000);
  energyTTot_tgated = new TH1F("energyTTot_tgated","Total Energy (TDC Channel Multiplicity > 0) with Time Gates", 2048, 0, 8192);
  energyT0Tot = new TH1F("energyT0Tot","Total Energy (TDC Channel Multiplicity = 0", 2048, 0, 8192);
  energyT0Tot->GetXaxis()->SetTitle("Caesar Energy (keVish)");
  energyT1Tot = new TH1F("energyT1Tot","Total Energy (TDC Channel Multiplicity = 1", 2048, 0, 8192);
  energyT1Tot->GetXaxis()->SetTitle("Caesar Energy (keVish)");
  energyT2Tot = new TH1F("energyT2Tot","Total Energy (TDC Channel Multiplicity > 1", 2048, 0, 8192);
  energyT2Tot->GetXaxis()->SetTitle("Caesar Energy (keVish)");
  enAddback = new TH1F("enAddback","Addback Energy",2048, 0, 8192);
  enAddbackvsT = new TH2F("enAddbackvsT","Addback Energy vs T",2048, 0, 8192,5000,0,5000);
  enAddback_tgated = new TH1F("enAddback_tgated","Addback Energy with Time Gates",2048, 0, 8192);

  energyTTot_Mg20 = new TH1F("energyTTot_Mg20","Total Energy Gated on ^{20}Mg (TDC Channel Multiplicity > 0", 2048, 0, 8192);
  energyTTot_Mg20->GetXaxis()->SetTitle("Caesar Energy (keVish)");

  energyTTot_Mg21 = new TH1F("energyTTot_Mg21","Total Energy Gated on ^{21}Mg (TDC Channel Multiplicity > 0", 2048, 0, 8192);
  energyTTot_Mg21->GetXaxis()->SetTitle("Caesar Energy (keVish)");

  enAddback_Si24 = new TH1F("enAddback_Si24","enAddback_Si24",2048, 0, 8192);
  enAddback_Si24vsT = new TH2F("enAddback_Si24vsT","",2048, 0, 8192,5000,0,5000);

  enAddback_Si23 = new TH1F("enAddback_Si23","enAddback_Si23",2048, 0, 8192);
  enAddback_Si23_mult1 = new TH1F("enAddback_Si23_mult1","enAddback_Si23_mult1",2048, 0, 8192);
  enAddback_Si23vsT = new TH2F("enAddback_Si23vsT","",2048, 0, 8192,5000,0,5000);

  enAddback_Al22 = new TH1F("enAddback_Al22","enAddback_Al22",2048, 0, 8192);
  enAddback_Al22_antiTG = new TH1F("enAddback_Al22_antiTG","enAddback_Al22",2048, 0, 8192);
  enAddback_Al22_mult1 = new TH1F("enAddback_Al22_mult1","enAddback_Al22",2048, 0, 8192);
  enAddback_Al22_mult1_antiTG = new TH1F("enAddback_Al22_mult1_antiTG","enAddback_Al22",2048, 0, 8192);
  enAddback_Al22vsT = new TH2F("enAddback_Al22vsT","",2048, 0, 8192,5000,0,5000);
  enAddback_Al22vsT_mult1 = new TH2F("enAddback_Al22vsT_mult1","",2048, 0, 8192,5000,0,5000);
  enAddback_Al22vsCh = new TH2F("enAddback_Al22vsCh","",2048, 0, 8192,192,0,191);
  enSelect_Al22 = new TH1F("enSelect_Al22","enSelect_Al22",2048, 0, 8192);

  caeDir->cd();

  multDir = caeDir->mkdir("Multiplicity");
  multDir->cd();
  energyM = new TH1F("energyM","energyM", Ncaesar, 0, Ncaesar);
  timeM = new TH1F("timeM","timeM", 192, 0, 192);
  adctdcM = new TH2F("adctdcM","ADC_V_TDC_Multiplicity", 32, 0, 32, 32, 0, 32);
  caeDir->cd();
  file->cd();

  for(int i = 0; i < 6; i++) {
    adcDir[i] = caeDir->mkdir(Form("ADC%i",i+1));
  }
  for(int ic = 0; ic < Ncaesar; ic++) {
    int aN = ic / 32;
    adcDir[aN]->cd();
    detN[ic] = adcDir[aN]->mkdir(Form("Det_%3.3i", ic + 1));
    detN[ic]->cd();
    char hname[64];
    sprintf(hname,"E%3.3i",ic);
    ECaesar[ic] = new TH1F(hname, hname, 8192, 0, 16384);
    ECaesar[ic]->GetXaxis()->SetTitle("Caesar Energy (keVish)");

    char htname[64];
    sprintf(htname,"ET%3.3i",ic);
    ETCaesar[ic] = new TH2F(htname, htname, 25, 0, 25, 8192, 0, 16384);
    ETCaesar[ic]->GetXaxis()->SetTitle("TDC Mult");
    ETCaesar[ic]->GetYaxis()->SetTitle("Caesar Energy (keVish)  ");

    char tname[64];
    sprintf(tname,"T%3.3i",ic);
    TCaesar[ic] = new TH1F(tname,tname, 21001, 0, 21000);
    TCaesar[ic]->GetXaxis()->SetTitle("Caesar Time [ns]");

    char htename[64];
    sprintf(htename,"TE%3.3i",ic);
    TECaesar[ic] = new TH2F(htename, htename, 5, 0, 5, 2101, 0, 2100);
    TECaesar[ic]->GetXaxis()->SetTitle("ADC Mult");
    TECaesar[ic]->GetYaxis()->SetTitle("Caesar Time [ns]");

    char hmult[64];
    sprintf(hmult,"ADCTDC_MULT%3.3i",ic);
    mult[ic] = new TH2I(hmult, Form("ADC_V_TDC_MultiplicityChannel_%3.3i", ic), 6, 0, 6, 6, 0, 6);
    mult[ic]->GetXaxis()->SetTitle("ADC Mult");
    mult[ic]->GetYaxis()->SetTitle("TDC Mult");

    adcDir[aN]->cd();
    file->cd();
  }

  //S800 Spectra

  dirS800Raw->cd();
  //for (int i=0;i<3;i++) //will never really have e2 and e3 so no loop needed
  //{
  name.str("");
  name << "e1up_R";
  e1up = new TH1I(name.str().c_str(),"",501,0,2040);
    
  name.str("");
  name << "e1down_R";
  e1down = new TH1I(name.str().c_str(),"",501,0,2040);
 
  e1upvsQDC = new TH2I("e1up_1 vs QDC 16","",1000,0,8000,1000,0,8000);

  Te1up = new TH1I("Te1up","",501,-500,500);
  Te1down = new TH1I("Te1down","",501,-500,500);
  Tobj = new TH1I("Tobj","",501,-500,500);
  Tobj_mult = new TH1I("Tobj_mult","",10,0,10);

  ICSummary = new TH2I("ICSummary","",16,-0.5,15.5,1024,0,4095);

  CRDC1Summary = new TH2I("CRDC1Summary","",480,-0.5,480,1200,0,1200);
  CRDC2Summary = new TH2I("CRDC2Summary","",480,-0.5,480,1200,0,1200);
  CRDC1Summary_cal = new TH2I("CRDC1Summary_cal","",480,-0.5,223.5,1200,0,1200);
  CRDC2Summary_cal = new TH2I("CRDC2Summary_cal","",480,-0.5,223.5,1200,0,1200);
  CRDC1raw = new TH1I("CRCD1raw","",2000,0,2000);
  CRDC2raw = new TH1I("CRCD2raw","",2000,0,2000);
  CRDC1PadMult = new TH1I("CRDC1PadMult","",480,-0.5,223.5);  
  CRDC2PadMult = new TH1I("CRDC2PadMult","",480,-0.5,223.5);  
  CRDC1Tac = new TH1I("CRDC1Tac","",400,30000,34000);
  CRDC2Tac = new TH1I("CRDC2Tac","",400,30000,34000);
  CRDC1AnodevsTac = new TH2I("CRDC1AnodevsTac","",512,0,4096,400,30000,34000);
  CRDC2AnodevsTac = new TH2I("CRDC2AnodevsTac","",512,0,4096,400,30000,34000);
  CRDC1XrawY = new TH2I("CRDC1XrawY","",1024,0,600,1200,27000,33000);
  CRDC2XrawY = new TH2I("CRDC2XrawY","",1024,0,600,1200,27000,33000);

  dirS800Cal->cd();
  CRDC1X = new TH1I("CRDC1X","",512,-300,300);
  CRDC1Y = new TH1I("CRDC1Y","",500,-7500,-7000);
  CRDC2X = new TH1I("CRDC2X","",512,-300,300);
  CRDC2Y = new TH1I("CRDC2Y","",500,8250,8750);
  CRDC1XCalY = new TH2I("CRDC1XCalY","",1024,-300,300,512,-100,100);
  CRDC2XCalY = new TH2I("CRDC2XCalY","",1024,-300,300,512,-100,100);  
  CRDC1XY = new TH2I("CRDC1XY","",512,-300,300,600,30000,33000);
  CRDC2XY = new TH2I("CRDC2XY","",512,-300,300,600,30000,33000);
  ata1D = new TH1I("ata1D","",100,-5,5);
  bta1D = new TH1I("bta1D","",100,-5,5);
  atavsbta = new TH2I("atavsbta","",100,-5,5,100,-5,5);
  ThetavsPhi = new TH2I("ThetavsPhi","",100,0,100,720,-180,180);
  ThetavsPhi_F17 = new TH2I("ThetavsPhi_F17","",100,0,100,720,-180,180);
  ThetavsPhi_F17_fiber = new TH2I("ThetavsPhi_F17_fiber","",100,0,100,720,-180,180);
  ThetavsPhi_fiber = new TH2I("ThetavsPhi_fiber","",50,0,5,720,-180,180);
  ThetaFibervsThetaS800 = new TH2I("ThetaFibervsThetaS800","",200,-5,5,200,-5,5);
  PhiFibervsPhiS800 = new TH2I("PhiFibervsPhiS800","",720,-180,180,720,-180,180);
  CRDC1X_K35 = new TH1I("CRDC1X_K35","",512,-300,300);
  CRDC1X_K36 = new TH1I("CRDC1X_K36","",512,-300,300);
  //needed for CRDC gain matching TODO change for Si23 fragments? Not sure since we don't use CRDCs anymore
  CRDC1_Ca36_Ca36 = new TH2I("CRDC1_Ca36_Ca36","",480,-0.5,223.5,1200,0,1200);
  CRDC2_Ca36_Ca36 = new TH2I("CRDC2_Ca36_Ca36","",480,-0.5,223.5,1200,0,1200);
  CRDC1_Ca36_K35 = new TH2I("CRDC1_Ca36_K35","",480,-0.5,223.5,1200,0,1200);
  CRDC2_Ca36_K35 = new TH2I("CRDC2_Ca36_K35","",480,-0.5,223.5,1200,0,1200);
  CRDC1_K35_K35 = new TH2I("CRDC1_K35_K35","",480,-0.5,223.5,1200,0,1200);
  CRDC2_K35_K35 = new TH2I("CRDC2_K35_K35","",480,-0.5,223.5,1200,0,1200);
  CRDC1_K35_K36 = new TH2I("CRDC1_K35_K36","",480,-0.5,223.5,1200,0,1200);
  CRDC2_K35_K36 = new TH2I("CRDC2_K35_K36","",480,-0.5,223.5,1200,0,1200);
  CRDC1_Ca36_Ar33 = new TH2I("CRDC1_Ca36_Ar33","",480,-0.5,223.5,1200,0,1200);
  CRDC2_Ca36_Ar33 = new TH2I("CRDC2_Ca36_Ar33","",480,-0.5,223.5,1200,0,1200);
  CRDC1_Ca36_Ar34 = new TH2I("CRDC1_Ca36_Ar34","",480,-0.5,223.5,1200,0,1200);
  CRDC2_Ca36_Ar34 = new TH2I("CRDC2_Ca36_Ar34","",480,-0.5,223.5,1200,0,1200);
  CRDC1_K35_Ar34 = new TH2I("CRDC1_K35_Ar34","",480,-0.5,223.5,1200,0,1200);
  CRDC2_K35_Ar34 = new TH2I("CRDC2_K35_Ar34","",480,-0.5,223.5,1200,0,1200);
  CRDC1_K35_Ar35 = new TH2I("CRDC1_K35_Ar35","",480,-0.5,223.5,1200,0,1200);
  CRDC2_K35_Ar35 = new TH2I("CRDC2_K35_Ar35","",480,-0.5,223.5,1200,0,1200);
  CRDC1_Ca36_Cl32 = new TH2I("CRDC1_Ca36_Cl32","",480,-0.5,223.5,1200,0,1200);
  CRDC2_Ca36_Cl32 = new TH2I("CRDC2_Ca36_Cl32","",480,-0.5,223.5,1200,0,1200);

  crdcx23Si_23Si = new TH1F("crdcx23Si_23Si","",600,-300,300); //To check unreacted beam centering
  crdcx23S_20Mg_2p = new TH1F("crdcx23S_20Mg_2p","",600,-300,300);
  crdcx23P_22Si_p = new TH1F("crdcx23P_22Si_p","",600,-300,300);
  dirS800PID->cd();  
  ObjvsXFP = new TH2I("ObjvsXFP","",600,-450,-150,500,0,250);
  ObjvsXFP->GetXaxis()->SetTitle("T Object");
  ObjvsXFP->GetYaxis()->SetTitle("T XFP");

  ObjvsXFPwithAlpha1 = new TH2I("ObjvsXFPwithAlpha1","",600,-300,0,600,0,300);
  ObjvsXFPwithAlpha1->GetXaxis()->SetTitle("T Object");
  ObjvsXFPwithAlpha1->GetYaxis()->SetTitle("T XFP");


  ObjvsXFPwithProton1 = new TH2I("ObjvsXFPwithProton1","",600,-300,0,600,0,300);
  ObjvsXFPwithProton1->GetXaxis()->SetTitle("T Object");
  ObjvsXFPwithProton1->GetYaxis()->SetTitle("T XFP");


  ObjvsXFPwithProton2 = new TH2I("ObjvsXFPwithProton2","",600,-300,0,600,0,300);
  ObjvsXFPwithProton2->GetXaxis()->SetTitle("T Object");
  ObjvsXFPwithProton2->GetYaxis()->SetTitle("T XFP"); 
  
  ObjvsICsum = new TH2I("ObjvsICsum","",1500,-500,0,2048,0,4095);
  ObjvsICsum->GetXaxis()->SetTitle("T Objct");
  ObjvsICsum->GetYaxis()->SetTitle("IC Sum");

  ObjvsICsum_nobeam = new TH2I("ObjvsICsum_nobeam","",1500,-500,0,2048,0,4095);

  ObjvsICsum_corr = new TH2I("ObjvsICsum_corr","",300,-300,-150,2000,0,2000);
  ObjvsICsum_corr->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_corr->GetYaxis()->SetTitle("IC Sum");

  ObjvsICsum_Si23 = new TH2I("ObjvsICsum_Si23","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Si23_1p = new TH2I("ObjvsICsum_Si23_1p","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Si23_2p = new TH2I("ObjvsICsum_Si23_2p","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Si23_3p = new TH2I("ObjvsICsum_Si23_3p","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Si23_alpha = new TH2I("ObjvsICsum_Si23_alpha","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Al22 = new TH2I("ObjvsICsum_Al22","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Al22_1p = new TH2I("ObjvsICsum_Al22_1p","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Al22_2p = new TH2I("ObjvsICsum_Al22_2p","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Al22_3p = new TH2I("ObjvsICsum_Al22_3p","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Al22_alpha = new TH2I("ObjvsICsum_Al22_alpha","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Mg21 = new TH2I("ObjvsICsum_Mg21","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Mg21_1p = new TH2I("ObjvsICsum_Mg21_1p","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Mg21_2p = new TH2I("ObjvsICsum_Mg21_2p","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Mg21_3p = new TH2I("ObjvsICsum_Mg21_3p","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Mg21_alpha = new TH2I("ObjvsICsum_Mg21_alpha","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Na20 = new TH2I("ObjvsICsum_Na20","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Na20_1p = new TH2I("ObjvsICsum_Na20_1p","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Na20_2p = new TH2I("ObjvsICsum_Na20_2p","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Na20_3p = new TH2I("ObjvsICsum_Na20_3p","",300,-300,-150,2000,0,2000);
  ObjvsICsum_Na20_alpha = new TH2I("ObjvsICsum_Na20_alpha","",300,-300,-150,2000,0,2000);

  Timing1vsICsum_Ca37 = new TH2I("Timing1vsICsum_Ca37","",1200,-200,-50,2048,0,4095);
  Timing1vsICsum_Ca37->GetXaxis()->SetTitle("T diff");
  Timing1vsICsum_Ca37->GetYaxis()->SetTitle("IC Sum");
  Timing2vsICsum_Ca37 = new TH2I("Timing2vsICsum_Ca37","",1200,-200,-50,2048,0,4095);
  Timing2vsICsum_Ca37->GetXaxis()->SetTitle("T diff");
  Timing2vsICsum_Ca37->GetYaxis()->SetTitle("IC Sum");
  XFPvsICsum_Ca37 = new TH2I("XFPvsICsum_Ca37","",600,0,300,2048,0,4095);
  XFPvsICsum_Ca37->GetXaxis()->SetTitle("XFP time");
  XFPvsICsum_Ca37->GetYaxis()->SetTitle("IC Sum");


  ObjvsICsum_wProt = new TH2I("ObjvsICsum_wProt","",300,-300,-150,2000,0,2000);
  ObjvsICsum_wProt->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_wProt->GetYaxis()->SetTitle("IC Sum");

  ObjvsICsum_wProt_wFiber = new TH2I("ObjvsICsum_wProt_wFiber","",1200,-200,-50,2048,0,4095);
  ObjvsICsum_wProt_wFiber->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_wProt_wFiber->GetYaxis()->SetTitle("IC Sum");


  ObjvsICsum_wProt2 = new TH2I("ObjvsICsum_wProt2","",300,-300,-150,2000,0,2000);
  ObjvsICsum_wProt2->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_wProt2->GetYaxis()->SetTitle("IC Sum");

  Objvsafp = new TH2I("Objvsafp","",600,-300,0,1000,-1,1);
  Objvsafp->GetXaxis()->SetTitle("Obj Time");
  Objvsafp->GetYaxis()->SetTitle("afp");
  
  ObjvsCRDC1X = new TH2I("ObjvsCRDC1X","",600,-300,0,512,-300,300);
  ObjvsCRDC1X->GetXaxis()->SetTitle("Obj Time");
  ObjvsCRDC1X->GetYaxis()->SetTitle("CRDC1 X");
  ObjvsCRDC2X = new TH2I("ObjvsCRDC2X","",600,-300,0,512,-300,300);
  ObjvsCRDC2X->GetXaxis()->SetTitle("Obj Time");
  ObjvsCRDC2X->GetYaxis()->SetTitle("CRDC2 X");
  
  XfpUncvsCRDC1X = new TH2I("XfpUncvsCRDC1X","",600,0,300,512,-300,300);
  XfpUncvsCRDC1X->GetXaxis()->SetTitle("FP Time");
  XfpUncvsCRDC1X->GetYaxis()->SetTitle("CRDC1 X");
  XfpUncvsCRDC2X = new TH2I("XfpUncvsCRDC2X","",600,-0,300,512,-300,300);
  XfpUncvsCRDC2X->GetXaxis()->SetTitle("FP Time");
  XfpUncvsCRDC2X->GetYaxis()->SetTitle("CRDC2 X");



  rigidityK35 = new TH1I("rigidityK35","",100,1.85,2.25);
  rigidityK35->GetXaxis()->SetTitle("Rigidity (Tm)");
  rigidityCa36 = new TH1I("rigidityCa36","",100,1.85,2.25);
  rigidityCa36->GetXaxis()->SetTitle("Rigidity (Tm)");
  rigidityCa37beam = new TH1I("rigidityCa37beam","",1000,1.85,2.25);
  rigidityCa37beam->GetXaxis()->SetTitle("Rigidity (Tm)");
  rigidityCa38beam = new TH1I("rigidityCa38beam","",1000,1.85,2.25);
  rigidityCa38beam->GetXaxis()->SetTitle("Rigidity (Tm)");
  rigidityK36beam = new TH1I("rigidityK36beam","",1000,1.85,2.25);
  rigidityK36beam->GetXaxis()->SetTitle("Rigidity (Tm)");
  rigidityAr35beam = new TH1I("rigidityAr35beam","",1000,1.85,2.25);
  rigidityAr35beam->GetXaxis()->SetTitle("Rigidity (Tm)");

  energySi23_Si23 = new TH1I("energySi23_Si23","",500,1500,2500);
  energySi23_Si23_Fib = new TH1I("energySi23_Si23_Fib","",500,1500,2500);
  energySi23_Si23_FibTarg = new TH1I("energySi23_Si23_FibTarg","",500,1500,2500);

  energyMg21_Mg21 = new TH1I("energyMg21_Mg21","",500,1500,2500);
  energyMg21_Mg21_Fib = new TH1I("energyMg21_Mg21_Fib","",500,1500,2500);
  energyMg21_Mg21_FibTarg = new TH1I("energyMg21_Mg21_FibTarg","",500,1500,2500);

  BeamCa37_energy = new TH1I("BeamCa37_energy","",350,1500,2200);
  BeamCa37_ptot = new TH1I("BeamCa37_ptot","",300,10000,13000);
  BeamCa37_ppar = new TH1I("BeamCa37_ppar","",300,10000,13000);
  BeamCa37_ptra = new TH1I("BeamCa37_ptra","",400,0,800);
  BeamCa37_TofAngle = new TH2I("BeamCa37_TofAngle","",100,-135,-125,100,-0.15,0.1);
  BeamCa37_TofCRDC1raw = new TH2I("BeamCa37_TofCRDC1raw","",100,-135,-125,300,0,1200);
  BeamCa37_TofCRDC1cal = new TH2I("BeamCa37_TofCRDC1cal","",100,-135,-125,300,0,1200);
  BeamCa37_ObjvsCRDC1X = new TH2I("BeamCa37_ObjvsCRDC1X","",600,-300,0,512,-300,300);
  BeamCa37_Fiber_XY = new TH2I("BeamCa37_Fiber_XY","",64,-8,8,64,-8,8);


  BeamCa38_energy = new TH1I("BeamCa38_energy","",350,1500,2200);
  BeamCa38_ptot = new TH1I("BeamCa38_ptot","",300,10000,13000);
  BeamCa38_ppar = new TH1I("BeamCa38_ppar","",300,10000,13000);
  BeamCa38_ptra = new TH1I("BeamCa38_ptra","",400,0,800);
  BeamCa38_TofAngle = new TH2I("BeamCa38_TofAngle","",100,-135,-125,100,-0.15,0.1);
  BeamCa38_TofCRDC1raw = new TH2I("BeamCa38_TofCRDC1raw","",100,-135,-125,300,0,1200);
  BeamCa38_TofCRDC1cal = new TH2I("BeamCa38_TofCRDC1cal","",100,-135,-125,300,0,1200);
  BeamCa38_ObjvsCRDC1X = new TH2I("BeamCa38_ObjvsCRDC1X","",600,-300,0,512,-300,300);
  BeamCa38_Fiber_XY = new TH2I("BeamCa38_Fiber_XY","",64,-8,8,64,-8,8);

  BeamK36_energy = new TH1I("BeamK36_energy","",350,1500,2200);
  BeamK36_ptot = new TH1I("BeamK36_ptot","",300,10000,13000);
  BeamK36_ppar = new TH1I("BeamK36_ppar","",300,10000,13000);
  BeamK36_ptra = new TH1I("BeamK36_ptra","",400,0,800);

  BeamAr35_energy = new TH1I("BeamAr35_energy","",350,1500,2200);
  BeamAr35_ptot = new TH1I("BeamAr35_ptot","",300,10000,13000);
  BeamAr35_ppar = new TH1I("BeamAr35_ppar","",300,10000,13000);
  BeamAr35_ptra = new TH1I("BeamAr35_ptra","",400,0,800);



  dirS800veldist->cd();
  residuetheta = new TH1I("residuetheta","",200,0,0.05);

  //********************************************************
  //Correlation Combs
  //********************************************************
  dirCorrCombs->cd();
  SiBeam_CorrComb = new TH2I("SiBeam_CorrComb","",12,0,12,7,0,7);
  SiBeam_CorrComb->SetMinimum(1.0);
  AlBeam_CorrComb = new TH2I("AlBeam_CorrComb","",8,0,8,7,0,7);
  AlBeam_CorrComb->SetMinimum(1.0);

  //********************************************************
  //Correlation hists
  //********************************************************
  //Proton Directory
  dir1H->cd();
  protonKE = new TH1I("protonKE","",300,0,200); //500 keV bins
  protonKE_Si23Beam = new TH1I("protonKE_Si23Beam","",300,0,200); //500 keV bins
  protonKE_Al22Beam = new TH1I("protonKE_Al22Beam","",300,0,200); //500 keV bins
  protonKE_Mg21Beam = new TH1I("protonKE_Mg21Beam","",300,0,200); //500 keV bins

  //F16 -> p + O15
  dir16F->cd();
  Erel_16F_p15O = new TH1I("Erel_16F_p15O","",500,-5,15);
  Erel_16F_p15O_trans = new TH1I("Erel_16F_p15O_trans","",500,-5,15);
  Erel_16F_p15O_long = new TH1I("Erel_16F_p15O_long","",500,-5,15);
  Ex_16F_p15O = new TH1I("Ex_16F_p15O","",500,-5,15);
  Ex_16F_p15O_trans = new TH1I("Ex_16F_p15O_trans","",500,-5,15);
  Ex_16F_p15O_long = new TH1I("Ex_16F_p15O_long","",500,-5,15);
  ThetaCM_16F_p15O = new TH1I("ThetaCM_16F_p15O","",200,0,10);
  VCM_16F_p15O = new TH1I("VCM_16F_p15O","",400,7,15);
  Erel_p15O_costhetaH = new TH2I("Erel_p15O_costhetaH","",800,0,8,50,-1,1);
  Erel_p15O_costhetaH->SetMinimum(1.0);
  p15O_VCMvsErel = new TH2I("p15O_VCMvsErel","",400,7,15,500,-5,15);
  p15O_gammasADD = new TH1I("p15O_gammasADD","",2048, 0, 8192);
  p15O_gammasADDvsErel = new TH2I("p15O_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p15O_gammasADDvsErel->SetMinimum(1.0);
  p15O_gammasADDvsgammasADD = new TH2I("p15O_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p15O_gammasADDvsgammasADD->SetMinimum(1.0);

  //F17 -> p + O16
  dir17F->cd();
  Erel_17F_p16O = new TH1I("Erel_17F_p16O","",500,-5,15);
  Erel_17F_p16O_trans = new TH1I("Erel_17F_p16O_trans","",500,-5,15);
  Erel_17F_p16O_long = new TH1I("Erel_17F_p16O_long","",500,-5,15);
  Ex_17F_p16O = new TH1I("Ex_17F_p16O","",500,-5,15);
  Ex_17F_p16O_trans = new TH1I("Ex_17F_p16O_trans","",500,-5,15);
  Ex_17F_p16O_long = new TH1I("Ex_17F_p16O_long","",500,-5,15);
  ThetaCM_17F_p16O = new TH1I("ThetaCM_17F_p16O","",200,0,10);
  VCM_17F_p16O = new TH1I("VCM_17F_p16O","",400,7,15);
  Erel_p16O_costhetaH = new TH2I("Erel_p16O_costhetaH","",800,0,8,50,-1,1);
  Erel_p16O_costhetaH->SetMinimum(1.0);
  p16O_VCMvsErel = new TH2I("p16O_VCMvsErel","",400,7,15,500,-5,15);
  p16O_gammasADD = new TH1I("p16O_gammasADD","",2048, 0, 8192);
  p16O_gammasADDvsErel = new TH2I("p16O_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p16O_gammasADDvsErel->SetMinimum(1.0);
  p16O_gammasADDvsgammasADD = new TH2I("p16O_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p16O_gammasADDvsgammasADD->SetMinimum(1.0);

  //Ne18 -> p + F17
  dir18Ne->cd();
  Erel_18Ne_p17F = new TH1I("Erel_18Ne_p17F","",500,-5,15);
  Erel_18Ne_p17F_trans = new TH1I("Erel_18Ne_p17F_trans","",500,-5,15);
  Erel_18Ne_p17F_long = new TH1I("Erel_18Ne_p17F_long","",500,-5,15);
  Ex_18Ne_p17F = new TH1I("Ex_18Ne_p17F","",500,-5,15);
  Ex_18Ne_p17F_trans = new TH1I("Ex_18Ne_p17F_trans","",500,-5,15);
  Ex_18Ne_p17F_long = new TH1I("Ex_18Ne_p17F_long","",500,-5,15);
  ThetaCM_18Ne_p17F = new TH1I("ThetaCM_18Ne_p17F","",200,0,10);
  VCM_18Ne_p17F = new TH1I("VCM_18Ne_p17F","",400,7,15);
  Erel_p17F_costhetaH = new TH2I("Erel_p17F_costhetaH","",800,0,8,50,-1,1);
  Erel_p17F_costhetaH->SetMinimum(1.0);
  p17F_VCMvsErel = new TH2I("p17F_VCMvsErel","",400,7,15,500,-5,15);
  p17F_gammasADD = new TH1I("p17F_gammasADD","",2048, 0, 8192);
  p17F_gammasADDvsErel = new TH2I("p17F_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p17F_gammasADDvsErel->SetMinimum(1.0);
  p17F_gammasADDvsgammasADD = new TH2I("p17F_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p17F_gammasADDvsgammasADD->SetMinimum(1.0);

  //Ne19 -> p + F18
  dir19Ne->cd();
  Erel_19Ne_p18F = new TH1I("Erel_19Ne_p18F","",500,-5,15);
  Erel_19Ne_p18F_trans = new TH1I("Erel_19Ne_p18F_trans","",500,-5,15);
  Erel_19Ne_p18F_long = new TH1I("Erel_19Ne_p18F_long","",500,-5,15);
  Ex_19Ne_p18F = new TH1I("Ex_19Ne_p18F","",500,-5,15);
  Ex_19Ne_p18F_trans = new TH1I("Ex_19Ne_p18F_trans","",500,-5,15);
  Ex_19Ne_p18F_long = new TH1I("Ex_19Ne_p18F_long","",500,-5,15);
  ThetaCM_19Ne_p18F = new TH1I("ThetaCM_19Ne_p18F","",200,0,10);
  VCM_19Ne_p18F = new TH1I("VCM_19Ne_p18F","",400,7,15);
  Erel_p18F_costhetaH = new TH2I("Erel_p18F_costhetaH","",800,0,8,50,-1,1);
  Erel_p18F_costhetaH->SetMinimum(1.0);
  p18F_VCMvsErel = new TH2I("p18F_VCMvsErel","",400,7,15,500,-5,15);
  p18F_gammasADD = new TH1I("p18F_gammasADD","",2048, 0, 8192);
  p18F_gammasADDvsErel = new TH2I("p18F_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p18F_gammasADDvsErel->SetMinimum(1.0);
  p18F_gammasADDvsgammasADD = new TH2I("p18F_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p18F_gammasADDvsgammasADD->SetMinimum(1.0);

  //Na19 -> p + Ne18
  dir19Na->cd();
  Erel_19Na_p18Ne = new TH1I("Erel_19Na_p18Ne","",500,-5,15);
  Erel_19Na_p18Ne_trans = new TH1I("Erel_19Na_p18Ne_trans","",500,-5,15);
  Erel_19Na_p18Ne_long = new TH1I("Erel_19Na_p18Ne_long","",500,-5,15);
  Ex_19Na_p18Ne = new TH1I("Ex_19Na_p18Ne","",500,-5,15);
  Ex_19Na_p18Ne_trans = new TH1I("Ex_19Na_p18Ne_trans","",500,-5,15);
  Ex_19Na_p18Ne_long = new TH1I("Ex_19Na_p18Ne_long","",500,-5,15);
  Erel_19Na_p18Ne_set1 = new TH1I("Erel_19Na_p18Ne_set1","",500,-5,15);
  Erel_19Na_p18Ne_set1a = new TH1I("Erel_19Na_p18Ne_set1a","",500,-5,15);
  ThetaCM_19Na_p18Ne = new TH1I("ThetaCM_19Na_p18Ne","",200,0,10);
  VCM_19Na_p18Ne = new TH1I("VCM_19Na_p18Ne","",400,7,15);
  Erel_p18Ne_costhetaH = new TH2I("Erel_p18Ne_costhetaH","",800,0,8,50,-1,1);
  Erel_p18Ne_costhetaH->SetMinimum(1.0);
  p18Ne_VCMvsErel = new TH2I("p18Ne_VCMvsErel","",400,7,15,500,-5,15);
  p18Ne_gammasADD = new TH1I("p18Ne_gammasADD","",2048, 0, 8192);
  p18Ne_gammasADD_tgate = new TH1I("p18Ne_gammasADD_tgate","",2048, 0, 8192);
  p18Ne_gammasADDvsErel = new TH2I("p18Ne_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p18Ne_gammasADDvsErel->SetMinimum(1.0);
  p18Ne_gammasADDvsgammasADD = new TH2I("p18Ne_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p18Ne_gammasADDvsgammasADD->SetMinimum(1.0);

  //Na20 -> p + Ne19
  dir20Na->cd();
  Erel_20Na_p19Ne = new TH1I("Erel_20Na_p19Ne","",500,-5,15);
  Erel_20Na_p19Ne_trans = new TH1I("Erel_20Na_p19Ne_trans","",500,-5,15);
  Ex_20Na_p19Ne = new TH1I("Ex_20Na_p19Ne","",500,0,20);
  ThetaCM_20Na_p19Ne = new TH1I("ThetaCM_20Na_p19Ne","",200,0,10);
  VCM_20Na_p19Ne = new TH1I("VCM_20Na_p19Ne","",400,7,15);
  Erel_p19Ne_costhetaH = new TH2I("Erel_p19Ne_costhetaH","",800,0,8,50,-1,1);
  Erel_p19Ne_costhetaH->SetMinimum(1.0);
  p19Ne_VCMvsErel = new TH2I("p19Ne_VCMvsErel","",400,7,15,500,-5,15);
  p19Ne_gammasADD = new TH1I("p19Ne_gammasADD","",2048, 0, 8192);
  p19Ne_gammasADD_tgate = new TH1I("p19Ne_gammasADD_tgate","",2048, 0, 8192);
  p19Ne_gammasADDvsErel = new TH2I("p19Ne_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p19Ne_gammasADDvsErel->SetMinimum(1.0);
  p19Ne_gammasADDvsgammasADD = new TH2I("p19Ne_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p19Ne_gammasADDvsgammasADD->SetMinimum(1.0);

  //Na21 -> p + Ne20
  dir21Na->cd();
  Erel_21Na_p20Ne = new TH1I("Erel_21Na_p20Ne","",500,-5,15);
  Erel_21Na_p20Ne_trans = new TH1I("Erel_21Na_p20Ne_trans","",500,-5,15);
  Ex_21Na_p20Ne = new TH1I("Ex_21Na_p20Ne","",500,0,20);
  ThetaCM_21Na_p20Ne = new TH1I("ThetaCM_21Na_p20Ne","",200,0,10);
  VCM_21Na_p20Ne = new TH1I("VCM_21Na_p20Ne","",400,7,15);
  Erel_p20Ne_costhetaH = new TH2I("Erel_p20Ne_costhetaH","",800,0,8,50,-1,1);
  Erel_p20Ne_costhetaH->SetMinimum(1.0);
  p20Ne_VCMvsErel = new TH2I("p20Ne_VCMvsErel","",400,7,15,500,-5,15);
  p20Ne_gammasADD = new TH1I("p20Ne_gammasADD","",2048, 0, 8192);
  p20Ne_gammasADD_tgate = new TH1I("p20Ne_gammasADD_tgate","",2048, 0, 8192);
  p20Ne_gammasADDvsErel = new TH2I("p20Ne_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p20Ne_gammasADDvsErel->SetMinimum(1.0);
  p20Ne_gammasADDvsgammasADD = new TH2I("p20Ne_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p20Ne_gammasADDvsgammasADD->SetMinimum(1.0);

  //Mg19 -> 2p + Ne17
  dir19Mg->cd();
  Erel_19Mg_2p17Ne = new TH1I("Erel_19Mg_2p17Ne","",500,-5,15);
  Erel_19Mg_2p17Ne_trans = new TH1I("Erel_19Mg_2p17Ne_trans","",500,-5,15);
  Ex_19Mg_2p17Ne = new TH1I("Ex_19Mg_2p17Ne","",500,0,20);
  ThetaCM_19Mg_2p17Ne = new TH1I("ThetaCM_19Mg_2p17Ne","",200,0,10);
  VCM_19Mg_2p17Ne = new TH1I("VCM_19Mg_2p17Ne","",400,7,15);
  Erel_2p17Ne_costhetaH = new TH2I("Erel_2p17Ne_costhetaH","",800,0,8,25,-1,1);
  Erel_2p17Ne_costhetaH->SetMinimum(1.0);
  pp17Ne_VCMvsErel = new TH2I("pp17Ne_VCMvsErel","",400,7,15,500,-5,15);
  pp17Ne_gammasADD = new TH1I("pp17Ne_gammasADD","",2048, 0, 8192);
  pp17Ne_gammasADD_tgate = new TH1I("pp17Ne_gammasADD_tgate","",2048, 0, 8192);
  pp17Ne_gammasADDvsErel = new TH2I("pp17Ne_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  pp17Ne_gammasADDvsErel->SetMinimum(1.0);
  pp17Ne_gammasADDvsgammasADD = new TH2I("pp17Ne_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  pp17Ne_gammasADDvsgammasADD->SetMinimum(1.0);

  //Mg20 -> 2p + Ne18
  dir20Mg->cd();
  Erel_20Mg_2p18Ne = new TH1I("Erel_20Mg_2p18Ne","",500,-5,15);
  Erel_20Mg_2p18Ne_trans = new TH1I("Erel_20Mg_2p18Ne_trans","",500,-5,15);
  Ex_20Mg_2p18Ne = new TH1I("Ex_20Mg_2p18Ne","",500,0,20);
  ThetaCM_20Mg_2p18Ne = new TH1I("ThetaCM_20Mg_2p18Ne","",200,0,10);
  VCM_20Mg_2p18Ne = new TH1I("VCM_20Mg_2p18Ne","",400,7,15);
  Erel_2p18Ne_costhetaH = new TH2I("Erel_2p18Ne_costhetaH","",800,0,8,50,-1,1);
  Erel_2p18Ne_costhetaH->SetMinimum(1.0);
  pp18Ne_VCMvsErel = new TH2I("pp18Ne_VCMvsErel","",400,7,15,500,-5,15);
  pp18Ne_gammasADD = new TH1I("pp18Ne_gammasADD","",2048, 0, 8192);
  pp18Ne_gammasADD_tgate = new TH1I("pp18Ne_gammasADD_tgate","",2048, 0, 8192);
  pp18Ne_gammasADDvsErel = new TH2I("pp18Ne_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  pp18Ne_gammasADDvsErel->SetMinimum(1.0);
  pp18Ne_gammasADDvsgammasADD = new TH2I("pp18Ne_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  pp18Ne_gammasADDvsgammasADD->SetMinimum(1.0);

  //Mg21 -> p + Na20
  dir21Mg->cd();
  Erel_21Mg_p20Na = new TH1I("Erel_21Mg_p20Na","",500,-5,15);
  Erel_21Mg_p20Na_trans = new TH1I("Erel_21Mg_p20Na_trans","",500,-5,15);
  Erel_21Mg_p20Na_long = new TH1I("Erel_21Mg_p20Na_long","",500,-5,15);
  Ex_21Mg_p20Na = new TH1I("Ex_21Mg_p20Na","",500,0,20);
  Ex_21Mg_p20Na_trans = new TH1I("Ex_21Mg_p20Na_trans","",500,0,20);
  Ex_21Mg_p20Na_long = new TH1I("Ex_21Mg_p20Na_long","",500,0,20);
  Erel_21Mg_p20Na_set1 = new TH1I("Erel_21Mg_p20Na_set1","",500,-5,15);
  Erel_21Mg_p20Na_set1a = new TH1I("Erel_21Mg_p20Na_set1a","",500,-5,15);
  ThetaCM_21Mg_p20Na = new TH1I("ThetaCM_21Mg_p20Na","",200,0,10);
  VCM_21Mg_p20Na = new TH1I("VCM_21Mg_p20Na","",400,7,15);
  Erel_p20Na_costhetaH = new TH2I("Erel_p20Na_costhetaH","",800,0,8,50,-1,1);
  Erel_p20Na_costhetaH->SetMinimum(1.0);
  p20Na_VCMvsErel = new TH2I("p20Na_VCMvsErel","",400,7,15,500,-5,15);
  p20Na_gammasADD = new TH1I("p20Na_gammasADD","",2048, 0, 8192);
  p20Na_gammasADD_tgate = new TH1I("p20Na_gammasADD_tgate","",2048, 0, 8192);
  p20Na_gammasADDvsErel = new TH2I("p20Na_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p20Na_gammasADDvsErel->SetMinimum(1.0);
  p20Na_gammasADDvsgammasADD = new TH2I("p20Na_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p20Na_gammasADDvsgammasADD->SetMinimum(1.0);
  Mg21_pKE = new TH1I("Mg21_pKE","",200, 0, 200);  

  //Mg21 -> 2p + 19Ne
  dir2p19Ne->cd();
  Erel_21Mg_2p19Ne = new TH1I("Erel_21Mg_2p19Ne","",500,-5,15);
  Erel_21Mg_2p19Ne_trans = new TH1I("Erel_21Mg_2p19Ne_trans","",500,-5,15);
  Ex_21Mg_2p19Ne = new TH1I("Ex_21Mg_2p19Ne","",500,0,20);
  ThetaCM_21Mg_2p19Ne = new TH1I("ThetaCM_21Mg_2p19Ne","",200,0,10);
  VCM_21Mg_2p19Ne = new TH1I("VCM_21Mg_2p19Ne","",400,7,15);
  Erel_2p19Ne_costhetaH = new TH2I("Erel_2p219Ne_costhetaH","",800,0,8,50,-1,1);
  Erel_2p19Ne_costhetaH->SetMinimum(1.0);
  pp19Ne_VCMvsErel = new TH2I("pp19Ne_VCMvsErel","",400,7,15,500,-5,15);
  pp19Ne_gammasADD = new TH1I("pp19Ne_gammasADD","",2048, 0, 8192);
  pp19Ne_gammasADD_tgate = new TH1I("pp19Ne_gammasADD_tgate","",2048, 0, 8192);
  pp19Ne_gammasADDvsErel = new TH2I("pp19Ne_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  pp19Ne_gammasADDvsErel->SetMinimum(1.0);
  pp19Ne_gammasADDvsgammasADD = new TH2I("pp19Ne_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  pp19Ne_gammasADDvsgammasADD->SetMinimum(1.0);
  Mg21_2p_pKE = new TH1I("Mg21_2p_pKE","",200, 0, 200);  

  //Mg22 -> p + Na21
  dir22Mg->cd();
  Erel_22Mg_p21Na = new TH1I("Erel_22Mg_p21Na","",500,-5,15);
  Erel_22Mg_p21Na_trans = new TH1I("Erel_22Mg_p21Na_trans","",500,-5,15);
  Ex_22Mg_p21Na = new TH1I("Ex_22Mg_p21Na","",500,0,20);
  ThetaCM_22Mg_p21Na = new TH1I("ThetaCM_22Mg_p21Na","",200,0,10);
  VCM_22Mg_p21Na = new TH1I("VCM_22Mg_p21Na","",400,7,15);
  Erel_p21Na_costhetaH = new TH2I("Erel_p21Na_costhetaH","",800,0,8,50,-1,1);
  Erel_p21Na_costhetaH->SetMinimum(1.0);
  p21Na_VCMvsErel = new TH2I("p21Na_VCMvsErel","",400,7,15,500,-5,15);
  p21Na_gammasADD = new TH1I("p21Na_gammasADD","",2048, 0, 8192);
  p21Na_gammasADD_tgate = new TH1I("p21Na_gammasADD_tgate","",2048, 0, 8192);
  p21Na_gammasADDvsErel = new TH2I("p21Na_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p21Na_gammasADDvsErel->SetMinimum(1.0);
  p21Na_gammasADDvsgammasADD = new TH2I("p21Na_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p21Na_gammasADDvsgammasADD->SetMinimum(1.0);

  //Al21 -> p + Mg20
  dir21Al->cd();
  Erel_21Al_p20Mg = new TH1I("Erel_21Al_p20Mg","",500,-5,15);
  Erel_21Al_p20Mg_trans = new TH1I("Erel_21Al_p20Mg_trans","",500,-5,15);
  Ex_21Al_p20Mg = new TH1I("Ex_21Al_p20Mg","",500,0,20);
  ThetaCM_21Al_p20Mg = new TH1I("ThetaCM_21Al_p20Mg","",200,0,10);
  VCM_21Al_p20Mg = new TH1I("VCM_21Al_p20Mg","",400,7,15);
  Erel_p20Mg_costhetaH = new TH2I("Erel_p20Mg_costhetaH","",800,0,8,50,-1,1);
  Erel_p20Mg_costhetaH->SetMinimum(1.0);
  p20Mg_VCMvsErel = new TH2I("p20Mg_VCMvsErel","",400,7,15,500,-5,15);
  p20Mg_gammasADD = new TH1I("p20Mg_gammasADD","",2048, 0, 8192);
  p20Mg_gammasADD_tgate = new TH1I("p20Mg_gammasADD_tgate","",2048, 0, 8192);
  p20Mg_gammasADDvsErel = new TH2I("p20Mg_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p20Mg_gammasADDvsErel->SetMinimum(1.0);
  p20Mg_gammasADDvsgammasADD = new TH2I("p20Mg_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p20Mg_gammasADDvsgammasADD->SetMinimum(1.0);
  Al21_pKE = new TH1I("Al21_pKE","",200, 0, 200);  

  //Al22 -> p + Mg21
  dir22Al->cd();
  Erel_22Al_p21Mg = new TH1I("Erel_22Al_p21Mg","",500,-5,15);
  Erel_22Al_p21Mg_trans = new TH1I("Erel_22Al_p21Mg_trans","",500,-5,15);
  Ex_22Al_p21Mg = new TH1I("Ex_22Al_p21Mg","",500,0,20);
  ThetaCM_22Al_p21Mg = new TH1I("ThetaCM_22Al_p21Mg","",200,0,10);
  VCM_22Al_p21Mg = new TH1I("VCM_22Al_p21Mg","",400,7,15);
  Erel_p21Mg_costhetaH = new TH2I("Erel_p21Mg_costhetaH","",800,0,8,50,-1,1);
  Erel_p21Mg_costhetaH->SetMinimum(1.0);
  p21Mg_VCMvsErel = new TH2I("p21Mg_VCMvsErel","",400,7,15,500,-5,15);
  p21Mg_gammasADD_nodopp = new TH1I("p21Mg_gammasADD_nodopp","",2048, 0, 8192);
  p21Mg_gammasADD = new TH1I("p21Mg_gammasADD","",2048, 0, 8192);
  p21Mg_gammasADD_tgate = new TH1I("p21Mg_gammasADD_tgate","",2048, 0, 8192);
  p21Mg_gammasADDvsErel = new TH2I("p21Mg_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p21Mg_gammasADDvsErel->SetMinimum(1.0);
  p21Mg_gammasADD_peak2 = new TH1I("p21Mg_gammasADD_peak2","",2048, 0, 8192);
  p21Mg_gammasADD_peak2_tgate = new TH1I("p21Mg_gammasADD_peak2_tgate","",2048, 0, 8192);
  p21Mg_gammasADD_peak3 = new TH1I("p21Mg_gammasADD_peak3","",2048, 0, 8192);
  p21Mg_gammasADD_peak3_tgate = new TH1I("p21Mg_gammasADD_peak3_tgate","",2048, 0, 8192);
  p21Mg_gammasADDvsgammasADD = new TH2I("p21Mg_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p21Mg_gammasADDvsgammasADD->SetMinimum(1.0);
  Al22_pKE = new TH1I("Al22_pKE","",200, 0, 200);  

  //Al23 -> p + Mg22
  dir23Al->cd();
  Erel_23Al_p22Mg = new TH1I("Erel_23Al_p22Mg","",500,-5,15);
  Erel_23Al_p22Mg_trans = new TH1I("Erel_23Al_p22Mg_trans","",500,-5,15);
  Erel_23Al_p22Mg_long = new TH1I("Erel_23Al_p22Mg_long","",500,-5,15);
  Ex_23Al_p22Mg = new TH1I("Ex_23Al_p22Mg","",500,0,20);
  Ex_23Al_p22Mg_trans = new TH1I("Ex_23Al_p22Mg_trans","",500,0,20);
  Ex_23Al_p22Mg_long = new TH1I("Ex_23Al_p22Mg_long","",500,0,20);
  Erel_23Al_p22Mg_set1 = new TH1I("Erel_23Al_p22Mg_set1","",500,-5,15);
  Erel_23Al_p22Mg_set1a = new TH1I("Erel_23Al_p22Mg_set1a","",500,-5,15);
  ThetaCM_23Al_p22Mg = new TH1I("ThetaCM_23Al_p22Mg","",200,0,10);
  VCM_23Al_p22Mg = new TH1I("VCM_23Al_p22Mg","",400,7,15);
  Erel_p22Mg_costhetaH = new TH2I("Erel_p22Mg_costhetaH","",800,0,8,50,-1,1);
  Erel_p22Mg_costhetaH->SetMinimum(1.0);
  p22Mg_VCMvsErel = new TH2I("p22Mg_VCMvsErel","",400,7,15,500,-5,15);
  p22Mg_gammasADD = new TH1I("p22Mg_gammasADD","",2048, 0, 8192);
  p22Mg_gammasADD_tgate = new TH1I("p22Mg_gammasADD_tgate","",2048, 0, 8192);
  p22Mg_gammasADDvsErel = new TH2I("p22Mg_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p22Mg_gammasADDvsErel->SetMinimum(1.0);
  p22Mg_gammasADDvsgammasADD = new TH2I("p22Mg_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p22Mg_gammasADDvsgammasADD->SetMinimum(1.0);

  //Si22 -> 2p + 20Mg
  dir22Si->cd();
  Erel_22Si_2p20Mg = new TH1I("Erel_22Si_2p20Mg","",500,-5,15);
  Erel_22Si_2p20Mg_trans = new TH1I("Erel_22Si_2p20Mg_trans","",500,-5,15);
  Ex_22Si_2p20Mg = new TH1I("Ex_22Si_2p20Mg","",500,0,20);
  ThetaCM_22Si_2p20Mg = new TH1I("ThetaCM_22Si_2p20Mg","",200,0,10);
  VCM_22Si_2p20Mg = new TH1I("VCM_22Si_2p20Mg","",400,7,15);
  Erel_2p20Mg_costhetaH = new TH2I("Erel_2p20Mg_costhetaH","",800,0,8,50,-1,1);
  Erel_2p20Mg_costhetaH->SetMinimum(1.0);
  pp20Mg_VCMvsErel = new TH2I("pp20Mg_VCMvsErel","",400,7,15,500,-5,15);
  pp20Mg_gammasADD = new TH1I("pp20Mg_gammasADD","",2048, 0, 8192);
  pp20Mg_gammasADD_tgate = new TH1I("pp20Mg_gammasADD_tgate","",2048, 0, 8192);
  pp20Mg_gammasADDvsErel = new TH2I("pp20Mg_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  pp20Mg_gammasADDvsErel->SetMinimum(1.0);
  pp20Mg_gammasADDvsgammasADD = new TH2I("pp20Mg_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  pp20Mg_gammasADDvsgammasADD->SetMinimum(1.0);
  Si22_pKE = new TH1I("Si22_pKE","",200, 0, 200);  
  //Jacobi
  Si22_JacobiT_xy_s = new TH2I("Si22_JacobiT_xy_s","",50,0,1,50,-1,1);
  Si22_JacobiY_xy_s = new TH2I("Si22_JacobiY_xy_s","",50,0,1,50,-1,1);

  //Si23 -> p + 22Al
  dirp22Al->cd();
  Erel_23Si_p22Al = new TH1I("Erel_23Si_p22Al","",500,-5,15);
  Erel_23Si_p22Al_trans = new TH1I("Erel_23Si_p22Al_trans","",500,-5,15);
  Ex_23Si_p22Al = new TH1I("Ex_23Si_p22Al","",500,0,20);
  ThetaCM_23Si_p22Al = new TH1I("ThetaCM_23Si_p22Al","",200,0,10);
  VCM_23Si_p22Al = new TH1I("VCM_23Si_p22Al","",400,7,15);
  Erel_p22Al_costhetaH = new TH2I("Erel_p22Al_costhetaH","",800,0,8,50,-1,1);
  Erel_p22Al_costhetaH->SetMinimum(1.0);
  p22Al_VCMvsErel = new TH2I("p22Al_VCMvsErel","",400,7,15,500,-5,15);
  p22Al_gammasADD = new TH1I("p22Al_gammasADD","",2048, 0, 8192);
  p22Al_gammasADD_tgate = new TH1I("p22Al_gammasADD_tgate","",2048, 0, 8192);
  p22Al_gammasADDvsErel = new TH2I("p22Al_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p22Al_gammasADDvsErel->SetMinimum(1.0);
  p22Al_gammasADDvsgammasADD = new TH2I("p22Al_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p22Al_gammasADDvsgammasADD->SetMinimum(1.0);
  Si23_pKE = new TH1I("Si23_pKE","",200, 0, 200);  

  //Si23 -> 2p + 21Mg
  dir2p21Mg->cd();
  Erel_23Si_2p21Mg = new TH1I("Erel_23Si_2p21Mg","",500,-5,15);
  Erel_23Si_2p21Mg_trans = new TH1I("Erel_23Si_2p21Mg_trans","",500,-5,15);
  Ex_23Si_2p21Mg = new TH1I("Ex_23Si_2p21Mg","",500,0,20);
  ThetaCM_23Si_2p21Mg = new TH1I("ThetaCM_23Si_2p221Mg","",200,0,10);
  VCM_23Si_2p21Mg = new TH1I("VCM_23Si_2p21Mg","",400,7,15);
  Erel_2p21Mg_costhetaH = new TH2I("Erel_2p21Mg_costhetaH","",800,0,8,50,-1,1);
  Erel_2p21Mg_costhetaH->SetMinimum(1.0);
  pp21Mg_VCMvsErel = new TH2I("pp21Mg_VCMvsErel","",400,7,15,500,-5,15);
  pp21Mg_gammasADD = new TH1I("pp21Mg_gammasADD","",2048, 0, 8192);
  pp21Mg_gammasADD_tgate = new TH1I("pp21Mg_gammasADD_tgate","",2048, 0, 8192);
  pp21Mg_gammasADDvsErel = new TH2I("pp21Mg_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  pp21Mg_gammasADDvsErel->SetMinimum(1.0);
  pp21Mg_gammasADDvsgammasADD = new TH2I("pp21Mg_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  pp21Mg_gammasADDvsgammasADD->SetMinimum(1.0);

  //Si24 -> p + 23Al
  dir24Si->cd();
  Erel_24Si_p23Al = new TH1I("Erel_24Si_p23Al","",500,-5,15);
  Erel_24Si_p23Al_trans = new TH1I("Erel_24Si_p23Al_trans","",500,-5,15);
  Ex_24Si_p23Al = new TH1I("Ex_24Si_p23Al","",500,0,20);
  ThetaCM_24Si_p23Al = new TH1I("ThetaCM_24Si_p23Al","",200,0,10);
  VCM_24Si_p23Al = new TH1I("VCM_24Si_p23Al","",400,7,15);
  Erel_p23Al_costhetaH = new TH2I("Erel_p23Al_costhetaH","",800,0,8,50,-1,1);
  Erel_p23Al_costhetaH->SetMinimum(1.0);
  p23Al_VCMvsErel = new TH2I("p23Al_VCMvsErel","",400,7,15,500,-5,15);
  p23Al_gammasADD = new TH1I("p23Al_gammasADD","",2048, 0, 8192);
  p23Al_gammasADD_tgate = new TH1I("p23Al_gammasADD_tgate","",2048, 0, 8192);
  p23Al_gammasADDvsErel = new TH2I("p23Al_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p23Al_gammasADDvsErel->SetMinimum(1.0);
  p23Al_gammasADDvsgammasADD = new TH2I("p23Al_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p23Al_gammasADDvsgammasADD->SetMinimum(1.0);

  //P23 -> 3p + 20Mg
  dir3p20Mg->cd();
  Erel_23P_3p20Mg = new TH1I("Erel_23P_3p20Mg","",500,-5,15);
  Erel_23P_3p20Mg_trans = new TH1I("Erel_23P_3p20Mg_trans","",500,-5,15);
  Ex_23P_3p20Mg = new TH1I("Ex_23P_3p20Mg","",500,0,20);
  ThetaCM_23P_3p20Mg = new TH1I("ThetaCM_23P_3p20Mg","",200,0,10);
  VCM_23P_3p20Mg = new TH1I("VCM_23P_3p20Mg","",400,7,15);
  Erel_3p20Mg_costhetaH = new TH2I("Erel_3p20Mg_costhetaH","",800,0,8,50,-1,1);
  Erel_p23Al_costhetaH->SetMinimum(1.0);
  ppp20Mg_VCMvsErel = new TH2I("ppp20Mg_VCMvsErel","",400,7,15,500,-5,15);
  ppp20Mg_gammasADD = new TH1I("ppp20Mg_gammasADD","",2048, 0, 8192);
  ppp20Mg_gammasADD_tgate = new TH1I("ppp20Mg_gammasADD_tgate","",2048, 0, 8192);
  ppp20Mg_gammasADDvsErel = new TH2I("ppp20Mg_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  ppp20Mg_gammasADDvsErel->SetMinimum(1.0);
  ppp20Mg_gammasADDvsgammasADD = new TH2I("ppp20Mg_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  ppp20Mg_gammasADDvsgammasADD->SetMinimum(1.0);

  //P23 -> p + 22Si
  dirp22Si->cd();
  Erel_23P_p22Si = new TH1I("Erel_23P_p22Si","",500,0,20);
  Erel_23P_p22Si_trans = new TH1I("Erel_23P_p22Si_trans","",500,0,20);
  Ex_23P_p22Si = new TH1I("Ex_23P_p22Si","",500,0,20);
  ThetaCM_23P_p22Si = new TH1I("ThetaCM_23P_p22Si","",200,0,10);
  VCM_23P_p22Si = new TH1I("VCM_23P_p22Si","",400,7,15);
  Erel_p22Si_costhetaH = new TH2I("Erel_p22Si_costhetaH","",800,0,8,50,-1,1);
  Erel_p22Si_costhetaH->SetMinimum(1.0);
  p22Si_VCMvsErel = new TH2I("p22Si_VCMvsErel","",400,7,15,500,-5,15);
  p22Si_gammasADD = new TH1I("p22Si_gammasADD","",2048, 0, 8192);
  p22Si_gammasADD_tgate = new TH1I("p22Si_gammasADD_tgate","",2048, 0, 8192);
  p22Si_gammasADDvsErel = new TH2I("p22Si_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p22Si_gammasADDvsErel->SetMinimum(1.0);
  p22Si_gammasADDvsgammasADD = new TH2I("p22Si_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p22Si_gammasADDvsgammasADD->SetMinimum(1.0);

  //P24 -> p + 23Si
  dir24P->cd();
  Erel_24P_p23Si = new TH1I("Erel_24P_p23Si","",500,-5,15);
  Erel_24P_p23Si_trans = new TH1I("Erel_24P_p23Si_trans","",500,-5,15);
  Ex_24P_p23Si = new TH1I("Ex_24P_p23Si","",500,0,20);
  ThetaCM_24P_p23Si = new TH1I("ThetaCM_24P_p23Si","",200,0,10);
  VCM_24P_p23Si = new TH1I("VCM_24P_p23Si","",400,7,15);
  Erel_p23Si_costhetaH = new TH2I("Erel_p23Si_costhetaH","",800,0,8,50,-1,1);
  Erel_p23Si_costhetaH->SetMinimum(1.0);
  p23Si_VCMvsErel = new TH2I("p23Si_VCMvsErel","",400,7,15,500,-5,15);
  p23Si_gammasADD = new TH1I("p23Si_gammasADD","",2048, 0, 8192);
  p23Si_gammasADD_tgate = new TH1I("p23Si_gammasADD_tgate","",2048, 0, 8192);
  p23Si_gammasADDvsErel = new TH2I("p23Si_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  p23Si_gammasADDvsErel->SetMinimum(1.0);
  p23Si_gammasADDvsgammasADD = new TH2I("p23Si_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p23Si_gammasADDvsgammasADD->SetMinimum(1.0);

  //Al21 -> 3p + 18Ne
  dir3p18Ne->cd();
  Erel_21Al_3p18Ne = new TH1I("Erel_21Al_3p18Ne","",500,-5,15);
  Erel_21Al_3p18Ne_trans = new TH1I("Erel_21Al_3p18Ne_trans","",500,-5,15);
  Ex_21Al_3p18Ne = new TH1I("Ex_21Al_3p18Ne","",500,0,20);
  ThetaCM_21Al_3p18Ne = new TH1I("ThetaCM_21Al_3p18Ne","",200,0,10);
  VCM_21Al_3p18Ne = new TH1I("VCM_21Al_3p18Ne","",400,7,15);
  Erel_3p18Ne_costhetaH = new TH2I("Erel_3p18Ne_costhetaH","",800,0,8,50,-1,1);
  Erel_3p18Ne_costhetaH->SetMinimum(1.0);
  ppp18Ne_VCMvsErel = new TH2I("ppp18Ne_VCMvsErel","",400,7,15,500,-5,15);
  ppp18Ne_gammasADD = new TH1I("ppp18Ne_gammasADD","",2048, 0, 8192);
  ppp18Ne_gammasADD_tgate = new TH1I("ppp18Ne_gammasADD_tgate","",2048, 0, 8192);
  ppp18Ne_gammasADDvsErel = new TH2I("ppp18Ne_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  ppp18Ne_gammasADDvsErel->SetMinimum(1.0);
  ppp18Ne_gammasADDvsgammasADD = new TH2I("ppp18Ne_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  ppp18Ne_gammasADDvsgammasADD->SetMinimum(1.0);

  //Al20 -> 3p + 17Ne
  dir3p17Ne->cd();
  Erel_20Al_3p17Ne = new TH1I("Erel_20Al_3p17Ne","",500,-5,15);
  Erel_20Al_3p17Ne_trans = new TH1I("Erel_20Al_3p17Ne_trans","",500,-5,15);
  Ex_20Al_3p17Ne = new TH1I("Ex_20Al_3p17Ne","",500,0,20);
  ThetaCM_20Al_3p17Ne = new TH1I("ThetaCM_20Al_3p17Ne","",200,0,10);
  VCM_20Al_3p17Ne = new TH1I("VCM_20Al_3p17Ne","",400,7,15);
  Erel_3p17Ne_costhetaH = new TH2I("Erel_3p17Ne_costhetaH","",800,0,8,50,-1,1);
  Erel_3p17Ne_costhetaH->SetMinimum(1.0);
  ppp17Ne_VCMvsErel = new TH2I("ppp17Ne_VCMvsErel","",400,7,15,500,-5,15);
  ppp17Ne_gammasADD = new TH1I("ppp17Ne_gammasADD","",2048, 0, 8192);
  ppp17Ne_gammasADD_tgate = new TH1I("ppp17Ne_gammasADD_tgate","",2048, 0, 8192);
  ppp17Ne_gammasADDvsErel = new TH2I("ppp17Ne_gammasADDvsErel","",500,-5,15,2048, 0, 8192);
  ppp17Ne_gammasADDvsErel->SetMinimum(1.0);
  ppp17Ne_gammasADDvsgammasADD = new TH2I("ppp17Ne_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  ppp17Ne_gammasADDvsgammasADD->SetMinimum(1.0);


}

//Needed to clear Janus instances
void histo_sort::clear() {
  event.clear();
	red.clear();
	blue.clear();
}

//*********************************************
void histo_sort::write()
{

  file->Write();

  cout << "histo written" << endl;
  file->Close();
  cout << "here!!!!" << endl;
  /*
    for (int i=0;i<Ntele;i++)
    {
    delete red[i];
    }
    delete [] red;
  */
}
