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

  //directory for DeltaE-E plots
  dirDEEplots = new TDirectoryFile("DEEplots","DEEplots");
  dirhitmaps = new TDirectoryFile("hitmaps","hitmaps");

  //directories for Ceasar
  dirC = new TDirectoryFile("ceasar","ceasar");
  dirEnergy = dirC->mkdir("raw","energyC");
  dirTime = dirC->mkdir("timeC","timeC");
  dirTimeS800 = dirC->mkdir("time&S800","time&S800");
  dirEcal = dirC->mkdir("cal","cal");
  dirRSum = dirC->mkdir("RSum","Rsum");
  dirCSum = dirC->mkdir("CSum","Csum");
  dirCa36 = dirC->mkdir("Ca36","Ca36");

  //fiber
  dirFiber = new TDirectoryFile("Fiber","Fiber");
  dir1dBoard0TOA_R = dirFiber->mkdir("Board0TOA_R","Board0TOA_R");
  dir1dBoard0TOT_R = dirFiber->mkdir("Board0TOT_R","Board0TOT_R");
  dir1dBoard1TOA_R = dirFiber->mkdir("Board1TOA_R","Board1TOA_R");
  dir1dBoard1TOT_R = dirFiber->mkdir("Board1TOT_R","Board0TOT_R");
  dirHitMap = dirFiber->mkdir("FiberHitMap","FiberHitMap");

  //Inv mass directories
  dirInvMass = new TDirectoryFile("InvMass", "InvMass");
  dir1H = dirInvMass->mkdir("1H","1H");
  dir21Al = dirInvMass->mkdir("21Al","21Al");
  dir22Si = dirInvMass->mkdir("22Si","22Si");
  dir23P = dirInvMass->mkdir("23P","23P");
  dir3p20Mg = dir23P->mkdir("3p20Mg","3p20Mg");
  dirp22Si = dir23P->mkdir("22Si","22Si");


  //for calibrated spectra
  int Nbin = 5000;
  float Ecal_Emax = 50.0;
  float CsI_Emax = 120;

  
  NCsI = 16;
  Ntele = 4;
  Nstrip = 32; //strips per silicon
  Nceasar = 192;
  
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
  sumCsIE_cal = new TH2I("sumCsIE_cal","",16,0,16,Nbin,0,CsI_Emax);
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
  sumCsITime_R = new TH2I("sumCsITime_R","",16,0,16,500,-1500,1000);
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
      FrontE_R[board_i][chan_i] = new TH1I(name.str().c_str(),"",2048,0,8192);

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
      BackE_R[board_i][chan_i] = new TH1I(name.str().c_str(),"",2048,0,8192);

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
    CsI_Energy_cal[ichan] = new TH1I(name.str().c_str(),"",600,0,60);

    dir1dCsI_Time->cd();
    name.str("");
    name << "CsI_Time_" << ichan << "R_unmatched";
    CsI_Time_R_um[ichan] = new TH1I(name.str().c_str(),"",500,-1500,1000);

    name.str("");
    name << "CsI_Time_" << ichan << "R_matched";
    CsI_Time_R[ichan] = new TH1I(name.str().c_str(),"",500,-1500,1000);

    name.str("");
    name << "CsI_Time_" << ichan << "cal";
    CsI_Time_cal[ichan] = new TH1I(name.str().c_str(),"",2000,-1000,1000);
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
      name.str("");
      name << "DEE_CsI" << tele << "_" << icsi;   
      DEE_CsI[tele][icsi] = new TH2I(name.str().c_str(),"",4096,0,4096,400,0,100);

      name.str("");
      name << "DEE_CsI_0deg" << tele << "_" << icsi;   
      DEE_CsI_0deg[tele][icsi] = new TH2I(name.str().c_str(),"",4096,0,4096,400,0,100);

      name.str("");
      name << "DEE_CsI_lowgain" << tele << "_" << icsi;   
      DEE_CsI_lowgain[tele][icsi] = new TH2I(name.str().c_str(),"",4096,0,4096,400,0,100);

      name.str("");
      name << "timediff" << tele << "_" << icsi;    
      timediff_CsI[tele][icsi] = new TH1I(name.str().c_str(),"",1000,-2000,2000);
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
     name << "GFrontStrip_CsI_" << i;
     GFrontStrip_CsI[i] = new TH2I(name.str().c_str(),"",4,-0.5,3.5,32,-0.5,31.5);

     name.str("");
     name << "GBackStrip_CsI_" << i;
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

	// singles histograms
  int fnum = 64;
  //Make 1d plots
  for (int i=0;i<fnum;i++)
  {
    dir1dBoard0TOA_R->cd();
    name.str("");
    name << "Board0_TOA_" << i << "_R";
    B0TOA_R[i] = new TH1I(name.str().c_str(),"",1024,0,4095);

    dir1dBoard0TOT_R->cd();
    name.str("");
    name << "Board0_TOT_" << i << "_R";
    B0TOT_R[i] = new TH1I(name.str().c_str(),"",1024,0,4095);

    dir1dBoard1TOA_R->cd();
    name.str("");
    name << "Board1_TOA_" << i << "_R";
    B1TOA_R[i] = new TH1I(name.str().c_str(),"",1024,0,4095);

    dir1dBoard1TOT_R->cd();
    name.str("");
    name << "Board1_TOT_" << i << "_R";
    B1TOT_R[i] = new TH1I(name.str().c_str(),"",1024,0,4095);
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

  dirS800 = new TDirectoryFile("S800","S800");
  dirS800Raw = dirS800->mkdir("Raw","Raw");
  dirS800Cal = dirS800->mkdir("Cal","Cal");
  dirS800PID = dirS800->mkdir("PID","PID");
  dirS800Trig = dirS800->mkdir("Trig","Trig");
  dirS800veldist = dirS800->mkdir("veldist","veldist");
  

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

  //CAESAR
  ECeasar = new TH1I*[Nceasar];
  TCeasar = new TH1I*[Nceasar];
  TCeasarS800 = new TH1I*[Nceasar];
  ECCeasar = new TH1I*[Nceasar];
  Ca36gamma_nodopp = new TH1I*[Nceasar];
  Ca36gamma = new TH1I*[Nceasar];

  dirCSum->cd();
  TECeasar = new TH1I("TECeasar","",512,0,10);
  TECeasar->GetXaxis()->SetTitle("E_{#gamma} [MeV]");

  TECeasar_difbin = new TH1I("TECeasar_difbin","",300,0,6);
  TECeasar_difbin->GetXaxis()->SetTitle("E_{#gamma} [MeV]");

  Egated = new TH1I("Egated","",512,0,10);
  Egated->GetXaxis()->SetTitle("E_{#gamma} [MeV]");



  TvsEchn42 = new TH2I("TvsEchn42","",1024,0,4095, 500,-1000,1000);
  TvsEchn42->GetXaxis()->SetTitle("E_{#gamma} [chn]");
  TvsEchn42->GetYaxis()->SetTitle("#Delta T [ns]");

  TvsEchn48 = new TH2I("TvsEchn48","",1024,0,4095, 500,-1000,1000);
  TvsEchn48->GetXaxis()->SetTitle("E_{#gamma} [chn]");
  TvsEchn48->GetYaxis()->SetTitle("#Delta T [ns]");

  TvsEgated = new TH2I("TvsEgated","",300,0,6, 500,-1000,1000);
  TvsEgated->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  TvsEgated->GetYaxis()->SetTitle("#Delta T [ns]");

  TEC_Dop = new TH1I("TEC_Dop","",4095,0,10);
  TEC_Dop->GetXaxis()->SetTitle("E_{#gamma} [MeV]");
  TEC_Dop->SetTitle("Dopler corrected");

  map_iring_iring = new TH2I("Map_iring_iring","",10,-0.5,9.5,10,-0.5,9.5);
  map_idet_idet = new TH2I("Map_idet_idet","",192,0,192,192,0,192);
  map_DataEC_id_id = new TH2I("map_DataEC_id_id","",192,0,192,192,0,192);
  map_DataTC_id_id = new TH2I("map_DataTC_id_id","",192,0,192,192,0,192);

  delta_phi_gamma = new TH1I("delta_phi_gamma","",180,0,180);
  
  CEMult = new TH1I("CEMult","",10,-0.5,9.5);
  CTMult = new TH1I("CTMult","",10,-0.5,9.5);
  CETMult = new TH1I("CETMult","",10,-0.5,9.5);
  ETOF_Ceasar = new TH2I("ETOF_Ceasar","",250,0,1000,1024,-0.5,9.5);
  TCeasarRaw1 = new TH1I("TCeasarRaw1","",1000,-3000,3000);
  TCeasarRaw2 = new TH1I("TCeasarRaw2","",1000,-3000,3000);
  TCeasarSum_gated = new TH2I("TCeasarSum_gated","",192,-0.5,191.5,1000,-3000,3000);


  TCeasarRawS800 = new TH1I("TCeasarRawS800","",1000,-3000,3000);
  TCeasarCal1 = new TH1I("TCeasarCal1","",1000,-3000,3000);
  TCeasarCal2 = new TH1I("TCeasarCal2","",1000,-3000,3000);

  dirTime->cd();
  TCeasarRawSum = new TH2I("TCeasarRawSum","",192,-0.5,191.5,1000,-3000,3000);
  TCeasarMatchedSum = new TH2I("TCeasarMatchedSum","",192,-0.5,191.5,1000,-3000,3000);

  dirEnergy->cd();
  ECeasarRawSum = new TH2I("ECeasarRawSum","",192,-0.5,191.5,1024,0,4095);
  ECeasarCalSum = new TH2I("ECeasarCalSum","",192,-0.5,191.5,300,0,6);
  ECeasarDopSum = new TH2I("ECeasarDopSum","",192,-0.5,191.5,300,0,6);
  ECeasarMatchedRawSum = new TH2I("ECeasarMatchedRawSum","",192,-0.5,191.5,1024,0,4095);
  ECeasarMatchedCalSum = new TH2I("ECeasarMatchedCalSum","",192,-0.5,191.5,300,0,6);

  ECeasarRawSum_idet = new TH2I("ECeasarRawSum_idet","",192,-0.5,191.5,1024,0,4095);
  ECeasarCalSum_idet = new TH2I("ECeasarCalSum_idet","",192,-0.5,191.5,300,0,6);
  ECeasarDopSum_idet = new TH2I("ECeasarDopSum_idet","",192,-0.5,191.5,300,0,6);
  ECeasarMatchedRawSum_idet = new TH2I("ECeasarMatchedRawSum_idet","",192,-0.5,191.5,1024,0,4095);
  ECeasarMatchedCalSum_idet = new TH2I("ECeasarMatchedCalSum_idet","",192,-0.5,191.5,300,0,6);
  TCeasarMatchedSum_idet = new TH2I("TCeasarMatchedSum_idet","",192,-0.5,191.5,1000,-3000,3000);


  for(int ic = 0;ic<Nceasar;ic++)
  {
    name.str("");
    if (ic < 10)      name << "EC00" << ic;
    else if (ic < 100) name << "EC0" << ic;
    else name << "EC" << ic;
      
    dirEnergy->cd();
    ECeasar[ic] = new TH1I(name.str().c_str(),"",1024,0,4095);
    ECeasar[ic]->GetXaxis()->SetTitle("Ceasar Energy [channel]");
      
    name.str("");
    if (ic < 10)      name << "ECcal00" << ic;
    else if (ic < 100) name << "ECcal0" << ic;
    else name << "ECcal" << ic;

    dirEcal->cd();
    ECCeasar[ic] = new TH1I(name.str().c_str(),"",500,0,10);
    ECCeasar[ic]->GetXaxis()->SetTitle("Ceasar Energy [MeV]");

    name.str("");
    if (ic < 10)      name << "TC00" << ic;
    else if (ic < 100) name << "TC0" << ic;
    else name << "TC" << ic;

    dirTime->cd();
    TCeasar[ic] = new TH1I(name.str().c_str(),"",1000,-3000,3000);
    TCeasar[ic]->GetXaxis()->SetTitle("Ceasar Time [ns]");

    name << "&s800";
    dirTimeS800->cd();
    TCeasarS800[ic] = new TH1I(name.str().c_str(),"",1000,-3000,3000);
    TCeasarS800[ic]->GetXaxis()->SetTitle("Ceasar Time [ns]");

    name.str("");
    if (ic < 10)      name << "Ca36gamma00" << ic;
    else if (ic < 100) name << "Ca36gamma0" << ic;
    else name << "Ca36gamma" << ic;
    dirCa36->cd();
    Ca36gamma[ic] = new TH1I(name.str().c_str(),"",300,0,6);
    Ca36gamma[ic]->GetXaxis()->SetTitle("Ceasar Energy [MeV]");

    name << "_nodopp";
    dirCa36->cd();
    Ca36gamma_nodopp[ic] = new TH1I(name.str().c_str(),"",300,0,6);
    Ca36gamma_nodopp[ic]->GetXaxis()->SetTitle("Ceasar Energy [MeV]");


  }

  //S800 Spectra

  dirS800Raw->cd();
  Te1up = new TH1I("Te1up","",501,-500,500);
  Te1down = new TH1I("Te1down","",501,-500,500);
  Tobj = new TH1I("Tobj","",501,-500,500);
  Tobj_mult = new TH1I("Tobj_mult","",10,0,10);

  ICSummary = new TH2I("ICSummary","",16,-0.5,15.5,1024,0,4095);

  CRDC1Summary = new TH2I("CRDC1Summary","",224,-0.5,223.5,1200,0,1200);
  CRDC2Summary = new TH2I("CRDC2Summary","",224,-0.5,223.5,1200,0,1200);
  CRDC1Summary_cal = new TH2I("CRDC1Summary_cal","",224,-0.5,223.5,1200,0,1200);
  CRDC2Summary_cal = new TH2I("CRDC2Summary_cal","",224,-0.5,223.5,1200,0,1200);
  CRDC1raw = new TH1I("CRCD1raw","",2000,0,2000);
  CRDC2raw = new TH1I("CRCD2raw","",2000,0,2000);
  CRDC1PadMult = new TH1I("CRDC1PadMult","",224,-0.5,223.5);  
  CRDC2PadMult = new TH1I("CRDC2PadMult","",224,-0.5,223.5);  
  CRDC1Tac = new TH1I("CRDC1Tac","",2000,0,2000);
  CRDC2Tac = new TH1I("CRDC2Tac","",2000,0,2000);
  CRDC1AnodevsTac = new TH2I("CRDC1AnodevsTac","",512,0,4096,512,0,4096);
  CRDC2AnodevsTac = new TH2I("CRDC2AnodevsTac","",512,0,4096,512,0,4096);
  CRDC1XrawY = new TH2I("CRDC1XrawY","",1024,0,200,512,0,1000);
  CRDC2XrawY = new TH2I("CRDC2XrawY","",1024,0,200,512,0,1000);

  dirS800Cal->cd();
  CRDC1X = new TH1I("CRDC1X","",512,-300,300);
  CRDC1Y = new TH1I("CRDC1Y","",512,-100,100);
  CRDC2X = new TH1I("CRDC2X","",512,-300,300);
  CRDC2Y = new TH1I("CRDC2Y","",512,-100,100);
  CRDC1XCalY = new TH2I("CRDC1XCalY","",1024,-300,300,512,-100,100);
  CRDC2XCalY = new TH2I("CRDC2XCalY","",1024,-300,300,512,-100,100);  
  CRDC1XY = new TH2I("CRDC1XY","",512,-300,300,512,-100,100);
  CRDC2XY = new TH2I("CRDC2XY","",512,-300,300,512,-100,100);
  atavsbta = new TH2I("atavsbta","",100,-5,5,100,-5,5);
  ThetavsPhi = new TH2I("ThetavsPhi","",100,0,100,720,-180,180);
  ThetavsPhi_fiber = new TH2I("ThetavsPhi_fiber","",50,0,5,720,-180,180);
  ThetaFibervsThetaS800 = new TH2I("ThetaFibervsThetaS800","",200,-5,5,200,-5,5);
  PhiFibervsPhiS800 = new TH2I("PhiFibervsPhiS800","",720,-180,180,720,-180,180);
  CRDC1X_K35 = new TH1I("CRDC1X_K35","",512,-300,300);
  CRDC1X_K36 = new TH1I("CRDC1X_K36","",512,-300,300);
  //needed for CRDC gain matching
  CRDC1_Ca36_Ca36 = new TH2I("CRDC1_Ca36_Ca36","",224,-0.5,223.5,1200,0,1200);
  CRDC2_Ca36_Ca36 = new TH2I("CRDC2_Ca36_Ca36","",224,-0.5,223.5,1200,0,1200);
  CRDC1_Ca36_K35 = new TH2I("CRDC1_Ca36_K35","",224,-0.5,223.5,1200,0,1200);
  CRDC2_Ca36_K35 = new TH2I("CRDC2_Ca36_K35","",224,-0.5,223.5,1200,0,1200);
  CRDC1_K35_K35 = new TH2I("CRDC1_K35_K35","",224,-0.5,223.5,1200,0,1200);
  CRDC2_K35_K35 = new TH2I("CRDC2_K35_K35","",224,-0.5,223.5,1200,0,1200);
  CRDC1_K35_K36 = new TH2I("CRDC1_K35_K36","",224,-0.5,223.5,1200,0,1200);
  CRDC2_K35_K36 = new TH2I("CRDC2_K35_K36","",224,-0.5,223.5,1200,0,1200);
  CRDC1_Ca36_Ar33 = new TH2I("CRDC1_Ca36_Ar33","",224,-0.5,223.5,1200,0,1200);
  CRDC2_Ca36_Ar33 = new TH2I("CRDC2_Ca36_Ar33","",224,-0.5,223.5,1200,0,1200);
  CRDC1_Ca36_Ar34 = new TH2I("CRDC1_Ca36_Ar34","",224,-0.5,223.5,1200,0,1200);
  CRDC2_Ca36_Ar34 = new TH2I("CRDC2_Ca36_Ar34","",224,-0.5,223.5,1200,0,1200);
  CRDC1_K35_Ar34 = new TH2I("CRDC1_K35_Ar34","",224,-0.5,223.5,1200,0,1200);
  CRDC2_K35_Ar34 = new TH2I("CRDC2_K35_Ar34","",224,-0.5,223.5,1200,0,1200);
  CRDC1_K35_Ar35 = new TH2I("CRDC1_K35_Ar35","",224,-0.5,223.5,1200,0,1200);
  CRDC2_K35_Ar35 = new TH2I("CRDC2_K35_Ar35","",224,-0.5,223.5,1200,0,1200);
  CRDC1_Ca36_Cl32 = new TH2I("CRDC1_Ca36_Cl32","",224,-0.5,223.5,1200,0,1200);
  CRDC2_Ca36_Cl32 = new TH2I("CRDC2_Ca36_Cl32","",224,-0.5,223.5,1200,0,1200);

  dirS800PID->cd();  
  ObjvsXFP = new TH2I("ObjvsXFP","",600,-300,0,600,0,300);
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
  
  ObjvsICsum = new TH2I("ObjvsICsum","",1200,-200,-50,2048,0,4095);
  ObjvsICsum->GetXaxis()->SetTitle("T Object");
  ObjvsICsum->GetYaxis()->SetTitle("IC Sum");

  ObjvsICsum_corr = new TH2I("ObjvsICsum_corr","",1200,-200,-50,2048,0,4095);
  ObjvsICsum_corr->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_corr->GetYaxis()->SetTitle("IC Sum");


  ObjvsICsum_Ca37 = new TH2I("ObjvsICsum_Ca37","",1200,-200,-50,2048,0,4095);
  ObjvsICsum_Ca37->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_Ca37->GetYaxis()->SetTitle("IC Sum");
  ObjUncvsICsum_Ca37 = new TH2I("ObjUncvsICsum_Ca37","",800,-200,-50,512,0,4095);
  ObjUncvsICsum_Ca37->GetXaxis()->SetTitle("T Object (Uncorrected)");
  ObjUncvsICsum_Ca37->GetYaxis()->SetTitle("IC Sum");

  Timing1vsICsum_Ca37 = new TH2I("Timing1vsICsum_Ca37","",1200,-200,-50,2048,0,4095);
  Timing1vsICsum_Ca37->GetXaxis()->SetTitle("T diff");
  Timing1vsICsum_Ca37->GetYaxis()->SetTitle("IC Sum");
  Timing2vsICsum_Ca37 = new TH2I("Timing2vsICsum_Ca37","",1200,-200,-50,2048,0,4095);
  Timing2vsICsum_Ca37->GetXaxis()->SetTitle("T diff");
  Timing2vsICsum_Ca37->GetYaxis()->SetTitle("IC Sum");
  XFPvsICsum_Ca37 = new TH2I("XFPvsICsum_Ca37","",600,0,300,2048,0,4095);
  XFPvsICsum_Ca37->GetXaxis()->SetTitle("XFP time");
  XFPvsICsum_Ca37->GetYaxis()->SetTitle("IC Sum");




  ObjvsICsum_K36 = new TH2I("ObjvsICsum_K36","",1200,-200,-50,2048,0,4095);
  ObjvsICsum_K36->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_K36->GetYaxis()->SetTitle("IC Sum");
  ObjvsICsum_Ar35 = new TH2I("ObjvsICsum_Ar35","",1200,-200,-50,2048,0,4095);
  ObjvsICsum_Ar35->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_Ar35->GetYaxis()->SetTitle("IC Sum");
  ObjvsICsum_Cl34 = new TH2I("ObjvsICsum_Cl34","",1200,-200,-50,2048,0,4095);
  ObjvsICsum_Cl34->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_Cl34->GetYaxis()->SetTitle("IC Sum");




  ObjvsICsum_wProt = new TH2I("ObjvsICsum_wProt","",1200,-200,-50,2048,0,4095);
  ObjvsICsum_wProt->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_wProt->GetYaxis()->SetTitle("IC Sum");

  ObjvsICsum_wProt_wFiber = new TH2I("ObjvsICsum_wProt_wFiber","",1200,-200,-50,2048,0,4095);
  ObjvsICsum_wProt_wFiber->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_wProt_wFiber->GetYaxis()->SetTitle("IC Sum");


  ObjvsICsum_wProt2 = new TH2I("ObjvsICsum_wProt2","",1200,-200,-50,2048,0,4095);
  ObjvsICsum_wProt2->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_wProt2->GetYaxis()->SetTitle("IC Sum");

  ObjvsICsum_K36_wProt = new TH2I("ObjvsICsum_K36_wProt","",1200,-200,-50,2048,0,4095);
  ObjvsICsum_K36_wProt->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_K36_wProt->GetYaxis()->SetTitle("IC Sum");

  ObjvsICsum_Ar35_wProt = new TH2I("ObjvsICsum_Ar35_wProt","",1200,-200,-50,2048,0,4095);
  ObjvsICsum_Ar35_wProt->GetXaxis()->SetTitle("T Object");
  ObjvsICsum_Ar35_wProt->GetYaxis()->SetTitle("IC Sum");

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
  Vlab_Ca38 = new TH1I("Vlab_Ca38","",100,0.1,0.5);
  Ppar_Ca38 = new TH1I("Ppar_Ca38","",300,12000,15000);
  Ptra_Ca38 = new TH1I("Ptra_Ca38","",400,0,800);
  Vlab_Ca37 = new TH1I("Vlab_Ca37","",100,0.1,0.5);
  Ppar_Ca37 = new TH1I("Ppar_Ca37","",300,12000,15000);
  Ptra_Ca37 = new TH1I("Ptra_Ca37","",400,0,800);
  Elab_Ca37 = new TH1I("Elab_Ca37","",1000,1000,4000);
  Etar_Ca37 = new TH1I("Etar_Ca37","",1000,1000,4000);
  Etarfib_Ca37 = new TH1I("Etarfib_Ca37","",1000,1000,4000);
  Efib_Ca37 = new TH1I("Efib_Ca37","",1000,1000,4000);

  Vlab_Ca35 = new TH1I("Vlab_Ca35","",100,0.1,0.5);
  Vlab_Ca36 = new TH1I("Vlab_Ca36","",100,0.1,0.5);
  Vlab_Ca36_wgamma = new TH1I("Vlab_Ca36_wgamma","",100,0.1,0.5);
  Vlab_K35 = new TH1I("Vlab_K35","",100,0.1,0.5);
  Vlab_K36 = new TH1I("Vlab_K36","",100,0.1,0.5);
  Vlab_Ar32 = new TH1I("Vlab_Ar32","",100,0.1,0.5);
  Vlab_Ar33 = new TH1I("Vlab_Ar33","",100,0.1,0.5);
  Vlab_Ar34 = new TH1I("Vlab_Ar34","",100,0.1,0.5);
  Vlab_Ar35 = new TH1I("Vlab_Ar35","",100,0.1,0.5);
  Vlab_Cl31 = new TH1I("Vlab_Cl31","",100,0.1,0.5);
  Vlab_Cl32 = new TH1I("Vlab_Cl32","",100,0.1,0.5);
  Vlab_Cl33 = new TH1I("Vlab_Cl33","",100,0.1,0.5);
  Vlab_Cl34 = new TH1I("Vlab_Cl34","",100,0.1,0.5);


  //********************************************************
  //Correlation hists
  //********************************************************
  //Proton Directory
  dir1H->cd();
  protonKE = new TH1I("protonKE","",300,0,150); //500 keV bins

  //Al21 -> p + Mg20
  dir21Al->cd();
  Erel_21Al_p20Mg = new TH1I("Erel_21Al_p20Mg","",100,0,20);
  Ex_21Al_p20Mg = new TH1I("Ex_21Al_p20Mg","",375,0,15);
  ThetaCM_21Al_p20Mg = new TH1I("ThetaCM_21Al_p20Mg","",200,0,10);
  VCM_21Al_p20Mg = new TH1I("VCM_21Al_p20Mg","",400,2,8);
  Erel_p20Mg_costhetaH = new TH2I("Erel_p20Mg_costhetaH","",100,0,8,25,-1,1);
  p20Mg_VCMvsErel = new TH2I("p20Mg_VCMvsErel","",400,2,8,100,-5,15);

  //Si22 -> 2p + 20Mg
  dir22Si->cd();
  Erel_22Si_2p20Mg = new TH1I("Erel_22Si_2p20Mg","",100,0,20);
  Ex_22Si_2p20Mg = new TH1I("Ex_22Si_2p20Mg","",375,0,15);
  ThetaCM_22Si_2p20Mg = new TH1I("ThetaCM_22Si_2p20Mg","",200,0,10);
  VCM_22Si_2p20Mg = new TH1I("VCM_22Si_2p20Mg","",400,2,8);
  Erel_2p20Mg_costhetaH = new TH2I("Erel_2p20Mg_costhetaH","",100,0,8,25,-1,1);
  pp20Mg_VCMvsErel = new TH2I("pp20Mg_VCMvsErel","",400,2,8,100,-5,15);

  //P23 -> 3p + 20Mg
  dir3p20Mg->cd();
  Erel_23P_3p20Mg = new TH1I("Erel_23P_3p20Mg","",100,0,20);
  Ex_23P_3p20Mg = new TH1I("Ex_23P_3p20Mg","",375,0,15);
  ThetaCM_23P_3p20Mg = new TH1I("ThetaCM_23P_3p20Mg","",200,0,10);
  VCM_23P_3p20Mg = new TH1I("VCM_23P_3p20Mg","",400,2,8);
  Erel_3p20Mg_costhetaH = new TH2I("Erel_3p20Mg_costhetaH","",100,0,8,25,-1,1);
  ppp20Mg_VCMvsErel = new TH2I("ppp20Mg_VCMvsErel","",400,2,8,100,-5,15);

  //P23 -> p + 22Si
  dirp22Si->cd();
  Erel_23P_p22Si = new TH1I("Erel_23P_p22Si","",100,0,20);
  Ex_23P_p22Si = new TH1I("Ex_23P_p22Si","",375,0,15);
  ThetaCM_23P_p22Si = new TH1I("ThetaCM_23P_p22Si","",200,0,10);
  VCM_23P_p22Si = new TH1I("VCM_23P_p22Si","",400,2,8);
  Erel_p22Si_costhetaH = new TH2I("Erel_p22Si_costhetaH","",100,0,8,25,-1,1);
  p22Si_VCMvsErel = new TH2I("p22Si_VCMvsErel","",400,2,8,100,-5,15);


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
  /*
    for (int i=0;i<Ntele;i++)
    {
    delete red[i];
    }
    delete [] red;
  */
}
