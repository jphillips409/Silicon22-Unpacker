#include "pixels_si22.h"
#include "loss.h"
#include "TFile.h"
#include "TH1I.h"
#include "TTree.h"
#include "kin.h"
#include "calibrate.h"
#include "TH2I.h"
#include "constants.h"
//#include "calCorrection.h"
//#include "fineCal.h"
#include <sstream>
#include <cmath>
#include <limits>

//from p+17F 
double suppress(double E)
{
   const double c0 = .71969;
   const double c1 = 1.067;
   return ((1.+exp(c0/c1))/(1.+exp(-(E-c0)/c1))-1.)*exp(-c0/c1);
}


int main()
{
  const int nn = 200;
  save savep[nn];
  save save19[nn];
  int nwrite = 0;
  int nsave = 0.;

  bool B17[nn];
  bool B15[nn];


  kin Kin;

  //fineCal FineCal;

  calibrate FrontEcal(4, 32, "/home/Silicon22/Silicon22sort/cal2/GobbiFrontEcal.dat", 1);
  calibrate BackEcal(4, 32, "/home/Silicon22/Silicon22sort/cal2/GobbiBackEcal.dat", 1);
  calibrate CsIEcal(4, 4, "/home/Silicon22/Silicon22sort/cal/calMacros/CsI_cal_out.dat", 1);
  CLoss lossp_Be("Hydrogen_Be.loss",Mass_p/m0);
  CLoss lossp_Al("Hydrogen_Al.loss",Mass_p/m0);

  double Bethick = 97.77768; // ~0.5 mm
  double Althick = 1351.; // 5 mm
  double gobbidist = 55.2;

  pixels_si22 Pix(gobbidist,Bethick);

  TFile file("p17F_save.root");
  //TFile file("p17F_pospivot.root");
  TTree * tree = (TTree*)file.Get("T");
  
  int N = tree->GetEntries();
  cout << "Number of ENtries = " << N << endl;


  vel ee1;
  vel ee2;
  vel ee3;
  vel CM;

  float MM1;
  float MM2;
  float M1[3];
  float M2[3];

  int id2;
  int id1;
  int ifront2;
  int iback2;
  int ifront1;
  int iback1;

  int itele1;

  float et1;
  float et2;
  float denergy1;
  float energy_p1;

  float theta1;
  float theta2;
  float theta_s800;
  float phi1;
  float phi2;
  float phi_s800;

  float denergy2;
  float energy_p2;
  float denergyR;
  float energyR;

  float timeCsI1;
  float timeCsI2;

  int Ngamma;
  float Egamma[15];
  float Tgamma[15];
  int Chgamma[15];

  float Erel;
  float Ex;
  float Vcm;
  float thetaCM;
  float cos_thetaH;

  int runnum;
  int beamZ;  



  tree->SetBranchAddress("M1",M1);
  tree->SetBranchAddress("M2",M2);
  tree->SetBranchAddress("et1",&et1);
  tree->SetBranchAddress("et2",&et2);
  tree->SetBranchAddress("id2",&id2);
  tree->SetBranchAddress("id1",&id1);
  tree->SetBranchAddress("id1",&id1);
  tree->SetBranchAddress("itele1",&itele1);
  tree->SetBranchAddress("ifront2",&ifront2);
  tree->SetBranchAddress("iback2",&iback2);
  tree->SetBranchAddress("ifront1",&ifront1);
  tree->SetBranchAddress("iback1",&iback1);


  //tree->SetBranchAddress("denergy1",&denergy1);
  tree->SetBranchAddress("energy_p1",&energy_p1);
  tree->SetBranchAddress("denergy1_R",&denergyR);
  tree->SetBranchAddress("energy_p1_R",&energyR);


  //tree->SetBranchAddress("denergy2",&denergy2);
  tree->SetBranchAddress("energy_p2",&energy_p2);

  tree->SetBranchAddress("theta1",&theta1);
  tree->SetBranchAddress("theta2",&theta2);
  tree->SetBranchAddress("theta_s800",&theta_s800);
  tree->SetBranchAddress("phi1",&phi1);
  tree->SetBranchAddress("phi2",&phi2);
  tree->SetBranchAddress("phi_s800",&phi_s800);

  tree->SetBranchAddress("Ngamma",&Ngamma);
  tree->SetBranchAddress("Egamma",Egamma);
  tree->SetBranchAddress("Tgamma",Tgamma);
  tree->SetBranchAddress("Chgamma",Chgamma);

  tree->SetBranchAddress("Erel",&Erel);
  tree->SetBranchAddress("Ex",&Ex);
  tree->SetBranchAddress("Vcm",&Vcm);
  tree->SetBranchAddress("thetaCM",&thetaCM);
  tree->SetBranchAddress("cos_thetaH",&cos_thetaH);

  tree->SetBranchAddress("runnum",&runnum);
  tree->SetBranchAddress("beamZ",&beamZ);

  tree->SetBranchAddress("time1",&timeCsI1);
  //tree->SetBranchAddress("timeCsI2",&timeCsI2);


  //name ="Carbon.loss";
  //CLoss loss_Core(name,Acore);
  //name="Hydrogen.loss";
  //CLoss loss_proton(name,1.);

  TFile fileOut("sort.root","RECREATE");

  TDirectoryFile * SortCode = new TDirectoryFile("SortCode","SortCode");
  TDirectoryFile * ReCalc = new TDirectoryFile("ReCalc","ReCalc");
  TDirectory * ReCalc_CsIShift = ReCalc->mkdir("ReCalc_CsIShift","ReCalc_CsIShift");

  SortCode->cd();
  TH1I * Erel_hist = new TH1I("Erel_hist","",250,0,10);
  TH1I * Erel_trans = new TH1I("Erel_trans","",250,0,10);
  TH1I * Erel_long = new TH1I("Erel_long","",250,0,10);
  TH1I * Erel_Sibeam = new TH1I("Erel_Sibeam","",250,0,10);
  TH1I * Erel_Albeam = new TH1I("Erel_Albeam","",250,0,10);
  TH1I * Erel_Mgbeam = new TH1I("Erel_Mgbeam","",250,0,10);
  TH1I * Erel_Nabeam = new TH1I("Erel_Nabeam","",250,0,10);
  TH1I * Ex_hist = new TH1I("Ex_hist","",250,0,10);
  TH1I * Ex_trans = new TH1I("Ex_trans","",500,-5,15);
  TH1I * Ex_trans_0_4 = new TH1I("Ex_trans_0_4","",250,0,10);
  TH1I * Ex_long = new TH1I("Ex_long","",500,-5,15);
  TH1I * Ex_narrowtime = new TH1I("Ex_narrowtime","",250,0,10);
  TH1I * Ex_0_5costhetaH = new TH1I("Ex_0_5costhetaH","",250,0,10);
  //TH1I * Erel_19Na_p18Ne_set1 = new TH1I("Erel_19Na_p18Ne_set1","",500,-5,15);
  //TH1I * Erel_19Na_p18Ne_set1a = new TH1I("Erel_19Na_p18Ne_set1a","",500,-5,15);
  TH1I * ThetaCM = new TH1I("ThetaCM","",200,0,10);
  TH1I * VCM = new TH1I("VCM","",400,7,15);
  TH2I * Erel_costhetaH = new TH2I("Erel_costhetaH","",200,0,8,25,-1,1);
  Erel_costhetaH->SetMinimum(1.0);
  TH2I * VCMvsErel = new TH2I("VCMvsErel","",400,7,15,250,0,10);
  VCMvsErel->SetMinimum(1.0);
  TH1I * gammasADD = new TH1I("gammasADD","",2048, 0, 8192);
  //TH1I * p18Ne_gammasADD_tgate = new TH1I("p18Ne_gammasADD_tgate","",2048, 0, 8192);
  TH2I * gammasADDvsErel = new TH2I("gammasADDvsErel","",250,0,10,2048/8, 0, 2048);
  gammasADDvsErel->SetMinimum(1.0);
  TH2I * gammasADDvsEx = new TH2I("gammasADDvsEx","",250,0,10,2048/8, 0, 2048);
  gammasADDvsEx->SetMinimum(1.0);
  TH2I * gammasADDvsgammasADD = new TH2I("gammasADDvsgammasADD","",2048/4, 0, 8192,2048/4, 0, 8192);
  gammasADDvsgammasADD->SetMinimum(1.0);

  //Gamma-Gamma gates
  TH2I * gammasADDvsErel_greater2 = new TH2I("gammasADDvsErel_greater2","",250,0,10,2048, 0, 8192);
  gammasADDvsErel_greater2->SetMinimum(1.0);
  TH2I * gammasADDvsErel_with200 = new TH2I("gammasADDvsErel_with200","",250,0,10,2048, 0, 8192);
  gammasADDvsErel_with200->SetMinimum(1.0);
  TH2I * gammasADDvsErel_with600 = new TH2I("gammasADDvsErel_with600","",250,0,10,2048, 0, 8192);
  gammasADDvsErel_with600->SetMinimum(1.0);

  TH2I * gammasADDvsgammasADD_with200 = new TH2I("gammasADDvsgammasADD_with200","",2048/4, 0, 8192,2048/4, 0, 8192);
  gammasADDvsgammasADD_with200->SetMinimum(1.0);
  TH2I * gammasADDvsgammasADD_with600 = new TH2I("gammasADDvsgammasADD_with600","",2048/4, 0, 8192,2048/4, 0, 8192);
  gammasADDvsgammasADD_with600->SetMinimum(1.0);

  //Erel gated on gammas
  TH1I * Erel_with500 = new TH1I("Erel_with500","",250,0,10);
  TH1I * Erel_no500 = new TH1I("Erel_no500","",250,0,10);

  //Ex gated on gammas
  TH1I * Ex_with500 = new TH1I("Ex_with500","",250,0,10);
  TH1I * Ex_no500 = new TH1I("Ex_no500","",250,0,10);
  TH1I * Ex_no500_0_5costhetaH = new TH1I("Ex_no500_0_5costhetaH","",250,0,10);

  TH1I * Ex_nogamma= new TH1I("Ex_nogamma","",250,0,10);
  TH1I * Ex_nogamma_0_5costhetaH = new TH1I("Ex_nogamma_0_5costhetaH","",250,0,10);

  //2D gamma t vs channel
  TH2I * gammaTvsChannel = new TH2I("gammaTvsChannel","",192,0,192,5000,-5000,5000);
  TH2I * gammaTvsE = new TH2I("gammaTvsE","",2048/4, 0, 8192,5000,-5000,5000);

  //2D csi time vs channel
  TH2I * CsITvsChannel = new TH2I("CsITvsChannel","",16,-0.5,15.5,2000,-2000,2000);
  //vs Ex
  TH2I * CsITvsEx = new TH2I("CsITvsEx","",250,0,10,2000,-2000,2000);

  ReCalc->cd();

  //S800 vs Janus angles
  TH2I * janus_thetavsphi = new TH2I("janus_thetavsphi","",20,0,10,720,-180,180);
  janus_thetavsphi->SetMinimum(1.0);
  TH2I * s800_thetavsphi = new TH2I("s800_thetavsphi","",20,0,10,720,-180,180);
  s800_thetavsphi->SetMinimum(1.0);

  TH2I * theta_janusvss800 = new TH2I("theta_janusvss800","",1000,0,1,1000,0,1);
  theta_janusvss800->SetMinimum(1.0);
  TH2I * phi_janusvss800 = new TH2I("phi_janusvss800","",3000,-3.5,3.5,3000,-3.5,3.5);
  phi_janusvss800->SetMinimum(1.0);


  //Recalculated Erel
  TH1I * CsIE = new TH1I("CsIE","",400,0,200);
  TH2I * ErelvsCsIE = new TH2I("ErelvsCsIE","",250,0,10,400,0,200);
  TH1I * CsI_theta = new TH1I("CsI_theta","",300,0,60);
  TH1I * HE = new TH1I("HE","",400,0,200);
  TH1I * H_theta = new TH1I("H_theta","",300,0,60);
  TH1I * Erel_re = new TH1I("Erel_re","",250,0,10);
  TH1I * Erel_re_set1 = new TH1I("Erel_re_set1","",250,0,10);
  TH1I * Erel_re_set1a = new TH1I("Erel_re_set1a","",250,0,10);

  TH1I * VCM_re = new TH1I("VCM_re","",400,7,15);

  TH2I * Erel_costhetaH_re = new TH2I("Erel_costhetaH_re","",200,0,8,25,-1,1);
  Erel_costhetaH_re->SetMinimum(1.0);
  TH2I * Erel_costhetaH_re_set1 = new TH2I("Erel_costhetaH_re_set1","",200,0,8,25,-1,1);
  Erel_costhetaH_re_set1->SetMinimum(1.0);
  TH2I * Erel_costhetaH_re_set1a = new TH2I("Erel_costhetaH_re_set1a","",200,0,8,25,-1,1);
  Erel_costhetaH_re_set1a->SetMinimum(1.0);

  TH1I * Erel_long_re = new TH1I("Erel_long_re","",250,0,10);
  TH1I * Erel_long_re_set1 = new TH1I("Erel_long_re_set1","",250,0,10);
  TH1I * Erel_long_re_set1a = new TH1I("Erel_long_re_set1a","",250,0,10);

  TH1I * Erel_trans_re = new TH1I("Erel_trans_re","",250,0,10);
  TH1I * Erel_trans_re_set1 = new TH1I("Erel_trans_re_set1","",250,0,10);
  TH1I * Erel_trans_re_set1a = new TH1I("Erel_trans_re_set1a","",250,0,10);

  TH1I * Ex_nogamma_re = new TH1I("Ex_nogamma_re","",250,0,10);
  TH1I * Ex_no500_re = new TH1I("Ex_no500_re","",250,0,10);
  TH1I * Ex_no500_trans_re = new TH1I("Ex_no500_trans_re","",250,0,10);
  TH1I * Ex_no500_s800_re = new TH1I("Ex_no500_s800_re","",250,0,10);

  TH2I * pXY_s = new TH2I("pXY_s","",100,-10,10,100,-10,10);
  TH2I * HXY_s = new TH2I("HXY_s","",100,-10,10,100,-10,10);

  //E* for different absolute calibration shifts
  ReCalc_CsIShift->cd(); 
  TH1I * Ex_re_CsIShift_min1 = new TH1I("Ex_re_CsIShift_min1","",250,0,10);
  TH1I * Ex_re_CsIShift_minp9 = new TH1I("Ex_re_CsIShift_minp9","",250,0,10);
  TH1I * Ex_re_CsIShift_minp8 = new TH1I("Ex_re_CsIShift_minp8","",250,0,10);
  TH1I * Ex_re_CsIShift_minp7 = new TH1I("Ex_re_CsIShift_minp7","",250,0,10);
  TH1I * Ex_re_CsIShift_minp6 = new TH1I("Ex_re_CsIShift_minp6","",250,0,10);
  TH1I * Ex_re_CsIShift_minp5 = new TH1I("Ex_re_CsIShift_minp5","",250,0,10);
  TH1I * Ex_re_CsIShift_minp4 = new TH1I("Ex_re_CsIShift_minp4","",250,0,10);
  TH1I * Ex_re_CsIShift_minp3 = new TH1I("Ex_re_CsIShift_minp3","",250,0,10);
  TH1I * Ex_re_CsIShift_minp2 = new TH1I("Ex_re_CsIShift_minp2","",250,0,10);
  TH1I * Ex_re_CsIShift_minp1 = new TH1I("Ex_re_CsIShift_minp1","",250,0,10);
  TH1I * Ex_re_CsIShift_min0 = new TH1I("Ex_re_CsIShift_min0","",250,0,10);
  TH1I * Ex_re_CsIShift_plsp1 = new TH1I("Ex_re_CsIShift_plsp1","",250,0,10);
  TH1I * Ex_re_CsIShift_plsp2 = new TH1I("Ex_re_CsIShift_plsp2","",250,0,10);
  TH1I * Ex_re_CsIShift_plsp3 = new TH1I("Ex_re_CsIShift_plsp3","",250,0,10);
  TH1I * Ex_re_CsIShift_plsp4 = new TH1I("Ex_re_CsIShift_plsp4","",250,0,10);
  TH1I * Ex_re_CsIShift_plsp5 = new TH1I("Ex_re_CsIShift_plsp5","",250,0,10);
  TH1I * Ex_re_CsIShift_plsp6 = new TH1I("Ex_re_CsIShift_plsp6","",250,0,10);
  TH1I * Ex_re_CsIShift_plsp7 = new TH1I("Ex_re_CsIShift_plsp7","",250,0,10);
  TH1I * Ex_re_CsIShift_plsp8 = new TH1I("Ex_re_CsIShift_plsp8","",250,0,10);
  TH1I * Ex_re_CsIShift_plsp9 = new TH1I("Ex_re_CsIShift_plsp9","",250,0,10);
  TH1I * Ex_re_CsIShift_pls1 = new TH1I("Ex_re_CsIShift_pls1","",250,0,10);

  file.cd();
  for (int i=0;i<N;i++)
  {
    if (i%3000 == 0)cout << i << endl;
    
    tree->GetEntry(i); //Pull data from tree
    //cout << timeCsI1 << endl;

    //Set runnum
    bool set1 = true; //true for setting 1, false for setting1a
    if (runnum < 216) gobbidist = 55.2;
    if (runnum >= 216)
    {
      gobbidist = 43.896;//45.08;
      set1 = false;
    }

    if (timeCsI1 < 960 || timeCsI1 > 1030) continue;

    Erel_hist->Fill(Erel);
    Ex_hist->Fill(Ex);
    ThetaCM->Fill(thetaCM*180./pi);
    VCM->Fill(Vcm);
    VCMvsErel->Fill(Vcm,Erel);
    Erel_costhetaH->Fill(Erel,cos_thetaH);

    //Fill time spectra 
    CsITvsChannel->Fill(itele1*4 + id1,timeCsI1);
    CsITvsEx->Fill(Ex,timeCsI1);

    if (beamZ == 14) Erel_Sibeam->Fill(Erel);
    if (beamZ == 13) Erel_Albeam->Fill(Erel);
    if (beamZ == 12) Erel_Mgbeam->Fill(Erel);
    if (beamZ == 11) Erel_Nabeam->Fill(Erel);

    //CsI time gate between 970 and 1030
    if (timeCsI1 >= 970 && timeCsI1 <= 1030) Ex_narrowtime->Fill(Ex);

    //ThetaH gates
    if (fabs(cos_thetaH) <= 0.4)
    {
      Erel_trans->Fill(Erel);
      Ex_trans->Fill(Ex);
    }
    if (fabs(cos_thetaH) >= 0.8)
    {
      Erel_long->Fill(Erel);
      Ex_long->Fill(Ex);
    }

    //ThetaH gates
    if (fabs(cos_thetaH) <= 0.4)
    {
      Ex_trans_0_4->Fill(Ex);
    }

    //Gamma coincidence gates
    int with500 = 0;
    int no500 = 1;

    for (int j=0;j<Ngamma;j++) {
      if (Tgamma[j] < 675 || Tgamma[j] > 715) continue;
      if (Egamma[j] >= 0.443 && Egamma[j] <= 0.522) with500 = 1;
      if (Egamma[j] >= 0.443 && Egamma[j] <= 0.522) no500 = 0;
    } 

    //Erel gated on gammas
    if (with500 == 1) Erel_with500->Fill(Erel);
    if (no500 == 1) Erel_no500->Fill(Erel);

    //Ex gated on gammas
    if (with500 == 1) Ex_with500->Fill(Ex);
    if (with500 == 0) Ex_no500->Fill(Ex);

    if (fabs(cos_thetaH) <= 0.5) Ex_0_5costhetaH->Fill(Ex);

    if (Ngamma == 0) Ex_nogamma->Fill(Ex);
    if (Ngamma == 0 && fabs(cos_thetaH) <= 0.5) Ex_nogamma_0_5costhetaH->Fill(Ex);
    if (with500 == 0 && fabs(cos_thetaH) <= 0.5) Ex_no500_0_5costhetaH->Fill(Ex);

    //Gammas - already doppler shifted
    for (int j=0;j<Ngamma;j++)
    {
      if (Ngamma >= 30) {
        cout << "Ngamma >= 30, skip " << endl;
        break;
      }

      //Gamma time gate
      if (Tgamma[j] < 675 || Tgamma[j] > 715) continue;

      gammasADD->Fill(Egamma[j]*1000.);
      gammasADDvsErel->Fill(Erel,Egamma[j]*1000.);
      gammasADDvsEx->Fill(Ex,Egamma[j]*1000.);

      gammaTvsChannel->Fill(Chgamma[j],Tgamma[j]);
      gammaTvsE->Fill(Egamma[j]*1000.,Tgamma[j]);
      
      //cout << Tgamma[j] << endl;

      if (Ngamma > 2) gammasADDvsErel_greater2->Fill(Erel,Egamma[j]*1000.);

      for (int k=j+1;k<Ngamma;k++)
      {
        gammasADDvsgammasADD->Fill(Egamma[j]*1000.,Egamma[k]*1000.);

      }
    }

    //Only look at Mg beam
    //if (beamZ != 12) continue;

    //Fill ECsI_R
    //ECsI_R[id1 + itele1*4]->Fill(energyR,denergyR);

    //Apply calibration
    double denergy = FrontEcal.getEnergy(itele1, ifront1, denergyR);
    double energy = CsIEcal.getEnergy(itele1, id1, energyR);

    //Do fine calibration here, modify CsI "energy" by some a + b(x-x_0) where x_0 is avg p E
    CsIE->Fill(energy);
    ErelvsCsIE->Fill(Erel,energy);
    CsI_theta->Fill(theta1*180./acos(-1.));
    double avgE_p = 77.;
    float a = 0.;//0.15;
    float b = 0.0;//-0.015;

    HE->Fill(energy_p2/17.);
    H_theta->Fill(theta2*180./acos(-1.));

    energy += a + b*(energy - avgE_p);

    //Fill ECsI_cal
    //ECsI_cal[id1 + itele1*4]->Fill(energy,denergy);

    //Shift the detector distance and recalculate angles
    if (set1 == true) gobbidist += 0.0;
    else gobbidist += 0.0;
    
    Pix.SetTarget(gobbidist,Bethick);
    //cout << theta1 << " ";
    theta1 = Pix.getAngle(itele1, ifront1, iback1);
    //cout << theta1 << endl;
    phi1 = Pix.phi;
    
    
    //Do Egain calculations
    double sumEnergy = denergy+energy;
    sumEnergy = lossp_Al.getEin(sumEnergy,Althick/cos(theta1));
    sumEnergy = lossp_Be.getEin(sumEnergy,(Bethick/2.)/cos(theta1));


    //Recalculating Erel. User sort angle and E for heavy fragment
    et1 = sumEnergy + Mass_p;
    double  pc = sqrt(pow(et1,2) - pow(Mass_p,2));

    M1[0] = pc*sin(theta1)*cos(phi1);
    M1[1] = pc*sin(theta1)*sin(phi1);
    M1[2] = pc*cos(theta1); 

    //Assume theta & phi == 0
    /*theta2 = 0.;
    phi2 = 0.;
    double pc_H = sqrt(pow(et2,2) - pow(Mass_17F,2));
    M2[0] = pc_H*sin(theta2)*cos(phi2);
    M2[1] = pc_H*sin(theta2)*sin(phi2);
    M2[2] = pc_H*cos(theta2);*/

    //Recalculate theta/phi based on more accurate fiber distance
    double fibdist_og = gobbidist + 10.151+2.664;
    double fibdist = gobbidist + 10.557;
    double x_fib = fibdist_og*tan(theta2)*cos(phi2);
    double y_fib = fibdist_og*tan(theta2)*sin(phi2);
    double r_fib = sqrt(pow(fibdist,2)+pow(x_fib,2)+pow(y_fib,2));
    theta2 = acos(fibdist/r_fib);
    phi2 = atan2(y_fib,x_fib);
    double pc_H = sqrt(pow(et2,2) - pow(Mass_17F,2));
    M2[0] = pc_H*sin(theta2)*cos(phi2);
    M2[1] = pc_H*sin(theta2)*sin(phi2);
    M2[2] = pc_H*cos(theta2);

    CM = Kin.findCM(M1,M2,et1,et2);
    VCM_re->Fill(CM.vv);
    ee1 = Kin.trans(CM,M1,et1);
    ee2 = Kin.trans(CM,M2,et2);
    Erel = ee1.vv + ee2.vv - Mass_17F - Mass_p;

    //Finding the heavy fragment cos(theta)
    //TODO heavy fragment always N-1 index?
    float mv = 0.;
    for (int j=0;j<3;j++)
    {
      mv += pow(ee2.v[j],2);
    }
    mv = sqrt(mv);
    cos_thetaH = ee2.v[2]/mv;

    Erel_re->Fill(Erel);
    Ex = Erel - (Mass_18Ne - (Mass_17F + Mass_p));

    if (set1 == true)
    {
      Erel_re_set1->Fill(Erel);
      Erel_costhetaH_re_set1->Fill(Erel,cos_thetaH);
      if (fabs(cos_thetaH) <= 0.4) Erel_trans_re_set1->Fill(Erel);
      if (fabs(cos_thetaH) >= 0.7) Erel_long_re_set1->Fill(Erel);
    }
    else
    {
      Erel_re_set1a->Fill(Erel);
      Erel_costhetaH_re_set1a->Fill(Erel,cos_thetaH);
      if (fabs(cos_thetaH) <= 0.4) Erel_trans_re_set1a->Fill(Erel);
      if (fabs(cos_thetaH) >= 0.7) Erel_long_re_set1a->Fill(Erel);
    }

    Erel_costhetaH_re->Fill(Erel,cos_thetaH);

    if (fabs(cos_thetaH) <= 0.4) Erel_trans_re->Fill(Erel);
    if (fabs(cos_thetaH) >= 0.7) Erel_long_re->Fill(Erel);

    if (Ngamma == 0) Ex_nogamma_re->Fill(Ex);
    if (with500 == 0)
    {
      Ex_no500_re->Fill(Ex);
      if (fabs(cos_thetaH) <= 0.4) Ex_no500_trans_re->Fill(Ex);
    }

    //Plots janus vs s800 angles for theta and phi

    //cout << "janus and s800 theta: " << theta2 << " " << theta_s800 << endl;
    //cout << "janus and s800 phi: " << phi2 << " " << phi_s800 << " sub " << phi2-phi_s800 << endl;

    janus_thetavsphi->Fill(theta2*180./acos(-1.),phi2*180./acos(-1.));
    s800_thetavsphi->Fill(theta_s800*180./acos(-1.),phi_s800*180./acos(-1.));

    theta_janusvss800->Fill(theta2,theta_s800);
    phi_janusvss800->Fill(phi2,phi_s800);

    //Hit map
    double x = theta1*180./acos(-1.)*cos(phi1);
    double y = theta1*180./acos(-1.)*sin(phi1);
    pXY_s->Fill(x,y);

    x = theta2*180./acos(-1.)*cos(phi2);
    y = theta2*180./acos(-1.)*sin(phi2);
    HXY_s->Fill(x,y);

    //Calculate the Erel using the S800 theta/phi
    pc_H = sqrt(pow(et2,2) - pow(Mass_17F,2));
    M2[0] = pc_H*sin(theta_s800)*cos(phi_s800);
    M2[1] = pc_H*sin(theta_s800)*sin(phi_s800);
    M2[2] = pc_H*cos(theta_s800);

    CM = Kin.findCM(M1,M2,et1,et2);
    VCM_re->Fill(CM.vv);
    ee1 = Kin.trans(CM,M1,et1);
    ee2 = Kin.trans(CM,M2,et2);
    double Erel_s800 = ee1.vv + ee2.vv - Mass_17F - Mass_p;
    double Ex_s800 = Erel_s800 - (Mass_18Ne - (Mass_17F + Mass_p));
    
    if (with500 == 0) Ex_no500_s800_re->Fill(Ex_s800); 

    //Iterate over absolute CsI calibration shifts, -1 to 1 MeV, 0.1 MeV increment.
    float calchange = -0.1;
    float calend = 0.1;
    float calit = 0.01;
    int calnum = (calend - calchange)/calit;
    /*for (int j=0;j<=calnum;j++)
    {
      //Only no 495 keV gamma
      if (no500 == 0) break;

      Pix.SetTarget(gobbidist,Bethick);
      theta1 = Pix.getAngle(itele1, ifront1, iback1);
      phi1 = Pix.phi;

      //Do fine calibration here, modify CsI "energy" by some a + b(x-x_0) where x_0 is avg p E
      double avgE_p = 77.;
      float a = 0;
      float b = calchange;//-0.015;

      energy = 0;
      energy = CsIEcal.getEnergy(itele1, id1, energyR);
      energy += a + b*(energy - avgE_p);

      //Do Egain calculations
      double sumEnergy_2 = denergy+energy;
      sumEnergy_2 = lossp_Al.getEin(sumEnergy_2,Althick/cos(theta1));
      sumEnergy_2 = lossp_Be.getEin(sumEnergy_2,(Bethick/2.)/cos(theta1));


      //Recalculating Erel. User sort angle and E for heavy fragment
      et1 = sumEnergy_2 + Mass_p;
      double  pc = sqrt(pow(et1,2) - pow(Mass_p,2));

      M1[0] = pc*sin(theta1)*cos(phi1);
      M1[1] = pc*sin(theta1)*sin(phi1);
      M1[2] = pc*cos(theta1); 

      CM = Kin.findCM(M1,M2,et1,et2);
      VCM_re->Fill(CM.vv);
      ee1 = Kin.trans(CM,M1,et1);
      ee2 = Kin.trans(CM,M2,et2);
      Erel = ee1.vv + ee2.vv - Mass_17F - Mass_p;
      Ex = Erel - (Mass_18Ne - (Mass_17F + Mass_p));

      mv = 0.;
      for (int j=0;j<3;j++)
      {
        mv += pow(ee2.v[j],2);
      }
      mv = sqrt(mv);
      cos_thetaH = ee2.v[2]/mv;

      if (fabs(calchange - (-.1)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_min1->Fill(Ex);
      if (fabs(calchange - (-0.09)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_minp9->Fill(Ex);
      if (fabs(calchange - (-0.08)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_minp8->Fill(Ex);
      if (fabs(calchange - (-0.07)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_minp7->Fill(Ex);
      if (fabs(calchange - (-0.06)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_minp6->Fill(Ex);
      if (fabs(calchange - (-0.05)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_minp5->Fill(Ex);
      if (fabs(calchange - (-0.04)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_minp4->Fill(Ex);
      if (fabs(calchange - (-0.03)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_minp3->Fill(Ex);
      if (fabs(calchange - (-0.02)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_minp2->Fill(Ex);
      if (fabs(calchange - (-0.01)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_minp1->Fill(Ex);
      if (fabs(calchange - (0)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_min0->Fill(Ex);
      if (fabs(calchange - (0.01)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_plsp1->Fill(Ex);
      if (fabs(calchange - (0.02)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_plsp2->Fill(Ex);
      if (fabs(calchange - (0.03)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_plsp3->Fill(Ex);
      if (fabs(calchange - (0.04)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_plsp4->Fill(Ex);
      if (fabs(calchange - (0.05)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_plsp5->Fill(Ex);
      if (fabs(calchange - (0.06)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_plsp6->Fill(Ex);
      if (fabs(calchange - (0.07)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_plsp7->Fill(Ex);
      if (fabs(calchange - (0.08)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_plsp8->Fill(Ex);
      if (fabs(calchange - (0.09)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_plsp9->Fill(Ex);
      if (fabs(calchange - (.1)) < std::numeric_limits<float>::epsilon()) Ex_re_CsIShift_pls1->Fill(Ex);
      //cout << calchange << endl;
      calchange += calit;
    }*/


/*
//et1 = FineCal.cal_p(energy_p1,denergy1,ifront1,iback1,id1,193.,M1);
//double theta1 = FineCal.theta;
//double phi1 =FineCal.phi;
//et2 = FineCal.cal_O14(energy_p2,denergy2,ifront2,iback2,id2,193.,M2);
//double theta2 = FineCal.theta;
//double phi2 =FineCal.phi;

  double theta1;
  double phi1;
  
  double theta2;
  double phi2;

/*
      double theta1 = Pix.getAngle(id1/4,ifront1,iback1);
      double phi1 = Pix.phi;
      double E = energy_p1;
      E += denergy1;
        double thick = 193./2./cos(theta1);
        E = loss_p.getEin(E,thick);

        et1 = E + Mass_p;
        double  pc = sqrt(pow(et1,2) - pow(Mass_p,2));

        M1[0] = pc*sin(theta1)*cos(phi1);
        M1[1] = pc*sin(theta1)*sin(phi1);
        M1[2] = pc*cos(theta1);



      double theta2 = Pix.getAngle(id2/4,ifront2,iback2);
      double phi2 = Pix.phi;
       E = energy_p2;
       E = cal.getEnergy(0,id2,E);
       E += 2.8;
       E += denergy2;
        thick = 193./2./cos(theta2);
        E = loss_O.getEin(E,thick);

        et2 = E + Mass_14O;
        pc = sqrt(pow(et2,2) - pow(Mass_14O,2));

        M2[0] = pc*sin(theta2)*cos(phi2);
        M2[1] = pc*sin(theta2)*sin(phi2);
        M2[2] = pc*cos(theta2);


      //hist_Ek_14O->Fill(et2-Mass_14O);

      
      //CM = Kin.findCM(M1,M2,et1,et2);
      //hist_Vcm->Fill(CM.vv);
      //ee1 = Kin.trans(CM,M1,et1);
      //ee2 = Kin.trans(CM,M2,et2);
      //float Erel = ee1.vv + ee2.vv - Mass_14O - Mass_p;





      bool O15beam = false;
      bool Ne17beam = false;
      if (CM.vv > 9.39)
      {
      Ne17beam = true;
      }
     if (CM.vv < 9.27)
      {
       O15beam = true;
      }



     bool Bgamma = false;
     for (int j=0;j<Ngamma;j++)
       {
	 if (Egamma[j] > 0.01)
	   {
            hist_Egamma->Fill(Egamma[j]);
            map_Erel_Egamma->Fill(Erel,Egamma[j]);
            if(O15beam)hist_Egamma_15O->Fill(Egamma[j]);
            if(Ne17beam)hist_Egamma_17Ne->Fill(Egamma[j]);
            if (Egamma[j] > 4.16 && Egamma[j] < 5.6) Bgamma = true;
	   }
       }
*/


//if (CM.vv < 9.36) continue;    //sort only 17Ne beam contribution
/*
      hist_Vcm->Fill(CM.vv);
      hist_Erel->Fill(Erel);


/if (Erel > 2.5 && Erel < 3. && CM.vv > 9.36)
{
 hist_Vcm_peak->Fill(CM.vv);
 hist_thetaCM_peak->Fill(acos(CM.v[2]/CM.vv)*180./acos(-1.));
}
      double cosTheta = ee2.v[2]/sqrt(pow(ee2.v[0],2)+pow(ee2.v[1],2)+pow(ee2.v[2],2));

      map_Erel_cos->Fill(Erel,cosTheta);
      map_Vcm_cos->Fill(CM.vv,cosTheta);
      if (fabs(cosTheta) < .5) hist_Erel_trans->Fill(Erel);
      if (fabs(cosTheta) < .2) hist_Erel_trans_narrow->Fill(Erel);
      if (cosTheta > .7) hist_Erel_for->Fill(Erel);
      if (cosTheta < -.4 && cosTheta > -.7) hist_Erel_mid_back->Fill(Erel);
      if (cosTheta < -.7) hist_Erel_back->Fill(Erel);

      bool BB17 = false;
      bool BB15 = false;
      if (Bgamma) hist_Erel_gamma->Fill(Erel);



      if (CM.vv > 9.39)
      {

      hist_Erel_17Ne->Fill(Erel);
      if (Bgamma) hist_Erel_17Ne_gamma->Fill(Erel);
      if (fabs(cosTheta) < .5) hist_Erel_trans_17Ne->Fill(Erel);
      if (fabs(cosTheta) < .2) hist_Erel_trans_narrow_17Ne->Fill(Erel);
      BB17 = true;
      map_Erel_cos_17Ne->Fill(Erel,cosTheta);
      }
     else if (CM.vv < 9.27)
      {

       hist_Erel_15O->Fill(Erel);
      if (Bgamma) hist_Erel_15O_gamma->Fill(Erel);
      if (fabs(cosTheta) < .5) hist_Erel_trans_15O->Fill(Erel);
      if (fabs(cosTheta) < .2) hist_Erel_trans_narrow_15O->Fill(Erel);
      BB15 = true;
      map_Erel_cos_15O->Fill(Erel,cosTheta);
      }




   map_proton->Fill(theta1*cos(phi1),theta1*sin(phi1));    
   map_core->Fill(theta2*cos(phi2),theta2*sin(phi2));    
      save saveP(id1,et1,M1);
      save saveO(id2,et2,M2);
*/
      /*for (int j=0;j<nsave;j++)
        {

           if(savep[j].id != saveP.id && savep[j].id != saveO.id && BB17 && B17[j])
             {
              CM = Kin.findCM(savep[j],saveP,saveO);
              ee1 = Kin.trans(CM,savep[j]);
              ee2 = Kin.trans(CM,saveO);
              ee3 = Kin.trans(CM,saveP);
              double Erel = ee1.vv + ee2.vv + ee3.vv - 2.*Mass_p - Mass_14O;

              CM = Kin.findCM(savep[j],saveO);
              ee1 = Kin.trans(CM,savep[j]);
              ee2 = Kin.trans(CM,saveO);
              double Erel_p_O = ee1.vv + ee2.vv;
              double supp = suppress(Erel_p_O);
               hist_Erel_16Ne_mix_weight->Fill(Erel,supp);
               double cosTheta = ee3.v[2]/sqrt(pow(ee3.v[0],2)+pow(ee3.v[1],2)+pow(ee3.v[2],2));
               if(fabs(cosTheta) < .2) hist_Erel_16Ne_mix_weight_trans_narrow->Fill(Erel,supp);
             }

           if(saveP.id != savep[j].id && saveP.id != save18[j].id && BB17 && B17[j])
             {
              CM = Kin.findCM(savep[j],saveP,save18[j]);
              ee1 = Kin.trans(CM,savep[j]);
              ee2 = Kin.trans(CM,save18[j]);
              ee3 = Kin.trans(CM,saveP);
              double Erel = ee1.vv + ee2.vv + ee3.vv - 2.*Mass_p - Mass_14O;

              CM = Kin.findCM(saveP,save18[j]);
              ee1 = Kin.trans(CM,savep[j]);
              ee2 = Kin.trans(CM,saveO);
              double Erel_p_O = ee1.vv + ee2.vv;
              double supp = suppress(Erel_p_O);
              hist_Erel_16Ne_mix_weight->Fill(Erel,supp);
               double cosTheta = ee3.v[2]/sqrt(pow(ee3.v[0],2)+pow(ee3.v[1],2)+pow(ee3.v[2],2));
               if(fabs(cosTheta) < .2) hist_Erel_16Ne_mix_weight_trans_narrow->Fill(Erel,supp);

             }


           if (savep[j].id != saveO.id)
             {
              CM = Kin.findCM(savep[j],saveO);
              ee1 = Kin.trans(CM,savep[j]);
              ee2 = Kin.trans(CM,saveO);
              double Erel = ee1.vv + ee2.vv - Mass_p - Mass_14O;
              double supp = suppress(Erel);
              //cout << Erel << " " << supp << endl;
              hist_Erel_mix->Fill(Erel,supp);
               if (BB17 && B17[j])hist_Erel_17Ne_mix->Fill(Erel,supp);
               if (BB15 && B15[j])hist_Erel_15O_mix->Fill(Erel,supp);
              double cosTheta = ee2.v[2]/sqrt(pow(ee2.v[0],2)+pow(ee2.v[1],2)+pow(ee2.v[2],2));
              if (fabs(cosTheta) < .2)
                {

                 hist_Erel_trans_narrow_mix->Fill(Erel,supp);
                 if (BB17 && B17[j])hist_Erel_trans_narrow_17Ne_mix->Fill(Erel,supp);
                 if (BB15 && B15[j])hist_Erel_trans_narrow_15O_mix->Fill(Erel,supp);
                }
             }

           if (saveP.id != save18[j].id)
             {
              CM = Kin.findCM(saveP,save18[j]);
              ee1 = Kin.trans(CM,saveP);
              ee2 = Kin.trans(CM,save18[j]);
              double Erel = ee1.vv + ee2.vv - Mass_p - Mass_14O;
              double supp = suppress(Erel);
              //cout << Erel << " " << supp << endl;
              hist_Erel_mix->Fill(Erel,supp);
                 if (BB17 && B17[j])hist_Erel_17Ne_mix->Fill(Erel,supp);
                 if (BB15 && B15[j])hist_Erel_15O_mix->Fill(Erel,supp);


              double cosTheta = ee2.v[2]/sqrt(pow(ee2.v[0],2)+pow(ee2.v[1],2)+pow(ee2.v[2],2));
              if (fabs(cosTheta) < .2)
                {

                hist_Erel_trans_narrow_mix->Fill(Erel,supp);
                 if (BB17 && B17[j])hist_Erel_trans_narrow_17Ne_mix->Fill(Erel,supp);
                 if (BB15 && B15[j])hist_Erel_trans_narrow_15O_mix->Fill(Erel,supp);
             }
           }
        }*/

      //B17[nwrite] = BB17;
      //B15[nwrite] = BB15;
      //savep[nwrite] = saveP;
      //save18[nwrite] = saveO;
      //nwrite++;
      //if(nwrite == nn) nwrite = 0;
      //nsave++;
      //if(nsave == nn) nsave = nn -1;
    }

  fileOut.Write();
}
