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

//from p+20Mg 
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
  save save18[nn];
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

  TFile file("p20Mg_save.root");
  //TFile file("p20Mg.root");
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
  float denergyR;
  float energyR;

  float theta1;
  float theta2;
  float phi1;
  float phi2;


  float denergy2;
  float energy_p2;


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
  tree->SetBranchAddress("phi1",&phi1);
  tree->SetBranchAddress("phi2",&phi2);

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



  //tree->SetBranchAddress("timeCsI1",&timeCsI1);
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
  TH1I * Erel_21Al_p20Mg = new TH1I("Erel_21Al_p20Mg","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_vfilt = new TH1I("Erel_21Al_p20Mg_vfilt","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_set1 = new TH1I("Erel_21Al_p20Mg_set1","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_set1a = new TH1I("Erel_21Al_p20Mg_set1a","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_tgate = new TH1I("Erel_21Al_p20Mg_tgate","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_trans = new TH1I("Erel_21Al_p20Mg_trans","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_trans0_2 = new TH1I("Erel_21Al_p20Mg_trans0_2","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_trans_set1 = new TH1I("Erel_21Al_p20Mg_trans_set1","",250,0,10);
  TH1I * Erel_21Al_p20Mg_trans_set1a = new TH1I("Erel_21Al_p20Mg_trans_set1a","",250,0,10);
  TH1I * Erel_21Al_p20Mg_trans0_2_set1a = new TH1I("Erel_21Al_p20Mg_trans0_2_set1a","",250,0,10);
  TH1I * Erel_21Al_p20Mg_trans_narr = new TH1I("Erel_21Al_p20Mg_trans_narr","",250,0,10);
  TH1I * Erel_21Al_p20Mg_long = new TH1I("Erel_21Al_p20Mg_long","",250,0,10);
  TH1I * Erel_21Al_p20Mg_Sibeam = new TH1I("Erel_21Al_p20Mg_Sibeam","",250,0,10);
  TH1I * Erel_21Al_p20Mg_Albeam = new TH1I("Erel_21Al_p20Mg_Albeam","",250,0,10);
  TH1I * Erel_21Al_p20Mg_Mgbeam = new TH1I("Erel_21Al_p20Mg_Mgbeam","",250,0,10);
  TH1I * Ex_21Al_p20Mg = new TH1I("Ex_21Al_p20Mg","",500,-5,15);
  TH1I * Ex_21Al_p20Mg_trans = new TH1I("Ex_21Al_p20Mg_trans","",500,-5,15);
  TH1I * Ex_21Al_p20Mg_long = new TH1I("Ex_21Al_p20Mg_long","",500,-5,15);
  //TH1I * Erel_19Na_p18Ne_set1 = new TH1I("Erel_19Na_p18Ne_set1","",500,-5,15);
  //TH1I * Erel_19Na_p18Ne_set1a = new TH1I("Erel_19Na_p18Ne_set1a","",500,-5,15);
  TH1I * ThetaCM_21Al_p20Mg = new TH1I("ThetaCM_21Al_p20Mg","",200,0,10);
  TH1I * VCM_21Al_p20Mg = new TH1I("VCM_21Al_p20Mg","",400,7,15);
  TH2I * Erel_p20Mg_costhetaH = new TH2I("Erel_p20Mg_costhetaH","",200,0,8,25,-1,1);
  Erel_p20Mg_costhetaH->SetMinimum(1.0);
  TH2I * p20Mg_VCMvsErel = new TH2I("p20Mg_VCMvsErel","",400,7,15,500,0,10);
  p20Mg_VCMvsErel->SetMinimum(1.0);
  TH1I * p20Mg_gammasADD = new TH1I("p20Mg_gammasADD","",2048, 0, 8192);
  TH1I * p20Mg_gammasADD_mult1 = new TH1I("p20Mg_gammasADD_mult1","",2048, 0, 8192);
  TH1I * gammasADD_mult = new TH1I("gammasADD_mult","",20, 0, 20);
  //TH1I * p18Ne_gammasADD_tgate = new TH1I("p18Ne_gammasADD_tgate","",2048, 0, 8192);
  TH2I * p20Mg_gammasADDvsErel = new TH2I("p20Mg_gammasADDvsErel","",250,0,10,2048, 0, 8192);
  p20Mg_gammasADDvsErel->SetMinimum(1.0);
  TH2I * p20Mg_gammasADDvsErel_trans = new TH2I("p20Mg_gammasADDvsErel_trans","",250,0,10,2048, 0, 8192);
  p20Mg_gammasADDvsErel_trans->SetMinimum(1.0);
  TH2I * p20Mg_gammasADDvsgammasADD = new TH2I("p20Mg_gammasADDvsgammasADD","",2048, 0, 8192,2048, 0, 8192);
  p20Mg_gammasADDvsgammasADD->SetMinimum(1.0);
  TH1I * p20Mg_gammasCh = new TH1I("p20Mg_gammasCh","",192, -0.5, 191.5);
  TH2I * p20Mg_gammasADDvsCh = new TH2I("p20Mg_gammasADDvsCh","",192, 0, 192,2048/4, 0, 8192);
  p20Mg_gammasADDvsCh->SetMinimum(1.0);
  TH2I * p20Mg_gammasADDvsCh_mult1 = new TH2I("p20Mg_gammasADDvsCh_mult1","",192, 0, 192,2048/4, 0, 8192);
  p20Mg_gammasADDvsCh_mult1->SetMinimum(1.0);


  //Erel gated on gammas
  TH1I * Erel_21Al_p20Mg_with1600 = new TH1I("Erel_21Al_p20Mg_with1600","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_with1600_mult1 = new TH1I("Erel_21Al_p20Mg_with1600_mult1","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_with1600_set1a = new TH1I("Erel_21Al_p20Mg_with1600_set1a","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_no1600 = new TH1I("Erel_21Al_p20Mg_no1600","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_no1600_set1a = new TH1I("Erel_21Al_p20Mg_no1600_set1a","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_with1600_trans = new TH1I("Erel_21Al_p20Mg_with1600_trans","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_no1600_trans = new TH1I("Erel_21Al_p20Mg_no1600_trans","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_nogamma = new TH1I("Erel_21Al_p20Mg_nogamma","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_nogamma_trans = new TH1I("Erel_21Al_p20Mg_nogamma_trans","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_withgamma = new TH1I("Erel_21Al_p20Mg_withgamma","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_left1600 = new TH1I("Erel_21Al_p20Mg_left1600","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_right1600 = new TH1I("Erel_21Al_p20Mg_right1600","",500,-5,15);
  TH1I * Erel_21Al_p20Mg_withgamma_trans = new TH1I("Erel_21Al_p20Mg_withgamma_trans","",500,-5,15);

  TH2I * gammaTvsE = new TH2I("gammaTvsE","",2048/4, 0, 8192,5000,-5000,5000);

  //vs Ex
  TH2I * CsITvsEx = new TH2I("CsITvsEx","",250,0,10,2000,-2000,2000);

  ReCalc->cd();
  //Recalculated Erel
  TH1I * CsIE = new TH1I("CsIE","",400,0,200);
  TH1I * Erel_re = new TH1I("Erel_re","",250,0,10);
  TH1I * Erel_re_set1 = new TH1I("Erel_re_set1","",250,0,10);
  TH1I * Erel_re_set1a = new TH1I("Erel_re_set1a","",250,0,10);
  TH1I * Ex_re = new TH1I("Ex_re","",250,0,10);
  TH1I * Ex_re_set1 = new TH1I("Ex_re_set1","",250,0,10);
  TH1I * Ex_re_set1a = new TH1I("Ex_re_set1a","",250,0,10);

  //Different distances for setting 1
  TH1I * Ex_re_set1_minp5 = new TH1I("Ex_re_set1_minp5","",250,0,10);
  TH1I * Ex_re_set1_min1 = new TH1I("Ex_re_set1_min1","",250,0,10);
  TH1I * Ex_re_set1_plsp5 = new TH1I("Ex_re_set1_plsp5","",250,0,10);
  TH1I * Ex_re_set1_pls1 = new TH1I("Ex_re_set1_pls1","",250,0,10);

  TH1I * VCM_re = new TH1I("VCM_re","",400,7,15);

  TH2I * Erel_costhetaH_re = new TH2I("Erel_costhetaH_re","",200,0,8,25,-1,1);
  Erel_costhetaH_re->SetMinimum(1.0);
  TH2I * Erel_costhetaH_re_set1 = new TH2I("Erel_costhetaH_re_set1","",100,0,8,25,-1,1);
  Erel_costhetaH_re_set1->SetMinimum(1.0);
  TH2I * Erel_costhetaH_re_set1a = new TH2I("Erel_costhetaH_re_set1a","",100,0,8,25,-1,1);
  Erel_costhetaH_re_set1a->SetMinimum(1.0);

  TH1I * Erel_long_re = new TH1I("Erel_long_re","",250,0,10);
  TH1I * Erel_long_re_set1 = new TH1I("Erel_long_re_set1","",250,0,10);
  TH1I * Erel_long_re_set1a = new TH1I("Erel_long_re_set1a","",250,0,10);

  TH1I * Erel_trans_re = new TH1I("Erel_trans_re","",250,0,10);
  TH1I * Erel_trans_re_set1 = new TH1I("Erel_trans_re_set1","",250,0,10);
  TH1I * Erel_trans_re_set1a = new TH1I("Erel_trans_re_set1a","",250,0,10);

  TH2I * gammasADDvsErel_re = new TH2I("gammasADDvsErel_re","",250,0,10,2048, 0, 8192);
  gammasADDvsErel_re->SetMinimum(1.0);

  TH1I * Erel_with1600_re = new TH1I("Erel_with1600_re","",250,0,10);
  TH1I * Ex_nogamma_re = new TH1I("Ex_nogamma_re","",250,0,10);

  file.cd();
  for (int i=0;i<N;i++)
  {
    if (i%3000 == 0)cout << i << endl;
      
    tree->GetEntry(i); //Pull data from tree

    //Uses CsI 0,0 except for runs with a gain of 90
    if (itele1 == 0 && id1 == 0 && runnum <= 131) continue;

    bool set1 = true; //true for setting 1, false for setting1a
    if (runnum < 216) gobbidist = 55.2;
    if (runnum >= 216)
    {
      gobbidist = 43.896;
      set1 = false;
    }

    //Fill time spectra 
    CsITvsEx->Fill(Erel,timeCsI1);
    if (timeCsI1 >= 975 && timeCsI1 <= 1025) Erel_21Al_p20Mg_tgate->Fill(Erel);

    gammasADD_mult->Fill(Ngamma);

    Erel_21Al_p20Mg->Fill(Erel);
    if (set1 == true) Erel_21Al_p20Mg_set1->Fill(Erel);
    else Erel_21Al_p20Mg_set1a->Fill(Erel);

    cout << beamZ << endl;
    if (beamZ == 14) Erel_21Al_p20Mg_Sibeam->Fill(Erel);
    if (beamZ == 13) Erel_21Al_p20Mg_Albeam->Fill(Erel);
    if (beamZ == 12) Erel_21Al_p20Mg_Mgbeam->Fill(Erel);

    Ex_21Al_p20Mg->Fill(Ex);
    ThetaCM_21Al_p20Mg->Fill(thetaCM*180./pi);
    VCM_21Al_p20Mg->Fill(Vcm);
    p20Mg_VCMvsErel->Fill(Vcm,Erel);
    Erel_p20Mg_costhetaH->Fill(Erel,cos_thetaH);

    if (Vcm >= 12.7 && Vcm <= 12.8) Erel_21Al_p20Mg_vfilt->Fill(Erel);

    //ThetaH gates
    if (fabs(cos_thetaH) <= 0.4)
    {
      Erel_21Al_p20Mg_trans->Fill(Erel);
      Ex_21Al_p20Mg_trans->Fill(Ex);

      if (set1 == true) Erel_21Al_p20Mg_trans_set1->Fill(Erel);
      else Erel_21Al_p20Mg_trans_set1a->Fill(Erel);
    }
    if (fabs(cos_thetaH) <= 0.1)
    {
      Erel_21Al_p20Mg_trans_narr->Fill(Erel);
    }
    if (fabs(cos_thetaH) >= 0.8)
    {
      Erel_21Al_p20Mg_long->Fill(Erel);
      Ex_21Al_p20Mg_long->Fill(Ex);
    }
    if (fabs(cos_thetaH) <= 0.2)
    {
      Erel_21Al_p20Mg_trans0_2->Fill(Erel);
      if (set1 == false) Erel_21Al_p20Mg_trans0_2_set1a->Fill(Erel);
    }

    //Gamma coincidence gates
    int with1600 = 0;
    int withleft = 0;
    int withright = 0;
    for (int j=0;j<Ngamma;j++) {
      if (Tgamma[j] < 685 || Tgamma[j] > 730) continue;
      if (Egamma[j] >= 1.45 && Egamma[j] <= 1.7) with1600 = 1;
      if (Egamma[j] >= 1. && Egamma[j] < 1.4) withleft = 1; //region left of peak
      if (Egamma[j] > 1.7 && Egamma[j] <= 2.2) withright = 1;
    } 

    if (with1600 == 1) Erel_21Al_p20Mg_with1600->Fill(Erel);
    if (with1600 == 1 && Ngamma == 1) Erel_21Al_p20Mg_with1600_mult1->Fill(Erel);
    if (with1600 == 1 && set1 == false) Erel_21Al_p20Mg_with1600_set1a->Fill(Erel);
    if (with1600 == 1 && fabs(cos_thetaH) <= 0.4) Erel_21Al_p20Mg_with1600_trans->Fill(Erel);
    if (with1600 == 0) Erel_21Al_p20Mg_no1600->Fill(Erel);
    if (with1600 == 0 && set1 == false) Erel_21Al_p20Mg_no1600_set1a->Fill(Erel);
    if (with1600 == 0 && fabs(cos_thetaH) <= 0.4) Erel_21Al_p20Mg_no1600_trans->Fill(Erel);
    if (Ngamma == 0) Erel_21Al_p20Mg_nogamma->Fill(Erel);
    if (Ngamma == 0 && fabs(cos_thetaH) <= 0.4) Erel_21Al_p20Mg_nogamma_trans->Fill(Erel);
    if (Ngamma > 0) Erel_21Al_p20Mg_withgamma->Fill(Erel);
    if (Ngamma > 0 && fabs(cos_thetaH) <= 0.4) Erel_21Al_p20Mg_withgamma_trans->Fill(Erel);

    if (withleft == 1) Erel_21Al_p20Mg_left1600->Fill(Erel);
    if (withright == 1) Erel_21Al_p20Mg_right1600->Fill(Erel);

    //Gammas - already doppler shifted
    for (int j=0;j<Ngamma;j++)
    {
      if (Tgamma[j] < 685 || Tgamma[j] > 730) continue;
      p20Mg_gammasADD->Fill(Egamma[j]*1000.);
      if (Ngamma == 1) p20Mg_gammasADD_mult1->Fill(Egamma[j]*1000.);
      p20Mg_gammasADDvsErel->Fill(Erel,Egamma[j]*1000.);
      if (fabs(cos_thetaH) <= 0.4) p20Mg_gammasADDvsErel_trans->Fill(Erel,Egamma[j]*1000.);
      gammaTvsE->Fill(Egamma[j]*1000.,Tgamma[j]);

      p20Mg_gammasCh->Fill(Chgamma[j]);
      p20Mg_gammasADDvsCh->Fill(Chgamma[j],Egamma[j]*1000.);
      if (Ngamma == 1) p20Mg_gammasADDvsCh_mult1->Fill(Chgamma[j],Egamma[j]*1000.);
     
      for (int k=j+1;k<Ngamma;k++)
      {
        p20Mg_gammasADDvsgammasADD->Fill(Egamma[j]*1000,Egamma[k]*1000.);
      }
    }

    //Apply calibration
    double denergy = FrontEcal.getEnergy(itele1, ifront1, denergyR);
    double energy = CsIEcal.getEnergy(itele1, id1, energyR);

    //Do fine calibration here, modify CsI "energy" by some a + b(x-x_0) where x_0 is avg p E
    CsIE->Fill(energy);
    double avgE_p = 77.;
    float a = 0.0;
    float b = 0.0;

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

    //Recalculate theta/phi based on closer fiber distance
    double fibdist_og = gobbidist + 10.151+2.664;
    double fibdist = gobbidist + 10.557;
    double x_fib = fibdist_og*tan(theta2)*cos(phi2);
    double y_fib = fibdist_og*tan(theta2)*sin(phi2);
    double r_fib = sqrt(pow(fibdist,2)+pow(x_fib,2)+pow(y_fib,2));
    theta2 = acos(fibdist/r_fib);
    phi2 = atan2(y_fib,x_fib);
    double pc_H = sqrt(pow(et2,2) - pow(Mass_20Mg,2));
    M2[0] = pc_H*sin(theta2)*cos(phi2);
    M2[1] = pc_H*sin(theta2)*sin(phi2);
    M2[2] = pc_H*cos(theta2);

    CM = Kin.findCM(M1,M2,et1,et2);
    VCM_re->Fill(CM.vv);
    ee1 = Kin.trans(CM,M1,et1);
    ee2 = Kin.trans(CM,M2,et2);
    Erel = ee1.vv + ee2.vv - Mass_20Mg - Mass_p;

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
    Ex = Erel - (Mass_21Al - (Mass_20Mg + Mass_p));
    Ex_re->Fill(Ex);

    if (set1 == true)
    {
      Erel_re_set1->Fill(Erel);
      Ex_re_set1->Fill(Ex);
      Erel_costhetaH_re_set1->Fill(Erel,cos_thetaH);
      if (fabs(cos_thetaH) <= 0.4) Erel_trans_re_set1->Fill(Erel);
      if (fabs(cos_thetaH) >= 0.7) Erel_long_re_set1->Fill(Erel);
    }
    else
    {
      Erel_re_set1a->Fill(Erel);
      Ex_re_set1a->Fill(Ex);
      Erel_costhetaH_re_set1a->Fill(Erel,cos_thetaH);
      if (fabs(cos_thetaH) <= 0.4) Erel_trans_re_set1a->Fill(Erel);
      if (fabs(cos_thetaH) >= 0.7) Erel_long_re_set1a->Fill(Erel);
    }

    Erel_costhetaH_re->Fill(Erel,cos_thetaH);

    if (fabs(cos_thetaH) <= 0.4) Erel_trans_re->Fill(Erel);
    if (fabs(cos_thetaH) >= 0.7) Erel_long_re->Fill(Erel);

    if (Ngamma == 0) Ex_nogamma_re->Fill(Ex);
    if (with1600 == 1) Erel_with1600_re->Fill(Erel);

    //Gammas - already doppler shifted
    for (int j=0;j<Ngamma;j++)
    {
      gammasADDvsErel_re->Fill(Erel,Egamma[j]*1000.);
    }


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
