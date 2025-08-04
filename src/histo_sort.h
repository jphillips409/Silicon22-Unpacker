#ifndef histo_sort_
#define histo_sort_
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "TH1I.h"
#include "TH2I.h"
#include "TH3I.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "CAENd5202.h" //Needs to have own instances of 5202 and fiber to work
#include "fiber.h"

using namespace std;

class histo_sort
{
 protected:
  
  TFile * file; //!< output root file

  //Fiber instances needed
	Event event;
	fiber fib;
	Event red;
	Event blue;

  //de
  TDirectoryFile * dirdee; //!< directory for dee maps
  

  //fiberArray
  TDirectoryFile * dirFiber; //!< directory for fiber array
  //Directories for 1d TOA and TOT plots for both boards
  TDirectory * dir1dBoard0TOA_R;
  TDirectory * dir1dBoard0TOT_R;

  TDirectory * dir1dBoard1TOA_R;
  TDirectory * dir1dBoard1TOT_R;

  TDirectory* dirHitMap; // directory for plotting xy-correlated event parameters

  //Gobbi
  TDirectoryFile * dirSummary;
  TDirectory * dir1dFrontE_R;
  TDirectory * dir1dBackE_R;
  TDirectory * dir1dNewBackE_R;

  TDirectory * dir1dFrontlowE_R;
  TDirectory * dir1dBacklowE_R;

  TDirectory * dir1dFrontE_cal;
  TDirectory * dir1dBackE_cal;
  TDirectory * dir1dFrontlowE_cal;
  TDirectory * dir1dBacklowE_cal;

  TDirectory * dir1dFrontTime_R;
  TDirectory * dir1dBackTime_R;

  TDirectory * dir1dCsI_Energy;
  TDirectory * dir1dCsI_Time;
  TDirectory * dir1dCsI_QDC;

  TDirectory * dirpTheta_Etot; //CsI theta dependence for proton

  TDirectoryFile * dirDEEplots; //!< directory for deltaE-E plots used in particle identificaiton
  TDirectory * dirPSD; //!< directory for PSD plots
  TDirectoryFile * dirhitmaps; //!< directory for all particle type hitmaps


  TDirectory * dirRSum;

  //Summary
  TDirectoryFile * dirSum; //!< directory for summary spectra

  
  //CAESAR
  TDirectory *caeDir;
  TDirectory *scalDir;
  TDirectory *summaryDir;
  TDirectory *sumDir;
  TDirectory *multDir;
  TDirectory *adcDir[6];
  TDirectory *detN[192];


  //S800 directory
  TDirectoryFile *dirS800;
  TDirectory * dirS800Raw;
  TDirectory * dirS800Cal;
  TDirectory * dirS800PID;
  TDirectory * dirS800Trig;
  TDirectory * dirS800veldist;

 public:
  histo_sort(string suffix);                  //!< constructor
  ~histo_sort(){};
  void clear(); //This and next 2 needed for Janus
	Event* GetRedEvent() { return &red; }
	Event* GetBlueEvent() { return &blue; }
  void write(); //!< write the root spectra to file
  

  //saved scalar values
  //TTree* tree;

  //Delta t between different rings
  TH1D * DeltaT_VME_Janus;
  TH1D * DeltaT_VME_S800;
  TH1D * DeltaT_S800_Janus;

  //T diff vs event number
  TH2D * DeltaT_VME_Janus_vsEvtCnt;
  TH2D * DeltaT_VME_S800_vsEvtCnt;
  TH2D * DeltaT_S800_Janus_vsEvtCnt;

  //Histogram the number of words in S800 events
  TH1D * S800_NumWords_R;
  TH1D * S800_NumWords_Accepted;

  //2D for e1up time and obj time vs gobbi copies
  TH2D * DB5T_vs_gCopy;
  TH2D * objT_vs_gCopy;
  TH2D * RFT_vs_gCopy;
  TH1I * S800_RFTime;
  TH1I * S800_RFTime_coin;
  TH1I * S800_RFTime_wbeam;
  TH1I * S800_RFTime_Sibeam;
  TH1I * S800_RFTime_Albeam;
  TH1I * S800_RFTime_Mgbeam;
  TH1I * S800_RFTime_Nabeam;
  TH1I * Gobbi_RFTime;
  
  int Ntele;
  int Nstrip;
  int Ncaesar;
  int NCsI;

  //Gobbi summaries
  //all Gobbi summaries
  //Energies, Raw+Calibrated
  TH2I * sumFrontE_R;
  TH2I * sumBackE_R;
  TH2I * sumFrontlowE_R;
  TH2I * sumBacklowE_R;

  TH2I * sumFrontE_cal;
  TH2I * sumBackE_cal;
  TH2I * sumFrontlowE_cal;
  TH2I * sumFrontlowE_cal_thetacorr;
  TH2I * sumBacklowE_cal;
  TH2I * sumBacklowE_cal_thetacorr;

  TH2I * sumCsIE_R;
  TH2I * sumCsIE_cal;

  //Comparing b1 strip 2 zline energies
  TH2I * b1s2Zline;
  TH2I * b1s2FvsB;

  //Gobbi times
  TH2I * sumFrontTime_R;
  TH2I * sumFrontTime_cal;
  TH2I * sumBackTime_R;
  TH2I * sumBackTime_cal;
  TH2I * sumCsITime_R;
  TH2I * sumCsITime_cal;


  //create all Gobbi 1d spectra
  //Energies  (quadrant, channum)
  TH1I * FrontE_R[4][32];
  TH1I * FrontElow_R[4][32];
  TH1I * FrontTime_R[4][32];
  TH1I * FrontE_cal[4][32];
  TH1I * FrontElow_V[4][32];

  TH1I * FrontlowE_cal[4][32]; //I only care about W low gain
  TH1I * FrontlowE_cal_thetacorr[4][32]; //I only care about W low gain
  TH1I * BacklowE_cal[4][32];
  TH1I * BackE_R[4][32];
  TH1I * BackElow_R[4][32];
  TH1I * BackTime_R[4][32];
  TH1I * BackE_cal[4][32];

  //create all the CsI 1d spectra
  TH1I * CsI_Energy_R_um[16];
  TH1I * CsI_Energy_R[16];
  TH1I * CsI_Energy_cal[16];
  TH1I * CsI_Energyp_cal[16];
  TH1I * CsI_Time_R_um[19]; //Plugged into the TDC, chans 1-3 are S800 DB5 (called XFP in code), object, and RF times
  TH1I * CsI_Time_R[19];    //chans 16-31 are CsI. Remapped when filling so that 0-15 are CsI.
  TH1I * CsI_Time_cal[19];
  TH1I * CsI_QDC_R[32];
  TH1I * CsI_QDC_matched[32];
  //Take the center of each CsI
  TH1I * CsI_Energy_R_center[16];
  TH1I * CsI_Energy_protonBE[16]; //proton beam energy


  //DEE
  //DeltaE-E plots
  TH1I * timediff[4];
  TH2I * DEE_CsI[4][4];
  TH2I * DEE_CsI_S800Coinc[4][4];
  TH2I * DEE_CsI_S800Coinc_Si23[4][4];
  TH2I * DEE_CsI_0deg[4][4];
  TH1I * timediff_CsI[4][4];
  TH2I * DEE_lowgain[4];
  TH2I * DEE_CsI_lowgain[4][4];

  TH2I * PSD_CsI_unmatched[4][4];
  TH2I * PSD_CsI_matched[4][4];

  TH2I * AnomalousStrips[4];

  TH2I * GFBE[4];
  TH2I * GFBE_Matched[4];
  TH2I * GFBE_strip[16];
  TH2I * GFEvsChan_strip[16];
  TH2I * FrontElow_V_Strip[16];
  TH2I * GdEvsChan_strip[16]; //Is this for dE or for CsI? 

  TH2I * GFrontStrip_CsI[4];
  TH2I * GBackStrip_CsI[4];

  //Hit maps
  TH2I * testinghitmap;
  TH2I * testWhit;
  TH2I * xyhitmap;
  TH2I * protonhitmap;
  TH2I * tphitmap;
  TH2I * CsIHitMap[4];

  
  //fiber
	TH1I* toa_hist;
	TH2I* tot_summary_blue;
	TH2I* tot_summary_red;
	TH2I* tot_summary_blue_matched;
	TH2I* tot_summary_red_matched;

  //1d spectra
  TH1I * B0TOA_R[64];
  TH1I * B0TOT_R[64];

  TH1I * B1TOA_R[64];
  TH1I * B1TOT_R[64];

	TH2I* Fiber_ixiy;
	TH2I* Fiber_xy;
	TH2I* Fiber_tot_summary_x;
	TH2I* Fiber_tot_summary_y;
	TH1F* Fiber_totx;
	TH1F* Fiber_toty;
  TH1F* Fiber_postotx;
  TH1F* Fiber_postoty;
  TH1F* Fiber_postoax;
  TH1F* Fiber_postoay;
  TH1I* Fiber_toax;
  TH1I* Fiber_toay;
  
  TH1I * Fiber_X;
  TH1I * Fiber_Y;
  TH2I * Fiber_XY;
  TH1I * Fiber_Xid;
  TH1I * Fiber_Yid;
  TH2I * Fiber_XYid;
  TH1I * Fiber_XBeam;
  TH1I * Fiber_YBeam;
  TH2I * Fiber_XYBeam;

  TH1I * Txfp;
  TH1I * TRF;
  
  TH2I * PCSum;
  TH2I * PCSum_AfterAddback;
  TH2I * PCLSum;
  TH2I * RCSum;
  TH2I * RCSum_AfterAddback;
  TH2I * PSum;
  TH2I * RSum;

  TH2I * PievsRing;
  TH2I * EdiffvsPie;
  TH1I * PiesMult;
  TH1I * RingsMult;

  TH1I* S800_Csi_time;
  TH1I* S800_Csi_time_with_proton;
  TH1I* singles_trig_time;

  TH1I * T_RFSCIN;
  TH1I * T_RFCYC;
  TH1I * T_A1900;
    
  TH2I **  dee;
  TH2I** dee_S800;

  //CAESAR
  TH1F *ECaesar[192];
  TH2F *ETCaesar[192];
  TH1F *TCaesar[192];
  TH2F *TECaesar[192];
  TH2I *mult[192];
  TH2F *energySum;
  TH2F *energyiSum;
  TH2F *energyTSum;
  TH2F *timeSum;
  TH1F *energyTot;
  TH2F *energyTTotvsT;
  TH1F *energyTTot_tgated;
  TH1F *energyTTot;
  TH1F *energyT0Tot;
  TH1F *energyT1Tot;
  TH1F *energyT2Tot;
  TH1F *enAddback_tgated;
  TH1F *enAddback;
  TH2F *enAddbackvsT;
  TH1F *energyM;
  TH1F *timeM;
  TH2F *adctdcM;

  TH1F *enAddback_Si24;
  TH2F *enAddback_Si24vsT;

  TH1F *enAddback_Si23;
  TH1F *enAddback_Si23_mult1;
  TH2F *enAddback_Si23vsT;

  TH1F *enAddback_Al22;
  TH1F *enAddback_Al22_antiTG;
  TH1F *enAddback_Al22_mult1;
  TH1F *enAddback_Al22_mult1_antiTG;
  TH2F *enAddback_Al22vsT;
  TH2F *enAddback_Al22vsT_mult1;
  TH2F *enAddback_Al22vsCh;
  TH1F *enSelect_Al22;

  TH1F *energyTTot_Mg20;
  TH1F *energyTTot_Mg21;


  //S800

  //Raw
  TH1I * e1up;
  TH1I * e1down;
  TH2I * e1upvsQDC;
  TH1I * Te1up;
  TH1I * Te1down;
  TH1I * Tobj;
  TH1I * Tobj_mult;

  TH2I * ICSummary;
  TH2I * CRDC1Summary;
  TH2I * CRDC2Summary;
  TH2I * CRDC1Summary_cal;
  TH2I * CRDC2Summary_cal;
  TH1I * CRDC1raw;
  TH1I * CRDC2raw;
  TH1I * CRDC1PadMult;
  TH1I * CRDC2PadMult;
  TH1I * CRDC1Tac;
  TH1I * CRDC2Tac;
  TH2I * CRDC1AnodevsTac;
  TH2I * CRDC2AnodevsTac;
  TH2I * CRDC1XrawY;
  TH2I * CRDC2XrawY;

  //Cal
  TH1I * CRDC1X;
  TH1I * CRDC1Y;
  TH1I * CRDC2X;
  TH1I * CRDC2Y;
  TH2I * CRDC1XY;
  TH2I * CRDC2XY;

  TH2I * CRDC1XCalY;
  TH2I * CRDC2XCalY;

  TH1I * ata1D;
  TH1I * bta1D;
  TH2I * atavsbta;
  TH2I * ThetavsPhi;
  TH2I * ThetavsPhi_F17;
  TH2I * ThetavsPhi_F17_fiber;
  TH2I * ThetavsPhi_fiber;
  TH2I * ThetaFibervsThetaS800;
  TH2I * PhiFibervsPhiS800;

  TH1I * CRDC1X_K35;
  TH1I * CRDC1X_K36;


  //Cal/CRDC
  TH2I * CRDC1_Ca36_Ca36;
  TH2I * CRDC2_Ca36_Ca36;
  TH2I * CRDC1_Ca36_K35;
  TH2I * CRDC2_Ca36_K35;
  TH2I * CRDC1_K35_K35;
  TH2I * CRDC2_K35_K35;
  TH2I * CRDC1_K35_K36;
  TH2I * CRDC2_K35_K36;
  TH2I * CRDC1_Ca36_Ar33;
  TH2I * CRDC2_Ca36_Ar33;
  TH2I * CRDC1_Ca36_Ar34;
  TH2I * CRDC2_Ca36_Ar34;
  TH2I * CRDC1_K35_Ar34;
  TH2I * CRDC2_K35_Ar34;
  TH2I * CRDC1_K35_Ar35;
  TH2I * CRDC2_K35_Ar35;
  TH2I * CRDC1_Ca36_Cl32;
  TH2I * CRDC2_Ca36_Cl32;
  
  TH1F * crdcx23Si_23Si; //To check unreacted beam centering
  TH1F * crdcx23S_20Mg_2p;
  TH1F * crdcx23P_22Si_p;

  //PID
  TH2I * ObjvsXFP;
  TH2I * ObjvsXFPwithAlpha1;
  TH2I * ObjvsXFPwithProton1;
  TH2I * ObjvsXFPwithProton2;  
  TH2I * ObjvsICsum_nobeam;
  TH2I * ObjvsICsum;
  TH2I * ObjvsICsum_corr;
  TH2I * ObjvsICsum_Si23; 
  TH2I * ObjvsICsum_Si23_1p; 
  TH2I * ObjvsICsum_Si23_2p; 
  TH2I * ObjvsICsum_Si23_3p; 
  TH2I * ObjvsICsum_Si23_alpha; 
  TH2I * ObjvsICsum_Al22; 
  TH2I * ObjvsICsum_Al22_1p; 
  TH2I * ObjvsICsum_Al22_2p; 
  TH2I * ObjvsICsum_Al22_3p; 
  TH2I * ObjvsICsum_Al22_alpha; 
  TH2I * ObjvsICsum_Mg21; 
  TH2I * ObjvsICsum_Mg21_1p; 
  TH2I * ObjvsICsum_Mg21_2p; 
  TH2I * ObjvsICsum_Mg21_3p; 
  TH2I * ObjvsICsum_Mg21_alpha; 
  TH2I * ObjvsICsum_Na20; 
  TH2I * ObjvsICsum_Na20_1p; 
  TH2I * ObjvsICsum_Na20_2p; 
  TH2I * ObjvsICsum_Na20_3p;
  TH2I * ObjvsICsum_Na20_alpha; 
  TH2I * ObjvsICsum_Ca37; 
  TH2I * ObjUncvsICsum_Ca37;
  TH2I * Timing1vsICsum_Ca37;
  TH2I * Timing2vsICsum_Ca37;
  TH2I * XFPvsICsum_Ca37;
  TH2I * ObjvsICsum_K36;
  TH2I * ObjvsICsum_Ar35;
  TH2I * ObjvsICsum_Cl34;
  TH2I * ObjvsICsum_wProt;
  TH2I * ObjvsICsum_wProt_wFiber;
  TH2I * ObjvsICsum_wProt2;
  TH2I * ObjvsICsum_K36_wProt;
  TH2I * ObjvsICsum_Ar35_wProt;

  TH2I * Objvsafp;
  TH2I * ObjvsCRDC1X;
  TH2I * ObjvsCRDC2X;
  TH2I * XfpUncvsCRDC1X;
  TH2I * XfpUncvsCRDC2X;

  TH1I * energySi23_Si23;
  TH1I * energySi23_Si23_Fib;
  TH1I * energySi23_Si23_FibTarg;

  TH1I * energyMg21_Mg21;
  TH1I * energyMg21_Mg21_Fib;
  TH1I * energyMg21_Mg21_FibTarg;

  TH1I * rigidityK35;
  TH1I * rigidityCa36;
  TH1I * rigidityCa37beam;
  TH1I * rigidityCa38beam;
  TH1I * rigidityK36beam;
  TH1I * rigidityAr35beam;

  TH1I * BeamCa37_energy;
  TH1I * BeamCa37_ptot;
  TH1I * BeamCa37_ppar;
  TH1I * BeamCa37_ptra;
  TH2I * BeamCa37_TofAngle;
  TH2I * BeamCa37_TofCRDC1raw;
  TH2I * BeamCa37_TofCRDC1cal;
  TH2I * BeamCa37_ObjvsCRDC1X;
  TH2I * BeamCa37_Fiber_XY;
  TH1I * BeamCa38_energy;
  TH1I * BeamCa38_ptot;
  TH1I * BeamCa38_ppar;
  TH1I * BeamCa38_ptra;
  TH2I * BeamCa38_TofAngle;
  TH2I * BeamCa38_TofCRDC1raw;
  TH2I * BeamCa38_TofCRDC1cal;
  TH2I * BeamCa38_ObjvsCRDC1X;
  TH2I * BeamCa38_Fiber_XY;

  TH1I * BeamK36_energy;
  TH1I * BeamK36_ptot;
  TH1I * BeamK36_ppar;
  TH1I * BeamK36_ptra;
  TH1I * BeamAr35_energy;
  TH1I * BeamAr35_ptot;
  TH1I * BeamAr35_ppar;
  TH1I * BeamAr35_ptra;

  //veldist
  TH1I * residuetheta;

  TDirectoryFile * dirInvMass; //!< directory for all correlations and inv-mass
  TDirectory * dirCorrCombs;
  TDirectory * dir1H;
  TDirectory * dir16F;
  TDirectory * dir17F;
  TDirectory * dir19Ne;
  TDirectory * dir18Ne;
  TDirectory * dir19Na;
  TDirectory * dir20Na;
  TDirectory * dir21Na;
  TDirectory * dir19Mg;
  TDirectory * dir20Mg;
  TDirectory * dir21Mg;
  TDirectory * dir2p19Ne;
  TDirectory * dir22Mg;
  TDirectory * dir21Al;
  TDirectory * dir22Al;
  TDirectory * dir23Al;
  TDirectory * dir22Si;
  TDirectory * dir23Si;
  TDirectory * dir24Si;
  TDirectory * dirp22Al;
  TDirectory * dir2p21Mg;
  TDirectory * dir23P;
  TDirectory * dir3p20Mg;
  TDirectory * dirp22Si;
  TDirectory * dir24P;
  TDirectory * dir3p18Ne;
  TDirectory * dir3p17Ne;

  //********************************************************
  //Correlation combs
  TH2I * SiBeam_CorrComb;
  TH2I * AlBeam_CorrComb;
  //********************************************************

  //********************************************************
  //Correlation hists
  //********************************************************
  //protons  
  TH1I * protonKE;
  TH1I * protonKE_Si23Beam;
  TH1I * protonKE_Al22Beam;
  TH1I * protonKE_Mg21Beam;

  //F16 -> p + 15O
  TH1I * Erel_16F_p15O;
  TH1I * Erel_16F_p15O_trans;
  TH1I * Erel_16F_p15O_long;
  TH1I * Ex_16F_p15O;
  TH1I * Ex_16F_p15O_trans;
  TH1I * Ex_16F_p15O_long;
  TH1I * ThetaCM_16F_p15O;
  TH1I * VCM_16F_p15O;
  TH2I * Erel_p15O_costhetaH;
  TH2I * p15O_VCMvsErel;
  TH1I * p15O_gammasADD;
  TH2I * p15O_gammasADDvsErel;
  TH2I * p15O_gammasADDvsgammasADD;

  //F17 -> p + 16O
  TH1I * Erel_17F_p16O;
  TH1I * Erel_17F_p16O_trans;
  TH1I * Erel_17F_p16O_long;
  TH1I * Ex_17F_p16O;
  TH1I * Ex_17F_p16O_trans;
  TH1I * Ex_17F_p16O_long;
  TH1I * ThetaCM_17F_p16O;
  TH1I * VCM_17F_p16O;
  TH2I * Erel_p16O_costhetaH;
  TH2I * p16O_VCMvsErel;
  TH1I * p16O_gammasADD;
  TH2I * p16O_gammasADDvsErel;
  TH2I * p16O_gammasADDvsgammasADD;

  //Ne19 -> p + 18F
  TH1I * Erel_19Ne_p18F;
  TH1I * Erel_19Ne_p18F_trans;
  TH1I * Erel_19Ne_p18F_long;
  TH1I * Ex_19Ne_p18F;
  TH1I * Ex_19Ne_p18F_trans;
  TH1I * Ex_19Ne_p18F_long;
  TH1I * ThetaCM_19Ne_p18F;
  TH1I * VCM_19Ne_p18F;
  TH2I * Erel_p18F_costhetaH;
  TH2I * p18F_VCMvsErel;
  TH1I * p18F_gammasADD;
  TH2I * p18F_gammasADDvsErel;
  TH2I * p18F_gammasADDvsgammasADD;

  //Ne18 -> p + 17F
  TH1I * Erel_18Ne_p17F;
  TH1I * Erel_18Ne_p17F_trans;
  TH1I * Erel_18Ne_p17F_long;
  TH1I * Ex_18Ne_p17F;
  TH1I * Ex_18Ne_p17F_trans;
  TH1I * Ex_18Ne_p17F_long;
  TH1I * Erel_18Ne_p17F_set1;
  TH1I * Erel_18Ne_p17F_set1a;
  TH1I * ThetaCM_18Ne_p17F;
  TH1I * VCM_18Ne_p17F;
  TH2I * Erel_p17F_costhetaH;
  TH2I * p17F_VCMvsErel;
  TH1I * p17F_gammasADD;
  TH1I * p17F_gammasADD_tgate;
  TH2I * p17F_gammasADDvsErel;
  TH2I * p17F_gammasADDvsgammasADD;

  //Na19 -> p + 18Ne
  TH1I * Erel_19Na_p18Ne;
  TH1I * Erel_19Na_p18Ne_trans;
  TH1I * Erel_19Na_p18Ne_long;
  TH1I * Ex_19Na_p18Ne;
  TH1I * Ex_19Na_p18Ne_trans;
  TH1I * Ex_19Na_p18Ne_long;
  TH1I * Erel_19Na_p18Ne_set1;
  TH1I * Erel_19Na_p18Ne_set1a;
  TH1I * ThetaCM_19Na_p18Ne;
  TH1I * VCM_19Na_p18Ne;
  TH2I * Erel_p18Ne_costhetaH;
  TH2I * p18Ne_VCMvsErel;
  TH1I * p18Ne_gammasADD;
  TH1I * p18Ne_gammasADD_tgate;
  TH2I * p18Ne_gammasADDvsErel;
  TH2I * p18Ne_gammasADDvsgammasADD;

  //Na20 -> p + 19Ne
  TH1I * Erel_20Na_p19Ne;
  TH1I * Erel_20Na_p19Ne_trans;
  TH1I * Ex_20Na_p19Ne;
  TH1I * ThetaCM_20Na_p19Ne;
  TH1I * VCM_20Na_p19Ne;
  TH2I * Erel_p19Ne_costhetaH;
  TH2I * p19Ne_VCMvsErel;
  TH1I * p19Ne_gammasADD;
  TH1I * p19Ne_gammasADD_tgate;
  TH2I * p19Ne_gammasADDvsErel;
  TH2I * p19Ne_gammasADDvsgammasADD;

  //Na21 -> p + 20Ne
  TH1I * Erel_21Na_p20Ne;
  TH1I * Erel_21Na_p20Ne_trans;
  TH1I * Ex_21Na_p20Ne;
  TH1I * ThetaCM_21Na_p20Ne;
  TH1I * VCM_21Na_p20Ne;
  TH2I * Erel_p20Ne_costhetaH;
  TH2I * p20Ne_VCMvsErel;
  TH1I * p20Ne_gammasADD;
  TH1I * p20Ne_gammasADD_tgate;
  TH2I * p20Ne_gammasADDvsErel;
  TH2I * p20Ne_gammasADDvsgammasADD;

  //Mg19 -> 2p + 17Ne
  TH1I * Erel_19Mg_2p17Ne;
  TH1I * Erel_19Mg_2p17Ne_trans;
  TH1I * Ex_19Mg_2p17Ne;
  TH1I * ThetaCM_19Mg_2p17Ne;
  TH1I * VCM_19Mg_2p17Ne;
  TH2I * Erel_2p17Ne_costhetaH;
  TH2I * pp17Ne_VCMvsErel;
  TH1I * pp17Ne_gammasADD;
  TH1I * pp17Ne_gammasADD_tgate;
  TH2I * pp17Ne_gammasADDvsErel;
  TH2I * pp17Ne_gammasADDvsgammasADD;

  //Mg20 -> 2p + 18Ne
  TH1I * Erel_20Mg_2p18Ne;
  TH1I * Erel_20Mg_2p18Ne_trans;
  TH1I * Ex_20Mg_2p18Ne;
  TH1I * ThetaCM_20Mg_2p18Ne;
  TH1I * VCM_20Mg_2p18Ne;
  TH2I * Erel_2p18Ne_costhetaH;
  TH2I * pp18Ne_VCMvsErel;
  TH1I * pp18Ne_gammasADD;
  TH1I * pp18Ne_gammasADD_tgate;
  TH2I * pp18Ne_gammasADDvsErel;
  TH2I * pp18Ne_gammasADDvsgammasADD;

  //Mg21 -> p + 20Na
  TH1I * Erel_21Mg_p20Na;
  TH1I * Erel_21Mg_p20Na_trans;
  TH1I * Erel_21Mg_p20Na_long;
  TH1I * Ex_21Mg_p20Na;
  TH1I * Ex_21Mg_p20Na_trans;
  TH1I * Ex_21Mg_p20Na_long;
  TH1I * Erel_21Mg_p20Na_set1;
  TH1I * Erel_21Mg_p20Na_set1a;
  TH1I * ThetaCM_21Mg_p20Na;
  TH1I * VCM_21Mg_p20Na;
  TH2I * Erel_p20Na_costhetaH;
  TH2I * p20Na_VCMvsErel;
  TH1I * p20Na_gammasADD;
  TH1I * p20Na_gammasADD_tgate;
  TH2I * p20Na_gammasADDvsErel;
  TH2I * p20Na_gammasADDvsgammasADD;
  TH1I * Mg21_pKE;

  //Mg21 -> 2p + 19Ne
  TH1I * Erel_21Mg_2p19Ne;
  TH1I * Erel_21Mg_2p19Ne_trans;
  TH1I * Ex_21Mg_2p19Ne;
  TH1I * ThetaCM_21Mg_2p19Ne;
  TH1I * VCM_21Mg_2p19Ne;
  TH2I * Erel_2p19Ne_costhetaH;
  TH2I * pp19Ne_VCMvsErel;
  TH1I * pp19Ne_gammasADD;
  TH1I * pp19Ne_gammasADD_tgate;
  TH2I * pp19Ne_gammasADDvsErel;
  TH2I * pp19Ne_gammasADDvsgammasADD;
  TH1I * Mg21_2p_pKE;

  //Mg22 -> p + 21Na
  TH1I * Erel_22Mg_p21Na;
  TH1I * Erel_22Mg_p21Na_trans;
  TH1I * Ex_22Mg_p21Na;
  TH1I * ThetaCM_22Mg_p21Na;
  TH1I * VCM_22Mg_p21Na;
  TH2I * Erel_p21Na_costhetaH;
  TH2I * p21Na_VCMvsErel;
  TH1I * p21Na_gammasADD;
  TH1I * p21Na_gammasADD_tgate;
  TH2I * p21Na_gammasADDvsErel;
  TH2I * p21Na_gammasADDvsgammasADD;

  //Al21 -> p + 20Mg
  TH1I * Erel_21Al_p20Mg;
  TH1I * Erel_21Al_p20Mg_trans;
  TH1I * Ex_21Al_p20Mg;
  TH1I * ThetaCM_21Al_p20Mg;
  TH1I * VCM_21Al_p20Mg;
  TH2I * Erel_p20Mg_costhetaH;
  TH2I * p20Mg_VCMvsErel;
  TH1I * p20Mg_gammasADD;
  TH1I * p20Mg_gammasADD_tgate;
  TH2I * p20Mg_gammasADDvsErel;
  TH2I * p20Mg_gammasADDvsgammasADD;
  TH1I * Al21_pKE;

  //Al22 -> p + 21Mg
  TH1I * Erel_22Al_p21Mg;
  TH1I * Erel_22Al_p21Mg_trans;
  TH1I * Ex_22Al_p21Mg;
  TH1I * ThetaCM_22Al_p21Mg;
  TH1I * VCM_22Al_p21Mg;
  TH2I * Erel_p21Mg_costhetaH;
  TH2I * p21Mg_VCMvsErel;
  TH1I * p21Mg_gammasADD_nodopp;
  TH1I * p21Mg_gammasADD;
  TH1I * p21Mg_gammasADD_tgate;
  TH2I * p21Mg_gammasADDvsErel;
  TH1I * p21Mg_gammasADD_peak2;
  TH1I * p21Mg_gammasADD_peak2_tgate;
  TH1I * p21Mg_gammasADD_peak3;
  TH1I * p21Mg_gammasADD_peak3_tgate;
  TH2I * p21Mg_gammasADDvsgammasADD;
  TH1I * Al22_pKE;

  //Al23 -> p + 22Mg
  TH1I * Erel_23Al_p22Mg;
  TH1I * Erel_23Al_p22Mg_trans;
  TH1I * Erel_23Al_p22Mg_long;
  TH1I * Ex_23Al_p22Mg;
  TH1I * Ex_23Al_p22Mg_trans;
  TH1I * Ex_23Al_p22Mg_long;
  TH1I * Erel_23Al_p22Mg_set1;
  TH1I * Erel_23Al_p22Mg_set1a;
  TH1I * ThetaCM_23Al_p22Mg;
  TH1I * VCM_23Al_p22Mg;
  TH2I * Erel_p22Mg_costhetaH;
  TH2I * p22Mg_VCMvsErel;
  TH1I * p22Mg_gammasADD;
  TH1I * p22Mg_gammasADD_tgate;
  TH2I * p22Mg_gammasADDvsErel;
  TH2I * p22Mg_gammasADDvsgammasADD;

  //Si22 -> 2p + 20Mg
  TH1I * Erel_22Si_2p20Mg;
  TH1I * Erel_22Si_2p20Mg_trans;
  TH1I * Ex_22Si_2p20Mg;
  TH1I * ThetaCM_22Si_2p20Mg;
  TH1I * VCM_22Si_2p20Mg;
  TH2I * Erel_2p20Mg_costhetaH;
  TH2I * pp20Mg_VCMvsErel;
  TH1I * pp20Mg_gammasADD;
  TH1I * pp20Mg_gammasADD_tgate;
  TH2I * pp20Mg_gammasADDvsErel;
  TH2I * pp20Mg_gammasADDvsgammasADD;
  TH1I * Si22_pKE;
  //Jacobi
  TH2I * Si22_JacobiT_xy_s;
  TH2I * Si22_JacobiY_xy_s;

  //Si23 -> p + 22Al
  TH1I * Erel_23Si_p22Al;
  TH1I * Erel_23Si_p22Al_trans;
  TH1I * Ex_23Si_p22Al;
  TH1I * ThetaCM_23Si_p22Al;
  TH1I * VCM_23Si_p22Al;
  TH2I * Erel_p22Al_costhetaH;
  TH2I * p22Al_VCMvsErel;
  TH1I * p22Al_gammasADD;
  TH1I * p22Al_gammasADD_tgate;
  TH2I * p22Al_gammasADDvsErel;
  TH2I * p22Al_gammasADDvsgammasADD;
  TH1I * Si23_pKE;

  //Si23 -> 2p + 21Mg
  TH1I * Erel_23Si_2p21Mg;
  TH1I * Erel_23Si_2p21Mg_trans;
  TH1I * Ex_23Si_2p21Mg;
  TH1I * ThetaCM_23Si_2p21Mg;
  TH1I * VCM_23Si_2p21Mg;
  TH2I * Erel_2p21Mg_costhetaH;
  TH2I * pp21Mg_VCMvsErel;
  TH1I * pp21Mg_gammasADD;
  TH1I * pp21Mg_gammasADD_tgate;
  TH2I * pp21Mg_gammasADDvsErel;
  TH2I * pp21Mg_gammasADDvsgammasADD;

  //Si24 -> p + 23Al
  TH1I * Erel_24Si_p23Al;
  TH1I * Erel_24Si_p23Al_trans;
  TH1I * Ex_24Si_p23Al;
  TH1I * ThetaCM_24Si_p23Al;
  TH1I * VCM_24Si_p23Al;
  TH2I * Erel_p23Al_costhetaH;
  TH2I * p23Al_VCMvsErel;
  TH1I * p23Al_gammasADD;
  TH1I * p23Al_gammasADD_tgate;
  TH2I * p23Al_gammasADDvsErel;
  TH2I * p23Al_gammasADDvsgammasADD;

  //P23 -> 3p + 20Mg
  TH1I * Erel_23P_3p20Mg;
  TH1I * Erel_23P_3p20Mg_trans;
  TH1I * Ex_23P_3p20Mg;
  TH1I * ThetaCM_23P_3p20Mg;
  TH1I * VCM_23P_3p20Mg;
  TH2I * Erel_3p20Mg_costhetaH;
  TH2I * ppp20Mg_VCMvsErel;
  TH1I * ppp20Mg_gammasADD;
  TH1I * ppp20Mg_gammasADD_tgate;
  TH2I * ppp20Mg_gammasADDvsErel;
  TH2I * ppp20Mg_gammasADDvsgammasADD;

  //P23 -> p + 22Si
  TH1I * Erel_23P_p22Si;
  TH1I * Erel_23P_p22Si_trans;
  TH1I * Ex_23P_p22Si;
  TH1I * ThetaCM_23P_p22Si;
  TH1I * VCM_23P_p22Si;
  TH2I * Erel_p22Si_costhetaH;
  TH2I * p22Si_VCMvsErel;
  TH1I * p22Si_gammasADD;
  TH1I * p22Si_gammasADD_tgate;
  TH2I * p22Si_gammasADDvsErel;
  TH2I * p22Si_gammasADDvsgammasADD;

  //P24 -> p + 23Si
  TH1I * Erel_24P_p23Si;
  TH1I * Erel_24P_p23Si_trans;
  TH1I * Ex_24P_p23Si;
  TH1I * ThetaCM_24P_p23Si;
  TH1I * VCM_24P_p23Si;
  TH2I * Erel_p23Si_costhetaH;
  TH2I * p23Si_VCMvsErel;
  TH1I * p23Si_gammasADD;
  TH1I * p23Si_gammasADD_tgate;
  TH2I * p23Si_gammasADDvsErel;
  TH2I * p23Si_gammasADDvsgammasADD;

  //Al21 -> 3p + 18Ne
  TH1I * Erel_21Al_3p18Ne;
  TH1I * Erel_21Al_3p18Ne_trans;
  TH1I * Ex_21Al_3p18Ne;
  TH1I * ThetaCM_21Al_3p18Ne;
  TH1I * VCM_21Al_3p18Ne;
  TH2I * Erel_3p18Ne_costhetaH;
  TH2I * ppp18Ne_VCMvsErel;
  TH1I * ppp18Ne_gammasADD;
  TH1I * ppp18Ne_gammasADD_tgate;
  TH2I * ppp18Ne_gammasADDvsErel;
  TH2I * ppp18Ne_gammasADDvsgammasADD;

  //Al20 -> 3p + 17Ne
  TH1I * Erel_20Al_3p17Ne;
  TH1I * Erel_20Al_3p17Ne_trans;
  TH1I * Ex_20Al_3p17Ne;
  TH1I * ThetaCM_20Al_3p17Ne;
  TH1I * VCM_20Al_3p17Ne;
  TH2I * Erel_3p17Ne_costhetaH;
  TH2I * ppp17Ne_VCMvsErel;
  TH1I * ppp17Ne_gammasADD;
  TH1I * ppp17Ne_gammasADD_tgate;
  TH2I * ppp17Ne_gammasADDvsErel;
  TH2I * ppp17Ne_gammasADDvsgammasADD;

};
#endif
