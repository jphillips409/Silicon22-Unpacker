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

  //correlations
  //de
  TDirectoryFile * dirdee; //!< directory for dee maps
  

  //fiberArray
  TDirectoryFile * dirFiber; //!< directory for fiber array
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

  TDirectory * dirpTheta_Etot; //CsI theta dependence for proton

  TDirectoryFile * dirDEEplots; //!< directory for deltaE-E plots used in particle identificaiton
  TDirectoryFile * dirhitmaps; //!< directory for all particle type hitmaps


  TDirectory * dirRSum;

  //Summary
  TDirectoryFile * dirSum; //!< directory for summary spectra

  
  //CEARAR
  TDirectoryFile *dirC; //!< directory for the CEASAR info
  TDirectory * dirEnergy; //!< CsI energies
  TDirectory * dirTime; //!< CsI time
  TDirectory * dirTimeS800;
  TDirectory * dirEcal; //!< CsI Sum
  TDirectory * dirRsum; //!< CsI gated on Z
  TDirectory * dirCSum;
  TDirectory * dirCa36;


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
  
  int Ntele;
  int Nstrip;
  int Nceasar;
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
  TH1I * CsI_Time_R_um[16];
  TH1I * CsI_Time_R[16];
  TH1I * CsI_Time_cal[16];


  //DEE
  //DeltaE-E plots
  TH1I * timediff[4];
  TH2I * DEE_CsI[4][4];
  TH2I * DEE_CsI_0deg[4][4];
  TH1I * timediff_CsI[4][4];
  TH2I * DEE_lowgain[4];
  TH2I * DEE_CsI_lowgain[4][4];

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

  //CEASAR
  TH2I* TCeasarRawSum;
  TH2I* ECeasarRawSum;
  TH2I* ECeasarCalSum;
  TH2I* ECeasarDopSum;
  TH2I* TCeasarMatchedSum;
  TH2I* ECeasarMatchedRawSum;
  TH2I* ECeasarMatchedCalSum;
  TH2I* ECeasarRawSum_idet;
  TH2I* ECeasarCalSum_idet;
  TH2I* ECeasarDopSum_idet;
  TH2I* TCeasarMatchedSum_idet;
  TH2I* ECeasarMatchedRawSum_idet;
  TH2I* ECeasarMatchedCalSum_idet;
  TH1I** ECeasar;
  TH1I** TCeasar;
  TH1I** TCeasarS800;
  TH1I** ECCeasar;
  TH1I** Ca36gamma;
  TH1I** Ca36gamma_nodopp;
  TH1I* TECeasar;
  TH1I* TECeasar_difbin;
  TH1I* Egated;
  TH2I* TvsEgated;
  TH2I* TvsEchn42;
  TH2I* TvsEchn48;
  TH1I* TEC_Dop;

  TH2I* map_iring_iring;
  TH2I* map_idet_idet;
  TH2I* map_DataEC_id_id;
  TH2I* map_DataTC_id_id;
  TH1I* delta_phi_gamma;

  TH1I* CEMult;
  TH1I* CTMult;
  TH1I* TCeasarRaw1;
  TH1I* TCeasarRaw2;
  TH2I* TCeasarSum_gated;
  TH1I* TCeasarRawS800;
  TH1I* TCeasarCal1;
  TH1I* TCeasarCal2;

  TH1I* CETMult;
  TH2I* ETOF_Ceasar;


  //S800

  //Raw
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

  TH2I * atavsbta;
  TH2I * ThetavsPhi;
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
  
  //PID
  TH2I * ObjvsXFP;
  TH2I * ObjvsXFPwithAlpha1;
  TH2I * ObjvsXFPwithProton1;
  TH2I * ObjvsXFPwithProton2;  
  TH2I * ObjvsICsum;
  TH2I * ObjvsICsum_corr;
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
  TH1I * Vlab_Ca38;
  TH1I * Ppar_Ca38;
  TH1I * Ptra_Ca38;
  TH1I * Vlab_Ca37;
  TH1I * Ppar_Ca37;
  TH1I * Ptra_Ca37;
  TH1I * Elab_Ca37;
  TH1I * Etar_Ca37;
  TH1I * Etarfib_Ca37;
  TH1I * Efib_Ca37;

  TH1I * Vlab_Ca35;
  TH1I * Vlab_Ca36;
  TH1I * Vlab_Ca36_wgamma;
  TH1I * Vlab_K35;
  TH1I * Vlab_K36;
  TH1I * Vlab_Ar32;
  TH1I * Vlab_Ar33;
  TH1I * Vlab_Ar34;
  TH1I * Vlab_Ar35;
  TH1I * Vlab_Cl31;
  TH1I * Vlab_Cl32;
  TH1I * Vlab_Cl33;
  TH1I * Vlab_Cl34;


};
#endif
