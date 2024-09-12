#ifndef histo_read_
#define histo_read_
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "TH1F.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TH3I.h"
#include "TFile.h"
#include "TDirectory.h"

using namespace std;

class histo_read
{
protected:

  TFile * file_read; //!< output root file

  //correlations
  TDirectoryFile * dirCorr; //!< directory for correlations

  //gammas
  TDirectoryFile * dirGamma;

  //CRDC Gates
  TDirectoryFile * dirCRDC_proj;



  TDirectoryFile *dirdEE; //!< directory for the dEE maps

  TDirectoryFile *dirProjections;
  TDirectoryFile *dirProjectionsRaw;
  

  TDirectoryFile *dirCsI;
  TDirectory *dirCsIGate;

  TDirectoryFile *dirSi;
  TDirectory *dirELoss;

  //Front & Back
  TDirectoryFile * dirFB; //!< directory for front/back mixed spectra


  
public:
    histo_read(string suffix);                  //!< constructor
    ~histo_read(){};
    void write(); //!< write the root spectra to file


    int Ntele;
    int Nstrip;
    int Nceasar;
    int NCsI;
    int Nring;

    TH2I *CorrelationTable;

    TH2I ** dEE;
    TH2I * HitMap;

    TH1I *** dEE_Projections;
    TH1I *** dEE_ProjectionsRaw;

    //CsI spectra
    TH1I ** Light;

    //ELoss spectra
    TH1I ** ELoss;
      
 
    TH2I ** FvB_Tele;
    

    //FB
    TH2I ** FBMult;
    TH1I ** FBDiff;
    TH1I ** FBDiffLG;
    TH2I ** FB;


    //34Ca
    TH1I * Erel_34Ca_2p32Ar;
    TH1I * Ex_34Ca_2p32Ar;
    TH1I * ThetaCM_34Ca_2p32Ar;
    TH1I * VCM_34Ca_2p32Ar;

    //35Ca
    TH1I * Erel_35Ca_2p33Ar;
    TH1I * Ex_35Ca_2p33Ar;
    TH1I * ThetaCM_35Ca_2p33Ar;
    TH1I * VCM_35Ca_2p33Ar;
    TH1I * Erel_fake34K_p33Ar;
    TH1I * Ex_fake34K_p33Ar;
    TH1I * ExY1_35Ca;
    TH1I * ExY2_35Ca;
    TH1I * ExT_35Ca;
    TH2I * Jacobi_T_35Ca;
    TH2I * Jacobi_Y_35Ca;
    TH1I * Erel_35Ca_EvntMix;
    TH1D * Erel_35Ca_EvntMixWeighted;

    //36Ca
    TH1I * Erel_36Ca_p35K;
    TH1I * Ex_36Ca_p35K;
    TH1I * Erel_36Ca_p35K_lowres;
    TH1I * Ex_36Ca_p35K_lowres;
    TH1I * ThetaCM_36Ca_p35K;
    TH1I * VCM_36Ca_p35K;
    TH1I * Vlab_HF_p35K;
    TH1I * Vlab_LF_p35K;
    TH1I * Vlab_HF_p35K_before;
    TH1I * Vlab_LF_p35K_before;
    TH2I * Vlab_LF_vs_Ex;
    TH1I * Timing_CsI_36Ca_p35K;
    TH2I * Timing_vs_Erel_36Ca_p35K;
    TH2I * VCM_vs_ThetaCM;
    TH1I * Ca36_ppar;
    TH1I * Ca36_deltappar;
    TH1I * Ca36_ptra;

    TH1I * Rigidity_protongated2plus;
    TH1I * Rigidity_gammagated2plus;
    TH1I * Rigidity_gammagated2plus_strict;
    TH2I * cosbeamCMtoHF_Ex_p35K;
    TH2I * cosbeamCMtoHF_Erel_p35K;
    TH1I * Energy_protongated2plus;
    TH1I * EnergyKin_protongated2plus;
    TH2I * Ex_ineachCsI;


    TH1I * Erel_36Ca_2p34Ar;
    TH1I * Ex_36Ca_2p34Ar;
    TH1I * ThetaCM_36Ca_2p34Ar;
    TH1I * VCM_36Ca_2p34Ar;

    //37Ca
    TH1I * Erel_37Ca_p36K;
    TH1I * Ex_37Ca_p36K;
    TH1I * Erel_37Ca_p36K_lowres;
    TH1I * Ex_37Ca_p36K_lowres;
    TH1I * ThetaCM_37Ca_p36K;
    TH1I * VCM_37Ca_p36K;
    TH1I * Vlab_HF_p36K;
    TH1I * Vlab_LF_p36K;
    TH1I * Vlab_HF_p36K_before;
    TH1I * Vlab_LF_p36K_before;
    TH1I * proton_KE_37Ca;
    TH2I * cosbeamCMtoHF_Ex_p36K;
    TH1I * Timing_CsI_37Ca_p36K;
    TH2I * Timing_vs_Erel_37Ca_p36K;

    TH1I * costheta_protongated2plus;
    TH1I * Erel_37Ca_p36Kd;
    TH1I * Ex_37Ca_p36Kd;
    TH1I * VCM_37Ca_p36Kd;
    TH1I * Erel_37Ca_p36Kt;
    TH1I * Ex_37Ca_p36Kt;
    TH1I * VCM_37Ca_p36Kt;
    TH1I * Erel_37Ca_EvntMix;
    TH1I * Ex_37Ca_EvntMix;
    TH2I * cosbeamCMtoHF_37Ca_EvntMix;

    TH1I * Erel_37Ca_a33Ar;
    TH1I * Ex_37Ca_a33Ar;
    TH1I * Erel_37Ca_2p35Ar;
    TH1I * VCM_37Ca_2p35Ar;

    //37Sc
    TH1I * Erel_37Sc_p36Ca;
    TH1I * Ex_37Sc_p36Ca;
    TH1I * ThetaCM_37Sc_p36Ca;
    TH1I * VCM_37Sc_p36Ca;

    TH1I * Erel_37Sc_2p35K;
    TH1I * Ex_37Sc_2p35K;
    TH1I * ThetaCM_37Sc_2p35K;
    TH1I * VCM_37Sc_2p35K;
    TH2I * cosbeamCMtoHF_Erel_p36Ca;

    TH1I * Ca36_rigidity;
    TH1I * Sc37_ppar;
    TH1I * Sc37_ptra;

    //38Sc
    TH1I * Erel_38Sc_p37Ca;
    TH1I * Ex_38Sc_p37Ca;
    TH1I * Gamma37Ca_allfrom38Sc;
    TH1I * Gamma37Ca_1strom38Sc;
    TH1I * ThetaCM_38Sc_p37Ca;
    TH1I * VCM_38Sc_p37Ca;
    TH2I * cosbeamCMtoHF_Erel_p37Ca;
    TH1I * Erel_38Sc_p37Ca_vgate;
    TH1I * Timing_CsI_38Sc_p37Ca;
    TH2I * Timing_vs_Erel_38Sc_p37Ca;

    TH1I * Ca37_ppar;
    TH1I * Ca37_ptra;
    TH1I * Ca37_rigidity;
    TH1I * Sc38_ppar;
    TH1I * Sc38_deltappar;
    TH1I * Sc38_ptra;
    TH1I * Sc38_ppar_gs;
    TH1I * Sc38_ptra_gs;
    TH1I * Ca37_rigidity_gs;
    TH1I * Vlab_LF_p37Ca;
    TH1I * Vlab_HF_p37Ca;
    TH2I * Vlab_LF_vs_Erel_p37Ca;

    TH1I * Erel_38Sc_2p36K;
    TH1I * VCM_38Sc_2p36K;
    TH1I * Erel_38Sc_sub37Ca;
    TH1I * Ex_38Sc_sub37Ca;
    TH2I * Jacobi_T_38Sc;
    TH2I * Jacobi_Y_38Sc;
    TH1I * Timing_CsI_38Sc_2p36K;

    //39Sc
    TH1I * Erel_39Sc_d37Ca;
    TH1I * Ex_39Sc_d37Ca;

    //39Ti
    TH1I * Erel_39Ti_2p37Ca;
    TH1I * VCM_39Ti_2p37Ca;

    //38Ti
    TH1I * Erel_38Ti_2p36Ca;
    TH1I * Ex_38Ti_2p36Ca;
    TH1I * ThetaCM_38Ti_2p36Ca;
    TH1I * VCM_38Ti_2p36Ca;

    //K34
    TH1I * Erel_33K_p32Ar;
    TH1I * Erel_34K_p33Ar;
    TH1I * Ex_34K_p33Ar;
    TH1I * Gamma33Ar_allfrom34K;
    TH1I * Gamma33Ar_1stfrom34K;
    TH1I * Gamma33Ar_2ndfrom34K;
    TH1I * Gamma33Ar_upper34K;
    TH2I * Gamma33ArvsErel;
    TH2I * cosbeamCMtoHF_Erel_p33Ar;
    TH1I * Erel_34K_p33Ar_37CaBeam;
    TH1I * Erel_34K_p33Ar_36KBeam;
    TH1I * ThetaCM_34K_p33Ar;
    TH1I * VCM_34K_p33Ar;
    TH1I * VCM_34K_1stfrom34K;
    TH1I * VCM_34K_2ndfrom34K;
    TH1I * VCM_34K_upper34K;
    TH2I * Vlab_34K_LF_vs_Erel;
    TH1I * proton_KE_34K;
    TH1I * Erel_34K_EvntMix;
    TH1D * Erel_34K_EvntMixWeighted;
    TH2I * cosbeamCMtoHF_EvntMix;
    TH1I * Erel_34K_EvntMix_gated;


    //K35
    TH1I * Erel_35K_p34Ar;
    TH1I * Ex_35K_p34Ar;
    TH2I * cosbeamCMtoHF_Erel_p34Ar;
    TH1I * VCM_35K_p34Ar;
    TH1I * ThetaCM_35K_p34Ar;
    TH1I * K35_ppar_1st;
    TH1I * K35_ptra_1st;
    TH1I * K35_rigidity_1st;
    TH1I * Gamma34Ar_allfrom35K;
    TH1I * Gamma34Ar_1stfrom35K;
    TH1I * Erel_35K_p34Ar_otherbeam;
    TH1I * Ex_35K_p34Ar_otherbeam;
    TH1I * Erel_35K_EvntMix;
    TH1I * Ex_35K_EvntMix;
    TH1I * Erel_35K_2p33Cl;
    TH1I * VCMl_35K_2p33Cl;

    //K36
    TH1I * Erel_36K_p35Ar;
    TH1I * Ex_36K_p35Ar;

    //Ar33,34,35
    TH1I * Erel_33Ar_p32Cl;
    TH1I * Ex_33Ar_p32Cl;
    TH1I * Erel_34Ar_p33Cl;
    TH1I * Ex_34Ar_p33Cl;
    TH1I * Erel_35Ar_p34Cl;
    TH1I * Ex_35Ar_p34Cl;

    //Cl32,33,34
    TH1I * Erel_32Cl_p31S;
    TH1I * Ex_32Cl_p31S;
    TH1I * Erel_33Cl_p32S;
    TH1I * Ex_33Cl_p32S;
    TH1I * Erel_34Cl_p33S;
    TH1I * Ex_34Cl_p33S;

    //S31,32,33
    TH1I * Erel_31S_p30P;
    TH1I * Ex_31S_p30P;
    TH1I * Erel_32S_p31P;
    TH1I * Ex_32S_p31P;
    TH1I * Erel_33S_p32P;
    TH1I * Ex_33S_p32P;

    //P29
    TH1I * Erel_29P_p28Si;
    TH1I * Ex_29P_p28Si;

    //LiBe
    TH1I * Erel_6Li_ad;
    TH1I * Ex_6Li_ad;
    TH1I * Erel_7Li_at;
    TH1I * Ex_7Li_at;
    TH1I * Erel_8Be_aa;
    TH1I * Ex_8Be_aa;


    //Gammas
    TH1I * TEC_Ca38;
    TH2I * dopp_timing_Ca38;
    TH1I * TEC_Ca37;
    TH2I * dopp_timing_Ca37;
    TH1I * TEnoC_Ca37;
    TH1I * TEC_Ca36;
    TH1I * TEC_Ca36_simpletheta;
    TH1I * TEC_Ca36_DataEC;
    TH1I * TEC_Ca36_NoTadded;
    TH1I * TEC_Ca36_mult0;
    TH1I * TEC_Ca36_mult1;
    TH1I * TEC_Ca36_mult2;
    TH1I * TEC_Ca36_mult3;
    TH1I * TEC_Ca36_early;
    TH1I * TEC_Ca36_late;
    TH2I * dopp_timing;
    TH1I * TEC_Ca36_timecutearly;
    TH1I * TEC_Ca36_timecutlate;
    TH1I * TEC_Ca36_forward;
    TH1I * TEC_Ca36_slightforward;
    TH1I * TEC_Ca36_slightbackward;
    TH1I * TEC_Ca36_backward;

    TH1I * TEC_Ca36_noaddback_nodoppler;
    TH1I * TEC_Ca36_noaddback;
    TH1I * TEC_K36_noaddback;
    TH1I * TEC_K35;
    TH1I * TEC_K35_DataEC;
    TH1I * TEC_K35_DataEC_nodopp;
    TH1I * TEC_K35_NoTadded;

    TH1I * TEC_K36;
    TH1I * TEC_K36_forward;
    TH1I * TEC_K36_slightforward;
    TH1I * TEC_K36_slightbackward;
    TH1I * TEC_K36_backward;
    TH1I * TEC_Ar32;
    TH1I * TEC_Ar33;
    TH1I * TEC_Ar34;
    TH1I * TEC_Ar35;


    TH1I * TEC_K36_K36beam;
    TH1I * TEnoC_K36_K36beam;
};

#endif
