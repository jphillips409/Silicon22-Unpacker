#ifndef s800_
#define s800_
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "TRandom.h"
#include "histo_sort.h"
#include "histo_read.h"
#include "calibrate.h"
#include "TObject.h"
#include "TString.h"
#include "S800defs.h"
#include "histo_sort.h"
#include "histo_read.h"
#include "pid.h"

using namespace std;

class S800;


struct tof
{
  Short_t rf;
  Short_t obj;
  Short_t xfp;
  Short_t si;
  Short_t tac_obj;
  Short_t tac_xfp;
};

struct trig
{
  int registr;
  int s800;
  int external1;
  int external2;
  int secondary;
};

struct scint
{
  int de_up;
  int time_up;
  int de_down;
  int time_down;
};


//we will not be using the hodoscope, but its here for completeness
struct hodo
{
  int id;
  int energy;
};

struct crdc
{
  Short_t anode;
  Short_t tac;
  Short_t samplewidth =0;
  int PadMult=0;
  Short_t cdata[224][32]; // [NPads][NSamples] - could be much less than 32 but leaving it 
  float raw[224]; //can be at max all of the pads
  float cal[224]; //can be at max all of the pads
  float x_gravity;
  float padsum;
  float mom;
  float calY;
  Short_t pad[224]; //can be at max all of the pads
  void Reset()
  {
    for(int i=0;i<224;i++)
    {
      for(int j=0;j<32;j++)
      {
        cdata[i][j] = 0;
      }
    }

    for(int i=0;i<min(224,PadMult);i++)
    {
      raw[i] =0;
      cal[i] =0;
      pad[i] =0;
    }

    PadMult=0;
    calY =0;
    x_gravity=-1;
    padsum = 0;
    samplewidth=0;
    anode =0;
    tac=0; 
  }
};

struct ic
{
  float sum;
  Short_t raw[16];
  float Scaledsum;
  float ScalingFactor;
  void Reset()
  {
    for(int i =0;i<16;i++)
      raw[i] = 0;
    sum =0;
    Scaledsum = 0;
  }
};



//I'm going to use vectors here since this is multihit
//We may want to revist this at some point
//KB 
struct s8mtdc
{
  vector <float> e1up;
  vector <float> e1down;
  vector <float> xfp;
  vector <float> obj;
  vector <float> galotte;
  vector <float> rf;
  vector <float> hodoscope;
  vector <float> objCorrected;
  vector <float> xfpCorrected;
  void Reset()
  {
    e1up.clear();
    e1down.clear();
    xfp.clear();
    obj.clear();
    galotte.clear();
    rf.clear();
    hodoscope.clear();
    objCorrected.clear();
    xfpCorrected.clear();
  }
};


/*! \class S800Map
    \brief Class for implementation of an S800 inverse map.

    Class for use of an inverse map based on information from the 
    S800 fp detectors.  Access is by s800->fp.track.map.<..>.
*/

class S800Map : public TObject {
 private:
  
 public:
  UShort_t maxcoefficient[S800_TRACK_PARAMETERS];
  UShort_t order[S800_TRACK_PARAMETERS][S800_TRACK_COEFFICIENTS];
  UShort_t exponent[S800_TRACK_PARAMETERS][S800_TRACK_PARAMETERS][S800_TRACK_COEFFICIENTS];
  Double_t coefficient[S800_TRACK_PARAMETERS][S800_TRACK_COEFFICIENTS];
  
 public:
  /* Variables */
  double maxcoefficients;
  double maxparameters;
  double maxorder;
  
  TString mapFilename;

 public:
  S800Map();
  ~S800Map();
  void LoadInverseMap(TString filename);
  void Reset();
  Bool_t WasLoaded() 
  { if (maxorder == 0) return false;
    return true;
  }
  double Calculate(int calcorder, int parameter, double *input);
  
  //  ClassDef(S800Map, 1);
};

/**************************************************************/
/* Tracking s800->fp.track.<..>                               */
/**************************************************************/

class S800FpTrack : public TObject {

 public:
  S800Map ***map;
  int S800Setting;

  int runno;

  /* Parameters */
  double xfp; /*!< Focal plane dispersive direction position, in m */
  double afp; /*!< FP dispersive angle, in radians */
  double yfp; /*!< FP non-dispersive direction, in m */
  double bfp; /*!< FP non-dispersive angle, in radians */
  double ata; /*!< Dispersive angle at the target, in radians */
  double yta; /*!< Non-dispersive position at the target, in m */
  double bta; /*!< Non-dispersive angle at the targ  pid *** residue_pid;et, in radians */
  double dta; /*!< Fractional energy at the target, (E-E0)/E0, in parts */
  double azita; /*!< Phi polar angle at target, in radians */
  double scatter; /*!< Theta polar angle at target, in mradians */
  double energy; /*!< Energy at target, in MeV */
  double ptot; /*!< Total momentum, in MeV/c */
  double ppar; /*!< Parallel momentum, in MeV/c */
  double ptra; /*!< Transverse momentum, in MeV/c */
  double theta; //radians
  double phi; //radians
  double thetadeg; //degrees
  double phideg; //degrees
  
  /* Corrected angles for yta/dta correlation */
  double ata_cor;
  double bta_cor;
  double azita_cor;
  double scatter_cor;
  
  /* Variables */
  double anglea; /*!< Dispersive angle offset, in radians */
  double angleb; /*!< Non-dispersive angle offset, in radians */
  double brho; /*!< Brho of the S800, i.e. rigidity of a particle 
		 following a central trajectory */
  double mass;
  double charge;
  double order;
  double zfp; /*!< Displacement of FP position from CRDC1, 
		in the beam direction, in m */
  double gecorr;
  double gap;
  
  double anglea_cor;
  double angleb_cor;
  double ata_dtacor;
  double bta_ytacor;

  double beta0, deltabeta;

  crdc CRDC[2];

  S800FpTrack(S800Map ***map);
  ~S800FpTrack();
  void Reset();
  void CalculateTracking(int A, int Z);

  void GetThetaPhi(double ata0, double bta0);


};


class S800
{

 public:
  S800(TRandom * ran0, histo_sort * Histo1, int setting);
  void Reset();
  ~S800();
  bool unpack(unsigned short*& point,int runno);

  bool analyze();
  void CRDCIntegrate(int id);
  void CalculateGravity(int id);
  void reset();

  //Data structures
  tof ToF;
  trig Trig;
  scint Scint[3];
  hodo Hodo[32];
  crdc CRDC[2];
  ic IC;
  s8mtdc mTDC;
  pid * beam_pid;
  pid *** residue_pid;

  int nSettings;
  int nBeams;
  //  pid* residue_pid_Ar30;
  //pid* residue_pid_Cl30;

  void Init(int setting);

  S800Map ***map;
  S800FpTrack *track;

  int S800Setting; //0 for Cl and 1 for Al
  int BeamID;
  void GetBeamId(int);
   
  //Add timestamp
  long double tstamp;

  
  float Pressure[200];
  float Scaling[200];

  //Flags for residue PID, don't need to send fragments to PID subroutine if I don't have masses/loss files yet
  bool Ne17_flag;
  bool Ne18_flag;
  

  static const int NPads = 480;
  int CRDCPeds[2][NPads];
  float Chargeslope[2][2][NPads];
  float Chargeintercept[2][2][NPads];
  float CRDCSlopeX[2]; //[CRDC]
  float CRDCSlopeY[2]; //[CRDC]
  float CRDCOffsetX[2];//[CRDC #][run54/56 = 0, run109/110 = 1]
  float CRDCOffsetY[2][6];
  float ObjCorr[2]; //1 - obj/mrad  2 - obj/mm
  int gravity_width;

  float objcorr3;
  float xfpcorr3;

  //  S800FpTrack track;


  int SRStot = 0;
  int SRSgood = 0;
  int SRSbad = 0;

  int S800Huge = 0;

 private:
  TRandom *ran;
  histo_sort * Histo;
  histo_read * Histo_read;

  
 
 protected:
  
  unsigned short* DecodeS800Crdc(unsigned short *pevent, int id, int runno);
  unsigned short* DecodeS800CrdcRaw(unsigned short *pevent, int id);
  unsigned short* DecodeS800Scintillator(unsigned short *pevent, unsigned short updown, int id);
  unsigned short* DecodeS800IonChamber(unsigned short *pevent);
  unsigned short* DecodeS800TimeOfFlight(unsigned short *pevent);
  unsigned short* DecodeS800Trigger(unsigned short *pevent);
  unsigned short* DecodeS800HodoScope(unsigned short *pevent);
  unsigned short* DecodeS800MPDC(unsigned short *pevent, int runno);
  unsigned short* DecodeS800NewMultiHitTDC(unsigned short *pevent);

};



#endif
