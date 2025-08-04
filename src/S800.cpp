#include "S800.h"

S800::S800(TRandom * ran0, histo_sort * Histo1, int setting)
{
  Histo = Histo1;
  ran = ran0;
  Init(setting);
}

//***********************************************************************************
//During run 173, MPDC 2 (CRDC 2) crapped out. For fragments in MPDC 2, no longer have momentum resolution
// i.e. cannot use inverse maps or get a pid correction for those fragments that hit MPDC 2.
// Need to pass a run number to S800 so I can account for funky detector stuff like this.
// A possible solution would be to use previous runs momentum resolution as a standard for crapped MPDC 2 runs.
// Might not need an angular correction for the pid to work on those fragments, Si22 and Si23 might be well separated enough.
// Still get energy from those fragments, just no good position.
//***********************************************************************************





void S800::Reset()
{
  mTDC.Reset();
  IC.Reset();
  CRDC[0].Reset();
  CRDC[1].Reset();
  beam_pid->Z = 0;
  beam_pid->A = 0;

  for(int i =0;i<nSettings;i++)
  {
    for(int j=0;j<nBeams;j++)
    {
      residue_pid[i][j]->Z = 0;
      residue_pid[i][j]->A = 0;
    }
  } 

  track->Reset();

  BeamID = -99;
}

S800::~S800()
{
  delete beam_pid;
  for (int i=0; i<2; i++) {delete map[0][i];}
  for (int i=0; i<22; i++) {delete map[1][i];}
  for (int i=0; i<2; i++) {delete map[2][i];}
  delete [] map;
  delete track;
  /*for(int i =0;i<nSettings;i++) //Loop over setttings
  {
    for(int j=0;j<nBeams;j++) //Loop over beams
    {
      delete residue_pid[i][j];
    }
  }*/
  delete [] residue_pid;
  // delete residue_pid_O14;
  // delete residue_pid_Ne17;
}

void S800::Init(int setting)
{
  S800Setting = setting;

  //The inverse map is different for each S800 setting and will
  //need to be chaged for each charge and rigidity
  //manually load in all the inverse maps
  map = new S800Map**[3]; //initialize the inverse map for number of settings TODO add more settings
  map[0] = new S800Map*[2]; //initialize space for 1 Si23 setting inverse maps TODO add more maps
  for (int i=0; i<2; i++) {map[0][i] = new S800Map();}
  //Setting 0 is used for uncreated beam runs, switch the inverse map for the desired run
  map[0][0]->LoadInverseMap("s800inputs/UnreactedNoTarg_Si23.inv"); //load Si23 for just unreacted
  map[0][1]->LoadInverseMap("s800inputs/UnreactedNoTarg_Mg21.inv"); //load Mg21 for just unreacted
 // map[0][0]->LoadInverseMap("s800inputs/UnreactedFibers_Si23.inv"); //load Si23 for unreacted + fibers
  //map[0][1]->LoadInverseMap("s800inputs/UnreactedFibers_Mg21.inv"); //load Mg21 for unreacted + fibers
  //map[0][0]->LoadInverseMap("s800inputs/UnreactedFibersThinTarg_Si23.inv"); //load Si23 for unreacted + fibers + thin target
  //map[0][1]->LoadInverseMap("s800inputs/UnreactedFibersThinTarg_Mg21.inv"); //load Mg21 for unreacted + fibers + thin target


  map[1] = new S800Map*[22]; //initialize space for 5 Si23 setting inverse maps TODO add more maps
  for (int i=0; i<22; i++) {map[1][i] = new S800Map();}
  map[1][0]->LoadInverseMap("s800inputs/Setting1_Si23.inv"); //load Si23 setting
  map[1][1]->LoadInverseMap("s800inputs/Setting1_Al22.inv"); //load Al22 setting
  map[1][2]->LoadInverseMap("s800inputs/Setting1_Mg21.inv"); //load Mg20 setting
  map[1][3]->LoadInverseMap("s800inputs/Setting1_Mg20.inv"); //load Mg20 setting
  map[1][4]->LoadInverseMap("s800inputs/Setting1_Ne17.inv"); //load Ne17 setting
  map[1][5]->LoadInverseMap("s800inputs/Setting1_Ne18.inv"); //load Ne18 setting
  map[1][6]->LoadInverseMap("s800inputs/Setting1_Al23.inv"); //load Al23 setting
  map[1][7]->LoadInverseMap("s800inputs/Setting1_Na20.inv"); //load Na20 setting
  map[1][8]->LoadInverseMap("s800inputs/Setting1_Na21.inv"); //load Na21 setting
  map[1][9]->LoadInverseMap("s800inputs/Setting1_Mg22.inv"); //load Mg22 setting
  map[1][10]->LoadInverseMap("s800inputs/Setting1_Ne19.inv"); //load Ne19 setting
  map[1][11]->LoadInverseMap("s800inputs/Setting1_F18.inv"); //load F18 setting
  map[1][12]->LoadInverseMap("s800inputs/Setting1_F17_0_6604.inv"); //load F17 setting
  map[1][13]->LoadInverseMap("s800inputs/Setting1_O15.inv"); //load O15 setting
  map[1][14]->LoadInverseMap("s800inputs/Setting1_O16.inv"); //load O16 setting
  map[1][15]->LoadInverseMap("s800inputs/Setting1_Ne20.inv"); //load Ne20 setting
  map[1][16]->LoadInverseMap("s800inputs/Setting1_O14.inv"); //load O14 setting
  map[1][17]->LoadInverseMap("s800inputs/Setting1_N14.inv"); //load N14 setting
  map[1][18]->LoadInverseMap("s800inputs/Setting1_N13.inv"); //load N13 setting
  map[1][19]->LoadInverseMap("s800inputs/Setting1_C12.inv"); //load C12 setting
  map[1][20]->LoadInverseMap("s800inputs/Setting1_C11.inv"); //load C11 setting
  map[1][21]->LoadInverseMap("s800inputs/Setting1_Si24.inv"); //load Si24 setting

  //TODO get real inverse maps for this setting
  map[2] = new S800Map*[2]; //initialize space for 5 Si23 setting inverse maps TODO add more maps
  for (int i=0; i<2; i++) {map[2][i] = new S800Map();}
  map[2][0]->LoadInverseMap("s800inputs/Setting2_Si22.inv"); //load Si22 setting
  map[2][1]->LoadInverseMap("s800inputs/Setting2_Si23.inv"); //load Si22 setting


  //////////////////////////////////////
  track = new S800FpTrack(map);

  float brho;
  string beam_pid_name;

  //settings for the s800 settings: inverse map, brho
  if (S800Setting == 0)
  {
    //Use setting 0 for different calibration runs
    //Uncomment the correct one. Beam PID is the same as setting 1
    brho = 2.41493; //run 48, just unreacted beam
    //beam_pid_name = "Setting1_beam";

    //brho = 2.31911; //run 45 + 46, unreacted beam + fibers
    //beam_pid_name = "Setting1_beam";

    //brho = 2.24917; //run 49 + 50, unreacted beam + fibers + 0.5 mm targ
    beam_pid_name = "Setting1_beam";

    //Normally uncomment this, don't want to accidently use this if not calibrating
    //cout << "There is no setting 0, choose 1-3" << endl;
    //abort();
  }
  else if (S800Setting == 1)
  {
    brho = 2.35714;
    beam_pid_name = "Setting1_beam";
  }
  else if (S800Setting == 2)
  {
    brho = 1.9900;
    beam_pid_name = "Setting2_beam";
  }
  else
  {
    cout << "Unknown setting #" << setting << endl;
    abort();
  } 

  cout << "Loading beam_pid zline/" << beam_pid_name << ".zline" << endl;
  beam_pid = new pid(beam_pid_name, true); //directory for S800 zlines

  cout << "s800 rigidity brho = " << brho << endl;
  track->brho = brho;
  track->S800Setting = S800Setting;

  /////////////////////////////////////////////////////////////////////////////////////////
  //Load in the residue PID gates
  // beam numbers
  // 0 - Si23
  // 1 - Al22
  // 2 - Mg21
  // setting numbers
  // 0 - Used for unreacted beam calibrations
  // 1 - Furthest position to get 2+ state, brho = 2.3571
  // 2 - Closest position to get p + 22Si from 23P, brho = 1.9900

  string name;
  nSettings = 3;
  nBeams = 4;

  //string beams[4] = {"Ca37","K36","Ar35","Cl34"};
  string beams[nBeams] = {"Si23","Al22","Mg21","Na20"};
  string settings[nSettings] = {"Setting0","Setting1","Setting2"};
  residue_pid = new pid**[nSettings]; //Two arrays of residue pids for Ca and K settings

  for(int i =0;i<nSettings;i++) //Loop over setttings
  {
    residue_pid[i] = new pid*[nBeams];

    for(int j=0;j<nBeams;j++) //Loop over beams
    {
      if(i==0)
        //use for different calibrations, uncomment the needed one
        name = Form("UnreactedBeam/Unreacted_residue",beams[j].c_str());
        //name = Form("UnreactedBeam/Fiber_residue",beams[j].c_str());
        //name = Form("UnreactedBeam/FiberThinTarg_residue",beams[j].c_str());
      else if(i==1)
        name = Form("Setting1/s800_residue%s",beams[j].c_str());
      else if(i==2)
        name = Form("Setting2/s800_residue%s",beams[j].c_str());

      cout << "for s800setting " << i << " load " << name << endl;
      residue_pid[i][j] = new pid(name, true); //use S800 directory
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////
  
  ifstream infile("datfiles/CRDCPedestals.dat");
  if(!infile.is_open())
  {
    cout << "Couldn't open CRDC Pedestal file" << endl;
    abort();
  }

  int CRDC,set,pad;
  float ped;
  getline(infile,name);
  for(;;)
  {
    infile >> CRDC >> pad >> ped;

    CRDCPeds[CRDC-1][pad] = (int)ped; //currently all set to 0
    if(infile.eof()) break;
    if(infile.bad()) break;
  }
  infile.close();
  infile.clear();

  ifstream infile2("datfiles/CRDCgain.dat");
  if(!infile2.is_open())
  {
    cout << "Couldn't open CRDCpars file" << endl;
    abort();
  }
  getline(infile2,name);

  
  float slope, intercept;
  for(;;)
  {
    infile2 >> CRDC >> set >> pad >> intercept >> slope;
    Chargeslope[CRDC][set][pad] = slope;
    Chargeintercept[CRDC][set][pad] = intercept;
    //cout << CRDC << "   " << pad << "   " << slope << "   " << intercept << endl;
    if(infile2.eof()) break;
    if(infile2.bad()) break;
  }
  infile2.close();
  infile2.clear();

  ObjCorr[0] = 100; //obj/mrad
  ObjCorr[1] = 0.011; //obj/mm

  //CRDC calibration parameters
  //[CRDC #]
  // e19017
  
  CRDCSlopeX[0] = -1.23083;// mm/pad
  CRDCSlopeX[1] = -1.23083;// mm/pad
  CRDCSlopeY[0] = -0.03305;//-0.05858;// mm/ch
  CRDCSlopeY[1] = 0.033325;//0.05858;// mm/ch

  //[CRDC #]
  CRDCOffsetX[0] = 299.971;//mm
  CRDCOffsetX[1] = 299.971;//mm
  //[CRDC #][run54/56 = 0, run109/110 = 1]
  CRDCOffsetY[0][0] = 1034.55;//398.46;//5600; //mm
  CRDCOffsetY[1][0] = -1039.92;//-407.59;//-5616; //mm
  CRDCOffsetY[0][1] = -29.52 - 1.31 - 9.0;//-5.36;//50.614; //mm //overestimate to get runs right, mask will be off
  CRDCOffsetY[1][1] = 28.1 - 9.0;//3.11;//-61.917; //mm
  CRDCOffsetY[0][2] = 3.85 - 30.0; //mm
  CRDCOffsetY[1][2] = -4.83 - 20.0; //mm
  CRDCOffsetY[0][3] = 7.0 + 16.0;//-48.85; //mm    every run past 162, constant shift assumption
  CRDCOffsetY[1][3] = 6.0 + 4.0; //mm    every run past 162, constant shift assumption
  CRDCOffsetY[0][4] = 29 + 10.0;//-48.85; //mm    every run past 200, constant shift assumption
  CRDCOffsetY[1][4] = 7.0 + 9.0; //mm    every run past 200, constant shift assumption
  CRDCOffsetY[0][5] = 35.2;//-48.85; //mm    every run past 220, constant shift assumption
  CRDCOffsetY[1][5] = 10.0; //mm    every run past 220, constant shift assumption

  gravity_width = 12; //pads

}


bool S800::unpack(unsigned short *&point,int runno)
{
  unsigned long long int n;
  unsigned short plength, ptag, ID, nwords, words;

  nwords = *point++;
  Histo->S800_NumWords_R->Fill(nwords);
  if (nwords > 3000)
  {
    //cout << "huge S800 event nwords = " << nwords << endl;
    S800Huge++;
    return false; //S800 events should not be this massive, current limit set at 3000
  }
  Histo->S800_NumWords_Accepted->Fill(nwords);
  nwords--;
  if(*point++ != S800_PACKET)
    {
      cout << "Didn't find the S800 marker" << endl;
      return false;
    }
  else
    {
      nwords--;
    }
  

  // Unpack S800 data in tree
  // cout << *point << " " << S800_VERSION << endl;
  if(*point++ != S800_VERSION){
    std::cout
      << "DecodedEvent: " 
      << "Wrong version of S800 sub-event. Aborting ..." 
      << std::endl;
    return 1;
  }

  // if(runno <= 158 && runno >=81)
  //   IC.ScalingFactor = Scaling[runno];
  // else
  //   IC.ScalingFactor = 1.;


  // Histo->ScalingFactorvsRun->Fill(runno,IC.ScalingFactor);
  
  nwords--;
  //cout << "nwords " << nwords << endl;
  int counter = 0;

  int nwords_check = nwords;

  while(nwords > 0)
  {
    //cout << "HERE " << nwords << endl;
    plength = *point; ++point; 
    ptag = *point; ++point;

    if (plength > nwords || nwords > nwords_check || plength == 0)
    {
      cout << "Likely a problem with the MPDC unpacker" << endl;
      //Likely a problem with the MPDC unpacker
      SRSbad++;
      return false;
    }
    //cout << plength << " " << hex << ptag << dec << endl;
    //return false;
    switch(ptag){
    case S800_TRIGGER_PACKET:
      //cout << "Trigger: " << hex << S800_TRIGGER_PACKET << dec << endl;
      point = DecodeS800Trigger(point);
      break;

    case S800_TOF_PACKET:
      //cout << "Tofr: " << hex << S800_TOF_PACKET << dec << endl;
      point = DecodeS800TimeOfFlight(point);
      break;

    case S800_FP_SCINT_PACKET:
      //cout << "FP SCINT: " << hex << S800_FP_SCINT_PACKET << dec << endl;
      words = plength-2;
      while (words > 0) 
      {
        ID = ((*point)&0xF000)>>12;
        point = DecodeS800Scintillator(point, ID, ID/2);
        words -= 2;
      }
      break;
      
    case S800_FP_CRDC_PACKET:
      //cout << "CRDC: " << hex << S800_FP_CRDC_PACKET << dec << endl;
      ID = *point; ++point;
      //std::cout << "ID " << ID << std::endl;
      point = DecodeS800Crdc(point, ID, runno);
      break;

    case S800_FP_IC_PACKET:
      //cout << "FP_IC: " << hex << S800_FP_IC_PACKET << dec << endl;
      point = DecodeS800IonChamber(point);
      break;
    
    case S800_TIMESTAMP_PACKET:
      //cout << "TIMESTAMP: " << hex << S800_TIMESTAMP_PACKET << dec << endl;
      unsigned long long Tfrag;
      Tfrag = *point++;
      n = Tfrag;
      Tfrag = *point++;
      n = ((Tfrag<<16)|n);
      Tfrag = *point++;
      n = ((Tfrag<<32)|n);
      Tfrag = *point++;
      n = ((Tfrag<<48)|n);
      tstamp = n;
      //      this->SetInternalTS(n);
      break;
      
    case S800_EVENT_NUMBER_PACKET:
      //cout << "EVENT_NUMBER: " << hex << S800_EVENT_NUMBER_PACKET << dec << endl;
      n = *point++;
      n = (*point++<<16|n);
      n = (*point++<<16|n);
      //      this->SetEvtNr(n);
      break;
 
    case S800_FP_HODO_PACKET:
      //cout << "HODO: " << hex << S800_FP_HODO_PACKET << dec << endl;
      point = DecodeS800HodoScope(point);
      break;
    case S800_DC_PACKET:
      SRStot++;
      point = DecodeS800MPDC(point,runno); 
      
      break;
    case S800_VME_TDC_PACKET:
      //cout << "VME_TDC: " << hex << S800_VME_TDC_PACKET << dec << endl;
      point = DecodeS800NewMultiHitTDC(point);
      break;
    default: // S800_II_CRDC_PACKET, S800_II_TRACK_PACKET...
      point += plength - 2;
      break;
    }

    counter++;

    //if (counter >100) abort();
    nwords -= plength;
    //cout << nwords << endl;
  }
  //this->SetTS(ts);
  //if(ffirst_ts<1){
  //  ffirst_ts = ts;
  //}
  SRSgood++;
  return true;
 
  
}
unsigned short* S800::DecodeS800TimeOfFlight(unsigned short *p){

  UShort_t words = (*(p-2))-2, ch, dum;
  Short_t rf = -1;
  Short_t obj = -1;
  Short_t xfp = -1;
  Short_t si = -1;
  Short_t tac_obj = -1;
  Short_t tac_xfp = -1;
  while (words > 0) {
    ch = ((*p)&0xf000)>>12;
    int tmp = *p; ++p;
    if     (ch == 12) rf      = (tmp)&0xfff;
    else if(ch == 13) obj     = (tmp)&0xfff;
    else if(ch == 14) xfp     = (tmp)&0xfff;
    else if(ch == 15) si      = (tmp)&0xfff;
    else if(ch ==  5) tac_obj = (tmp)&0xfff;
    else if(ch ==  4) tac_xfp = (tmp)&0xfff;
    else if(ch>0 && ch<8) dum = (tmp)&0xfff;
    words--;
  }
  //  this->GetTimeOfFlight()->Set(rf, obj, xfp, si);
  //this->GetTimeOfFlight()->SetTAC(tac_obj, tac_xfp);

  ToF.rf = rf;
  ToF.obj = obj;
  ToF.xfp = xfp;
  ToF.si = si;
  ToF.tac_obj = tac_obj;
  ToF.tac_xfp = tac_xfp;

  return p;
}

unsigned short* S800::DecodeS800Trigger(unsigned short *p){
  UShort_t words = (*(p-2))-2, ch;
  int registr = -1;
  int s800 = -1;
  int external1 = -1;
  int external2 = -1;
  int secondary = -1;
  registr = *p++;
  words--;
  //cout << "words " << words << endl;
  while(words > 0){
    ch = ((*p)&0xf000)>>12;
    //cout << "CH " << ch << endl;
    if (ch ==  8) s800 = (*p++)&0xfff;
    if (ch ==  9) external1 = (*p++)&0xfff;
    if (ch == 10) external2 = (*p++)&0xfff;
    if (ch == 11) secondary = (*p++)&0xfff;
    words--;
  }
  //cout << registr << " " << s800 << " " << external1 << " " << external2 << " " << secondary << endl;
  //  this->GetTrigger()->Set(registr, s800, external1, external2, secondary);

  Trig.registr = registr;
  Trig.s800 = s800;
  Trig.external1 = external1;
  Trig.external2 = external2;
  Trig.secondary = secondary;
  
  return p;
}

unsigned short* S800::DecodeS800Scintillator(unsigned short *p, unsigned short updown, int id){
  int de_up     = -1;
  int time_up   = -1;
  int de_down   = -1;
  int time_down = -1;

  // Even updown: UP channels.  Odd ID: DOWN channels
  if (updown%2==0) {
    de_up     = (*p++)&0xfff;
    time_up   = (*p++)&0xfff;
    Scint[id].de_up = de_up;
    Scint[id].time_up = time_up;
  }
  else {
    de_down   = (*p++)&0xfff;
    time_down = (*p++)&0xfff;
    Scint[id].de_down = de_down;
    Scint[id].time_down = time_down; 
  }
  //  this->GetScintillator(id)->SetID(id);
  //this->GetScintillator(id)->Set(de_up, time_up, de_down, time_down);

  return p;
}

//Adding MPDC unpacker for SRS data - SG Oct29 2024
unsigned short* S800::DecodeS800MPDC(unsigned short *p, int runno){

  UShort_t dcId = *(p++);
  UShort_t length = *(p++);
  short i = length - 1; //the last bit isn't useful
  UShort_t nGoodHits = 0;

  int nHits = i/14;
  double tmpT = 0.0;
  double tmpTa = 0.0;
  uint64_t preT = 0;
  for(int iHit = 0; iHit < nHits; iHit++) {
//std::cout << std::hex << *(p++) << "\t" << *(p++) <<  "\t" << *(p++) << "\t" << *(p++) << std::endl;

    //unsigned short *pig = p;
    //for(int i = 0; i < 10; i++) std::cout << std::hex << *(pig++) << std::dec << " ";
    //std::cout << std::endl;

    int fecIdBuff = *(p++);
    int vmmIdBuff = *(p++);
    int chnoBuff = *(p++);
    int chnoMappedBuff = *(p++);
    int adcBuff = *(p++);
    int tdcBuff = *(p++);

    uint64_t recoTs = 0;
    for(int i = 0; i < 4; i++) {
      recoTs |= (static_cast<uint64_t>(*(p++))) << (i * 16);
    }

    short bcid = *(p++);
    short triggerOffset = *(p++);

    uint64_t recoMarker = 0;
    for(int i = 0; i < 2; i++) {
      recoMarker |= (static_cast<uint64_t>(*(p++))) << (i * 16);
    }

    double timeFec = triggerOffset * 4096.0 * 22.5;
    double timeChip = bcid * 22.5 + (1.5 * 22.5 - tdcBuff * 60/255.0);
    double driftTime = (timeChip + timeFec - recoMarker * 22.5);

    tmpT += adcBuff*driftTime;
    tmpTa += adcBuff;
    CRDC[dcId].raw[iHit] = adcBuff;
    CRDC[dcId].cal[iHit] = adcBuff;
    CRDC[dcId].pad[iHit] = chnoMappedBuff;
  }
  CRDC[dcId].PadMult = nHits;
  CRDC[dcId].tac = tmpT/tmpTa; //These are drift times weighted by pulse height
  //CRDC[dcId].calY = CRDC[dcId].tac * CRDCSlopeY[dcId] + CRDCOffsetY[dcId][0];

  //offset drifts slightly over course of experiment, interpolate linearly between
  //
  float offset;
  if (runno < 51) {offset = CRDCOffsetY[dcId][0];}
  else if (runno >= 51 && runno <= 116)
  {
    offset = CRDCOffsetY[dcId][0] + ((float)runno -50)*(CRDCOffsetY[dcId][1])/(113-50);
  }
  else if (runno > 116 && runno <= 162)
  {
    offset = CRDCOffsetY[dcId][0] + CRDCOffsetY[dcId][1] + ((float)runno -116)*(CRDCOffsetY[dcId][2])/(161-116);
  }
  else if (runno > 162 && runno < 200) //Assume the offset is double that from Mask runs 183, 184. Go up to run 267
  {
    //Start assuming constant offsets for runs
    offset = CRDCOffsetY[dcId][0] + CRDCOffsetY[dcId][1] + CRDCOffsetY[dcId][2] + (CRDCOffsetY[dcId][3]);
  }
  else if (runno >= 200 && runno <= 220) //Assume the offset is double that from Mask runs 183, 184. Go up to run 267
  {
    //Start assuming constant offsets for runs
    offset = CRDCOffsetY[dcId][0] + CRDCOffsetY[dcId][1] + CRDCOffsetY[dcId][2] + (CRDCOffsetY[dcId][4]);
  }
  else if (runno > 220) //Assume the offset is double that from Mask runs 183, 184. Go up to run 267
  {
    offset = CRDCOffsetY[dcId][0] + CRDCOffsetY[dcId][1] + CRDCOffsetY[dcId][2] + (CRDCOffsetY[dcId][5]); 
  }
  else offset = CRDCOffsetY[dcId][0];

  CRDC[dcId].calY = CRDC[dcId].tac*CRDCSlopeY[dcId] + offset;

  //CRDC[dcId].calY = tmpT/tmpTa * CRDCSlopeY[dcId] + CRDCOffsetY[dcId][0];
  //cout << dcId << " " << tmpT/tmpTa << " " << CRDC[dcId].tac << " " << CRDC[dcId].calY << " " << CRDCSlopeY[dcId] << " " << CRDCOffsetY[dcId][0] <<  endl;
  
  //cout << "dcId " << dcId << endl;
  //cout << "CRDC[dcId].PadMult " << CRDC[dcId].PadMult << "   CRDC[dcId].tac " << CRDC[dcId].tac << "    CRDC[dcId].calY " << CRDC[dcId].calY << endl;

  return p;
}

unsigned short* S800::DecodeS800HodoScope(unsigned short *p){
  UShort_t words = (*(p-2))-2;
  UShort_t id;
  UShort_t ch;
  UShort_t energy;
  int nHodo = 0;
  while (words > 0)
  {
    id = *p;
    if (id == 0)
    {
      p++;
      words--;
      while (words > 0)
      {
	      ch = (((*p)&0xF000)>>12);
	      energy = ((*p)&0x0FFF);
      	//	this->GetHodoscope(ch)->SetEnergy((Int_t)energy);
      	Hodo[nHodo].id = ch;
      	Hodo[nHodo].energy = (Int_t)energy;
      	nHodo++;
      	p++;
      	words--;
      }
    }   
    else if (id == 1) 
    {
      p++;
      words--;
      while (words > 0)
      {
      	ch = (((*p)&0xF000)>>12) + 16;
      	energy = ((*p)&0x0FFF);
      	//	this->GetHodoscope(ch)->SetEnergy((Int_t)energy);
      	Hodo[nHodo].id = ch;
      	Hodo[nHodo].energy = (Int_t)energy;
      	nHodo++;
      	p++;
      	words--;
      }
    } 
    else if (id == 2) 
    {
      p++;
      words--;
      while (words > 0) 
      {
      	// coincidence register A (for the first  16 channels)
      	p++; words--;
      	// coincidence register B (for the second 16 channels)
      	p++; words--;
      	// TAC time
      	p++; words--;
      }
    } 
    else
    {
      p++; words--;
    }
  }string name = "s800_beam_Caset";

  //cout << "hodo here " << energy << endl;

  return p;
}

unsigned short* S800::DecodeS800Crdc(unsigned short *p, int id, int runno){
  UShort_t anode=-1;
  UShort_t tac  =-1;  
  //  this->GetCrdc(id)->SetID(id);
  
  Int_t tag;
  tag = S800_FP_CRDC_PACKET;
  if (*(p+1) == tag+1) {
    p = DecodeS800CrdcRaw(p,id);
  }
  if (*(p+1) == tag+5) {
    anode = *(p+2);
    tac = *(p+3);
    p += 4;
  }

  //    this->GetCrdc(id)->SetAnodeTAC(anode, tac);
  CRDC[id].anode = anode;
  CRDC[id].tac = tac;
  //offset drifts slightly over course of experiment, interpolate linearly between
  //
  float offset;
 /* if (runno < 55) {offset = CRDCOffsetY[id][0]};
  else if (runno >= 55 && runno <= 109)
  {
    offset = CRDCOffsetY[id][0] + ((float)runno - (109-55))*(CRDCOffsetY[id][1]-CRDCOffsetY[id][0])/(109-55);
  }
  else if (runno >= 117 && runno <= 159)
  {
    offset = CRDCOffsetY[id][0] + ((float)runno - (159-117))*(CRDCOffsetY[id][2]-CRDCOffsetY[id][0])/(159-117);
  }
  else if (runno >= 255)
  {
    offset = CRDCOffsetY[id][0] + CRDCOffsetY[id][2]*2; //Assume the offset has roughly doubled
  }
  else offset = 0;*/


/*  if (runno < 56) { offset = CRDCOffsetY[id][0]; }
  else if (runno > 109) { offset = CRDCOffsetY[id][1]; }
  else
  {
    offset = CRDCOffsetY[id][0] + ((float)runno - 56.0)*
            (CRDCOffsetY[id][1]-CRDCOffsetY[id][0])/(55.0);
  }

  CRDC[id].calY = tac*CRDCSlopeY[id] + offset;
*/
  CRDC[id].calY = tac*CRDCSlopeY[id] + CRDCOffsetY[id][0];
  return p;
}

unsigned short* S800::DecodeS800IonChamber(unsigned short *p){
  UShort_t ch=-1;
  UShort_t raw=-1;
  if (*(p+1) == S800_FP_IC_ENERGY_PACKET) 
  {
    // IC packet with times
    UShort_t length;
    length = *p++;
    p++;
    length -= 2;
    while (length > 0) 
    {
      ch  = ((*p)&0xf000)>>12;
      raw = (*p++)&0xfff;
      length--;

      //      this->GetIonChamber()->Set(ch,raw);
      IC.raw[ch] = raw;
      IC.sum += raw;

      // cout << ch << " " << raw << " " << IC.sum <<endl;

      Histo->ICSummary->Fill(ch,raw);
      
    }
    IC.sum /= 16.;
  } 
  else 
  {
    // Old style IC packet
    UShort_t words = (*(p-2))-2;
    while(words > 0)
    {
      ch = ((*p)&0xf000)>>12;
      raw = (*p++)&0xfff;
      words--;
      //this->G      cout << ch << " " << raw << endl;etIonChamber()->Set(ch,raw);
      IC.raw[ch] = raw;
      IC.sum += raw;
    }
  }

  //IC.Scaledsum = IC.sum * IC.ScalingFactor;
  return p;
}

unsigned short* S800::DecodeS800CrdcRaw(unsigned short *p, int id){
  // cout << "Begin CRDC "<< id << " decoding" << endl;
  static ULong_t total=0, failed=0;
  Short_t sampleBegin = 0, sampleWidth =0, isample, ichannel, cdata[4], connector, previous_sample = 0, ch;
  Short_t sindex=0;
  Short_t previous_channel=0;
  Short_t maxwidth = S800_CRDC_MAXWIDTH;
  Short_t channels;
  channels = S800_FP_CRDC_CHANNELS;
    
  unsigned short *pStore = p;
  bool mes1=true, mes2=true, mes3=true, mes4=true;
  bool debug = S800_DEBUG;
  UShort_t length = *p++;
  short i = length-3;
  p++;  // skip packet id
  UShort_t threshold = *p++;

  //  cout << "Decoding CRDC " << id << " " << i << endl;
  while(i > 0)
  {
    if ((*p)>>15 != 1) 
    {
      std::cout << "DecodedEvent: " << "CRDC data is corrupted!" << std::endl;
      p++; i--;
      continue;
    }
    else
    {
      isample  = ((*p)&0x7FC0)>>6;
      ichannel =  (*p)&0x003F;
      if (i == length-3)
      {
        sampleBegin     = isample; 
        previous_sample = isample;
      }
//      cout << previous_channel << " " << ichannel << endl;
      if(previous_channel > ichannel) sindex++;
      previous_channel = ichannel;
    }
    p++; i--;
    memset(cdata, 0, sizeof(cdata));
    while ((*p)>>15 == 0) 
    {
      connector = ((*p)&0x0C00)>>10;
      cdata[connector] = (*p)&0x3FF;
      p++; i--;
      if (i == 0) break;
    }
    if(isample < sampleBegin || isample > sampleBegin+maxwidth)
    {
      if(debug)
        printf("Warning in Crdc Unpack: inconsistent sample number: %d (first: %d)\n", isample, sampleBegin);
      mes1 = false;
      //  continue;
    }
    if(isample < previous_sample)
    {
      if(debug)
        printf("Warning in Crdc Unpack: sample number lower than previous: %d (previous: %d)\n", isample, previous_sample);
      mes2 = false;
      //      continue;
    }
    previous_sample = isample;
    for(int j=0; j<4; j++)
    {
      ch = ichannel+j*64;
      if (cdata[j] != 0 && ch < channels)
      {
        if (cdata[j] < threshold)
        {
          if (debug)
            printf("Warning in Crdc Unpack: data lower than threshold: %d (threshold: %d)\n", cdata[j], threshold);
      	  mes3 = false;
	      }
        else
        {
	        //std::cout << "ch " << ch << " cdata[j]" << cdata[j] << " isample " << isample << std::endl;
	        //	  this->GetCrdc(id)->Set(cdata[j], isample, ch);
	        CRDC[id].cdata[ch][sindex] = cdata[j];
	      }
      } 
      else if (cdata[j] != 0 && ch >= channels)
      {
        if (debug) {
	        printf("Warning in Crdc Unpack: channel greater than limit: %d (limit: %d)\n", ch, channels);
	      }
	      mes4 = false;
      }
    }
    //    sampleWidth = isample - sampleBegin + 1;
//    cout << sindex << endl;
    sampleWidth = sindex + 1;

//cout << "sample width = " << sampleWidth << " " << isample - sampleBegin+1  << endl;
    //    this->GetCrdc(id)->SetSampleWidth(sampleWidth);
    CRDC[id].samplewidth = sampleWidth;
  }
 
  CRDCIntegrate(id);
  
  if (!mes1 || !mes2 || !mes3 || !mes4) failed++;
  total++;
  if (failed == 1000) {
    if (debug)
      printf ("Errors in Crdc Unpackings: %g%%\n", 1.0*failed/total*100);
    total = 0;
    failed = 0;
  }
  return (pStore+length);
}

void S800::CRDCIntegrate(int id)
{
  int PadMult=0;
  for(int i=0;i<NPads;i++)
  {
    float raw = -1;
    float cal = -1;
    int nsamples=0;
    for(int s=0;s<CRDC[id].samplewidth;s++)
    {
	    cout << "crdc " << id << " pad " << i << " sample "<< CRDC[id].cdata[i][s] << endl;
  	  if(CRDC[id].cdata[i][s] !=0)
      {
        nsamples++;
        raw += (CRDC[id].cdata[i][s] - CRDCPeds[id][i]);
        raw = (raw +1)/nsamples;
        cal = raw * Chargeslope[id][S800Setting][i] 
                  + Chargeintercept[id][S800Setting][i];

        //cout << "raw " << raw << " id " << id << " setting " << S800Setting << " i " << i << endl;
        //cout << "  slope " << Chargeslope[id][S800Setting][i] << "  int " << Chargeintercept[id][S800Setting][i] << endl;
        cal = raw;
      }
	  }
    if(nsamples>0)
  	{
  	  CRDC[id].raw[PadMult] = raw;
  	  CRDC[id].cal[PadMult] = cal;
  	  CRDC[id].pad[PadMult] = i;
  	  if(CRDC[id].raw[PadMult] !=0) PadMult++;
  	  
  	}
  }
  //  cout << "Padmult = "<< PadMult << endl;
  CRDC[id].PadMult=PadMult;
}



unsigned short* S800::DecodeS800NewMultiHitTDC(unsigned short *p){
  //Data should be interpreted in 16-bit words 

  //Declare temporary arrays to hold the raw data
  unsigned short data[32][32];
  signed short hits[32];
  unsigned short raw[32];
  for (int i=0;i<32;i++){
    hits[i]=-1;
    raw[i]=0;
    for (int j=0;j<32;j++){
      data[i][j]=0;
    }
  }

  UShort_t length, ch, hit;
  length = *(p-2);
  length -= 2;
  while (length > 0) {
    ch = (*p)&0xFF;
    hit = (*p++)>>8;
    if (hit < 32)
      data[ch][hit] = *p++;
    else
      p++;
    if (hit == 0) raw[ch] = data[ch][0];
    if (hit > hits[ch]) hits[ch] = hit;
    length -= 2;
  }
//MTDC channel assignment table
//////////////////////////////////////////////////////////////
//E1up        	0 	LeCroy Var. ampl. -> Mesytec MCFD ch #0 -> ECL-NIM -> Fan in/out -> 
//                  Fan in/out -> Fan in/out -> Gate Generator -> NIM-ECL
//E1down	      1 	LeCroy Var. ampl. -> Mesytec MCFD ch #1 -> ECL-NIM -> Fan in/out -> 
//                  NIM-ECL
//XFP 	        2 	Patch #1 (dU6) -> CANBERRA CFD (dU6) -> Patch #70 -> Fan in/out ->
//                  Mesytec MCFD ch #2 -> ECL-NIM -> Fan in/out -> NIM-ECL
//OBJ 	        3 	Patch #94 -> LeCroy Var. ampl. -> Mesytec MCFD ch #3 -> ECL-NIM -> 
//                  Fan in/out -> NIM-ECL
//Free 	        4 	
//RF 	          5 	Patch #69 -> Fan in/out -> Logic Unit -> ECL-NIM -> Fan in/out -> 
//                  NIM-ECL
//CRDC1 Anode 	6 	Tennelec Ampl. -> CANBERRA CFD -> ECL-NIM -> Fan in/out -> NIM-ECL
//CRDC2 Anode 	7 	Tennelec Ampl. -> CANBERRA CFD -> ECL-NIM -> Fan in/out -> NIM-ECL
//XFP 	        8 	Patch #1 (dU6) -> CANBERRA CFD (dU6) -> Patch #70 -> Fan in/out -> 
//                  NIM-ECL
//OBJ 	        9 	Patch #54 (dU6) -> CANBERRA CFD (dU6) -> Patch #62 -> Fan in/out -> 
//                  NIM-ECL
//Free 	        10-11 	
//Hodosc. OR 	  12 	Splitter att. -> CANBERRA CFD -> Fan in/out -> NIM-ECL
//Free  	      13 	
//E1 up 	      14 	LeCroy Var. ampl. -> Mesytec MCFD ch #0 -> ECL-NIM -> Fan in/out ->
//                  Fan in/out -> Fan in/out -> NIM-ECL
//E1 up 	      15 	LeCroy Var. ampl. -> Mesytec MCFD ch #0 -> ECL-NIM -> Fan in/out -> 
//                  Fan in/out -> NIM-ECL
//II PPAC2 anode 	16 	Mesytec MCFD2 ch #0
//II PPAC2 down 	17 	Mesytec MCFD2 ch #1
//II PPAC2 up 	  18 	Mesytec MCFD2 ch #2
//II PPAC2 right 	19 	Mesytec MCFD2 ch #3
//II PPAC2 left 	20 	Mesytec MCFD2 ch #4
//II PPAC1 anode 	21 	Mesytec MCFD2 ch #5
//II PPAC1 down 	22 	Mesytec MCFD2 ch #6
//II PPAC1 up 	  23 	Mesytec MCFD2 ch #7
//II PPAC1 right 	24 	Mesytec MCFD2 ch #8

  //raw[15] is ch15 which is E1 up

  //cout << "e1 up on 15: " << raw[15] << endl;
  //raw[15] = raw[14];

  if (raw[15] != 0) {
    for (int i=0; i<13; i++) {
      switch(i) {
        case 0: // e1up
	        if (hits[0] >= 0){   
	          //	  fMultiHitTOF.fE1Up.push_back( (data[0][0] - raw[15]) * 0.0625);
	          mTDC.e1up.push_back( (data[0][0] - raw[15]) * 0.0625);
            //cout << "e1up = (" << data[0][0] << " - " << raw[15] << ")*0.0625 = " << mTDC.e1up.at(0) << endl;
	        }
	        break;
        case 1: // e1down
	        if (hits[1] >= 0){
	          //	  fMultiHitTOF.fE1Down.push_back((data[1][0] - raw[15]) * 0.0625);
	          mTDC.e1down.push_back((data[1][0] - raw[15]) * 0.0625);
            //cout << "e1down = (" << data[1][0] << " - " << raw[15] << ")*0.0625 = " << mTDC.e1down.at(0) << endl;
	        }
	        break;
        case 2: // xf
	        if (hits[2] >= 0) {
	          for (int j=0; j<=hits[2]; j++){
	            //	    fMultiHitTOF.fXf.push_back((data[2][j] - raw[15]) * 0.0625);
	            mTDC.xfp.push_back((data[2][j] - raw[15]) * 0.0625);
              //cout << "xfp" << j << " = (" << data[2][j] << " - " << raw[15] << ")*0.0625 = " << mTDC.xfp.at(j) << endl;
	          }
	        }
	        break;
        case 3: // obj
	        if (hits[3] >= 0) {
	          for (int j=0; j<=hits[3]; j++){
	            //	    fMultiHitTOF.fObj.push_back( (data[3][j] - raw[15]) * 0.0625);
	            mTDC.obj.push_back( (data[3][j] - raw[15]) * 0.0625);
              //cout << "obj" << j << " = (" << data[3][j] << " - " << raw[15] << ")*0.0625 = " << mTDC.obj.at(j) << endl;
	          }
	        }
	        break;
        case 4: // galotte
	        if (hits[4] >= 0) {
	          for (int j=0; j<=hits[4]; j++){
	            // fMultiHitTOF.fGalotte.push_back( (data[4][j] - raw[15]) * 0.0625);
	            mTDC.galotte.push_back( (data[4][j] - raw[15]) * 0.0625);
	          }
	        }
	        break;
        case 5: // rf
	        if (hits[5] >= 0) {
	          for (int j=0; j<=hits[5]; j++){
	            //	    fMultiHitTOF.fRf.push_back( (data[5][j] - raw[15]) * 0.0625);
	            mTDC.rf.push_back((data[5][j] - raw[15]) * 0.0625);
	          }
	        }
	        break;
        case 12: // hodoscope
	        if (hits[12] >= 0) {
	          for (int j=0; j<=hits[12]; j++){
	            //	    fMultiHitTOF.fHodoscope.push_back( (data[12][j] - raw[15]) * 0.0625);
	            mTDC.hodoscope.push_back((data[12][j] - raw[15]) * 0.0625);
	          }
	        }
	        break;
        default:
	        break;
      }//end Switch i
    }//end for i
  }//end if raw[15]!=0

  // cout<<"fNewTOF.fE1Up.size() "<<fNewTOF.fE1Up.size()<<endl;
  // cout<<"fNewTOF.fE1Down.size() "<<fNewTOF.fE1Down.size()<<endl;
  // cout<<"fNewTOF.obj.size() "<<fNewTOF.fObj.size()<<endl;

  return p;
}


void S800::CalculateGravity(int id)
{
  //cout << "CalculateGravity" << endl;
  int maxpad = 0;
  int padmax = 0;
  //cout << "CRDC[" << id << "].PadMult " << CRDC[id].PadMult << endl;


  for(int i=0;i<CRDC[id].PadMult-1;i++)
  {
    //cout << "CRDC[" << id << "].cal[" << i << "] = " << CRDC[id].cal[i] << endl;
    //cout << "CRDC[" << id << "].pad[" << i << "] " << CRDC[id].pad[i] << endl;
    //if(CRDC[id].cal[i] > 0 && CRDC[id].cal[i+1] > 0 && CRDC[id].cal[i] > maxpad)
    //need to change the above IF statement, the mpdc doesn't always come in increasing pad order
    if(CRDC[id].cal[i] > maxpad)
    {
      maxpad = CRDC[id].cal[i];
      padmax = CRDC[id].pad[i];
    }
  }
  //cout << "1) maxpad " << maxpad << endl;
  //cout << "   PadMult " << CRDC[id].PadMult << endl;
  if(CRDC[id].PadMult%2 != 0) padmax = CRDC[id].pad[CRDC[id].PadMult/2];
  else  padmax = CRDC[id].pad[CRDC[id].PadMult/2+1];

  //cout << "2) padmax " << padmax << endl;

  if(padmax == 0)
  {
    CRDC[id].padsum = 0;
    CRDC[id].x_gravity = -500;
    //cout << "padmax = 0 for crdc " << id << " mult = " << CRDC[id].PadMult <<endl;
    return;
  }

  int lowpad = (int)padmax - (int)gravity_width/2;
  int highpad = lowpad + (int)gravity_width -1;
  if(lowpad <0)
  {
    lowpad =0;
  }
  if(highpad >= 480)
  {
    highpad = 480;
  }

  float sum=0.0, mom=0.0;
  for(int i=0;i<CRDC[id].PadMult;i++)
  {
    if(CRDC[id].pad[i] >=lowpad && CRDC[id].pad[i] <=highpad)
    {
      sum += CRDC[id].cal[i];
      mom += CRDC[id].pad[i]*CRDC[id].cal[i];
    }
  }
  float pad = mom/sum;
  CRDC[id].padsum = sum;
  CRDC[id].mom = mom;
  CRDC[id].x_gravity = mom/sum * CRDCSlopeX[id] + CRDCOffsetX[id];
  //cout << "CRDC[" << id << "].x_gravity " << CRDC[id].x_gravity << endl;
}


bool S800::analyze()
{
  //for(int i=0;i<3;i++) //will never really have e2 and e3 so no loop needed
  //{
  //cout << "e1up[" << i << "] = " << Scint[i].de_up << endl;
  Histo->e1up->Fill(Scint[0].de_up);
  Histo->e1down->Fill(Scint[0].de_down);
  //}
  for(int i=0;i<(int)mTDC.e1up.size();i++)
  {
    Histo->Te1up->Fill(mTDC.e1up.at(i));
  }
  for(int i=0;i<(int)mTDC.e1down.size();i++)
  {
    Histo->Te1down->Fill(mTDC.e1down.at(i));
  }
  for(int i=0;i<(int)mTDC.obj.size();i++)
  {
    Histo->Tobj->Fill(mTDC.obj.at(i));
  }
  //cout << "filling Tobj_mult " << mTDC.obj.size() << endl;
  Histo->Tobj_mult->Fill((int)mTDC.obj.size());
  
  Histo->CRDC1PadMult->Fill(CRDC[0].PadMult);
  Histo->CRDC2PadMult->Fill(CRDC[1].PadMult);
  Histo->CRDC1Tac->Fill(CRDC[0].tac);
  Histo->CRDC2Tac->Fill(CRDC[1].tac);
  Histo->CRDC1Y->Fill(CRDC[0].calY);
  Histo->CRDC2Y->Fill(CRDC[1].calY);
  Histo->CRDC1AnodevsTac->Fill(CRDC[0].anode,CRDC[0].tac);
  Histo->CRDC2AnodevsTac->Fill(CRDC[1].anode,CRDC[1].tac);
  
  for(int i=0;i<2;i++)
  {
    CalculateGravity(i);
    //      cout << "CRDC " << i << " Mult = " << CRDC[i].PadMult << endl;
    for(int p=0;p<CRDC[i].PadMult;p++)
    {
      if(i==0)
      {
        //cout << "CRDC " << i << " " << p << " raw = " << CRDC[i].raw[p] << endl;
        Histo->CRDC1raw->Fill(CRDC[i].raw[p]);
        Histo->CRDC1Summary->Fill(CRDC[i].pad[p],CRDC[i].raw[p]);
        Histo->CRDC1Summary_cal->Fill(CRDC[i].pad[p],CRDC[i].cal[p]);
      }
      else
      {
        Histo->CRDC2raw->Fill(CRDC[i].raw[p]);
        Histo->CRDC2Summary->Fill(CRDC[i].pad[p],CRDC[i].raw[p]);
	      Histo->CRDC2Summary_cal->Fill(CRDC[i].pad[p],CRDC[i].cal[p]);
      }
    }
  }

  if(CRDC[0].x_gravity == -500 || CRDC[1].x_gravity ==-500)
  {
    //cout << "0 x_gravity: " << CRDC[0].x_gravity << "    1 x_gravity " << CRDC[1].x_gravity << endl;
    return false;
  }

  Histo->CRDC1X->Fill(CRDC[0].x_gravity);
  Histo->CRDC2X->Fill(CRDC[1].x_gravity);
  //cout << "0 x_gravity " << CRDC[0].x_gravity << " " << CRDC[0].tac << endl;
  Histo->CRDC1XY->Fill(CRDC[0].x_gravity,CRDC[0].tac);
  Histo->CRDC2XY->Fill(CRDC[1].x_gravity,CRDC[1].tac);
  Histo->CRDC1XrawY->Fill(CRDC[0].mom/CRDC[0].padsum, CRDC[0].tac);
  Histo->CRDC2XrawY->Fill(CRDC[1].mom/CRDC[1].padsum, CRDC[1].tac);
  Histo->CRDC1XCalY->Fill(CRDC[0].x_gravity,CRDC[0].calY);
  Histo->CRDC2XCalY->Fill(CRDC[1].x_gravity,CRDC[1].calY);

  if((int)mTDC.obj.size()>0)
  {
    if((int)mTDC.xfp.size()>0)
    {
      Histo->ObjvsXFP->Fill(mTDC.obj.at(0),mTDC.xfp.at(0));
      Histo->ObjvsICsum_nobeam->Fill(mTDC.obj.at(0),IC.sum); //ungated pid
      bool stat = beam_pid->getPID(mTDC.obj.at(0),mTDC.xfp.at(0));
      //if (beam_pid->Z > 0) cout << beam_pid->Z << " " << beam_pid->A << endl;
      //cout << stat << endl;
      //cout << "did i get a beam particle? " << stat << endl;
      if(!stat)
      {
      //  cout << "   no :(" << endl;
        return false;
      }
      GetBeamId(beam_pid->Z); //assigns beamID 0-Ca37; 1-K36; 2-Ar35; 3-Cl34;

      track->CRDC[0] = CRDC[0];
      track->CRDC[1] = CRDC[1];

      //Tracking through the focal plane
      //this is for the initial tracking, will switch A, Z later and tack again
      //TODO tracking for initial Si-23?
      track->CalculateTracking(23, 14); //A,Z for K35
      float correctedObj = mTDC.obj.at(0) + ObjCorr[0]*track->afp +
                           ObjCorr[1]*CRDC[0].x_gravity;
      float correctedXFP = mTDC.xfp.at(0) + ObjCorr[0]*track->afp +
                           ObjCorr[1]*CRDC[0].x_gravity;

      mTDC.objCorrected.push_back(correctedObj);

      //Si23 beam
      if (beam_pid->Z == 14 && beam_pid->A == 23)
      {
        Histo->ObjvsICsum_Si23->Fill(mTDC.objCorrected.at(0),IC.sum);
      }
      //Al22 beam
      if (beam_pid->Z == 13 && beam_pid->A == 22)
      {
        Histo->ObjvsICsum_Al22->Fill(mTDC.objCorrected.at(0),IC.sum);
      }
      //Mg21 beam
      if (beam_pid->Z == 12 && beam_pid->A == 21)
      {
        Histo->ObjvsICsum_Mg21->Fill(mTDC.objCorrected.at(0),IC.sum);
      }
      //Na20 beam
      if (beam_pid->Z == 11 && beam_pid->A == 20)
      {
        Histo->ObjvsICsum_Na20->Fill(mTDC.objCorrected.at(0),IC.sum);
      }

      if(beam_pid->Z == 20 && beam_pid->A == 37)
      {
        Histo->ObjvsICsum_Ca37->Fill(mTDC.objCorrected.at(0),IC.sum);
        Histo->ObjUncvsICsum_Ca37->Fill(mTDC.obj.at(0),IC.sum);
        if (mTDC.e1up.size() >0)
          Histo->Timing1vsICsum_Ca37->Fill(mTDC.objCorrected.at(0)-mTDC.e1up.at(0),IC.sum);
        if (mTDC.e1down.size() >0)
          Histo->Timing2vsICsum_Ca37->Fill(mTDC.objCorrected.at(0)-mTDC.e1down.at(0),IC.sum);
        if (mTDC.xfp.size() >0)
        {     
          Histo->XFPvsICsum_Ca37->Fill(correctedXFP,IC.sum);
          //Histo->XFPUncvsICsum_Ca37->Fill(mTDC.xfp.at(0),IC.sum);
        }
      }
      if(beam_pid->Z == 19 && beam_pid->A == 36)
      {
        Histo->ObjvsICsum_K36->Fill(mTDC.objCorrected.at(0),IC.sum);
      }
      if(beam_pid->Z == 18 && beam_pid->A == 35)
      {
        Histo->ObjvsICsum_Ar35->Fill(mTDC.objCorrected.at(0),IC.sum);
      }
      if(beam_pid->Z == 17 && beam_pid->A == 34)
      {
        Histo->ObjvsICsum_Cl34->Fill(mTDC.objCorrected.at(0),IC.sum);
      }
//cout << mTDC.obj.at(0) << " " << mTDC.objCorrected.at(0) << endl;
      Histo->ObjvsICsum->Fill(mTDC.obj.at(0),IC.sum);
      Histo->ObjvsICsum_corr->Fill(mTDC.objCorrected.at(0),IC.sum);
  
      //Reset PID flags
      Ne17_flag = false;
      Ne18_flag = false;

      if(beam_pid->Z >0 && beam_pid->A >0)
      {
        //cout << S800Setting << " " << BeamID << endl;
        bool stat2 = residue_pid[S800Setting][BeamID]->getPID(mTDC.objCorrected.at(0),IC.sum);
        //if (residue_pid[S800Setting][BeamID]->A == 22) cout << "here" << endl;
        //cout << "here " << stat2 << endl;
        //cout << "here0" << endl;
        if(!stat2) return false;

        //Set PID flags
        if (residue_pid[S800Setting][BeamID]->Z == 10 && residue_pid[S800Setting][BeamID]->A == 17) Ne17_flag = true;
        if (residue_pid[S800Setting][BeamID]->Z == 10 && residue_pid[S800Setting][BeamID]->A == 18) Ne18_flag = true;

        //TODO need more loss files
       // if (residue_pid[S800Setting][BeamID]->Z == 10) return false;

        track->CRDC[0] = CRDC[0];
        track->CRDC[1] = CRDC[1];
        //cout << residue_pid[S800Setting][BeamID]->Z << " " << residue_pid[S800Setting][BeamID]->A << endl;
        //if (residue_pid[S800Setting][BeamID]->Z == 10) cout << residue_pid[S800Setting][BeamID]->A << endl;
        //second time tracking, this time calculating using the correct A,Z
        //cout << "here" << endl;
        track->CalculateTracking(residue_pid[S800Setting][BeamID]->A, residue_pid[S800Setting][BeamID]->Z);
        track->GetThetaPhi(track->ata,track->bta);
        Histo->atavsbta->Fill(track->ata*180./acos(-1),track->bta*180./acos(-1));
        Histo->ata1D->Fill(track->ata*180./acos(-1));
        Histo->bta1D->Fill(track->bta*180./acos(-1));
        //Looking at corrections
        Histo->Objvsafp->Fill(mTDC.objCorrected.at(0),track->afp);
        Histo->ObjvsCRDC1X->Fill(mTDC.objCorrected.at(0),CRDC[0].x_gravity);
        Histo->ObjvsCRDC2X->Fill(mTDC.objCorrected.at(0),CRDC[1].x_gravity);
        //cout << "XfpUncvsCRDC1X " << mTDC.xfp.at(0) << " " << CRDC[0].x_gravity << endl;
        Histo->XfpUncvsCRDC1X->Fill(mTDC.xfp.at(0),CRDC[0].x_gravity);
        Histo->XfpUncvsCRDC2X->Fill(mTDC.xfp.at(0),CRDC[1].x_gravity);

        Histo->ThetavsPhi->Fill(track->thetadeg,track->phideg);
      }
      else 
      {
        cout << "\nBAD PID" << endl;
        return false;
      }
    }
    else
    {
	    return false;
    }
  }
  else
  {
    //cout << "\n" << (int)mTDC.obj.size() << "\t" << (int)mTDC.obj.size() << endl;
    return false;
  }

  return true;
}

/************************************************************/
/* S800Map Class - Functions                                */
/************************************************************/

S800Map::S800Map() {
  /* Parameters */
  //  m_top = s800;

  /* Variables */
  maxcoefficients = S800_TRACK_COEFFICIENTS;
  maxparameters = S800_TRACK_PARAMETERS;
  maxorder = 0.0;

  /* Internal use values */
  for (int i=0; i<S800_TRACK_PARAMETERS; i++) {
    maxcoefficient[i] = 0;
    for (int j=0; j<S800_TRACK_COEFFICIENTS; j++) {
      order[i][j] = 0;
      coefficient[i][j] = 0;
      for (int k=0; k<S800_TRACK_PARAMETERS; k++) {
        exponent[i][k][j] = 0;
      }
    }
  }

  // LoadInverseMap(filename);

}

/* destructor */
S800Map::~S800Map()
{
}

void S800Map::LoadInverseMap(TString filename) {
  cout << "Loading S800 inverse map " << filename.Data() << "...";

  FILE *file;
  if ( (file = fopen(filename.Data(), "r")) == NULL ) {
    printf("\n Inverse map file %s was not found!\n", filename.Data());
    return;
  }
  char line[400];

  unsigned short index;
  double icoefficient;
  unsigned short iorder;
  unsigned short iexp0, iexp1, iexp2, iexp3, iexp4, iexp5;
  int parameter = 0;
  char* retVal;

  retVal = fgets(line, 400, file);
  while (strncmp(line, "     I", 6) != 0) { 
    retVal = fgets(line, 400, file);
  }
  parameter = 0;
  
  while (!feof(file)) {
    retVal = fgets(line, 400, file);
    while (strncmp(line, "    ---", 7) != 0) {
      sscanf(line, "%hd %lf %hd %hd %hd %hd %hd %hd %hd", &index, 
             &icoefficient, &iorder, &iexp0, &iexp1, &iexp2, &iexp3, 
             &iexp4, &iexp5);
      if (iorder > maxorder) { maxorder = (double)iorder; }
      if (index > maxcoefficients) {
        cout << "Error: too many coefficients in map" << endl;
        cout << " Please increase S800_TRACK_COEFFICIENTS and "<< endl
             << " recompile. " << endl;
        break;
      }

      if (parameter >= S800_TRACK_PARAMETERS) {
        cout << "Invalid parameter number, must be between 0 and " 
             << S800_TRACK_PARAMETERS-1 << endl;
      }
      if (index >= S800_TRACK_COEFFICIENTS) {
        cout << "Too many coefficients - maximum is " 
             << S800_TRACK_COEFFICIENTS << endl;
      }
      if (index >= maxcoefficient[parameter]) { 
        maxcoefficient[parameter] = index+1;
      }
      order[parameter][index] = iorder;
      exponent[parameter][0][index] = iexp0;
      exponent[parameter][1][index] = iexp1;
      exponent[parameter][2][index] = iexp2;
      exponent[parameter][3][index] = iexp3;
      exponent[parameter][4][index] = iexp4;
      exponent[parameter][5][index] = iexp5;
      coefficient[parameter][index] = icoefficient;
      retVal = fgets(line, 400, file);
    }
    retVal = fgets(line, 400, file);
    parameter++;
    if (parameter > maxparameters) {
      cout << "Error: too many parameters in map" << endl;
      cout << " Please increase S800_TRACK_PARAMETERS and " << endl
           << " recompile. " << endl;
      break;
    }
  }
  fclose(file);  
  cout << "Done!" << endl;  
}

void S800Map::Reset() {
  /* Internal use values */
  for (int i=0; i<S800_TRACK_PARAMETERS; i++) {
    //maxcoefficient[i] = 0;
    for (int j=0; j<S800_TRACK_COEFFICIENTS; j++) {
      //order[i][j] = 0;
      //coefficient[i][j] = 0;
      for (int k=0; k<S800_TRACK_PARAMETERS; k++) {
        //exponent[i][k][j] = 0;
      }
    }
  }
}

double S800Map::Calculate(int calcorder, int parameter, double *input) {
  double cumul = 0;
  double multiplicator;

  for (int index=0; index<maxcoefficient[parameter]; index++) {
    if (calcorder < order[parameter][index]) { break; }
    multiplicator = 1;

    for (int nex=0; nex<S800_TRACK_PARAMETERS; nex++) {
      if (exponent[parameter][nex][index] != 0) {
        multiplicator *= pow(input[nex], exponent[parameter][nex][index]);
      }
    }

    cumul += multiplicator * coefficient[parameter][index];
  }
  return cumul;
}


/************************************************************/
/* S800FpTrack Class - Functions                            */
/************************************************************/

S800FpTrack::S800FpTrack(S800Map ***map0) {

  map = map0;
  
  /* Parameters */
  xfp = -65000.;
  afp = -65000.;
  yfp = -65000.;
  bfp = -65000.;
  ata = -65000.;
  yta = -65000.;
  bta = -65000.;
  dta = -65000.;
  azita = -65000.;
  scatter = -65000.;
  energy = -65000.;
  ptot = -65000.;
  ppar = -65000.;
  ptra = -65000.;

  ata_cor = -65000.;
  bta_cor = -65000.;
  azita_cor = -65000.;
  scatter_cor = -65000.;

  /* Variables */
  anglea = 0.0133169; /* radians */ //starting value of zero
  angleb = 0.0061982; /* radians */

  anglea_cor = 0; /* degree */
  angleb_cor = 0; /* degree */
  bta_ytacor = 0; /* degree/yta */
  ata_dtacor = 0; /* degree/dta */

  //brho = 2.426; /* Tm */
  //mass = 35; /* amu */
  //charge = 19; /* q */
  order = 5; /* order */
  zfp = 0; /* m */
  gecorr = 0; /* keV/dta */
  gap = 1073.0; /* mm */

}
S800FpTrack::~S800FpTrack()
{
}

void S800FpTrack::Reset() {
  /* Parameters */
  xfp = -65000.;
  afp = -65000.;
  yfp = -65000.;
  bfp = -65000.;
  ata = -65000.;
  yta = -65000.;
  bta = -65000.;
  dta = -65000.;
  azita = -65000.;
  scatter = -65000.;
  energy = -65000.;
  ptot = -65000.;
  ppar = -65000.;
  ptra = -65000.;

  ata_cor = -65000.;
  bta_cor = -65000.;
  azita_cor = -65000.;
  scatter_cor = -65000.;
  
  CRDC[0].Reset();
  CRDC[1].Reset();

  theta = -999;
  phi = -999;
  thetadeg = -999;
  phideg = -999;


  //  map.Reset();
}

void S800FpTrack::CalculateTracking(int A, int Z) {
  double input[S800_TRACK_PARAMETERS];
  double pi = 3.14159265;
  double amu = 931.5016;
  double beta, betagamma0, gamma0, gamma, energy0, ptot0;
  
  //pick the right map for tracking
  int mapID = 0;
  if (S800Setting == 1 && Z==14 && A==23)
    mapID = 0;
  else if (S800Setting == 1 && Z==13 && A==22)
    mapID = 1;
  else if (S800Setting == 1 && Z==12 && A==21)
    mapID = 2;
  else if (S800Setting == 1 && Z==12 && A==20)
    mapID = 3;
  else if (S800Setting == 1 && Z==10 && A==17)
    mapID = 4;
  else if (S800Setting == 1 && Z==10 && A==18)
    mapID = 5;
  else if (S800Setting == 1 && Z==13 && A==23)
    mapID = 6;
  else if (S800Setting == 1 && Z==11 && A==20)
    mapID = 7;
  else if (S800Setting == 1 && Z==11 && A==21)
    mapID = 8;
  else if (S800Setting == 1 && Z==12 && A==22)
    mapID = 9;
  else if (S800Setting == 1 && Z==10 && A==19)
    mapID = 10;
  else if (S800Setting == 1 && Z==9 && A==18)
    mapID = 11;
  else if (S800Setting == 1 && Z==9 && A==17)
    mapID = 12;
  else if (S800Setting == 1 && Z==8 && A==15)
    mapID = 13;
  else if (S800Setting == 1 && Z==8 && A==16)
    mapID = 14;
  else if (S800Setting == 1 && Z==10 && A==20)
    mapID = 15;
  else if (S800Setting == 1 && Z==8 && A==14)
    mapID = 16;
  else if (S800Setting == 1 && Z==7 && A==14)
    mapID = 17;
  else if (S800Setting == 1 && Z==7 && A==13)
    mapID = 18;
  else if (S800Setting == 1 && Z==6 && A==12)
    mapID = 19;
  else if (S800Setting == 1 && Z==6 && A==11)
    mapID = 20;
  else if (S800Setting == 1 && Z==14 && A==24)
    mapID = 21;
  
  //TODO Get inverse maps for setting 2
  else if (S800Setting == 2 && Z==14 && A==22)
    mapID = 0;
  else if (S800Setting == 2 && Z==14 && A==23)
    mapID = 1;

  //TODO Get inverse maps for setting 0
  else if (S800Setting == 0 && Z==14 && A==23)
    mapID = 0;
  else if (S800Setting == 0 && Z==12 && A==21)
    mapID = 1;


  else
  {
    cout << "no map setting for Z=" << Z << " A=" << A << " on S800 Setting " << S800Setting << endl;
    abort(); //Don't accidently use the wrong inv map
  }
  
  if (!map[S800Setting][mapID]->WasLoaded()) {cout << "error!!!" << endl; return; }

  afp = atan((CRDC[1].x_gravity - CRDC[0].x_gravity) / gap);
  bfp = atan((CRDC[1].calY - CRDC[0].calY) / gap);
   
  xfp = CRDC[0].x_gravity / 1000. + zfp * tan(afp);
  yfp = CRDC[0].calY / 1000. + zfp * tan(bfp);
   
  /* x (dispersive) axis is reversed between COSY, and the
     local coordinate system, where x>0 is down. */
  if (afp == -65000.) 
  {
    cout << "\nBAD afp " << (CRDC[1].x_gravity - CRDC[0].x_gravity) / gap
         << "\t" << CRDC[1].x_gravity << "\t" << CRDC[0].x_gravity << endl;
    return;
  }
 
  input[0] = -xfp;
  input[1] = -afp;
  input[2] = yfp;
  input[3] = bfp;
  if ((double)order > (double)map[S800Setting][mapID]->maxorder)
  {
    order = map[S800Setting][mapID]->maxorder;
  }
  ata = map[S800Setting][mapID]->Calculate((int)order, 0, input);
  yta = map[S800Setting][mapID]->Calculate((int)order, 1, input);
  bta = map[S800Setting][mapID]->Calculate((int)order, 2, input);
  dta = map[S800Setting][mapID]->Calculate((int)order, 3, input);
  ata += anglea; /* Add theta offset in radians */
  bta += angleb; /* Add phi offset in radians */
  
  ata_cor = ata + dta*ata_dtacor;
  bta_cor = bta + yta*bta_ytacor;
  
  ata_cor += anglea_cor/180.*pi;
  bta_cor += angleb_cor/180.*pi;
  
  double xsin = sin(ata);
  double ysin = sin(bta);

  
  if ( (xsin > 0) && (ysin > 0) ) {
    azita = atan(ysin/xsin);
  } else if ( (xsin < 0) && (ysin > 0) ) {
    azita = pi - atan(ysin/fabs(xsin));
  } else if ( (xsin < 0) && (ysin < 0) ) {
    azita = pi + atan(fabs(ysin)/fabs(xsin));
  } else if ( (xsin > 0) && (ysin < 0) ) {
    azita = 2*pi - atan(fabs(ysin)/xsin);
  } else {
    azita = 0.0;
  }
  
  scatter = asin(sqrt( (xsin*xsin) + (ysin*ysin) ))*1000; /* mrad */
  /* Same for corrected ata_cor and bta_cor */
  xsin = sin(ata_cor);
  ysin = sin(bta_cor);
  
  if ( (xsin > 0) && (ysin > 0) ) {
    azita_cor = atan(ysin/xsin);
  } else if ( (xsin < 0) && (ysin > 0) ) {
    azita_cor = pi - atan(ysin/fabs(xsin));
  } else if ( (xsin < 0) && (ysin < 0) ) {
    azita_cor = pi + atan(fabs(ysin)/fabs(xsin));
  } else if ( (xsin > 0) && (ysin < 0) ) {
    azita_cor = 2*pi - atan(fabs(ysin)/xsin);
  } else {
    azita_cor = 0.0;
  }
  scatter_cor = asin(sqrt( (xsin*xsin) + (ysin*ysin) ))*1000; /* mrad */

  mass = A;
  charge = Z;
  betagamma0 = brho / 3.107 * charge / mass;
  gamma0 = sqrt( (betagamma0 * betagamma0) + 1);
  beta0 = betagamma0 / gamma0;
  energy0 = mass * amu * (gamma0 - 1); /* MeV */
  ptot0 = energy0 * sqrt(1 + 2 * mass * amu / energy0); /* MeV/c */
  energy = (1 + dta) * energy0; /* MeV */
 
  gamma = 1 + energy / mass / amu;
  beta = sqrt(1 - 1 / gamma / gamma);
  deltabeta = beta;// - beta0;
  // ptot = energy * sqrt(1 + 2 * mass * amu / energy); /* MeV/c */
  ptot = sqrt(energy*energy+2*mass*amu*energy);
  ppar = ptot * cos(scatter/1000.);
  ptra = ptot * sin(scatter/1000.);
}

void S800FpTrack::GetThetaPhi(double ata0, double bta0)
{
  double sineata = sin(ata0);
  double sinebta = sin(bta0);

  //  cout << "ata0 = " << ata0 << " bta0" << bta0 << endl;
  //cout << sineata << " " << sinebta << endl;

  double result = sqrt(sineata*sineata + sinebta*sinebta);

  if(result >1)
  {
    result =1;
  }
  else if(result < -1)
  {
    result =-1;
  }
  phi = atan2(sinebta,sineata);
  //  theta = asin(sqrt(sineata*sineata + sinebta*sinebta));
  theta = asin(result);
  
  phideg = phi*180./acos(-1);
  thetadeg = theta*180./acos(-1);
  
}

void S800::GetBeamId(int Z)
 {
  // beam numbers
  // 0 - Si23
  // 1 - Al22
  // 2 - Mg21
  // 3 - Na20

  if(Z==14)
    BeamID = 0;
  else if(Z==13)
    BeamID = 1;
  else if(Z==12)
    BeamID = 2;
  else if(Z==11)
    BeamID = 3;
  else
  {
    cout << "That beam is not associated with a good resiude PID, Z = " << Z << endl;
    cout << "Need a new pid in the array" << endl;
    abort();
  }
   
 }
 
