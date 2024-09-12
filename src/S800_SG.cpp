#include "S800.h"


S800::S800(TRandom * ran0, histo_sort * Histo0)
{
  Histo = Histo0;
  ran = ran0;

  Init();
}
S800::S800(TRandom * ran0, histo_read * Histo1)
{
  Histo_read = Histo1;
  ran = ran0;

  Init();
}

S800::S800(TRandom * ran0, histo_sort * Histo0, histo_read * Histo1)
{
  Histo = Histo0;
  Histo_read = Histo1;
  ran = ran0;

  Init();
}

void S800::Reset()
{
  mTDC.Reset();
  IC.Reset();
  CRDC[0].Reset();
  CRDC[1].Reset();
  beam_pid->Z = 0;
  beam_pid->A = 0;

  for(int i =0;i<2;i++)
    {
      for(int j=0;j<5;j++)
	{
	  residue_pid[i][j]->Z = 0;
	  residue_pid[i][j]->A = 0;
	}
    }

  for(int i = 0; i < 2; i++) track[i]->Reset();

  BeamID = -99;

  
}

S800::~S800()
{
  delete beam_pid;
  //  delete residue_pid;
  //delete residue_pid_Ar30;
  //delete residue_pid_Cl30;
  // delete residue_pid_O14;
  // delete residue_pid_Ne17;
}


void S800::Init()
{

  S800Setting =-99;
  
  string name = "s800_beam";
  beam_pid = new pid(name);



  /////////////////////////////////////////////////////////////////////////////////////////
  //Load in the residue PID gates
  //Setting 0 will be the Cl setting and setting 1 is the Ar settting
  //
  // beam numbers
  // 0 - Ar31
  // 1 - S29
  // 2 - P28
  // 3 - Si27
  // 4 - Al26


  string beams[5] = {"Ar31","S29","P28","Si27","Al26"};
		     
  // cout << "here" << endl;

  residue_pid = new pid**[2]; //Two arrays of residue pids for Ar and Cl settings

  for(int i =0;i<2;i++) //Loop over setttings
    {
      residue_pid[i] = new pid*[5];

      for(int j=0;j<5;j++) //Loop over beams
	{
	  if(i==0)
	    name = Form("ClSetting/s800_residue_Cl_%s",beams[j].c_str());
	  else
	    name = Form("ArSetting/s800_residue_Ar_%s",beams[j].c_str());

	  residue_pid[i][j] = new pid(name);
	}
    }

  //cout << "here 2" << endl;
  
  
  // //Cl30 rigidity
  // name = "s800_residue_Cl";
  // residue_pid_Cl30 = new pid(name);
    
  
  // name = "s800_residue_Ar";
  // residue_pid = new pid(name); //default to Ar30 setting
  // residue_pid_Ar30 = new pid(name);


  /////////////////////////////////////////////////////////////////////////////////////////
  
    ifstream infile("datfiles/CRDCPedestals.dat");
  if(!infile.is_open())
    {
      cout << "Couldn't open CRDC Pedestal file" << endl;
      abort();
    }

  int CRDC,pad;
  float ped;
  getline(infile,name);
  for(;;)
    {
      infile >> CRDC >> pad >> ped;

      CRDCPeds[CRDC-1][pad] = (int)ped;
      if(infile.eof()) break;
      if(infile.bad()) break;
    }
  infile.close();
  infile.clear();


  ifstream infile2("datfiles/CRDCpars.dat");
  if(!infile2.is_open())
    {
      cout << "Couldn't open CRDCpars file" << endl;
      abort();
    }
  getline(infile2,name);
  
  float slope, intercept;
  for(;;)
  {
    infile2 >> CRDC >> pad >> slope >> intercept;
    Chargeslope[CRDC][pad] = slope;
    Chargeintercept[CRDC][pad] = intercept;
    //cout << CRDC << "   " << pad << "   " << slope << "   " << intercept << endl;
    if(infile2.eof()) break;
    if(infile2.bad()) break;
  }
  infile2.close();
  infile2.clear();
  
  

  // //Read in the P10 gas pressure for the IC
  // //This will define the scaling factorfor the ICSum energy
  // for(int i=0;i<200;i++)
  //   {
  //     Pressure[i] = -1;
  //     Scaling[i] = 1;
  //   }
  
  // ifstream Gasin("gaslog_simple.txt");
  // if(!Gasin.is_open())
  //   {
  //     cout << "No gas pressure file found for S800" << endl;
  //     abort();
  //   }
  // char dum[200];
  // float in2;
  // int in1;
  // Gasin.getline(dum,200);
  // for(;;)
  //   {
  //     Gasin >> in1 >> in2;
  //     Pressure[in1] = in2;
  //     if(in2 !=0)
  // 	Scaling[in1] = 300./in2;

  //     if(Gasin.eof())break;
  //     if(Gasin.bad())break;
  //   }
  // Gasin.close();
  // Gasin.clear();


  //The inverse map is different for each S800 setting and will
  //need to be chaged.
  char mapname[64];
  sprintf(mapname,"s800inputs/Ar30_S28.inv");
//    sprintf(mapname,"s800inputs/Cl30_%s.inv");
  map = new S800Map(); //initialize the inverse map
  map->LoadInverseMap(mapname);
  track = new S800FpTrack(map);


  //CRDC calibration parameters
  //[CRDC #][X = 0 or Y =1]
  // from run 122 123
  // CRDCSlope[0][0] = 2.54;// mm/pad
  // CRDCSlope[1][0] = 2.54;// mm/pad
  // CRDCSlope[0][1] = -0.1847;// mm/ch
  // CRDCSlope[1][1] = 0.1869;// mm/ch

  // CRDCOffset[0][0] = -279.577;//mm
  // CRDCOffset[1][0] = -277.52;//mm
  // CRDCOffset[0][1] = 127.662; //mm
  // CRDCOffset[1][1] = -134.905; //mm

  ObjCorr[0] = 90; //obj/mrad
  ObjCorr[1] = 0.01; //obj/mm

  //   //CRDC calibration parameters
  // //[CRDC #][X = 0 or Y =1]
  // e16001
  CRDCSlope[0][0][0] = 2.54;// mm/pad
  CRDCSlope[1][0][0] = 2.54;// mm/pad

  CRDCSlope[0][1][0] = -0.151712;// mm/ch
  CRDCSlope[0][1][1] = -0.153206;// mm/ch
  CRDCSlope[0][1][2] = -0.183109;// mm/ch
  CRDCSlope[0][1][3] = -0.18868;// mm/ch
  CRDCSlope[0][1][4] = -0.190429;// mm/ch

  CRDCSlope[1][1][0] = 0.147407;// mm/ch
  CRDCSlope[1][1][1] = 0.141679;// mm/ch
  CRDCSlope[1][1][2] = 0.171754;// mm/ch
  CRDCSlope[1][1][3] = 0.184127;// mm/ch
  CRDCSlope[1][1][4] = 0.187908;// mm/ch


  CRDCOffset[0][0][0] = -281.94;//mm
  CRDCOffset[1][0][0] = -281.94;//mm

  CRDCOffset[0][1][0] = 138.394; //mm
  CRDCOffset[0][1][1] = 144.154; //mm
  CRDCOffset[0][1][2] = 143.1794; //mm
  CRDCOffset[0][1][3] = 137.103; //mm
  CRDCOffset[0][1][4] = 138.426; //mm

  CRDCOffset[1][1][0] = -137.915; //mm
  CRDCOffset[1][1][1] = -138.273; //mm
  CRDCOffset[1][1][2] = -137.843; //mm
  CRDCOffset[1][1][3] = -136.84; //mm
  CRDCOffset[1][1][4] = -136.685; //mm

  gravity_width = 12; //pads

  
}


bool S800::unpack(unsigned short *&point,int runno)
{
  long long int n;
  unsigned short plength, ptag, ID, nwords, words;

  nwords = *point++;
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
  while(nwords > 0){
    plength = *point; ++point; 
    ptag = *point; ++point;
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
      while (words > 0) {
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
      n = *point++;
      n = (*point++<<16|n);
      n = (*point++<<16|n);
      n = (*point++<<16|n);
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
    case S800_VME_TDC_PACKET:
      //cout << "VME_TDC: " << hex << S800_VME_TDC_PACKET << dec << endl;
      point = DecodeS800NewMultiHitTDC(point);
      break;
    default: // S800_II_CRDC_PACKET, S800_II_TRACK_PACKET...
      point += plength - 2;
      break;
    }
    nwords -= plength;
  }
  //this->SetTS(ts);
  //if(ffirst_ts<1){
  //  ffirst_ts = ts;
  //}

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
  while(words > 0){
    ch = ((*p)&0xf000)>>12;
    if (ch ==  8) s800 = (*p++)&0xfff;
    if (ch ==  9) external1 = (*p++)&0xfff;
    if (ch == 10) external2 = (*p++)&0xfff;
    if (ch == 11) secondary = (*p++)&0xfff;
    words--;
  }
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
  }
  else {
    de_down   = (*p++)&0xfff;
    time_down = (*p++)&0xfff;
  }
  //  this->GetScintillator(id)->SetID(id);
  //this->GetScintillator(id)->Set(de_up, time_up, de_down, time_down);

  Scint[id].de_up = de_up;
  Scint[id].time_up = time_up;
  Scint[id].de_down = de_down;
  Scint[id].time_down = time_down; 
  
  
  return p;
}

unsigned short* S800::DecodeS800HodoScope(unsigned short *p){
  UShort_t words = (*(p-2))-2;
  UShort_t id;
  UShort_t ch;
  UShort_t energy;
  int nHodo = 0;
  while (words > 0) {
    id = *p;
    if (id == 0) {
      p++;
      words--;
      while (words > 0) {
	ch = (((*p)&0xF000)>>12);
	energy = ((*p)&0x0FFF);
	//	this->GetHodoscope(ch)->SetEnergy((Int_t)energy);
	Hodo[nHodo].id = ch;
	Hodo[nHodo].energy = (Int_t)energy;
	nHodo++;
	p++;
	words--;
      }
    } else if (id == 1) {
      p++;
      words--;
      while (words > 0) {
	ch = (((*p)&0xF000)>>12) + 16;
	energy = ((*p)&0x0FFF);
	//	this->GetHodoscope(ch)->SetEnergy((Int_t)energy);
	Hodo[nHodo].id = ch;
	Hodo[nHodo].energy = (Int_t)energy;
	nHodo++;
	p++;
	words--;
      }
    } else if (id == 2) {
      p++;
      words--;
      while (words > 0) {
	// coincidence register A (for the first  16 channels)
	p++; words--;
	// coincidence register B (for the second 16 channels)
	p++; words--;
	// TAC time
	p++; words--;
      }
    } else {
      p++; words--;
    }
  }
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

  int r;
  int arraysize = sizeof(CRDC0Runs)/sizeof(CRDC0Runs[0]);
  for(int i = 0; i < arraysize; i++) {
    if(runno < CRDC0Runs[i]) {
      r = i-1;
      break;
    }
  }
//cout << "R == " << r << "\t" << runno << "\t" << CRDC0Runs[r+1] << "\t" << CRDC0Runs[r] << endl;
  int rdiff = CRDC0Runs[r+1] - CRDC0Runs[r];
  double mix1 = (double)(CRDC0Runs[r+1]-runno)/rdiff;
  double mix2 = (double)(runno - CRDC0Runs[r])/rdiff;
  double mix3 = (double)(CRDC1Runs[r+1]-runno)/rdiff;
  double mix4 = (double)(runno - CRDC1Runs[r])/rdiff;
//cout << "Rdiff == " << rdiff << endl;
//cout << "mix == " << mix1  <<"\t" << mix2 << endl;

  //    this->GetCrdc(id)->SetAnodeTAC(anode, tac);1
  CRDC[id].anode = anode;
  CRDC[id].tac = tac;
//  CRDC[id].calY = tac*CRDCSlope[id][1][0] + CRDCOffset[id][1][0];
  if(id == 0) CRDC[id].calY = (tac*CRDCSlope[id][1][r] + CRDCOffset[id][1][r])*mix1 + (tac*CRDCSlope[id][1][r+1] + CRDCOffset[id][1][r+1])*mix2;
  if(id == 1) CRDC[id].calY = (tac*CRDCSlope[id][1][r] + CRDCOffset[id][1][r])*mix3 + (tac*CRDCSlope[id][1][r+1] + CRDCOffset[id][1][r+1])*mix4;
//cout <<  CRDC[id].calY << endl;
  return p;
}

unsigned short* S800::DecodeS800IonChamber(unsigned short *p){
  UShort_t ch=-1;
  UShort_t raw=-1;
  if (*(p+1) == S800_FP_IC_ENERGY_PACKET) {
    // IC packet with times
    UShort_t length;
    length = *p++;
    p++;
    length -= 2;
    while (length > 0) {
      ch  = ((*p)&0xf000)>>12;
      raw = (*p++)&0xfff;
      length--;

      //      this->GetIonChamber()->Set(ch,raw);
      IC.raw[ch] = raw;
      IC.sum += raw;

//       cout << ch << " " << raw << " " << IC.sum <<endl;

      Histo->ICSummary->Fill(ch,raw);
//      cout << ch << " " << raw << " " << IC.sum <<endl;
    }
    IC.sum /= 16.;
  } else {
    // Old style IC packet
    UShort_t words = (*(p-2))-2;
    while(words > 0){
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
  Short_t sampleBegin = 0, sampleWidth=0, isample, ichannel, cdata[4], connector, previous_sample = 0, ch;
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
  p++;	// skip packet id
  UShort_t threshold = *p++;

  //  cout << "Decoding CRDC " << id << " " << i << endl;
  while(i > 0){
    if ((*p)>>15 != 1) {
      std::cout 
	<< "DecodedEvent: " 
	<< "CRDC data is corrupted!" 
	<< std::endl;
      p++; i--;
      continue;
    } else{
      isample  = ((*p)&0x7FC0)>>6;
      ichannel =  (*p)&0x003F;
      if (i == length-3){
	sampleBegin     = isample; 
	previous_sample = isample;
      }
      if(previous_channel > ichannel) sindex++;
      previous_channel = ichannel;
    }
    p++; i--;
    memset(cdata, 0, sizeof(cdata));
    while ((*p)>>15 == 0) {
      connector = ((*p)&0x0C00)>>10;
      cdata[connector] = (*p)&0x3FF;
      p++; i--;
      if (i == 0) break;
    }
    if(isample < sampleBegin || isample > sampleBegin+maxwidth){
      if(debug)
	printf("Warning in Crdc Unpack: inconsistent sample number: %d (first: %d)\n", 
	       isample, sampleBegin);
      mes1 = false;
      //  continue;
    }
    if(isample < previous_sample){
      if(debug)
	printf("Warning in Crdc Unpack: sample number lower than previous: %d (previous: %d)\n", 
	       isample, previous_sample);
      mes2 = false;
      //      continue;
    }
    previous_sample = isample;
    for(int j=0; j<4; j++){
      ch = ichannel+j*64;
      if (cdata[j] != 0 && ch < channels){
	if (cdata[j] < threshold) {
	  if (debug)
	    printf("Warning in Crdc Unpack: data lower than threshold: %d (threshold: %d)\n", cdata[j], threshold);
	  mes3 = false;
	} else {
	  //std::cout << "ch " << ch << " cdata[j]" << cdata[j] << " isample " << isample << std::endl;
	  //	  this->GetCrdc(id)->Set(cdata[j], isample, ch);
	  CRDC[id].cdata[ch][sindex] = cdata[j];
	  
	}
      } 
      else if (cdata[j] != 0 && ch >= channels){
	if (debug) {
	  printf("Warning in Crdc Unpack: channel greater than limit: %d (limit: %d)\n", ch, channels);
	}
	mes4 = false;
      }
    }
    //    sampleWidth = isample - sampleBegin + 1;
    sampleWidth = sindex + 1;
    //    cout << "sample width = " << sampleWidth << " " << isample - sampleBegin+1  << endl;
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

void S800::CRDCIntegrate(int id) //////////////// CRDC slope + offset 
 {

  int PadMult=0;
  for(int i=0;i<NPads;i++)
    {
      float raw = -1;
      float cal = -1;
      int nsamples=0;
      for(int s=0;s<CRDC[id].samplewidth;s++)
	{
//	  cout << "crdc " << id << " pad " << i << " sample "<< CRDC[id].cdata[i][s] << endl;
	  if(CRDC[id].cdata[i][s] !=0)
	    {
	      nsamples++;
	      raw += (CRDC[id].cdata[i][s] - CRDCPeds[id][i]);
	      raw = (raw +1)/nsamples;
	      cal = raw;
//	      cal = raw * Chargeslope[id][i] + Chargeintercept[id][i];
	      //cout << cal << endl;
	      //	      cal = raw;
	      //    cout << raw << " "<< cal << endl;
	    }
	}
      
      if(nsamples>0)
	{
          if(id == 0 && (i==16 || i==1 || i==56)) {
            if(raw < 222) {
  	      CRDC[id].raw[PadMult] = raw;
	      CRDC[id].cal[PadMult] = cal;
	      CRDC[id].pad[PadMult] = i;
	      if(CRDC[id].raw[PadMult] !=-1) PadMult++;
	    } else {
              CRDC[id].raw[PadMult] = -1;
              CRDC[id].cal[PadMult] = -1;
              CRDC[id].pad[PadMult] = i;
              if(CRDC[id].raw[PadMult] !=-1) PadMult++;
	    }
	  } else {
	    if(raw < 900) {
              CRDC[id].raw[PadMult] = raw;
              CRDC[id].cal[PadMult] = cal;
              CRDC[id].pad[PadMult] = i;
              if(CRDC[id].raw[PadMult] !=-1) PadMult++;
	    } else {
              CRDC[id].raw[PadMult] = -1;
              CRDC[id].cal[PadMult] = -1;
              CRDC[id].pad[PadMult] = i;
              if(CRDC[id].raw[PadMult] !=-1) PadMult++;
	    }
	  }
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



  if (raw[15] != 0) {
    for (int i=0; i<13; i++) {
      switch(i) {
      case 0: // e1up
	if (hits[0] >= 0){   
	  //	  fMultiHitTOF.fE1Up.push_back( (data[0][0] - raw[15]) * 0.0625);
	  mTDC.e1up.push_back( (data[0][0] - raw[15]) * 0.0625);
	}
	break;
      case 1: // e1down
	if (hits[1] >= 0){
	  //	  fMultiHitTOF.fE1Down.push_back((data[1][0] - raw[15]) * 0.0625);
	  mTDC.e1down.push_back((data[1][0] - raw[15]) * 0.0625);
	}
	break;
      case 2: // xf
	if (hits[2] >= 0) {
	  for (int j=0; j<=hits[2]; j++){
	    //	    fMultiHitTOF.fXf.push_back((data[2][j] - raw[15]) * 0.0625);
	    mTDC.xfp.push_back((data[2][j] - raw[15]) * 0.0625);
	  }
	}
	break;
      case 3: // obj
	if (hits[3] >= 0) {
	  for (int j=0; j<=hits[3]; j++){
	    //	    fMultiHitTOF.fObj.push_back( (data[3][j] - raw[15]) * 0.0625);
	    mTDC.obj.push_back( (data[3][j] - raw[15]) * 0.0625);
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

  int maxpad = 0;
  int padmax = 0;
  for(int i=0;i<CRDC[id].PadMult-1;i++)
  {
    if(CRDC[id].cal[i] > 0 && CRDC[id].cal[i+1] > 0 && CRDC[id].cal[i] > maxpad)
	  {
  	  maxpad = CRDC[id].cal[i];
  	  padmax = CRDC[id].pad[i];
	  }
  }
  if(CRDC[id].PadMult%2 != 0) padmax = CRDC[id].pad[CRDC[id].PadMult/2];
  else  padmax = CRDC[id].pad[CRDC[id].PadMult/2+1];
  if(padmax ==0)
    {
      CRDC[id].padsum = 0;
      CRDC[id].x_gravity = -500;
      //      cout << "padmax = 0 for crdc " << id << " mult = " << CRDC[id].PadMult <<endl;
      return;
    }

  int lowpad = (int)padmax - (int)gravity_width/2;
  int highpad = lowpad + (int)gravity_width -1;
  if(lowpad <0)
    {
      lowpad =0;
    }
  if(highpad >= 224)
    {
      highpad = 223;
    }
  float sum=0.0, mom =0.0;
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
  CRDC[id].x_gravity = mom/sum * CRDCSlope[id][0][0] + CRDCOffset[id][0][0];

}


bool S800::analyze()
{
  //  cout << "in analyze " << endl;
  //  cout << " size of vector = " << mTDC.e1up.size() << endl;

  for(int i=0;i<(int)mTDC.e1up.size();i++)
    {
      Histo->Te1up->Fill(mTDC.e1up.at(i));
    }
  for(int i=0;i<(int)mTDC.e1down.size();i++)
    {
      Histo->Te1down->Fill(mTDC.e1down.at(i));
    }


  
  
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
  // cout << " x_gravit " << CRDC[0].x_gravity << " " << CRDC[1].x_gravity << endl;
  // cout << " tac = " << CRDC[0].tac << " " << CRDC[1].tac << endl; 

  if(CRDC[0].x_gravity == -500 || CRDC[1].x_gravity ==-500)
  {
    return false;
  }
  Histo->CRDC1X->Fill(CRDC[0].x_gravity);
  Histo->CRDC2X->Fill(CRDC[1].x_gravity);
  Histo->CRDC1XY->Fill(CRDC[0].x_gravity,CRDC[0].tac);
  Histo->CRDC2XY->Fill(CRDC[1].x_gravity,CRDC[1].tac);
  Histo->CRDC1XCalY->Fill(CRDC[0].x_gravity,CRDC[0].calY);
  Histo->CRDC2XCalY->Fill(CRDC[1].x_gravity,CRDC[1].calY);

  if((int)mTDC.obj.size()>0)
  {
    if((int)mTDC.xfp.size()>0)
    {
      bool stat = beam_pid->getPID(mTDC.obj.at(0),mTDC.xfp.at(0));
      if(!stat) return false;
      GetBeamId(beam_pid->Z);
      if(beam_pid->Z == 18 && beam_pid->A == 31)
      {
        Histo->ObjvsICsum_Ar31_uncor->Fill(mTDC.obj.at(0),IC.sum);
      }

      track->CRDC[0] = CRDC[0];
      track->CRDC[1] = CRDC[1];

  	  Histo->ObjvsXFP->Fill(mTDC.obj.at(0),mTDC.xfp.at(0));
  	  //Tracking through the focal plane
  	  track->CalculateTracking(28, 16);
  	  float correctedObj = mTDC.obj.at(0) + ObjCorr[0]*track->afp + ObjCorr[1]*CRDC[0].x_gravity;

  	  mTDC.objCorrected.push_back(correctedObj);

  	  if(beam_pid->Z == 18 && beam_pid->A == 31)
      {
	      Histo->ObjvsICsum_Ar31->Fill(mTDC.objCorrected.at(0),IC.sum);
	    }
	    else if (beam_pid->Z == 16 && beam_pid->A == 29)
	    {
        Histo->ObjvsICsum_S29->Fill(mTDC.objCorrected.at(0),IC.sum);
        //	      cout << mTDC.objCorrected.at(0) << endl;
      }
  	  else if (beam_pid->Z == 15 && beam_pid->A == 28)Histo->ObjvsICsum_P28->Fill(mTDC.objCorrected.at(0),IC.sum);
  	  else if (beam_pid->Z == 14 && beam_pid->A == 27)Histo->ObjvsICsum_Si27->Fill(mTDC.objCorrected.at(0),IC.sum);

  	  //      cout << IC.sum << endl;
  	  Histo->ObjvsICsum->Fill(mTDC.obj.at(0),IC.sum);

  	  //	  if(beam_pid->Z == 18 && beam_pid->A ==31)
  	  if(beam_pid->Z >0 && beam_pid->A >0)
	    {
        bool stat2 = residue_pid[S800Setting][BeamID]->getPID(mTDC.objCorrected.at(0),IC.sum);
	      if(!stat2) return false;

        track->CRDC[0] = CRDC[0];
        track->CRDC[1] = CRDC[1];

	      track->CalculateTracking(residue_pid[S800Setting][BeamID]->A, residue_pid[S800Setting][BeamID]->Z);
	      track->GetThetaPhi(track->ata,track->bta, beam_pid->Z,  beam_pid->A);
		  //Looking at corrections
	      Histo->Objvsafp->Fill(mTDC.objCorrected.at(0),track->afp);
	      Histo->ObjvsCRDC1X->Fill(mTDC.objCorrected.at(0),CRDC[0].x_gravity);
//	      Histo->ThetavsPhi->Fill(track->thetadeg,track->phideg);
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
  cout << __LINE__ << endl;
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
    cout << __LINE__ << endl;
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



S800FpTrack::S800FpTrack(S800Map *map0) {

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
  anglea = 0; /* degree */
  angleb = 0; /* degree */

  anglea_cor = 0; /* degree */
  angleb_cor = 0; /* degree */
  bta_ytacor = 0; /* degree/yta */
  ata_dtacor = 0; /* degree/dta */

  brho = 2.426; /* Tm */
  mass = 17; /* amu */
  charge = 10; /* q */
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
  if (!map->WasLoaded()) { cout << "\nNo Map Loaded" << endl;return; }

  afp = atan((CRDC[1].x_gravity - CRDC[0].x_gravity) / gap);
  bfp = atan((CRDC[1].calY - CRDC[0].calY) / gap);

  xfp = CRDC[0].x_gravity / 1000. + zfp * tan(afp);
  yfp = CRDC[0].calY / 1000. + zfp * tan(bfp);

  if (afp == -65000.) { cout << "\nBAD afp " << (CRDC[1].x_gravity - CRDC[0].x_gravity) / gap << "\t" << CRDC[1].x_gravity << "\t" << CRDC[0].x_gravity << endl; return; }

  input[0] = -xfp;
  input[1] = -afp;
  input[2] = yfp;
  input[3] = bfp;

  if ((double)order > (double)map->maxorder) { order = map->maxorder; }
  ata = map->Calculate((int)order, 0, input);
  yta = map->Calculate((int)order, 1, input);
  bta = map->Calculate((int)order, 2, input);
  dta = map->Calculate((int)order, 3, input);
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
//  cout << "energy = " << energy << " mass = " << mass << " " << charge << endl;
  gamma = 1 + energy / mass / amu;
  beta = sqrt(1 - 1 / gamma / gamma);
  deltabeta = beta;// - beta0;
  // ptot = energy * sqrt(1 + 2 * mass * amu / energy); /* MeV/c */
  ptot = sqrt(energy*energy+2*mass*amu*energy);
  ppar = ptot * cos(scatter/1000.);
  ptra = ptot * sin(scatter/1000.);
}

void S800FpTrack::GetThetaPhi(double ata0, double bta0, int z, int a)
{
  double sineata = sin(ata0);
  double sinebta = sin(bta0);


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
  theta = asin(result);
  phideg = phi*180./acos(-1);
  thetadeg = theta*180./acos(-1);
}


void S800::GetBeamId(int Z)
 {
  // beam numbers
  // 0 - Ar31
  // 1 - S29
  // 2 - P28
  // 3 - Si27
  // 4 - Al26

   if(Z==18)
     BeamID =0;
   else if(Z==16)
     BeamID=1;
   else if(Z==15)
     BeamID=2;
   else if(Z==14)
     BeamID=3;
   else if(Z==13)
     BeamID=4;
   else
     {
       cout << "That beam is not associated with a good resiude PID, Z = " << Z << endl;
       cout << "Need a new pid in the array" << endl;
       abort();
     }
 
 }
