 #include "ceasar.h"

float r2d = 180.0/TMath::Pi();
ceasar::ceasar(TRandom *ran0, histo_sort * Histo0, histo_read * Histo1, float shift0)
{
  ran = ran0;
  Histo = Histo0;
  Histo_read = Histo1;
  chipmap = "ceasar.map";
  calfile = "cal/Caesar_new.cal";
  posmap = "cal/DetPos.txt";
  init(); //shift in mm down beam
}

/***********************************************************************/
/* Destructor **********************************************************/
/***********************************************************************/
ceasar::~ceasar() {
  cout << "start caesar destr" << endl;
  cout << " stop caesar destr" << endl;
}
void ceasar::init()
{
  ifstream ifile(chipmap.c_str());
  if (!ifile.is_open()) {
    cout << "caesar map not found" << endl;
    abort();
  }

  int iring,iloc,iadc,ichan;
  for(int iadc = 0; iadc < 6; iadc++) {
    for(int ich = 0; ich< 32; ich++) {
      MapC[iadc][ich].iRing = -1;
      MapC[iadc][ich].iLoc = -1;
    }
  }

  for (;;) {
    ifile >> iring >> iloc >> iadc >> ichan;
    if (ifile.eof()) break;
    if (ifile.bad()) break;
    MapC[iadc][ichan].iRing = iring;
    MapC[iadc][ichan].iLoc = iloc;
    MapC[iadc][ichan].iTDC = iadc;
    MapC[iadc][ichan].iTDCChan = ichan;
  }
  ifile.close();
  ifile.clear();

  // read in detector angles
  float x,y,z;
  for(int ir = 0; ir < 10; ir++) {
    for(int il = 0; il < 24; il++) {
      angle[ir][il] = -1000.;
      angle2[ir][il] = -1000.;
    }
  }
  ifile.open(posmap.c_str());
  if(!ifile.is_open()) {
    cout << "Couldn't open CAESAR position files" << endl;
    abort();
  }

  for (;;) {
    ifile >> iring >>  iloc >>  x >> y >> z;
    if (ifile.eof()) break;
    if (ifile.bad()) break;
    float r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    float theta = acos(z/r);
    float phi = atan2(y,-x);
    mag[iring-1][iloc-1] = r;
    angle[iring-1][iloc-1] = theta;
    angle2[iring-1][iloc-1] = phi;
  }

  ifile.close();
  ifile.clear();

  //read in calibrations
  int Ntele = 1;
  int Nstrip = 192;

  calCeasar = new calibrate(Ntele,Nstrip,calfile.c_str(),1, 0, false); //no weaving

}
/***********************************************************************/
/* Unpacker for CAESAR *************************************************/
/* DAQ consists of 6 MADCs and 6 MTDCS *********************************/
/***********************************************************************/
bool ceasar::unpack(unsigned short * point, int runno) {
  int NE = 0;
  int NT = 0;
  int NEM[192];
  int NTM[192];

  for(int i = 0; i < 192; i++) {
    NEM[i] = 0;
    NTM[i] = 0;
  }

  for (int iadc  = 0; iadc < 6 ; iadc++) {
    //check for ffff's
    unsigned short f3 = *point;
    unsigned short f4 = *(point+1);
    if (f3 == 0xffff && f4 == 0xffff) {
      //cout << "continue " << hex << f3 << " " << f4 << endl;
      point+=2;
      continue;
    }
    //cout << hex << f3 << " " << f4 << endl;
    Caen[iadc].number = 0;
    point = Caen[iadc].read(point);  // suck out the info in the adc
    for (int i = 0; i < Caen[iadc].number; i++) {
      bool overflow = 0.;
      bool underflow = 0.;
      if (Caen[iadc].underflow[i]) {
        underflow = 1;
      }
      if (Caen[iadc].overflow[i]) {
        overflow = 1;
      }

      int id = Caen[iadc].channel[i] + 32*iadc;
      int ienergy = Caen[iadc].data[i];
      if(overflow) {
//        cout << "overflow " << id << endl;
        ienergy = 5000;
      }
      if(underflow) {
//        cout << "underflow " << id << endl;
        ienergy = 5100;
      }
      float energy = calCeasar->getEnergy(0, id, ienergy + ran->Rndm());
      //cout << dec << "ADC ID " << id << hex << endl;
      //Save first value for multihit
      if (NEM[id] == 0)
      {
        DataEC[id].id = id;
        //if (id == 192) cout << id << endl;
        DataEC[id].ienergy = ienergy + ran->Rndm();
        DataEC[id].energy = energy;
        DataEC[id].iadc = iadc;
        DataEC[id].ichan = Caen[iadc].channel[i];
        int Ring = MapC[iadc][Caen[iadc].channel[i]].iRing;
        int Loc = MapC[iadc][Caen[iadc].channel[i]].iLoc;
        DataEC[id].iRing = Ring;
        DataEC[id].iLoc = Loc;
        DataEC[id].pos.SetMagThetaPhi(mag[Ring][Loc], angle[Ring][Loc], angle2[Ring][Loc]);
        DataEC[id].theta = angle[Ring][Loc];
        DataEC[id].phi = angle2[Ring][Loc];

        //int temp_id = id;
        //float temp_Mag = mag[Ring][Loc];
        //int temp_ring = Ring;
        //float temp_phi = angle2[Ring][Loc];
        //cout << "New i " << temp_id << " " << temp_ring << " " << temp_Mag << " " << temp_phi << endl;

      }
      
      Histo->ECaesar[id]->Fill(energy);
      Histo->energyM->Fill(id);
      Histo->energyiSum->Fill(id, ienergy + ran->Rndm());
      Histo->energySum->Fill(id, energy);
      Histo->energyTot->Fill(energy);
      NE++;
      NEM[id]++;

    }
    //check for ffff's
    unsigned short f1 = *point;
    point++;
    unsigned short f2 = *point;
    point++;
    
    //cout << f1 << " " << f2 << endl;
    if(f1 != 0xffff && f2 != 0xffff) return false;
  }

  for (int itdc = 0; itdc < 6; itdc++) {
    unsigned short f3 = *point;
    unsigned short f4 = *(point+1);
    if (f3 == 0xffff && f4 == 0xffff) {
      point+=2;
      continue;
    }
    Mtdc[itdc].trigger = false;
    point = Mtdc[itdc].read(point);  // suck out the info in the tdc
    if(Mtdc[itdc].number > 100) continue;
    uint64_t ts = (static_cast<uint64_t>(Mtdc[itdc].tstamplow)) | (static_cast<uint64_t>(Mtdc[itdc].tstamphigh) << 16) | (static_cast<uint64_t>(Mtdc[itdc].tstampxtend) << 32);
    for(int j = 0; j < Mtdc[itdc].number; j++) {
      if(Mtdc[itdc].channel[j] == 34) continue; //Extended timestamp
      if(Mtdc[itdc].channel[j] == 32 || Mtdc[itdc].channel[j] == 33) { //Trigger Channel
	      Mtdc[itdc].trigger = true;
	      continue;
      }
      int id = Mtdc[itdc].channel[j] + 32*itdc;
      int itime = Mtdc[itdc].data[j];
      //cout << Mtdc[itdc].channel[j] << endl;
      //cout << "Mtdc number " << Mtdc[itdc].number << " no itdc id " << Mtdc[itdc].channel[j] << " id " << id  << " itdc " << itdc << endl;
      //Save first value for multihit
      if (NTM[id] == 0)
      {
        DataTC[id].id  = id;
        DataTC[id].itime = itime;
        DataTC[id].time = itime * Mtdc[itdc].res / 1000.0; // Convert to ns
        DataTC[id].itdc = itdc;
        DataTC[id].ichan = Mtdc[itdc].channel[j];
      }

      Histo->TCaesar[id]->Fill(itime);
      Histo->timeSum->Fill(id, itime);
      Histo->timeM->Fill(id);

      NT++;
      NTM[id]++;

    }
    //check for ffff's
    unsigned short f1 = *point;
    point++;
    unsigned short f2 = *point;
    point++;
    if(f1 != 0xffff && f2 != 0xffff) return false;
  }

  Histo->adctdcM->Fill(NE, NT);

  if(NE > 192 || NT > 192) return false;
  Nselect = 0;
  int NET = 0;
  for(int i = 0; i < 192; i++) {
    if(NEM[i] > 0) {
      Histo->ETCaesar[i]->Fill(NTM[i], DataEC[i].energy);
      if(NTM[i] == 0) {
	Histo->energyT0Tot->Fill(DataEC[i].energy);
      }
      Histo->mult[i]->Fill(NEM[i], NTM[i]);
    }
    if(NTM[i] > 0) {
      Histo->TECaesar[i]->Fill(NEM[i], DataTC[i].itime);
      if(NEM[i] ==  0) Histo->mult[i]->Fill(NEM[i], NTM[i]); //TODO why == 0???
      if(NEM[i] >  0) {
	      Histo->energyTSum->Fill(DataEC[i].id, DataEC[i].energy);
	      Histo->energyTTot->Fill(DataEC[i].energy);
	      Histo->energyTTotvsT->Fill(DataEC[i].energy,DataTC[i].itime);
	      if (DataTC[i].itime >= 2000 && DataTC[i].itime <= 3700) Histo->energyTTot_tgated->Fill(DataEC[i].energy);
      	if(NEM[i] == 1 && NTM[i] == 1) Histo->energyT1Tot->Fill(DataEC[i].energy);
      	if(NEM[i] == 1 && NTM[i] > 1) Histo->energyT2Tot->Fill(DataEC[i].energy);
        DataEC[i].itime = DataTC[i].itime;
        DataEC[i].time = DataTC[i].time;
	      select[NET] = DataEC[i];
	      NET++;
        Nselect++;
       
      }
    }
  }

  NTSelect = NET;

  Nadded = 0;
  std:vector<bool> addback;
  addback.resize(NET, false);
  vector<int> addback_int;
  addback_int.resize(NET, 0);
  vector<int> jadded;
  for(int i = 0; i < NET; i++) {
    float sum = 0.0;
    sum = select[i].energy; 
    tempAdd = select[i]; //Pass select to temp
    if(addback.at(i)) continue;
    //if(DataEC[i].itime < 2000 || DataEC[i].itime > 3700) continue;
    int temp_id = select[i].id;
    //if (temp_id == 114) cout << "start" << endl;
    //jadded.resize(NET,0);
    for(int j = i+1; j < NET; j++) {
      //skip if added already
      if(addback.at(j)) continue;

      //Add 2 pi to negative values
      float phii = 0;
      float phij = 0;

      if (select[j].pos.Phi() < 0) phij = select[j].pos.Phi() + 2*acos(-1.);
      else phij = select[j].pos.Phi();      
  
      if (select[i].pos.Phi() < 0) phii = select[i].pos.Phi() + 2*acos(-1.);
      else phii = select[i].pos.Phi();

      //Old version. Doesn't properly take negative phi into account
      //Also didn't group the <= 1 and < 0.5 statements together
      //if((abs(select[i].iRing - select[j].iRing) <=1) && (abs(select[i].pos.Phi() - select[j].pos.Phi()) < 0.5)
	      //|| (abs(select[i].iRing - select[j].iRing) == 0 && abs(select[i].pos.Phi() - select[j].pos.Phi()) < 0.7))

      //New version, adds 2pi to negative angles
      if((abs(select[i].iRing - select[j].iRing) <=1 && (abs(phii - phij) < 0.5))
	      || (abs(select[i].iRing - select[j].iRing) == 0 && abs(phii - phij) < 0.7)) {


        //jadded.at(j) = 1;
        //if (temp_id == 114 || temp_id == 116) cout << select[i].id << " " << select[j].id << endl;

	      sum += select[j].energy;
	      if(select[i].energy < select[j].energy) 
        {
          //Don't change the select array, pass new variables to temp array 
          tempAdd = select[j]; //Passes all properties of j to tempADD, energy gets rewritten later by added = sum
          //select[i].id = select[j].id; //Old method
         // addback_int.at(i) = select[j].id;
          //addback_int.at(j) = select[j].id;
          //for(int k = i+1; k < NET; k++) {
            //if (jadded.at(k) == 1) addback_int.at(k) = select[j].id;
          //}

        }
       // else {
         // addback_int.at(j) = select[i].id;
          //addback_int.at(i) = select[i].id;
          //for(int k = i+1; k < NET; k++) {
            //if (added.at(k) == 1) addback_int.at(k) = select[i].id;
          //}
        //}


	      addback.at(j) = true;
      }
    }
    added[Nadded] = tempAdd; //select[i]; //Added is the temp array now
    //Correl Code uses energy in MeV
    added[Nadded].energy = sum/1000.;
    Nadded++;
  }

  //Print out the select gammas and the added gammas
  /*bool in_select116 = false;
  bool in_select114 = false;
  bool in_select115 = false;
  for (int i=0;i<NET;i++) {
    if (select[i].id == 116) in_select116 = true;
    if (select[i].id == 115) in_select115 = true;
    if (select[i].id == 114) in_select114 = true;
  }
  //See what combinations are present
  if (in_select116 == true && in_select114 == true) {
    for (int i=0;i<NET;i++) cout << "Select " << select[i].id << " Add array " << addback_int.at(i) << " ring " << select[i].iRing << " Phi " << select[i].phi << endl;
    for (int i=0;i<Nadded;i++) cout << "Added " << added[i].id << endl;
  }*/


  for(int i = 0; i < Nadded; i++) {
    Histo->enAddback->Fill(added[i].energy * 1000.);
    Histo->enAddbackvsT->Fill(added[i].energy * 1000.,added[i].itime);
    if (added[i].itime >= 2000 && added[i].itime <= 3700) Histo->enAddback_tgated->Fill(added[i].energy * 1000.); //Time gated based on first hit
  }
  return true;
}

/***********************************************************************/
/* Reset ***************************************************************/
/***********************************************************************/
void ceasar::Reset() {
  Nselect = 0;
  Nadded  = 0;
  for(int i = 0; i < 192; i++) {
    DataEC[i].energy = -1;
    DataEC[i].ienergy = -1;
    DataEC[i].itime = -1;
    DataEC[i].time = -1;

    DataTC[i].itime = -1;
    DataTC[i].time = -1;
    DataTC[i].itime = -1;
    DataTC[i].ichan = -1;
  }
}
