 #include "ceasar.h"

ceasar::ceasar(TRandom *ran0, histo_sort * Histo0, histo_read * Histo1, float shift0)
{
  ran = ran0;
  Histo = Histo0;
  Histo_read = Histo1;
  init(shift0); //shift in mm down beam
}


void ceasar::init(float shift)
{
  //  Doppler = new doppler(0.362); // beta for 65 Mev/A C-9
  //  Doppler = new doppler(0.3323); // v/c for Ca-36 = 0.3323
//  Doppler = new doppler(0.326477); // beta for 54 MeV/A Ne-17
  tdc = new TDC1190*[2];
  //tdc[0] = new TDC1190(3,20,128);
  //tdc[1] = new TDC1190(3,1,128);
  tdc[0] = new TDC1190(1,20,128);
  tdc[1] = new TDC1190(1,1,128);

  //make map of chips
  ifstream ifile("ceasar.map");

  if (!ifile.is_open())
  {
    cout << "ceasar map not found" << endl;
    abort();
  }


  //The CAESAR mapping has been modified since e10001 to go off of the actual qdc number (starting from
  //only the CAESAR qdcs (0-5)) rather than bank number
  //the map file has the ring# detector_position# qdc# channel#
  
  //  getline(ifile,name);
  int iring,iloc,iqdc,ichan;
  for(int iqdc = 0;iqdc<6;iqdc++)
  {
    for(int ich = 0;ich<32;ich++)
    {
      MapC[iqdc][ich].iRing = -1;
      MapC[iqdc][ich].iLoc = -1;
    }
  }

  for (;;)
  {
    ifile >> iring >> iloc >> iqdc >> ichan;
    if (ifile.eof()) break;
    if (ifile.bad()) break;
    iloc-=1; //The map starts at 1 rather than 0
    MapC[iqdc][ichan].iRing = iring;
    MapC[iqdc][ichan].iLoc = iloc;
    //      cout << iring << " " << iloc << " " << iqdc << " " << ichan << endl;

  }

  ifile.close();
  ifile.clear();
  //Read in QDC to TDC mapping
  int it,itc;
  ifile.open("CaesarQDCtoTDC.map");
  if(!ifile.is_open())
  {
    cout << "Couldn't open QDC to TDC map" << endl;
    abort();
  }

  for(;;)
  {
    ifile >> it >> itc >> iqdc >> ichan;
    if(ifile.eof()) break;
    if(ifile.bad()) break;

    MapC[iqdc][ichan].iTDC = it;
    MapC[iqdc][ichan].iTDCChan = itc;
    //      if(it !=0)
    //cout << it << " " <<  itc << " " << iqdc << " " << ichan << endl;
  }

  ifile.close();
  ifile.clear();
    
  
  //read in calibrations
  int Ntele = 1;
  int Nstrip = 192;
  string name;

  name = "cal/Caesar_new.cal";
  calCeasar = new calibrate(Ntele,Nstrip,name,1, 0, false); //no weaving
  // name = "cal/ceasar_time.cal";
  // calCeasarT = new calibrate(Ntele,Nstrip,name,1);

  // read in detector angles
  float x,y,z;

  for(int ir=0;ir<10;ir++)
  {
    for(int il = 0;il<24;il++)
    {
      //      cout << ir << " " << il << endl;
      angle[ir][il] = -1000.;
      angle2[ir][il] = -1000.;
    }
  }

  ifile.open("cal/DetPos.txt");
  for (;;)
  {
    ifile >> iring >>  iloc >>  x >> y >> z;
    if (ifile.eof()) break;
    if (ifile.bad()) break;
    
    z = z - shift; // shift in mm

    float r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    float theta = acos(z/r);
    float phi = atan2(y,-x);
    angle[iring-1][iloc-1] = theta; 
    angle2[iring-1][iloc-1] = phi; 
    // cout << iring-1 << " " << " " << iloc-1 << " " << phi*180/acos(-1) << endl;
  }
  //  abort();
  //   for(int i = 0;i<9;i++)
  //     {
  //       cout << abs(angle2[0][i] - angle2[0][i+2]) << endl;
  //     }
  //    abort();
  //  cout << "End init" << endl;
}


//bool ceasar::unpack(unsigned short * point)
bool ceasar::unpack(unsigned short * point, unsigned short *tdcpoint, int runno)
{
  Nxfp = 0;
  NE = 0;
  int Ring = 0;
  int QDC =0;
  int Loc =0;
  int chan = 0;

  NT = 0;
  // S Gillespie - 2020/11/02 
  // Changing point to tdcpoint in TDC unpacker
  for (int itdc = 0;itdc<2;itdc++)
  {
    //the itdc 0/1 are different in the use of tdcpointer vs point, itdc0 also has CsI on it which was read
    //earlier in the data stream
    if(itdc == 0) 
    { 
      //check for ffff's
      unsigned short f3 = *tdcpoint;
      unsigned short f4 = *(tdcpoint+1);
      if (f3 == 0xffff && f4 == 0xffff)
      {
        tdcpoint += 2;
        continue;
      } 
      tdcpoint = tdc[itdc]->read(tdcpoint);
      //     cout << "Tdc" << itdc << endl;
      //     cout << "Ncaesar  " << tdc[itdc]->Ndata << endl;
      for (int i =0;i< tdc[itdc]->Ndata;i++)
      {
        int id = -1;
        //Throwing away all but the CAESAR data
        if(itdc ==0 && tdc[itdc]->dataOut[i].channel <48)
          continue;
        if(itdc==1 && tdc[itdc]->dataOut[i].channel <16)
          continue;

        //mapping into unique ID number
        if(itdc==0)
          id = tdc[itdc]->dataOut[i].channel -48;
        else
          id = tdc[itdc]->dataOut[i].channel+64;

        // skipping channels that have issue with time
        //if (id<32 || (id >47 && id<64) || (id>95 && id<161) || id>175)
        //{
        //  continue;
        //}


        //cout <<"id = " << id << "   " << tdc[itdc]->dataOut[i].time << endl;
        int itime = tdc[itdc]->dataOut[i].time;

        if (tdc[itdc]->dataOut[i].order > 0) continue;

        DataTC[NT].id  = id;
        DataTC[NT].itime = itime;
        DataTC[NT].itdc = itdc;
        DataTC[NT].ichan = tdc[itdc]->dataOut[i].channel;
        
        float time;
        if (runno > 111)
          time = itime/10.;
        else //time was shifted slightly between early and later runs.
          time = itime/10. - 100;

        DataTC[NT].time = time;
        Histo->TCeasarCal1->Fill(time);
        Histo->TCeasarRaw1->Fill(itime/10.);
        Histo->TCeasar[id]->Fill(itime/10.);
        Histo->TCeasarRawSum->Fill(id, itime/10.);
        //cout << "id " << id << endl;
        //      Histo->TCSum->Fill(id,itime/10.);
        NT++; 
      }

      unsigned short f1 = *tdcpoint;
      tdcpoint++;
      unsigned short f2 = *tdcpoint;
      tdcpoint++;
      if(f1 != 0xffff && f2 != 0xffff) return false;
    }

    else if(itdc == 1)
    { 
      //check for ffff's
      unsigned short f3 = *point;
      unsigned short f4 = *(point+1);
      if (f3 == 0xffff && f4 == 0xffff)
      {
        point += 2;
        continue;
      } 
      point = tdc[itdc]->read(point);
      //     cout << "Tdc" << itdc << endl;
      //     cout << "Ncaesar  " << tdc[itdc]->Ndata << endl;
      for (int i =0;i< tdc[itdc]->Ndata;i++)
      {
        int id = -1;
        //Throwing away all but the CAESAR data
        if(itdc ==0 && tdc[itdc]->dataOut[i].channel <48)
          continue;
        if(itdc==1 && tdc[itdc]->dataOut[i].channel <16)
          continue;

        //mapping into unique ID number
        if(itdc==0)
          id = tdc[itdc]->dataOut[i].channel -48;
        else 
          id = tdc[itdc]->dataOut[i].channel+64;

        // skipping channels that have issue with time
        //if (id<32 || (id >47 && id<64) || (id>95 && id<161) || id>175)
        //{
        //  continue;
        //}



        // cout <<"id = " << id << endl;
        int itime = tdc[itdc]->dataOut[i].time;

        if (tdc[itdc]->dataOut[i].order > 0) continue;

        DataTC[NT].id  = id;
        DataTC[NT].itime = itime;
        DataTC[NT].itdc = itdc;
        DataTC[NT].ichan = tdc[itdc]->dataOut[i].channel;

        float time;
        if (runno > 111)
          time = itime/10.;
        else //time was shifted slightly between early and later runs.
          time = itime/10. - 100;

        DataTC[NT].time = time;
        Histo->TCeasarCal2->Fill(time);
        Histo->TCeasarRaw2->Fill(itime/10.);
        Histo->TCeasar[id]->Fill(itime/10.);
        Histo->TCeasarRawSum->Fill(id, itime/10.);
        NT++; 
      }

      unsigned short f1 = *point;
      point++;
      unsigned short f2 = *point;
      point++;
      if(f1 != 0xffff && f2 != 0xffff) return false;
    }
  }
  Histo->CTMult->Fill(NT);

  //These are actually ADCs for 36Ca exp in nov2020
  for (int iqdc = 0;iqdc<6;iqdc++)
  {

    //check for ffff's
    unsigned short f3 = *point;
    unsigned short f4 = *(point+1);
    if (f3 == 0xffff && f4 == 0xffff) 
    {
      point+=2;
      continue; 
    }

    Caen.number = 0;
    point = Caen.read(point);  // suck out the info in the qdc
    for (int i=0;i<Caen.number;i++)
    {
      bool overflow = 0.;
      bool underflow = 0.;


      if (Caen.underflow[i])
        underflow  =1;
      if (Caen.overflow[i])
        overflow =1;
        
      int id = Caen.channel[i] + 32*iqdc;
      int ienergy = Caen.data[i];
      chan = Caen.channel[i]; 
      if(id <0 || id > 192) continue;

      if(overflow)
      {
        //cout << "overflow " << id << endl;
        ienergy = 5000;
      }
      if(underflow)
      {
        //cout << "underflow " << id << endl;
        ienergy = 5100;
      }

      DataEC[NE].id = id;
      DataEC[NE].ienergy = ienergy;
      DataEC[NE].iqdc = iqdc;
      DataEC[NE].ichan = chan;
      
      // cout << ienergy << " "<< Bank << " " << chan << " ";
      Ring = MapC[iqdc][chan].iRing;
      Loc = MapC[iqdc][chan].iLoc;

      DataEC[NE].iRing = Ring;
      DataEC[NE].iLoc = Loc;
      DataEC[NE].theta = angle[Ring][Loc];
      DataEC[NE].phi= angle2[Ring][Loc];
        
      //cout << "ring " << Ring << "  loc " << Loc << "  id " << id << endl;
      if(Ring == -1 || Loc == -1) continue;

      //get calibrated ceasar energies, linear calibration used
      float energy = calCeasar->getEnergy(0,id,ienergy+ran->Rndm());
      DataEC[NE].energy = energy;
      float dop_energy = (1./sqrt(1.-pow(0.379,2)))*energy*(1.- 0.379*cos(DataEC[NE].theta));
      DataEC[NE].dop_energy = dop_energy;

      //Histo->RingSum[Ring]->Fill(Loc,ienergy);
      Histo->ECeasarRawSum->Fill(id,ienergy);
      Histo->ECeasarCalSum->Fill(id,energy);
      Histo->ECeasarDopSum->Fill(id,dop_energy);

      //Histo->RingSum_Cal[Ring]->Fill(Loc,energy);
      Histo->ECeasar[id]->Fill(ienergy);
      //      Histo->ECSum->Fill(id,ienergy);
      Histo->ECCeasar[id]->Fill(energy);
      Histo->TECeasar->Fill(energy);
      Histo->TECeasar_difbin->Fill(energy);
      DataEC[NE].Total += energy;

      NE++;
    }

    //check for ffff's
    unsigned short f1 = *point;
    point++;
    unsigned short f2 = *point;
    point++;
    if(f1 != 0xffff && f2 != 0xffff) return false;

  }
  Histo->CEMult->Fill(NE);

  //if (NE == 1 && NT == 1)
  //{
  //  Histo->ECeasar[DataEC[0].id]->Fill(DataEC[0].ienergy);
  //  Histo->ECeasarRawSum->Fill(DataEC[0].id,DataEC[0].ienergy);
  //}


  //addback on the energies without matching for time info
  N_NoTaddback = 0;
  int NoTaddbackmult = 0;
  std::vector<int> NoTmyvector;
  for(int i =0;i< NE;i++)
  {
    bool addedback = 0;
    float sum = 0.;
    sum = DataEC[i].energy;
    float maxE = sum;
    int maxE_indx = i;
    //check to see if any energies, at index values stored in myvector,
    //were already added back already. Skip these.
    if(i != 0)
    {
      for(int v = 0;v<(int)NoTmyvector.size();v++)
      {
        if(i == NoTmyvector.at(v))
          addedback = 1;
      }
    }

    if(addedback)
      continue;
    else
    {
      for(int j = i+1;j<NE;j++)
      {
        if(abs(DataEC[i].iRing - DataEC[j].iRing) <= 2)
        {
          if(abs(DataEC[i].phi - DataEC[j].phi) < 3.14/2)
          {
            sum += DataEC[j].energy;
            if(maxE < DataEC[j].energy)
            {
              maxE = DataEC[j].energy;
              maxE_indx = j;
            }

            NoTmyvector.push_back(j);
            NoTaddbackmult++;
        
          }
        }
      }
      NoTadded[N_NoTaddback]=DataEC[maxE_indx];
      NoTadded[N_NoTaddback].energy = sum;
      NoTadded[N_NoTaddback].addbackmult = NoTaddbackmult;
      N_NoTaddback++;
    }
  }




  //code to double check that each caesar det is mapped correctly
  int shift[10] = {0,10,24,48,72,96,120,144,168,182};

  if (NT >1)
  {
    for(int i1=0; i1<NT-1; i1++)
    {
      for(int i2=i1+2; i2<NT; i2++)
      {
        Histo->map_DataTC_id_id->Fill(DataTC[i1].id,DataTC[i2].id);
        Histo->map_DataTC_id_id->Fill(DataTC[i2].id,DataTC[i1].id);
      }
    }
  }

  if (NE >1)
  {
    for(int i1=0; i1<NE-1; i1++)
    {
      for(int i2=i1+2; i2<NE; i2++)
      {

        int idet1 = shift[DataEC[i1].iRing] + DataEC[i1].iLoc;
        int idet2 = shift[DataEC[i2].iRing] + DataEC[i2].iLoc;

        Histo->map_iring_iring->Fill(DataEC[i1].iRing,DataEC[i2].iRing);
        Histo->map_iring_iring->Fill(DataEC[i2].iRing,DataEC[i1].iRing);

        Histo->map_DataEC_id_id->Fill(DataEC[i1].id,DataEC[i2].id);
        Histo->map_DataEC_id_id->Fill(DataEC[i2].id,DataEC[i1].id);

        Histo->map_idet_idet->Fill(idet1,idet2);
        Histo->map_idet_idet->Fill(idet2,idet1);

        //double delta_phi = fabs(DataEC[i1].phi - DataEC[i2].phi)*180./acos(-1);
        //Histo->delta_phi_gamma->Fill(delta_phi);

        if (fabs(idet1 - idet2)-25 < 3)
        {
          double delta_phi = fabs(DataEC[i1].phi - DataEC[i2].phi)*180./acos(-1);
          Histo->delta_phi_gamma->Fill(delta_phi); 
        }

      }    
    }
  }

  for(int i1=0; i1<NE; i1++)
  { 
    int idet1 = shift[DataEC[i1].iRing] + DataEC[i1].iLoc;
    
    Histo->ECeasarRawSum_idet->Fill(idet1,DataEC[i1].ienergy);
    Histo->ECeasarCalSum_idet->Fill(idet1,DataEC[i1].energy);
    Histo->ECeasarDopSum_idet->Fill(idet1,DataEC[i1].dop_energy);
  }

  //  cout << "NE = " << NE << " NT = " << NT << endl;
  Nselect = 0;
  // match up energies to times
  for (int ie=0;ie<NE;ie++)
  {
    DataEC[ie].itime = -1;
    DataEC[ie].time = -1;
    for (int it=0;it<NT;it++)
    {
      //cout << "timeID " << DataTC[it].id << " energyID " << DataEC[ie].id << endl;
      //cout << "time Ring " << DataTC[ie].iRing << " time loc " << DataTC[i1].iLoc << endl;
      //cout << "ener Ring " << DataEC[ie].iRing << " ener loc " << DataEC[i1].iLoc << endl;

      if (DataEC[ie].id == DataTC[it].id) //we have matched
      {
        int iQDC = DataEC[ie].iqdc;
        int iChan = DataEC[ie].ichan;

        //cout << "match!!!! " << endl;
          
        if(DataTC[it].itdc != MapC[iQDC][iChan].iTDC ||
           DataTC[it].ichan != MapC[iQDC][iChan].iTDCChan)
        {
          cout << "You've d    int idet1 = shift[DataEC[i1].iRing] + DataEC[i1].iLoc;one a bad job mapping" << endl;
          abort();
        }

        //cout << DataEC[ie].id << endl;
        DataEC[ie].itime = DataTC[it].itime;
        DataEC[ie].time = DataTC[it].time;
            // if (DataEC[ie].itime > 3000 && DataEC[ie].itime < 6000 
        //       && DataEC[ie].iRing != -1)
        if(DataEC[ie].iRing !=-1)
        {
          //               Histo->ECMSum->Fill(DataEC[ie].id,DataEC[ie].ienergy);
          //float dop_energy = 0.;//Doppler->correct(DataEC[ie].energy,DataEC[ie].theta);
          //Histo->TEC_Dop->Fill(dop_energy);
          
          //cout << "DataEC[ie].itime " << DataEC[ie].itime/10 << " DataEC[ie].time " << DataEC[ie].time << endl;
          
          if (1) // (DataEC[ie].itime/10 > -750 && DataEC[ie].itime/10 < 250)
          {

            select[Nselect] = DataEC[ie];
            //select[Nselect].dop_energy = dop_energy;
            Nselect++;
          }
        }
        break;
      }//end if (DataEC[ie].id == DataTC[it].id)
      //else if (DataEC[ie].id < DataTC[it].id) break; // no match found
    } //end loop over NT
  }//end loop over NE


  //Used for Y-88, gate on the 1836keV peak, try to observe the 898keV.
  float gateL = 1.8;
  float gateR = 1.88;
  for (int i=0; i<Nselect;i++)
  {
    if (select[i].energy > gateL && select[i].energy < gateR)
    {
      for (int j=0; j<Nselect;j++)
      {
        if (i!=j)
        {
          Histo->Egated->Fill(select[j].energy);

          float Timedif = select[i].time - select[j].time;
          //cout << added[i].time << "   " << added[j].time << endl;
          Histo->TvsEgated->Fill(select[j].energy, Timedif);
          Histo->TCeasarSum_gated->Fill(select[j].id, select[j].time);

/*
            if (select[j].energy > 0.85 && select[j].energy < 0.95 && Timedif < -200)
            {
              cout << "NE " << NE << endl;
              for (int k=0; k<NE;k++)
              {
                cout << "k " << k << "  id " << DataEC[k].id << "  E=" << DataEC[k].energy << "  iqdc,ichan " << DataEC[k].iqdc << "," << DataEC[k].ichan << endl;
              }
              cout << " NT " << NT << endl;
              for (int k=0; k<NT;k++)
              {
                cout << " k " << k << "  id " << DataTC[k].id << "  T=" << DataTC[k].itime << "  iqdc,ichan " << DataTC[k].itdc << "," << DataTC[k].ichan << endl;
              }
              //abort();

            }
*/

        }
      }  
      continue;
    }
  }

  for (int i=0; i<Nselect;i++)
  {
    int idet = shift[select[i].iRing] + select[i].iLoc;
    Histo->TCeasarMatchedSum_idet->Fill(idet, select[i].itime/10);
    Histo->ECeasarMatchedRawSum_idet->Fill(idet, select[i].ienergy);
    Histo->ECeasarMatchedCalSum_idet->Fill(idet, select[i].energy);
  }


  Nadded = 0;
  int addbackmult = 0;
  std::vector<int> myvector;
  for(int i =0;i< Nselect;i++)
  {
    bool addedback = 0;
    float sum = 0.;
    sum = select[i].energy;
    float maxE = sum;
    int maxE_indx = i;
    //check to see if any energies, at index values stored in myvector,
    //were already added back already. Skip these.
    if(i != 0)
    {
      for(int v = 0;v<(int)myvector.size();v++)
      {
        if(i == myvector.at(v))
          addedback = 1;
      }
    }

    if(addedback)
      continue;
    else
    {
      for(int j = i+1;j<Nselect;j++)
      {
        //cout << select[i].iRing << ", " << select[j].iRing << "   " << select[i].phi << ", " << select[j].phi << endl;
        if(abs(select[i].iRing - select[j].iRing) <= 1)
        {
          if(abs(select[i].phi - select[j].phi) < 0.5)
          {
            sum += select[j].energy;
            if(maxE < select[j].energy)
            {
              maxE = select[j].energy;
              maxE_indx = j;
            }

            myvector.push_back(j);
            addbackmult++;
        
          }
        }
      }
      added[Nadded]=select[maxE_indx];
      added[Nadded].energy = sum;
      added[Nadded].addbackmult = addbackmult;
      Nadded++;
    }
  }

  //cout << "Nselect " << Nselect << " Nadded " << Nadded << endl;
  Histo->CETMult->Fill(Nselect);

  return true;

}
//************************************
ceasar::~ceasar()
{
  cout << "start ceasar destr" << endl;

  delete tdc[0];
  delete tdc[1];
  delete [] tdc;
  //delete Doppler;
  cout << " stop ceasar destr" << endl;
}
//****************************
void ceasar::reset()
{
  Nselect = 0;
}
