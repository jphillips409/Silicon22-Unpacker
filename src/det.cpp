#include "det.h"

#include <fstream>

#include <iostream>





//*********************************************************
/**
 * Constructor
 */
/*
  det::det(histo_sort * Histo0, forest * Forest0)
  {
  type = 0;
  Histo_sort = Histo0;
  Forest = Forest0;
  ran = new TRandom;
  Hira = new hira(ran,Histo_sort);
  s800 = new S800(ran,Histo_sort);
  string d_alpha("d_alpha");
  Ceasar = new ceasar(ran,Histo_sort);
  // Doppler = new doppler(0.326477); // beta for 54 MeV/A Ne-17


  }
*/

det::det(histo_sort * Histo0, histo_read * Histo1, int setting)
{
  Histo_read = Histo1;
  Histo_sort = Histo0;
  ran = new TRandom;
  Gobbi = new gobbi(ran,Histo_sort); //TODO do I need the ran?
  //last value is the s800 setting, 0=Ca36, 1=K35
  Janus = new janus(Histo_sort, 40); //TODO need to give the fibers the correct distance
  s800 = new S800(ran,Histo_sort,setting); 
  Ceasar = new ceasar(ran,Histo_sort,Histo_read, 13);
  Doppler = new doppler(0.3323); //beta for Ca-36

  losses_fiber = new CLosses(16,"_fiber.loss",true);
  losses_target = new CLosses(16,"_Be.loss",true);
  losses_fiberAl = new CLosses(16,"_Al.loss",true); // Losses in aluminum foil between fiber layers

  Nresidue = 0;
  Nbadresidue = 0;
  Ar33_37Cabeam = 0;
  Ar33_36Kbeam = 0;
  //TODO replace ring with gobbi equivalent
  Egamma0 = 0.;
  thetarel0 = 0.;
  thetagamma = 0.;
  thetares = 0.;
  res_Vel0 = 0.;
  detID = -1;
  Ring = -1;
  Loc = -1;
  Edop = 0.;

  //TODO turning off trees for now
  //Histo_sort->tree->Branch("Egamma0", &Egamma0);
  //Histo_sort->tree->Branch("thetarel0", &thetarel0);
  //Histo_sort->tree->Branch("thetagamma", &thetagamma);
  //Histo_sort->tree->Branch("thetares", &thetares);
  //Histo_sort->tree->Branch("res_Vel0", &res_Vel0);
  //Histo_sort->tree->Branch("detID", &detID);
  //Histo_sort->tree->Branch("Ring", &Ring);
  //Histo_sort->tree->Branch("Loc", &Loc);
  //Histo_sort->tree->Branch("Edop", &Edop);

  Nmaxmix = 100;
  NmaxmixSingles = 3000;

  //string name_a8B("Wood_a8B.root");
  //Wood_a8B = new wood(2,false,&name_a8B);
}

//************************************************
/**
 * destructor
 */
det::~det()
{
  for (int i=0; i<min(Nmaxmix,EventMixCounter); i++)
  {
    delete EventMixerP[i];
    delete EventMixer33Ar[i];
  }

  for (int i=0; i<min(NmaxmixSingles,EventSingleMixCounter); i++)
  {
    delete EventMixerSingleP[i];
  }
  delete losses_fiber;
  delete losses_target;
  
  //delete Hira;
  delete Gobbi;
  delete Janus;
  //delete ran;
  delete s800;
  delete Ceasar;
  cout << "You made it!" << endl;
}

void det::Reset()
{
  Gobbi->reset();
  s800->Reset();
}

//*************************************************************
/**
 * unpacks a physics event from the data stream
  point is the ifstream pointer
//TODO Make sure that you can skip detectors. Seems to work
//     It seems like the point will skip a fragment after it's been fed to det unpack
/      This means we can just comment out the unpack line for a detector we don't care about
*/
bool det::unpack(ifstream *point,int runno,int sourceID, int fragmentsize)
{
  //TODO unpack needs to be rewritten so that it accepts an ifstream pointer
  //Janus must be passed the ifstream pointer, Si and S800 should have their buffers read out and passed
  int buffersize = fragmentsize - 28; //should always be 28
  int bufferwords = 1;
  if (sourceID == SiID || sourceID == S800ID) bufferwords = buffersize/2;
  unsigned short buffer_arr[bufferwords];
  unsigned short *pointbuf; //pointer to the buffer, needed for Si, ADC, TDC, S800, CAESAR
  //cout << hex << fragmentsize << " " << buffersize << " " << bufferwords << dec << endl;
  //Read out buffer if silicon or S800
  if (sourceID == SiID || sourceID == S800ID)
  {
    point->read((char*)buffer_arr,buffersize);
    pointbuf = buffer_arr;
  }
  //TODO don't know if words does anything
  unsigned short  words;
  // S Gillespie 2020/11/2
  // Defining a new value to store tdc position to pass to CEASAR
  //TODO tdcpoint may be needed? Unclear
  unsigned short *tdcpoint;

  bool stat = false;

  //cout << "sourceID = " << sourceID << endl;
  if(sourceID == SiID)
  {
    pointbuf++; //skips ADC size, don't actually need the byte size.
    pointbuf += 6; //This skips the SIS3820 information which is already in the timestamp

    //To skip Gobbi, comment out the unpacker. The buffer has been read out so you've already skipped that info
    stat = Gobbi->unpack(pointbuf, tdcpoint, runno); // TODO make a gobbi unpack subroutine

    if(!stat)
    {
      //  cout << "Didn't read hira right" <<endl;
      return stat;
    }
    //To skip CAESAR, comment out the unpacker. Same as Gobbi
    //CAESAR should not have its own source ID
    //stat = Ceasar->unpack(pointbuf, tdcpoint, runno); // SG 2020/11/02
  }
  else if (sourceID == JanusID)
  { 
    //call janus unpacker here
    stat = Janus->unpack(point);//something something janus unpacker

    //Turn skip on if you want to not unpack Janus
    //point->ignore(buffersize);
    stat = true;
  }
  else if (sourceID == S800ID)
  {
    //To skip the S800, comment out the unpacker. Same as Gobbi
    stat = s800->unpack(pointbuf,runno);
    //if(stat)
      //NS800++;
  }
  else
  {
    cout << "found unexpected sourceID = " << sourceID << endl;
    return stat;
  }
  
  return stat;
}
//*********************************

void det::analyze(int event, int run)
{
  Correl.reset();
  bool foundresidue = false;

  //Events are already matched within boards. Call analyze() instead of MatchEvents()
  Janus->analyze();

  //TODO make sure to save janus info where we need it before clearing, clear() should really happen at the end
  Janus->janusevts.clear();
	Janus->Histo->clear();
  Janus->nevts++;
  

  //return; //TODO temp so we don't go far

  s800_results S800_results;
  S800_results.Reset();
  // S800_results.trig_coin=false;
  // S800_results.trig_singles=false;
  // S800_results.trig_s800_singles=false;
  
  if (s800->Trig.registr & 1) S800_results.trig_s800_singles =true;
  if (s800->Trig.registr & 2) S800_results.trig_coin = true;
  if (s800->Trig.registr & 16) S800_results.trig_singles = true;
  
  //ND adding in a counter for these events
  if (s800->Trig.registr & 1) N_s800_singles++;
  if (s800->Trig.registr & 2) N_coin++;
  if (s800->Trig.registr & 16) N_singles++;

  S800_results.Zbeam = -1;
  S800_results.Abeam = -1;
  S800_results.Zresidue = -1;
  S800_results.Aresidue = -1;

  //TODO keep this for runs directly into the beam?
  /*
  //ND, looking at run#29 with Ca37 beam directly into S800
  //Also needed for run#30-32 with S800_run30.cpp
  //This is mostly for mask runs and unreacted beam
  if(S800_results.trig_s800_singles) 
  {
    s800->track->runno = run;

    bool goods800 = s800->analyze();
    if(!goods800){ return; }
    if(s800->BeamID <0){ return; }
    S800_results.Zbeam = s800->beam_pid->Z;
    S800_results.Abeam = s800->beam_pid->A;
    S800_results.Zresidue = s800->residue_pid[s800->S800Setting][s800->BeamID]->Z;
    S800_results.Aresidue = s800->residue_pid[s800->S800Setting][s800->BeamID]->A;

    float ekin = s800->track->energy;

    //fiber array thickness
    float thickness = 51.1*1.30; // NB: 0.5mm -> 51.1 mg/cm2
    float ekinstep1 = losses_target->getEin(ekin,thickness,20,37);
    //target thickness
    thickness = Hira->RingCounter->TargetThickness/2*1.8;
    float ekinstep2 = losses_target->getEin(ekinstep1,thickness,20,37);
    float ekintarget = losses_target->getEin(ekin,thickness,20,37);
    
    
    //cout << "ekin " << ekin << endl;
    float mass = (37. - 13.1361/m0)*m0;
    float etot = ekin + mass;
    float pc = sqrt(pow(ekin+mass,2) - pow(mass,2));

    float velocity = pc/etot;
    //cout << "vel " << velocity << "     vel(cm/ns)" << velocity*30 << endl;

    //Histo_read->Vlab_HF_p35K->Fill(velocity);

    if(S800_results.Zbeam == 20 && S800_results.Abeam == 37)
    {
      Histo_sort->Vlab_Ca37->Fill(velocity);
      Histo_sort->Elab_Ca37->Fill(ekin);
      Histo_sort->Etarfib_Ca37->Fill(ekinstep2);
      Histo_sort->Etar_Ca37->Fill(ekintarget);
      Histo_sort->Efib_Ca37->Fill(ekinstep1);

    }

    return;
  }
  else { return;}
  */
  // s800->analyze();
  bool goods800 = s800->analyze();
  // if(!goods800)
  //   {
  //     return;
  //   }
  S800_results.Zbeam = s800->beam_pid->Z;
  S800_results.Abeam = s800->beam_pid->A;

  int release;

  //this method is doing a lot here.
  //  1. Applies S800 time gate to Si and CsI
  //  2. Adds neighboring Si strips within time gate
  //  3. matches up either E-CsI or dE-E events
  //  4. plots dE-E plots, and determines PID
  //  5. Does energy corrections for each telescope
  release = Gobbi->analyze(S800_results);  // TODO Gobbi needs to accept S800 results to apply time gates


  //TODO rewrite this for Gobbi
  /*if (Hira->RingCounter->proton_present)
  {
    for (int i=0;i<Hira->RingCounter->Nsolution; i++)
    {
      if (Hira->RingCounter->Solution[i].ipid == 1)
      {
        //cout << "denergy " << Hira->RingCounter->Solution[i].denergy << endl; 
        //cout << "energy " << Hira->RingCounter->Solution[i].energy << endl;
        //cout << "energyTot " << Hira->RingCounter->Solution[i].energyTot << endl;

        //cout << "Ekin " << Hira->RingCounter->Solution[i].Ekin << endl;
        //cout << "icsi " << Hira->RingCounter->Solution[i].icsi << endl;
        
        if (Hira->RingCounter->Solution[i].icsi>=0)
        {
          Histo_sort->p_Etot_theta[Hira->RingCounter->Solution[i].icsi]->Fill(Hira->RingCounter->Solution[i].theta, Hira->RingCounter->Solution[i].Ekin);

        Histo_sort->ProtonsAllEnergy->Fill(Hira->RingCounter->Solution[i].energyTot);
        Histo_sort->crossjointprotons->Fill(Hira->RingCounter->Solution[i].iRing);
        if (Hira->RingCounter->Solution[i].energyTot > 990 && Hira->RingCounter->Solution[i].energyTot < 1025)
        {
          Histo_sort->crossjointprotons_gated->Fill(Hira->RingCounter->Solution[i].iRing);
        }
      }
    }
  }*/

  
  if(s800->BeamID <0) return; //Need a good beam id
  

  S800_results.Zresidue = s800->residue_pid[s800->S800Setting][s800->BeamID]->Z;
  S800_results.Aresidue = s800->residue_pid[s800->S800Setting][s800->BeamID]->A;

  //cout << "s800->S800Setting " << s800->S800Setting << "  beamID " << s800->BeamID << endl;
  //cout << "S800_results.Zresidue " << S800_results.Zresidue << "  S800_results.Aresidue " << S800_results.Aresidue << endl;

  //ND, investigative counting
  //TODO rewrite this for Si22 experiment
  if (S800_results.Zresidue==18 && S800_results.Aresidue==33)
  {
    if (S800_results.Zbeam == 20 && S800_results.Abeam==37)
    {
      Ar33_37Cabeam++;
    }
    if (S800_results.Zbeam == 19 && S800_results.Abeam==36)
    {
      Ar33_36Kbeam++;
    }
  }



  //TODO Not sure what this is
/*
  if (S800_results.Zbeam == 20 && S800_results.Abeam==37)
  {
    NCa37beam++;
    if (S800_results.Zresidue==19 && S800_results.Aresidue==35)
    {
      beam_residue++;
      //want fiber data for good tracking
      if (Hira->XY_mon->has_data)
        beam_residue_fiber++;
      else
        beam_residue_nofiber++;

      //want protons in coincidence
      if (Hira->RingCounter->Nsolution > 0)
      {
        beam_residue_RCsoln++;
        if (Hira->XY_mon->has_data)
          beam_residue_RCsoln_fiber++; //everything needed for p+35K
        else
          beam_residue_RCsoln_nofiber++;
      }
      else
      {
        beam_residue_noRCsoln++;
      }
    }
  }
*/
  //TODO going to want gated position maps for Si22 exp fragments
  //fiber position maps for only K35 fragments
  if (S800_results.Zresidue==19 && S800_results.Aresidue==35)
  {
    //NK35residue++;t
    //if(Hira->XY_mon->has_data)
    //{
      //NK35_withfiber++;
      //Histo_sort->VertHitMap_35Kgate->Fill(Hira->XY_mon->vert->pmtx,Hira->XY_mon->vert->pmty);
      //Histo_sort->HorzHitMap_35Kgate->Fill(Hira->XY_mon->horz->pmtx,Hira->XY_mon->horz->pmty);
    //}
  }

  //if (S800_results.trig_coin) Histo_sort->S800_Csi_time->Fill(Hira->T_CsiTrig/10.);
  //if (S800_results.trig_singles) Histo_sort->singles_trig_time->Fill(Hira->T_CsiTrig/10.);

  //  cout << "before fiber" << endl;
  //demand the fibre array gave information of both x and y position
  //if(!Hira->XY_mon->has_data)return;
  //cout << "after fiber" << endl;

  //cout << S800_results.Zbeam << " " << S800_results.Abeam << endl;
 
  //if (S800_results.trig_coin && Hira->RingCounter->proton_present)
    //Histo_sort->S800_Csi_time_with_proton->Fill(Hira->T_CsiTrig/10.);


  //if (Hira->RingCounter->multAlpha == 1)
  //{
  //  if((int)s800->mTDC.objCorrected.size()>0)
  //  {
  //    if((int)s800->mTDC.xfp.size()>0)
  //    {
  //      Histo_sort->ObjvsXFPwithAlpha1->Fill(s800->mTDC.objCorrected.at(0),
  //                                           s800->mTDC.xfp.at(0));
  //    }
  //  }
  //}

  //TODO This seems to apply a proton gate to several hists, rewrite for Si22
  //big proton gate
  /*
  if (Hira->RingCounter->proton_present)
  {
    if ((int)s800->mTDC.objCorrected.size()>0)
    {
      if ((int)s800->mTDC.xfp.size()>0)
      {
        Histo_sort->ObjvsXFPwithProton1->Fill(s800->mTDC.objCorrected.at(0),
                                              s800->mTDC.xfp.at(0));
        Histo_sort->ObjvsXFPwithProton1->Fill(s800->mTDC.objCorrected.at(0),
                                              s800->mTDC.xfp.at(0));
      }
    
      if (S800_results.Zbeam == 20 && S800_results.Abeam==37)
      { 
        Histo_sort->ObjvsICsum_wProt->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);

        if (Hira->XY_mon->has_data)
        {
          Histo_sort->ObjvsICsum_wProt_wFiber->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);
        }
      }
      if (S800_results.Zbeam == 19 && S800_results.Abeam==36)
      { 
        Histo_sort->ObjvsICsum_K36_wProt->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);
      }
      if (S800_results.Zbeam == 18 && S800_results.Abeam==35)
      { 
        Histo_sort->ObjvsICsum_Ar35_wProt->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);
      }
    }
  }
  if (Hira->RingCounter->multProton == 2)
  {
    if((int)s800->mTDC.objCorrected.size()>0)
    {
      if((int)s800->mTDC.xfp.size()>0)
      {     
        Histo_sort->ObjvsXFPwithProton2->Fill(s800->mTDC.objCorrected.at(0),
                                              s800->mTDC.xfp.at(0));
      }
      Histo_sort->ObjvsICsum_wProt2->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);
    }       
  }*/

  //S800 beam properties
  int Zb = S800_results.Zbeam;
  int Ab = S800_results.Abeam;
  int Z = S800_results.Zresidue;
  int A = S800_results.Aresidue;
  
  //TODO don't need this for code to work, but could be useful to rewrite
  if (Zb==20 && Ab==37 && Z==20 && A==37)
  { 
    Histo_sort->BeamCa37_energy->Fill(s800->track->energy); //Energy at target, in MeV
    Histo_sort->BeamCa37_ptot->Fill(s800->track->ptot);     //Total momentum, in MeV/c
    Histo_sort->BeamCa37_ppar->Fill(s800->track->ppar);     //Parallel momentum, in MeV/c
    Histo_sort->BeamCa37_ptra->Fill(s800->track->ptra);     //Transverse momentum, in MeV/c
    Histo_sort->BeamCa37_TofAngle->Fill(s800->mTDC.objCorrected.at(0), s800->track->afp);
    Histo_sort->BeamCa37_ObjvsCRDC1X->Fill(s800->mTDC.objCorrected.at(0),s800->CRDC[0].x_gravity);

    for(int p=0;p<s800->CRDC[0].PadMult;p++)
    {
      Histo_sort->BeamCa37_TofCRDC1raw->Fill(s800->mTDC.objCorrected.at(0), s800->CRDC[0].raw[p]);
      Histo_sort->BeamCa37_TofCRDC1cal->Fill(s800->mTDC.objCorrected.at(0), s800->CRDC[0].cal[p]);
    }

    //if (Hira->XY_mon->has_data)
    //{
     // Histo_sort->BeamCa37_Fiber_XY->Fill(Hira->XY_mon->x,Hira->XY_mon->y);
    //}
    



/////////////////////////////////////////////////////////////////////

  }
  if (Zb==20 && Ab==37 && Z==20 && A==38)
  { 
    Histo_sort->BeamCa38_energy->Fill(s800->track->energy); //Energy at target, in MeV
    Histo_sort->BeamCa38_ptot->Fill(s800->track->ptot);     //Total momentum, in MeV/c
    Histo_sort->BeamCa38_ppar->Fill(s800->track->ppar);     //Parallel momentum, in MeV/c
    Histo_sort->BeamCa38_ptra->Fill(s800->track->ptra);     //Transverse momentum, in MeV/c
    Histo_sort->BeamCa38_TofAngle->Fill(s800->mTDC.objCorrected.at(0), s800->track->afp);

    Histo_sort->BeamCa38_ObjvsCRDC1X->Fill(s800->mTDC.objCorrected.at(0),s800->CRDC[0].x_gravity);


    for(int p=0;p<s800->CRDC[0].PadMult;p++)
    {
      Histo_sort->BeamCa38_TofCRDC1raw->Fill(s800->mTDC.objCorrected.at(0), s800->CRDC[0].raw[p]);
      Histo_sort->BeamCa38_TofCRDC1cal->Fill(s800->mTDC.objCorrected.at(0), s800->CRDC[0].cal[p]);

    }
    //if (Hira->XY_mon->has_data)
    //{
      //Histo_sort->BeamCa38_Fiber_XY->Fill(Hira->XY_mon->x,Hira->XY_mon->y);
    //}

/*
    cntCa37++;
    cout << "Ca 38 event #" << cntCa37 << endl;
    for(int p=0;p<s800->CRDC[0].PadMult;p++)
    {
      cout << s800->CRDC[0].pad[p] << " " << s800->CRDC[0].raw[p] << " " << s800->CRDC[0].cal[p] << endl;
    }
 
    cout << "x_gravity " << s800->CRDC[0].x_gravity << endl;

    if (s800->CRDC[0].x_gravity > -100 && s800->CRDC[0].x_gravity < -99)
      abort();
*/




  }
  if (Zb==19 && Ab==36 && Z==19 && A==36)
  { 
    Histo_sort->BeamK36_energy->Fill(s800->track->energy); //Energy at target, in MeV
    Histo_sort->BeamK36_ptot->Fill(s800->track->ptot);     //Total momentum, in MeV/c
    Histo_sort->BeamK36_ppar->Fill(s800->track->ppar);     //Parallel momentum, in MeV/c
    Histo_sort->BeamK36_ptra->Fill(s800->track->ptra);     //Transverse momentum, in MeV/c
  }
  if (Zb==18 && Ab==35 && Z==18 && A==35)
  { 
    Histo_sort->BeamAr35_energy->Fill(s800->track->energy); //Energy at target, in MeV
    Histo_sort->BeamAr35_ptot->Fill(s800->track->ptot);     //Total momentum, in MeV/c
    Histo_sort->BeamAr35_ppar->Fill(s800->track->ppar);     //Parallel momentum, in MeV/c
    Histo_sort->BeamAr35_ptra->Fill(s800->track->ptra);     //Transverse momentum, in MeV/c
  }


  //relative angles for multProton == 1
  //if (Hira->RingCounter->multProton == 1 && Hira->XY_mon->has_data &&
  //    S800_results.trig_coin)
  //{
  //  float xp,yp,zp;
  //  bool good = false;
  //  for (int i=0;i<Hira->RingCounter->Nsolution;i++)
  //  {
  //    if (Hira->RingCounter->Solution[i].ipid == 1)
  //    {
  //      zp = cos(Hira->RingCounter->Solution[i].theta);
  //      xp = sin(Hira->RingCounter->Solution[i].theta)*
  //                cos(Hira->RingCounter->Solution[i].phi);
  //      yp = sin (Hira->RingCounter->Solution[i].theta)*
  //                sin(Hira->RingCounter->Solution[i].phi);
  //      good = true;
  //      cout << "I exit here???" << endl;
  //      exit; //WHAT IS HAPPENING WITH THIS EXIT??????? -ND 9/12/2023
  //    }
  //  }
  //  if (good)
  //  {
  //    float xr,yr,zr;
  //    zr = cos(Hira->XY_mon->theta);
  //    xr = sin(Hira->XY_mon->theta)*cos(Hira->XY_mon->phi);
  //    yr = sin(Hira->XY_mon->theta)*sin(Hira->XY_mon->phi);

  //    float dot = xp*xr + yp*yr + zp*zr;
  //    float deltaTheta = acos(dot)*180./acos(-1);

      //if (S800_results.Zbeam == 16 && S800_results.Abeam == 29
      //    && S800_results.Zresidue == 15 && S800_results.Aresidue == 27)
      //{
      //  Histo_sort->Si28_deltaTheta->Fill(deltaTheta);
      //}
  //  }
  //}

  //TODO a lot of the following could be necessary, don't fully understand it yet.
  if(S800_results.Zresidue == 19 && S800_results.Aresidue == 35)
  {
    //Histo_sort->CornerMult->Fill(Hira->XY_mon->vert->Mult + Hira->XY_mon->horz->Mult);
    //Histo_sort->CornerMult_Vert->Fill(Hira->XY_mon->vert->Mult);
    //Histo_sort->CornerMult_Hori->Fill(Hira->XY_mon->horz->Mult);
    //Histo_sort->CornerMult_VvsH->Fill(Hira->XY_mon->horz->Mult,Hira->XY_mon->vert->Mult);
    //     cout << "hor = "<< Hira->XY_mon->horz->Mult << " vert = " << Hira->XY_mon->vert->Mult << endl;
  }

  ///////////////////////////////////////////// CRDC gain analysis
  int pad_max[2] = {0,0}; //max energy for a given pad
  int max_pad[2] = {-1,-1}; //pad with the max energy
  
  for(int i=0;i<2;i++)
  {
    for(int j=0;j< s800->CRDC[i].PadMult;j++)
    {
      //          cout << i << " " << j << "  " << s800->CRDC[i].raw[j] << endl;
      if(s800->CRDC[i].raw[j]>pad_max[i])
      {
        pad_max[i] = s800->CRDC[i].cal[j];
        max_pad[i] = s800->CRDC[i].pad[j];
      }
    }
  }
  if (s800->S800Setting == 0)
  {
    //Ca36
    if (S800_results.Zresidue==20 && S800_results.Aresidue==36)
    {          //Setting_residue
	    Histo_sort->CRDC1_Ca36_Ca36->Fill(max_pad[0],pad_max[0]);
	    Histo_sort->CRDC2_Ca36_Ca36->Fill(max_pad[1],pad_max[1]);
    }
    //K35
    if (S800_results.Zresidue==19 && S800_results.Aresidue==35)
    {          //Setting_residue
      Histo_sort->CRDC1_Ca36_K35->Fill(max_pad[0],pad_max[0]);
	    Histo_sort->CRDC2_Ca36_K35->Fill(max_pad[1],pad_max[1]);
    }
    //Ar33
    if (S800_results.Zresidue==18 && S800_results.Aresidue==33)
    {
      Histo_sort->CRDC1_Ca36_Ar33->Fill(max_pad[0],pad_max[0]);
      Histo_sort->CRDC2_Ca36_Ar33->Fill(max_pad[1],pad_max[1]);
    }
    //Ar34
    if (S800_results.Zresidue==18 && S800_results.Aresidue==34)
    {
      Histo_sort->CRDC1_Ca36_Ar34->Fill(max_pad[0],pad_max[0]);
	    Histo_sort->CRDC2_Ca36_Ar34->Fill(max_pad[1],pad_max[1]);
    }
    //Cl32
    if (S800_results.Zresidue==17 && S800_results.Aresidue==32)
    {
      Histo_sort->CRDC1_Ca36_Cl32->Fill(max_pad[0],pad_max[0]);
	    Histo_sort->CRDC2_Ca36_Cl32->Fill(max_pad[1],pad_max[1]);
    }
  }
  if (s800->S800Setting == 1)
  {
    //K35
    if (S800_results.Zresidue==19 && S800_results.Aresidue==35)
    {          //Setting_residue
      Histo_sort->CRDC1_K35_K35->Fill(max_pad[0],pad_max[0]);
	    Histo_sort->CRDC2_K35_K35->Fill(max_pad[1],pad_max[1]);

      //if (max_pad[1] < 20)
      //  abort();

    }
    //K36
    if (S800_results.Zresidue==19 && S800_results.Aresidue==36)
    {
      Histo_sort->CRDC1_K35_K36->Fill(max_pad[0],pad_max[0]);
	    Histo_sort->CRDC2_K35_K36->Fill(max_pad[1],pad_max[1]);
    }
    //Ar34
    if (S800_results.Zresidue==18 && S800_results.Aresidue==34)
    {                   //Setting_residue
	    Histo_sort->CRDC1_K35_Ar34->Fill(max_pad[0],pad_max[0]);
	    Histo_sort->CRDC2_K35_Ar34->Fill(max_pad[1],pad_max[1]);
    }
    //Ar35
    if (S800_results.Zresidue==18 && S800_results.Aresidue==35)
    {
      Histo_sort->CRDC1_K35_Ar35->Fill(max_pad[0],pad_max[0]);
	    Histo_sort->CRDC2_K35_Ar35->Fill(max_pad[1],pad_max[1]);
    }
  }
  ///////////////////////////////////////////// end CRDC gain analysis

  float Nsol = 0;
  int mult = 0;
  CsImult = 0;
  
  //TODO this looks like it's filling some fiber hists
  //TODO important section. Look over the number of Gobbi solutions and then create a new one for heavy S800 frag
  //Fibers will not be contained inside Hira for Si22. Whole new class for Janus
  //  if(S800_results.Zbeam ==18 && S800_results.Abeam ==31)
  if (S800_results.Zbeam >0 && S800_results.Abeam >0)
  {
    if (S800_results.Zresidue==20 && S800_results.Aresidue ==37)
    {
      //Histo_sort->Fiber_XBeam->Fill(Hira->XY_mon->x);
      //Histo_sort->Fiber_YBeam->Fill(Hira->XY_mon->y);
      //Histo_sort->Fiber_XYBeam->Fill(Hira->XY_mon->x,Hira->XY_mon->y);
    }
    if (S800_results.Zresidue==19 && S800_results.Aresidue ==35)
    {
      //Histo_sort->ThetavsPhi_fiber->Fill(Hira->XY_mon->thetadeg,Hira->XY_mon->phideg);
      //Histo_sort->ThetaFibervsThetaS800->Fill(Hira->XY_mon->thetadeg,s800->track->thetadeg);
      //Histo_sort->PhiFibervsPhiS800->Fill(Hira->XY_mon->phideg,s800->track->phideg);
    }
    //cout << "looking for residue" << endl;
    //Lets add the S800 to the Solution Class
    //TODO This is where I'll need to change stuff for Gobbi, maybe create 5th "Telescope"?
    if (S800_results.Zresidue >0 && S800_results.Aresidue >0)
    {
      foundresidue = true;
      foundresidue = LoadS800toSolution();


      if (foundresidue)
      {
        if (S800_results.Zresidue==20 && S800_results.Aresidue==37)
        {
          //if (Hira->RingCounter->multProton == 1)
          //{
            //Nresidue++;
          //}
        }
      }
      else
      {
        if (S800_results.Zresidue==20 && S800_results.Aresidue==37)
        {
          //if (Hira->RingCounter->multProton == 1)
          //{
            //Nbadresidue++;
          //}
        }
      }
    }
  }
  

  //Replace Hira with Gobbi and RingCounter with Telescope
  //IMPORTANT DIFFERENCE - each telescope has its own solution, so make a 5th "telescope" for S800 and add i=0 sol?
  //TODO will not compile while telescope does not have an index
  if(foundresidue)
  {
    int Nsol = Gobbi->Telescope[0]->Nsolution;

    for (int n=0;n<Gobbi->Telescope[0]->CsI.Nstore;n++)
    {
      //cout << "n " << n << endl;
      //cout << " energy " << Gobbi->RingCounter->Csi.Order[n].energyR << "    time " <<  Hira->RingCounter->Csi.Order[n].time << endl;
      //Histo_sort->ET_csi_res->Fill(Gobbi->Telescope->Csi.Order[n].energyR,Hira->RingCounter->Csi.Order[n].time/10.);
    }


    float res_Vel = Gobbi->Telescope[0]->Solution[Nsol-1].velocity;
    float phiR = Gobbi->Telescope[0]->Solution[Nsol-1].phi;//s800->track->phi;//Hira->RingCounter->Solution[Nsol-1].phi;
    float thetaR = Gobbi->Telescope[0]->Solution[Nsol-1].theta;//s800->track->theta;//Hira->RingCounter->Solution[Nsol-1].theta;
    float xR = sin(thetaR)*cos(phiR);
    float yR = sin(thetaR)*sin(phiR);
    float zR = cos(thetaR);

    // TODO would be good to track Si22 heavy frag velocities as well
    //track velocity distributions of all heavy fragments
    //velocities should be similar to beam velocity
    /*if(S800_results.Zbeam == 20 && S800_results.Abeam == 37)
    {
      if(S800_results.Zresidue == 20 && S800_results.Aresidue == 38)
      {
        Histo_sort->Vlab_Ca38->Fill(res_Vel);
        Histo_sort->Ppar_Ca38->Fill(Hira->RingCounter->Solution[Nsol-1].Mvect[2]); 
        float Ptra = sqrt(pow(Hira->RingCounter->Solution[Nsol-1].Mvect[0],2) + 
                          pow(Hira->RingCounter->Solution[Nsol-1].Mvect[1],2) );
        Histo_sort->Ptra_Ca38->Fill(Ptra);
      }
      if(S800_results.Zresidue == 20 && S800_results.Aresidue == 37)
      {
        Histo_sort->Vlab_Ca37->Fill(res_Vel); 
        Histo_sort->Ppar_Ca37->Fill(Hira->RingCounter->Solution[Nsol-1].Mvect[2]); 
        float Ptra = sqrt(pow(Hira->RingCounter->Solution[Nsol-1].Mvect[0],2) + 
                          pow(Hira->RingCounter->Solution[Nsol-1].Mvect[1],2) );
        Histo_sort->Ptra_Ca37->Fill(Ptra);
    //cout << "37Ca Ppar " << Hira->RingCounter->Solution[Nsol-1].Mvect[2] << " Ptra " << Ptra << endl;

      }
      if(S800_results.Zresidue == 20 && S800_results.Aresidue == 35)
        Histo_sort->Vlab_Ca35->Fill(res_Vel);
      if(S800_results.Zresidue == 20 && S800_results.Aresidue == 36)
        Histo_sort->Vlab_Ca36->Fill(res_Vel);
      if(S800_results.Zresidue == 19 && S800_results.Aresidue == 35)
        Histo_sort->Vlab_K35->Fill(res_Vel);
      if(S800_results.Zresidue == 19 && S800_results.Aresidue == 36)
        Histo_sort->Vlab_K36->Fill(res_Vel);
      if(S800_results.Zresidue == 18 && S800_results.Aresidue == 32)
        Histo_sort->Vlab_Ar32->Fill(res_Vel);
      if(S800_results.Zresidue == 18 && S800_results.Aresidue == 33)
        Histo_sort->Vlab_Ar33->Fill(res_Vel);
      if(S800_results.Zresidue == 18 && S800_results.Aresidue == 34)
        Histo_sort->Vlab_Ar34->Fill(res_Vel);
      if(S800_results.Zresidue == 18 && S800_results.Aresidue == 35)
        Histo_sort->Vlab_Ar35->Fill(res_Vel);
      if(S800_results.Zresidue == 17 && S800_results.Aresidue == 31)
        Histo_sort->Vlab_Cl31->Fill(res_Vel);
      if(S800_results.Zresidue == 17 && S800_results.Aresidue == 32)
        Histo_sort->Vlab_Cl32->Fill(res_Vel);
      if(S800_results.Zresidue == 17 && S800_results.Aresidue == 33)
        Histo_sort->Vlab_Cl33->Fill(res_Vel);
      if(S800_results.Zresidue == 17 && S800_results.Aresidue == 34)
        Histo_sort->Vlab_Cl34->Fill(res_Vel);
    }*/

    //this section has gamma spectra with NO add back and energies MAY OR MAY NOT be matched to time
    for(int m=0;m<Ceasar->NE;m++)
    {
      //Ceasar->DataEC[m].energy for no addback
      float phig = Ceasar->DataEC[m].phi;
      float thetag = Ceasar->DataEC[m].theta;
      float xg = sin(thetag)*cos(phig);
      float yg = sin(thetag)*sin(phig);
      float zg = cos(thetag);
      float dot = (xg*xR)+(yg*yR)+(zg*zR);
      float mag = sqrt(pow(xg,2)+pow(yg,2)+pow(zg,2))*sqrt(pow(xR,2)+pow(yR,2)+pow(zR,2));
      float thetarel = acos(dot/mag);   
      float dopp = Doppler->correct(Ceasar->DataEC[m].energy,thetarel, res_Vel);

      //TODO this seems to look at rigidity, also want that for Si22 heavy frags
      /*if(S800_results.Zbeam == 20 && S800_results.Abeam ==37)
      {
        //Ca36
        if(S800_results.Zresidue == 20 && S800_results.Aresidue == 36)
          if(dopp > 0)
          {
            Histo_read->TEC_Ca36_DataEC->Fill(dopp);
            if (dopp > 2.5 && dopp < 3.3)
              Histo_read->Rigidity_gammagated2plus->Fill(Hira->RingCounter->Solution[Nsol-1].rigidity);
            if (dopp > 2.8 && dopp < 3.3)
              Histo_read->Rigidity_gammagated2plus_strict->Fill(Hira->RingCounter->Solution[Nsol-1].rigidity);
          }
        //K35
        if(S800_results.Zresidue == 19 && S800_results.Aresidue == 35)
        {
          if(dopp > 0)
            Histo_read->TEC_K35_DataEC->Fill(dopp);
          Histo_read->TEC_K35_DataEC_nodopp->Fill(Ceasar->DataEC[m].energy);
        }

      }*/

    }
    
    //TODO more CAESAR stuff, comment out for now
    /*
    //this section has gamma spectra with add back and NO time matching
    for(int l=0;l<Ceasar->N_NoTaddback;l++)
    {
      //Ceasar->[l].energy for no addback
      float phig = Ceasar->NoTadded[l].phi;
      float thetag = Ceasar->NoTadded[l].theta;
      float xg = sin(thetag)*cos(phig);
      float yg = sin(thetag)*sin(phig);
      float zg = cos(thetag);
      float dot = (xg*xR)+(yg*yR)+(zg*zR);
      float mag = sqrt(pow(xg,2)+pow(yg,2)+pow(zg,2))*sqrt(pow(xR,2)+pow(yR,2)+pow(zR,2));
      float thetarel = acos(dot/mag);   
      float dopp = Doppler->correct(Ceasar->NoTadded[l].energy,thetarel, res_Vel);

      if(S800_results.Zbeam == 20 && S800_results.Abeam ==37)
      {
        //Ca36
        if(S800_results.Zresidue == 20 && S800_results.Aresidue == 36)
          if(dopp > 0)
            Histo_read->TEC_Ca36_NoTadded->Fill(dopp);
        //K35
        if(S800_results.Zresidue == 19 && S800_results.Aresidue == 35)
          if(dopp > 0)
            Histo_read->TEC_K35_NoTadded->Fill(dopp);
      }

    }

    //This section has gamma spectra with no add back but energies are matched to time
    for(int n=0;n<Ceasar->Nselect;n++)
    {
      float phig = Ceasar->select[n].phi;
      float thetag = Ceasar->select[n].theta;
      float xg = sin(thetag)*cos(phig);
      float yg = sin(thetag)*sin(phig);
      float zg = cos(thetag);
      float dot = (xg*xR)+(yg*yR)+(zg*zR);
      float mag = sqrt(pow(xg,2)+pow(yg,2)+pow(zg,2))*sqrt(pow(xR,2)+pow(yR,2)+pow(zR,2));
      float thetarel = acos(dot/mag); 
      float dopp = Doppler->correct(Ceasar->select[n].energy,thetarel,res_Vel);
   
      if(S800_results.Zbeam == 20 && S800_results.Abeam ==37)
      {
        if(S800_results.Zresidue == 20 && S800_results.Aresidue == 36)            
        {
          int id = Ceasar->select[n].id;
          Histo_read->TEC_Ca36_noaddback_nodoppler->Fill(Ceasar->select[n].energy);
          Histo_sort->Ca36gamma_nodopp[id]->Fill(Ceasar->select[n].energy);
          Histo_sort->Ca36gamma[id]->Fill(dopp);
          if(dopp > 0)            
            Histo_read->TEC_Ca36_noaddback->Fill(dopp);

        if(S800_results.Zresidue == 19 && S800_results.Aresidue == 36)
          if(dopp > 0)                        
            Histo_read->TEC_K36_noaddback->Fill(dopp);  
        }
      }
    }

    //This section has the add back gamma spectra
    for(int i=0;i<Ceasar->Nadded;i++)
    {
      //Ceasar->added[i] this hass addback
      float phig = Ceasar->added[i].phi;
      float thetag = Ceasar->added[i].theta;
      float xg = sin(thetag)*cos(phig);
      float yg = sin(thetag)*sin(phig);
      float zg = cos(thetag);
      float dot = (xg*xR)+(yg*yR)+(zg*zR);
      float mag = sqrt(pow(xg,2)+pow(yg,2)+pow(zg,2))*sqrt(pow(xR,2)+pow(yR,2)+pow(zR,2));
      float thetarel = acos(dot/mag);   
      float dopp = Doppler->correct(Ceasar->added[i].energy,thetarel, res_Vel);
      float dopp_simple = Doppler->correct(Ceasar->added[i].energy,thetag, res_Vel);


      Egamma0 = Ceasar->added[i].energy;
      thetarel0 = thetarel;
      thetagamma = thetag;
      thetares = thetaR;
      res_Vel0 = res_Vel;
      detID = Ceasar->added[i].id;
      Ring = Ceasar->added[i].iRing;
      Loc = Ceasar->added[i].iLoc;
      Edop = dopp;


                 
      //on any fragment hit into the S800, plot the time distributions
      Histo_sort->TCeasarS800[Ceasar->added[i].id]->Fill(Ceasar->added[i].itime/10);
      Histo_sort->TCeasarRawS800->Fill(Ceasar->added[i].itime/10);

      if(S800_results.Zbeam == 20 && S800_results.Abeam ==37)
      {
        //Ca37        
        if(S800_results.Zresidue == 20 && S800_results.Aresidue == 37)
          if(dopp > 0)
          {                   
            Histo_read->TEC_Ca37->Fill(dopp);
            Histo_read->dopp_timing_Ca37->Fill(dopp, Ceasar->added[i].time);
            Histo_read->TEnoC_Ca37->Fill(Ceasar->added[i].energy);
          }
        if(S800_results.Zresidue == 20 && S800_results.Aresidue == 38)
          if(dopp > 0)                        
          {                   
            Histo_read->TEC_Ca38->Fill(dopp);
            Histo_read->dopp_timing_Ca38->Fill(dopp, Ceasar->added[i].time);
          }

        //gamma spec for Ca36 s800 setting
        //Ca36
        if(S800_results.Zresidue == 20 && S800_results.Aresidue == 36)
          if(dopp > 0)
          {
            Histo_sort->tree->Fill();

            //cout << thetag << " " << thetarel << " " << thetaR << endl;
            Histo_sort->residuetheta->Fill(thetaR);
            Histo_sort->Vlab_Ca36_wgamma->Fill(res_Vel);
            Histo_read->TEC_Ca36->Fill(dopp);

            //Histo_read->TEC_Ca36_simpletheta->Fill(dopp_simple);
            if (Ceasar->added[i].addbackmult == 0)
              Histo_read->TEC_Ca36_mult0->Fill(dopp);
            if (Ceasar->added[i].addbackmult == 1)
              Histo_read->TEC_Ca36_mult1->Fill(dopp);
            if (Ceasar->added[i].addbackmult == 2)
              Histo_read->TEC_Ca36_mult2->Fill(dopp);
            if (Ceasar->added[i].addbackmult == 3)
              Histo_read->TEC_Ca36_mult3->Fill(dopp);

            if (cos(thetarel) > 0.5)
              Histo_read->TEC_Ca36_forward->Fill(dopp);
            if (cos(thetarel) > 0 && cos(thetarel) < 0.5)
              Histo_read->TEC_Ca36_slightforward->Fill(dopp);
            if (cos(thetarel) > -0.5 && cos(thetarel) < 0)
              Histo_read->TEC_Ca36_slightbackward->Fill(dopp);        
            if (cos(thetarel) < -0.5)
              Histo_read->TEC_Ca36_backward->Fill(dopp);

            if (run<54)
              Histo_read->TEC_Ca36_early->Fill(dopp);
            else
              Histo_read->TEC_Ca36_late->Fill(dopp);

            //timing plots
            Histo_read->dopp_timing->Fill(dopp, Ceasar->added[i].time);

            if (Ceasar->added[i].itime/10 < -500)
              Histo_read->TEC_Ca36_timecutearly->Fill(dopp);
            if (Ceasar->added[i].itime/10 > -500)
              Histo_read->TEC_Ca36_timecutlate->Fill(dopp);

          }
        //gamma spec for K35 s800 setting
        //K35        
        if(S800_results.Zresidue == 19 && S800_results.Aresidue == 35)
          if(dopp > 0)                        
            Histo_read->TEC_K35->Fill(dopp);
        //K36
        if(S800_results.Zresidue == 19 && S800_results.Aresidue == 36)
        {
          if(dopp > 0)                        
            Histo_read->TEC_K36->Fill(dopp);
          if (cos(thetarel) > 0.5)
            Histo_read->TEC_K36_forward->Fill(dopp);
          if (cos(thetarel) > 0 && cos(thetarel) < 0.5)
            Histo_read->TEC_K36_slightforward->Fill(dopp);
          if (cos(thetarel) > -0.5 && cos(thetarel) < 0)
            Histo_read->TEC_K36_slightbackward->Fill(dopp);        
          if (cos(thetarel) < -0.5)
            Histo_read->TEC_K36_backward->Fill(dopp);
        }

        //Ar32
        if(S800_results.Zresidue == 18 && S800_results.Aresidue == 32)
          if(dopp > 0)                        
            Histo_read->TEC_Ar32->Fill(dopp);
        //Ar33
        if(S800_results.Zresidue == 18 && S800_results.Aresidue == 33)
          if(dopp > 0)                        
            Histo_read->TEC_Ar33->Fill(dopp);
        //Ar34
        if(S800_results.Zresidue == 18 && S800_results.Aresidue == 34)
          if(dopp > 0)                        
            Histo_read->TEC_Ar34->Fill(dopp);
        //Ar35
        if(S800_results.Zresidue == 18 && S800_results.Aresidue == 35)
          if(dopp > 0)                        
            Histo_read->TEC_Ar35->Fill(dopp);
      }
      if(S800_results.Zbeam == 19 && S800_results.Abeam ==36)
      {
        if(S800_results.Zresidue == 19 && S800_results.Aresidue == 36)
        {
          if(dopp > 0)
          {             
            Histo_read->TEC_K36_K36beam->Fill(dopp);
            Histo_read->TEnoC_K36_K36beam->Fill(Ceasar->added[i].energy);
          }
        }
      }

    }*/

  }

  //////// end of gamma code

 
  //cout << "found residue = " << foundresidue << " " << S800_results.Zbeam << endl;
  //Load solutions into Correl
  //TODO telescope must be indexed
  int Mult = 0;
  for (int tele=0;tele<4;tele++)
  {
    for(int i=0;i<Gobbi->Telescope[tele]->Nsolution;i++)
    {
      if(Gobbi->Telescope[tele]->Solution[i].iZ >0)
      {
        Correl.load(&Gobbi->Telescope[tele]->Solution[i]);
        Mult++;
        solnZ++;
      }
      else
      {
        solnnoZ++;
      }
    }
  }
  
  //Should only be one S800 solution
  //Load S800 solution to Correl
  if (Gobbi->S800->Nsolution > 0 && Gobbi->S800->Solution[0].iZ > 0)
  {
    Correl.load(&Gobbi->S800->Solution[0]);
    Mult++;
    solnZ++;
  }
  
  //use this for searching for particle coincidences
  int CorrMultiplicity = 0;
  for (int p=0; p < Correl.Nparticles; p++)
  {
    CorrMultiplicity += Correl.particle[p]->mult;
  }
  //making a correlation table here to check if I'm missing any interesting channels
  if (CorrMultiplicity == 2)
  {
    int pos = 0;
    int particlenum[2] = {0,0};
    for (int p=0; p < Correl.Nparticles; p++)
    {
      if (Correl.particle[p]->mult == 1)
      {
        particlenum[pos] = p;
        pos++;
      }
      if (Correl.particle[p]->mult == 2)
      {
        particlenum[0] = p;
        particlenum[1] = p;
      }
    }
    Histo_read->CorrelationTable->Fill(particlenum[0],particlenum[1]);
  }

/*
  if (CorrMultiplicity == 3)
  {
    cout << "multiplicity 3 event" << endl;
    for (int p=0; p < Correl.Nparticles; p++)
    {
      if (Correl.particle[p]->mult >= 1)
      {
        for (int m=0; m < Correl.particle[p]->mult; m++)
        {
          cout << "z=" << Correl.particle[p]->Z <<",A=" << Correl.particle[p]->A << endl;
        }
      }
    }
  }
*/

  //if (CorrMultiplicity >= 2)
  //{
  //  if (Correl.H2.mult >=1)
  //  {
  //    cout << endl;
  //    cout << "deuteron with";
  //    for (int p=0; p < Correl.Nparticles; p++)
  //    {
  //      if (Correl.particle[p]->mult >= 1)
  //      {
  //        cout << " Z= " << Correl.particle[p]->Z << " A=" << Correl.particle[p]->A << "; ";
  //      }
  //    }
  //    cout << endl;
  //  }
  //}
 
  //if (Correl.H2.mult >=1)
  //{
  //  for (int p=0; p < Correl.Nparticles; p++)
  //  {
  //    if (Correl.particle[p]->mult >= 1)
  //    {
  //      cout << "a deutreron with particle " << Correl.particle[p]->Z << " A=" << Correl.particle[p]->A << endl;
  //    }
  //  }
  //}
    

  //TODO not sure what happens here
  /*
  bool K35frag = false;
  bool noZinfo = false;
  int loc = 0;
  //require K35
  if (S800_results.Zresidue==19 && S800_results.Aresidue==35)
  {
    //want fiber data for good tracking
    if (Hira->XY_mon->has_data)
    {
      //what solutions are saved for this?
      for (int i=0;i<Hira->RingCounter->Nsolution;i++)
      {
        if (Hira->RingCounter->Solution[i].iA == 35 && Hira->RingCounter->Solution[i].iZ == 19)
        {
          K35frag = true;
        }
        if (Hira->RingCounter->Solution[i].iZ <= 0)
        {
          loc = i;
          noZinfo = true;
        }
      }
      if (K35frag && noZinfo) 
      {
        Histo_sort->missingdee->Fill(Hira->RingCounter->Solution[loc].energyR,Hira->RingCounter->Solution[loc].denergy);
      }
    }
  }*/

  if(foundresidue)
  {
    //TODO make correlation functions and call here
    corr_21Al();
    corr_22Si();
    corr_23P();
  }
}


//TODO much of this should remain the same
//If I make a 5th "Telescope" for the S800, I can load heavy frag solution into Telescope[4]->Solution[0]
//The S800 frag will also get x-y from Janus now. Janus is not included in Gobbi
//Frag eloss will need to have a variable target thickness and 1 mm of BC400
bool det::LoadS800toSolution()
{

  int Nsol = Gobbi->S800->Nsolution;
  Gobbi->S800->Nsolution++;

  int Z = s800->residue_pid[s800->S800Setting][s800->BeamID]->Z;
  int A = s800->residue_pid[s800->S800Setting][s800->BeamID]->A;

  Gobbi->S800->Solution[Nsol].iZ = Z;
  Gobbi->S800->Solution[Nsol].iA = A;

  Gobbi->S800->Solution[Nsol].beamZ = s800->beam_pid->Z;
  Gobbi->S800->Solution[Nsol].beamA = s800->beam_pid->A;
  //cout << "theta fiber = " << Hira->XY_mon->theta << " theta S800 = " << s800->track->theta << " ratio = " << Hira->XY_mon->theta/s800->track->theta << endl;

  //using the theta phi from the fiber
  float thetaf = 0;
  float phif = 0;

  //TODO need Janus working for this
  if(Janus->Fiber->theta !=-999)
  {
    thetaf =  Janus->Fiber->theta;
    Gobbi->S800->Solution[Nsol].theta = thetaf;
    phif = Janus->Fiber->phi;
    Gobbi->S800->Solution[Nsol].phi = phif;
  }
  else
  {
    return false; //throw out events that use S800
    //cout << "do I get a lot of these?" << endl;
    //thetaf = s800->track->theta;
    //phif = s800->track->phi;
    //Hira->RingCounter->Solution[Nsol].theta = thetaf;
    //Hira->RingCounter->Solution[Nsol].phi = phif;
  }


  //thetaf = s800->track->theta;
  //phif = s800->track->phi;
  //TODO TEMP uses S800 angle instead of fiber
  //Gobbi->S800->Solution[Nsol].theta = thetaf;
  //Gobbi->S800->Solution[Nsol].phi = phif;
  
  double mass = Gobbi->S800->Pid->getMass(Z,A);

  float ekin = s800->track->energy;
  float pc = sqrt(pow(ekin+mass*m0,2) - pow(mass*m0,2));
  float rigidity = pc/1000.*3.3356/Z;
  Gobbi->S800->Solution[Nsol].rigidity = rigidity;

  if (Z == 19 && A == 35)
  {
    Histo_sort->rigidityK35->Fill(rigidity);
    Histo_read->Vlab_HF_p35K_before->Fill(pc/(ekin + mass*m0));
  }
  if (Z == 19 && A == 36)
  {
    Histo_read->Vlab_HF_p36K_before->Fill(pc/(ekin + mass*m0));
  }
  if (Z == 20 && A == 36)
  {
    Histo_sort->rigidityCa36->Fill(rigidity);
  }

  //collection of beam rigidities
  if (s800->beam_pid->Z == 20 && s800->beam_pid->A == 37)
  {
    if (Z == 20 && A == 37)
    {
      Histo_sort->rigidityCa37beam->Fill(rigidity);
    }
    if (Z == 20 && A == 38)
    {
      Histo_sort->rigidityCa38beam->Fill(rigidity);
    }
  }
  if (s800->beam_pid->Z == 19 && s800->beam_pid->A == 36)
    if (Z == 19 && A == 36)
    {
      Histo_sort->rigidityK36beam->Fill(rigidity);
    }
  if (s800->beam_pid->Z == 18 && s800->beam_pid->A == 35)
    if (Z == 18 && A == 35)
    {
      Histo_sort->rigidityAr35beam->Fill(rigidity);
    }

  float thickness = 103.2/cos(thetaf);// mg/cm^2    1 mm of BC400
  ekin = losses_fiber->getEin(ekin,thickness,Z,A);
  //changed to TargetThickness/2 -ND 8/17/21
  thickness = Gobbi->S800->TargetThickness/2/cos(thetaf);//*1.80;

  ekin = losses_target->getEin(ekin,thickness,Z,A);
  Gobbi->S800->Solution[Nsol].Ekin = ekin;
  
  Gobbi->S800->Solution[Nsol].mass = mass*m0; //mass*(amu->MeV)

  float momentum = Gobbi->S800->Solution[Nsol].Kinematics.
    getMomentum(Gobbi->S800->Solution[Nsol].Ekin,Gobbi->S800->Solution[Nsol].mass);

  Gobbi->S800->Solution[Nsol].momentum = momentum;

  Gobbi->S800->Solution[Nsol].Mvect[0] = momentum*sin(thetaf)*cos(phif);
  Gobbi->S800->Solution[Nsol].Mvect[1] = momentum*sin(thetaf)*sin(phif);
  Gobbi->S800->Solution[Nsol].Mvect[2] = momentum*cos(thetaf);



  Gobbi->S800->Solution[Nsol].energyTot = Gobbi->S800->Solution[Nsol].mass
                                            + Gobbi->S800->Solution[Nsol].Ekin;

  Gobbi->S800->Solution[Nsol].velocity = momentum/Gobbi->S800->Solution[Nsol].energyTot;

  return true;
}


/////////////////////////////////////////////////
//Correlation functions
/////////////////////////////////////////////////
void det::corr_1H()
{
  if(Correl.proton.mult > 0)
  {
    for (int i=0;i<Correl.proton.mult;i++)
    {
      Histo_sort->protonKE->Fill(Correl.proton.Sol[i]->Ekin);
    }
  }
}

void det::corr_21Al()
{
  //21Al -> 20Mg + p
  if(Correl.Mg20.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q21Al = mass_21Al - (mass_20Mg + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Mg20.mask[0]=1;
    Correl.makeArray(1);

    float Erel_21Al = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_21Al - Q21Al;


    Histo_sort->Erel_21Al_p20Mg->Fill(Erel_21Al);
    Histo_sort->Ex_21Al_p20Mg->Fill(Ex);

    Histo_sort->ThetaCM_21Al_p20Mg->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_21Al_p20Mg->Fill(Correl.velocityCM);
    Histo_sort->Erel_p20Mg_costhetaH->Fill(Erel_21Al,Correl.cos_thetaH);
    Histo_sort->p20Mg_VCMvsErel->Fill(Correl.velocityCM,Erel_21Al);
  }
}

void det::corr_22Si()
{
  //22Si -> 20Mg + 2p
  if(Correl.Mg20.mult == 1 && Correl.proton.mult == 2)
  {
    //TODO Si22 and P23 masses are not known
    //     Use calculations?
    //float const Q22Si = mass_22Si - (mass_20Mg + 2.*mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.Mg20.mask[0]=1;
    Correl.makeArray(1);

    float Erel_22Si = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    //float Ex = Erel_22Si - Q22Si;


    Histo_sort->Erel_22Si_2p20Mg->Fill(Erel_22Si);
    //Histo_sort->Ex_22Si_2p20Mg->Fill(Ex);

    Histo_sort->ThetaCM_22Si_2p20Mg->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_22Si_2p20Mg->Fill(Correl.velocityCM);
    Histo_sort->Erel_2p20Mg_costhetaH->Fill(Erel_22Si,Correl.cos_thetaH);
    Histo_sort->pp20Mg_VCMvsErel->Fill(Correl.velocityCM,Erel_22Si);
  }
}

void det::corr_23P()
{
  //23P -> 20Mg + 3p
  if(Correl.Mg20.mult == 1 && Correl.proton.mult == 3)
  {
    //float const Q23P = mass_23P - (mass_20Mg + 3.*mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.proton.mask[2]=1;
    Correl.Mg20.mask[0]=1;
    Correl.makeArray(1);

    float Erel_23P = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    //float Ex = Erel_22Si - Q22Si;


    Histo_sort->Erel_23P_3p20Mg->Fill(Erel_23P);
    //Histo_sort->Ex_23P_3p20Mg->Fill(Ex);

    Histo_sort->ThetaCM_23P_3p20Mg->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_23P_3p20Mg->Fill(Correl.velocityCM);
    Histo_sort->Erel_3p20Mg_costhetaH->Fill(Erel_23P,Correl.cos_thetaH);
    Histo_sort->ppp20Mg_VCMvsErel->Fill(Correl.velocityCM,Erel_23P);
  }

  //23P -> 22Si + p
  if(Correl.Si22.mult == 1 && Correl.proton.mult == 1)
  {
    //float const Q23P = mass_23P - (mass_22Si + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Si22.mask[0]=1;
    Correl.makeArray(1);

    float Erel_23P = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    //float Ex = Erel_22Si - Q22Si;


    Histo_sort->Erel_23P_p22Si->Fill(Erel_23P);
    //Histo_sort->Ex_23P_p22Si->Fill(Ex);

    Histo_sort->ThetaCM_23P_p22Si->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_23P_p22Si->Fill(Correl.velocityCM);
    Histo_sort->Erel_p22Si_costhetaH->Fill(Erel_23P,Correl.cos_thetaH);
    Histo_sort->p22Si_VCMvsErel->Fill(Correl.velocityCM,Erel_23P);
  }
}
