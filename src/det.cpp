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
  Gobbi = new gobbi(ran,Histo_sort);
  //last value is the s800 setting, 0=Ca36, 1=K35
  Janus = new janus(Histo_sort, -1); //Janus target is set through Janus subroutine, initialize as -1 to check
  s800 = new S800(ran,Histo_sort,setting); 
  Ceasar = new ceasar(ran,Histo_sort,Histo_read, 13);
  Doppler = new doppler(0.3323); //beta for Ca-36

  losses_fiber = new CLosses(16,"_fiber.loss",true);
  losses_target = new CLosses(16,"_Be.loss",true);
  //losses_fiberAl = new CLosses(16,"_Al.loss",true); // Losses in aluminum foil between fiber layers

  Corrcomb = new corrcomb(setting);

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

  s800coinc_cnt = 0;

  Nmaxmix = 100;
  NmaxmixSingles = 3000;

  //string name_a8B("Wood_a8B.root");
  //Wood_a8B = new wood(2,false,&name_a8B);

  //Initialize tree for gammas
  tgammas = new wood_gammas("tree/tgammas/treegammas.root");

  //Initialize correlation trees
  string corr_name;
  bool gamma = true; //Change to false before a channel with no gammas

  corr_name = "tree/p11C/p11C.root";
  p11C = new wood(2,gamma,corr_name);

  corr_name = "tree/p12C/p12C.root";
  p12C = new wood(2,gamma,corr_name);

  corr_name = "tree/p13N/p13N.root";
  p13N = new wood(2,gamma,corr_name);

  corr_name = "tree/p14N/p14N.root";
  p14N = new wood(2,gamma,corr_name);

  corr_name = "tree/p14O/p14O.root";
  p14O = new wood(2,gamma,corr_name);

  corr_name = "tree/pp14O/pp14O.root";
  pp14O = new wood(3,gamma,corr_name);

  corr_name = "tree/p15O/p15O.root";
  p15O = new wood(2,gamma,corr_name);

  corr_name = "tree/pp15O/pp15O.root";
  pp15O = new wood(3,gamma,corr_name);

  corr_name = "tree/p16O/p16O.root";
  p16O = new wood(2,gamma,corr_name);

  corr_name = "tree/pp16O/pp16O.root";
  pp16O = new wood(3,gamma,corr_name);

  corr_name = "tree/p18F/p18F.root";
  p18F = new wood(2,gamma,corr_name);

  corr_name = "tree/pp18F/pp18F.root";
  pp18F = new wood(3,gamma,corr_name);

  corr_name = "tree/p17F/p17F.root";
  p17F = new wood(2,gamma,corr_name);

  corr_name = "tree/pp17F/pp17F.root";
  pp17F = new wood(3,gamma,corr_name);

  corr_name = "tree/p18Ne/p18Ne.root";
  p18Ne = new wood(2,gamma,corr_name);

  corr_name = "tree/p19Ne/p19Ne.root";
  p19Ne = new wood(2,gamma,corr_name);

  corr_name = "tree/p20Ne/p20Ne.root";
  p20Ne = new wood(2,gamma,corr_name);

  corr_name = "tree/pp17Ne/pp17Ne.root";
  pp17Ne = new wood(3,gamma,corr_name);

  corr_name = "tree/pp18Ne/pp18Ne.root";
  pp18Ne = new wood(3,gamma,corr_name);

  corr_name = "tree/p20Na/p20Na.root";
  p20Na = new wood(2,gamma,corr_name);

  corr_name = "tree/pp19Ne/pp19Ne.root";
  pp19Ne = new wood(3,gamma,corr_name);

  corr_name = "tree/p21Na/p21Na.root";
  p21Na = new wood(2,gamma,corr_name);

  corr_name = "tree/ppp17Ne/ppp17Ne.root";
  ppp17Ne = new wood(4,gamma,corr_name);

  corr_name = "tree/p20Mg/p20Mg.root";
  p20Mg = new wood(2,gamma,corr_name);

  corr_name = "tree/ppp18Ne/ppp18Ne.root";
  ppp18Ne = new wood(4,gamma,corr_name);

  corr_name = "tree/p21Mg/p21Mg.root";
  p21Mg = new wood(2,gamma,corr_name);

  corr_name = "tree/p22Mg/p22Mg.root";
  p22Mg = new wood(2,gamma,corr_name);

  corr_name = "tree/pp20Mg/pp20Mg.root";
  pp20Mg = new wood(3,gamma,corr_name);

  corr_name = "tree/p22Al/p22Al.root";
  p22Al = new wood(2,gamma,corr_name);

  corr_name = "tree/pp21Mg/pp21Mg.root";
  pp21Mg = new wood(3,gamma,corr_name);

  corr_name = "tree/p23Al/p23Al.root";
  p23Al = new wood(2,gamma,corr_name);

  corr_name = "tree/ppp20Mg/ppp20Mg.root";
  ppp20Mg = new wood(4,gamma,corr_name);

  corr_name = "tree/p22Si/p22Si.root";
  p22Si = new wood(2,gamma,corr_name);

  corr_name = "tree/p23Si/p23Si.root";
  p23Si = new wood(2,gamma,corr_name);

  corr_name = "tree/pppp18Ne/pppp18Ne.root";
  pppp18Ne = new wood(5,gamma,corr_name);

  //alphas
  corr_name = "tree/a15O/a15O.root";
  a15O = new wood(2,gamma,corr_name);


  //delete corr_name;

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
 // delete losses_fiberAl;
  delete Gobbi;
  delete Janus;
  //delete ran;
  delete s800;
  delete Ceasar;
  delete Corrcomb;

  delete tgammas;
  
  //Delete correlation trees
  //cout << "tree here!!!" << endl;
  delete p11C;
  delete p12C;
  delete p13N;
  delete p14N;
  delete p14O;
  delete pp14O;
  delete p15O;
  delete pp15O;
  delete p16O;
  delete pp16O;
  delete p18F;
  delete pp18F;
  delete p17F;
  delete pp17F;
  delete p18Ne;
  delete p19Ne;
  delete p20Ne;
  delete pp17Ne;
  delete pp18Ne;
  delete p20Na;
  delete pp19Ne;
  delete p21Na;
  delete ppp17Ne;
  delete p20Mg;
  delete ppp18Ne;
  delete p21Mg;
  delete p22Mg;
  delete pp20Mg;
  delete p22Al;
  delete pp21Mg;
  delete p23Al;
  delete ppp20Mg;
  delete p22Si;
  delete p23Si;
  delete pppp18Ne;

  //Alphas
  delete a15O;

  cout << "You made it!" << endl;
}

void det::Reset()
{
  Gobbi->reset();
  s800->Reset();
  Janus->janusevts.clear();
	Janus->Histo->clear();
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

  runnum = runno; //Assign run number

  //TODO unpack needs to be rewritten so that it accepts an ifstream pointer
  //Janus must be passed the ifstream pointer, Si and S800 should have their buffers read out and passed
  int buffersize = fragmentsize - 28; //should always be 28
  int bufferwords = 1;
  if (sourceID == SiID || sourceID == S800ID) bufferwords = buffersize/2;
  unsigned short buffer_arr[bufferwords];
  unsigned short *pointbuf; //pointer to the buffer, needed for Si, ADC, TDC, S800, CAESAR
  //cout << hex << fragmentsize << " " << buffersize << " " << bufferwords << dec << endl;
  //Read out buffer if silicon or S800
  //cout << "runno " << runno << endl;
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
    //cout << "SiID SourceID " << sourceID << endl;
    //pointbuf++; //skips ADC size, don't actually need the byte size.
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
    stat = Ceasar->unpack(pointbuf, runno); // SG 2020/11/02
    
    if(!stat)
    {
      //  cout << "Didn't read hira right" <<endl;
      return stat;
    }
  }
  else if (sourceID == JanusID)
  { 
    //cout << "JanusID SourceID " << sourceID << endl;
    //call janus unpacker here
    stat = Janus->unpack(point);//something something janus unpacker

    //Turn skip on if you want to not unpack Janus
    //point->ignore(buffersize);
    stat = true;
  }
  else if (sourceID == S800ID)
  {
    //cout << "S800ID SourceID " << sourceID << endl;
    //To skip the S800, comment out the unpacker. Same as Gobbi
    pointbuf++; //skips the extra fragment size.
    stat = s800->unpack(pointbuf,runno);
    //if(stat)
      //NS800++;
  }
  else
  {
    cout << "found unexpected sourceID = " << sourceID << endl;
    return false;
  }
  
  return stat;
}
//*********************************

void det::analyze(int event, int run)
{
  Correl.reset();
  bool foundresidue = false;

  Si23_Mg20_2p_flag = false;
  Si23_Si22_p_flag = false;
  Si23_Si23_p_flag = false;

  Ne17_3p_flag = false;
  Ne18_3p_flag = false;
  
   //This is for doing S800 vs VME QDC 16 for Na22 testing
  //Check if QDC 16 is there in elist class in Gobbi
  double QDC_16 = 0;
  for (int i =0;i<Gobbi->NQ;i++)
  {
    if (Gobbi->DataQ[i].id == 16) QDC_16 = Gobbi->DataQ[i].ienergy;
  }
//cout << "" << endl;
  if (QDC_16 != 0) Histo_sort->e1upvsQDC->Fill(s800->Scint[0].de_up,QDC_16);
  
  //return;

  //return; //TODO temp so we don't go far

  s800_results S800_results;
  S800_results.Reset();
	S800_results.run = run;
  // S800_results.trig_coin=false;
  // S800_results.trig_singles=false;
  // S800_results.trig_s800_singles=false;
  
  if (s800->Trig.registr & 1) S800_results.trig_s800_singles =true;
  if (s800->Trig.registr & 2)
  {
    //cout << "coinc! " << endl;
    S800_results.trig_coin = true;
    s800coinc_cnt++;
    //if (s800coinc_cnt%3000 == 0) 
          //cout << '\xd'<< s800coinc_cnt << flush;
  }
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
  S800_results.tstamp = (double)s800->tstamp;

  //Events are already matched within boards. Call analyze() instead of MatchEvents()
  Janus->analyze(S800_results);

  //2d plots for DB5 time and obj time vs their Gobbi copies. Just take index 0
  //DB5 called XFP in code
  double TDC_16 = 0;
  int TDC16_count = 0;
  double TDC_17 = 0;
  int TDC17_count = 0;
  double TDC_18 = 0;
  int TDC18_count = 0;
  for (int i =0;i<Gobbi->NT;i++)
  {
    if (Gobbi->DataT[i].id == 16 && TDC_16 == 0)
    {
      TDC_16 = Gobbi->DataT[i].time;
      TDC16_count++;
    }
    if (Gobbi->DataT[i].id == 17 && TDC_17 == 0)
    {
      TDC_17 = Gobbi->DataT[i].time;
      TDC17_count++;
    }
    if (Gobbi->DataT[i].id == 18 && TDC_18 == 0)
    {
      TDC_18 = Gobbi->DataT[i].time;
      TDC18_count++;
    }
  }

  if (s800->mTDC.xfp.size() > 0 && TDC_16 > 0) Histo_sort->DB5T_vs_gCopy->Fill(s800->mTDC.xfp.at(0),TDC_16);
  if (s800->mTDC.obj.size() > 0 && TDC_17 > 0) Histo_sort->objT_vs_gCopy->Fill(s800->mTDC.obj.at(0),TDC_17);
  if (s800->mTDC.rf.size() > 0 && TDC_18 > 0) Histo_sort->RFT_vs_gCopy->Fill(s800->mTDC.rf.at(0),TDC_18);
  if (s800->mTDC.rf.size() > 0) Histo_sort->S800_RFTime->Fill(s800->mTDC.rf.at(0));
  if (s800->mTDC.rf.size() > 0) Histo_sort->Gobbi_RFTime->Fill(TDC_18);
  if (s800->mTDC.rf.size() > 0 && S800_results.trig_coin == true) Histo_sort->S800_RFTime_coin->Fill(s800->mTDC.rf.at(0));
  if (s800->mTDC.rf.size() > 0 && s800->beam_pid->Z == 14) Histo_sort->S800_RFTime_Sibeam->Fill(s800->mTDC.rf.at(0));
  if (s800->mTDC.rf.size() > 0 && s800->beam_pid->Z == 13) Histo_sort->S800_RFTime_Albeam->Fill(s800->mTDC.rf.at(0));
  if (s800->mTDC.rf.size() > 0 && s800->beam_pid->Z == 12) Histo_sort->S800_RFTime_Mgbeam->Fill(s800->mTDC.rf.at(0));
  if (s800->mTDC.rf.size() > 0 && s800->beam_pid->Z == 11) Histo_sort->S800_RFTime_Nabeam->Fill(s800->mTDC.rf.at(0));
  if (s800->mTDC.rf.size() > 0 && s800->beam_pid->Z > 0) Histo_sort->S800_RFTime_wbeam->Fill(s800->mTDC.rf.at(0));



  // if(!goods800)
  //   {
  //     return;
  //   }
  S800_results.Zbeam = s800->beam_pid->Z;
  S800_results.Abeam = s800->beam_pid->A;

  if (S800_results.Zbeam == 14 && S800_results.Abeam == 23)
  {
    //cout << "BEAM " << S800_results.Zbeam << " " << S800_results.Abeam << endl;
  }
  int release;

  //this method is doing a lot here.
  //  1. Applies S800 time gate to Si and CsI
  //  2. Adds neighboring Si strips within time gate
  //  3. matches up either E-CsI or dE-E events
  //  4. plots dE-E plots, and determines PID
  //  5. Does energy corrections for each telescope
  release = Gobbi->analyze(S800_results, run);  // TODO Gobbi needs to accept S800 results to apply time gates

  //Count the number of events with 2p and S800 Si23 beam
  if (S800_results.Zbeam == 14 && S800_results.Abeam == 23)
//  if (S800_results.Zbeam == 14 && S800_results.Abeam == 23 && S800_results.trig_coin == true)
  {
    //cout << "here" << endl;
    N_coin_Si23++;
    if (Gobbi->flag1p == true)
    {
      Histo_sort->ObjvsICsum_Si23_1p->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);
      N_coin_Si23_1p++;
    }
    if (Gobbi->flag2p == true)
    {
      Histo_sort->ObjvsICsum_Si23_2p->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);      
      N_coin_Si23_2p++;
    }
    if (Gobbi->flag3p == true)
    {
      Histo_sort->ObjvsICsum_Si23_3p->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);    
      N_coin_Si23_3p++;
    }
    if (Gobbi->flagalpha == true)
    {
      Histo_sort->ObjvsICsum_Si23_alpha->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);    
    }


  }

  //Count the number of events with 2p and S800 Al22 beam
  if (S800_results.Zbeam == 13 && S800_results.Abeam == 22)
  {
    //cout << "here" << endl;
    if (Gobbi->flag1p == true)
    {
      Histo_sort->ObjvsICsum_Al22_1p->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);
    }
    if (Gobbi->flag2p == true)
    {
      Histo_sort->ObjvsICsum_Al22_2p->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);      
    }
    if (Gobbi->flag3p == true)
    {
      Histo_sort->ObjvsICsum_Al22_3p->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);    
    }
    if (Gobbi->flagalpha == true)
    {
      Histo_sort->ObjvsICsum_Al22_alpha->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);    
    }

  }

  //Count the number of events with 2p and S800 Mg21 beam
  if (S800_results.Zbeam == 12 && S800_results.Abeam == 21)
  {
    //cout << "here" << endl;
    if (Gobbi->flag1p == true)
    {
      Histo_sort->ObjvsICsum_Mg21_1p->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);
    }
    if (Gobbi->flag2p == true)
    {
      Histo_sort->ObjvsICsum_Mg21_2p->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);      
    }
    if (Gobbi->flag3p == true)
    {
      Histo_sort->ObjvsICsum_Mg21_3p->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);    
    }
    if (Gobbi->flagalpha == true)
    {
      Histo_sort->ObjvsICsum_Mg21_alpha->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);    
    }

  }

  //Count the number of events with p and S800 Na20 beam
  if (S800_results.Zbeam == 11 && S800_results.Abeam == 20)
  {
    //cout << "here" << endl;
    if (Gobbi->flag1p == true)
    {
      Histo_sort->ObjvsICsum_Na20_1p->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);
    }
    if (Gobbi->flag2p == true)
    {
      Histo_sort->ObjvsICsum_Na20_2p->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);      
    }
    if (Gobbi->flag3p == true)
    {
      Histo_sort->ObjvsICsum_Na20_3p->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);    
    }
    if (Gobbi->flagalpha == true)
    {
      Histo_sort->ObjvsICsum_Na20_alpha->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);    
    }

  }

  if(s800->BeamID <0)
  {
    //TODO make sure to save janus info where we need it before clearing, clear() should really happen at the end
    Janus->janusevts.clear();
	  Janus->Histo->clear();
    Janus->nevts++;  
    return; //Need a good beam id
  }

  if (Gobbi->flag1p == true)
    Histo_sort->ObjvsICsum_wProt->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);
  if (Gobbi->flag2p == true) 
    Histo_sort->ObjvsICsum_wProt2->Fill(s800->mTDC.objCorrected.at(0),s800->IC.sum);
  

  S800_results.Zresidue = s800->residue_pid[s800->S800Setting][s800->BeamID]->Z;
  S800_results.Aresidue = s800->residue_pid[s800->S800Setting][s800->BeamID]->A;

  //Look at the Si23 unreacted beam centering
  if (S800_results.Zbeam == 14 && S800_results.Abeam == 23 && S800_results.Zresidue == 14 && S800_results.Aresidue == 23)
  {
    Histo_sort->crdcx23Si_23Si->Fill(s800->CRDC[0].x_gravity);
  }

  //Count the number of events with 2p, Mg20 fragment, and S800 Si23 beam
  if (S800_results.Zbeam == 14 && S800_results.Abeam == 23 && S800_results.Zresidue == 12 && S800_results.Aresidue == 20)
  //if (S800_results.Zbeam == 14 && S800_results.Abeam == 23 && S800_results.trig_coin == true && S800_results.Zresidue == 12 && S800_results.Aresidue == 20)
  {
    //cout << "here" << endl;
    N_coin_Si23_Mg20++;
    if (Gobbi->flag2p == true) Histo_sort->crdcx23S_20Mg_2p->Fill(s800->CRDC[0].x_gravity);
    if (Gobbi->flag1p == true) N_coin_Si23_Mg20_1p++;
    if (Gobbi->flag2p == true) N_coin_Si23_Mg20_2p++;
    if (Gobbi->flag2p == true) Si23_Mg20_2p_flag = true;

    //Addback gammas from Si23 beam / Mg20 res
    for(int i = 0; i < Ceasar->Nadded; i++) {
	    Histo_sort->energyTTot_Mg20->Fill(Ceasar->added[i].energy * 1000.); //E stored as MeV
    }
  }
  //Count the number of events with Mg21 fragment and S800 Si23 beam
  if (S800_results.Zbeam == 14 && S800_results.Abeam == 23 && S800_results.Zresidue == 12 && S800_results.Aresidue == 21)
  //if (S800_results.Zbeam == 14 && S800_results.Abeam == 23 && S800_results.trig_coin == true && S800_results.Zresidue == 12 && S800_results.Aresidue == 20)
  {
    //Addback gammas from Si23 beam / Mg21 res
    for(int i = 0; i < Ceasar->Nadded; i++) {
	    Histo_sort->energyTTot_Mg21->Fill(Ceasar->added[i].energy * 1000.); //E stored as MeV
    }
  }
  Si23_Si22_p_flag = false;
  //Count the number of events with p+22Si and S800 Si23 beam
  if (S800_results.Zbeam == 14 && S800_results.Abeam == 23 && S800_results.Zresidue == 14 && S800_results.Aresidue == 22)
  {
    if (Gobbi->flag1p == true) Histo_sort->crdcx23P_22Si_p->Fill(s800->CRDC[0].x_gravity);
    if (Gobbi->flag1p == true) N_coin_Si22_1p++;
    //if (Gobbi->flag1p == true) cout << " Si22 + p! " << endl;
    if (Gobbi->flag1p == true) Si23_Si22_p_flag = true;
  }

  //Count the number of events with p+23Si and S800 Si23 beam
  if (S800_results.Zbeam == 14 && S800_results.Abeam == 23 && S800_results.Zresidue == 14 && S800_results.Aresidue == 23)
  {
    if (Gobbi->flag1p == true) N_coin_Si23_1p++;
    if (Gobbi->flag1p == true) Si23_Si23_p_flag = true;
  }
  

  //TODO Skip Neon isotopes because they don't have a loss file.
  if (S800_results.Zresidue == 10)
  {

    //Count up using S800 PID flags
    if (Gobbi->flag3p == true && S800_results.Aresidue == 17) N_coin_Ne17_3p++;
    if (Gobbi->flag3p == true && S800_results.Aresidue == 17) Ne17_3p_flag = true;
    if (Gobbi->flag3p == true && S800_results.Aresidue == 18) N_coin_Ne18_3p++;
    if (Gobbi->flag3p == true && S800_results.Aresidue == 18) Ne18_3p_flag = true;


    /*Janus->janusevts.clear();
	  Janus->Histo->clear();
    Janus->nevts++;
    return;*/
  }

  //Plot the energy distribution for Si23 residues
  //Needed for S800 calibration
  if (S800_results.Zbeam == 14 && S800_results.Abeam == 23)
  {
    if (S800_results.Zresidue == 14 && S800_results.Aresidue == 23)
    {
      //Ekin of Si23 residue
      Histo_sort->energySi23_Si23->Fill(s800->track->energy);
    }
  }

  //Plot the energy distribution for Mg21 residues
  //Needed for S800 calibration
  if (S800_results.Zbeam == 12 && S800_results.Abeam == 21)
  {
    if (S800_results.Zresidue == 12 && S800_results.Aresidue == 21)
    {
      //Ekin of Mg21 residue
      Histo_sort->energyMg21_Mg21->Fill(s800->track->energy);
    }
  }

  //return;

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


  //Add gammas and coincident residue to the gamma-ray tree
  //If I have select gammas, I have "added" gammas
  if (S800_results.Zresidue > 0 && Ceasar->Nadded > 0 && Janus->Fiber->theta !=-999) {

    //reset tree variables
    //Must be done or the vectors go crazy
    tgammas->reset();

    //Fill residue information
    tgammas->Zres = S800_results.Zresidue;
    tgammas->Ares = S800_results.Aresidue;
    tgammas->theta_res = Janus->Fiber->theta;
    tgammas->phi_res = Janus->Fiber->phi;
    tgammas->ekin_res = s800->track->energy;

    //beamZ
    tgammas->beamZ = S800_results.Zbeam;

    //Gamma rays
    tgammas->Ngamma = Ceasar->Nadded;
    tgammas->Ngamma_Select = Ceasar->Nselect;

    for (int i=0;i<Ceasar->Nadded;i++) {
      tgammas->Egamma.push_back(Ceasar->added[i].energy);
      tgammas->Thetagamma.push_back(Ceasar->added[i].theta);
      tgammas->Phigamma.push_back(Ceasar->added[i].phi);
      tgammas->Tgamma.push_back(Ceasar->added[i].time);
      tgammas->Chgamma.push_back(Ceasar->added[i].id);
    }
    
    for (int i=0;i<Ceasar->Nselect;i++) {
      tgammas->Egamma_Select.push_back(Ceasar->select[i].energy);
      tgammas->Thetagamma_Select.push_back(Ceasar->select[i].theta);
      tgammas->Phigamma_Select.push_back(Ceasar->select[i].phi);
      tgammas->Tgamma_Select.push_back(Ceasar->select[i].time);
      tgammas->Chgamma_Select.push_back(Ceasar->select[i].id);
    }

    tgammas->T->Fill();
  }



  //For Si24
  if (S800_results.Zresidue==14 && S800_results.Aresidue==24 && Janus->Fiber->theta !=-999) {

    //Need to go ahead and do Egain and doppler corrections
    //Si24 does not come with protons so it won't matter.
    //Get fiber angle
    float si24theta = 0;
    float si24phi = 0;

    si24theta =  Janus->Fiber->theta;
    si24phi = Janus->Fiber->phi;

    //Get some necessary values
    double si24mass = Gobbi->S800->Pid->getMass(14,24);
    //cout << Z << " " << A << " " << mass << endl;
    float si24ekin = s800->track->energy;
    float si24pc = sqrt(pow(si24ekin+si24mass,2) - pow(si24mass,2));
    float si24rigidity = si24pc/1000.*3.3356/14;

    //Egain time
    float si24thickness = 112.66344/cos(si24theta);// mg/cm^2   103.2 is 1 mm of BC400, S800 gives 1.2 mm (123.84)
    si24ekin = losses_fiber->getEin(si24ekin,si24thickness,14,24);

    si24thickness = Gobbi->S800->TargetThickness/2./cos(si24theta);//*1.80; 
      
    si24ekin = losses_target->getEin(si24ekin,si24thickness,14,24);

    float si24momentum = sqrt(pow(si24ekin+si24mass,2) - pow(si24mass,2));

    float si24energyTot = si24mass + si24ekin;

    float si24velocity = si24momentum/si24energyTot;

    //Get Si24 unit vector
    xH = sin(si24theta)*cos(si24phi);
    yH = sin(si24theta)*sin(si24phi);
    zH = cos(si24theta);

    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Time gates
      if (Ceasar->added[i].time < 600 || Ceasar->added[i].time > 650) continue;

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,si24velocity);

     // cout << doppE << endl;

	   Histo_sort->enAddback_Si24->Fill(doppE * 1000.); //E stored as MeV
     Histo_sort->enAddback_Si24vsT->Fill(doppE * 1000.,Ceasar->added[i].time);
    }
  }

  //For Si23
  if (S800_results.Zresidue==14 && S800_results.Aresidue==23 && Janus->Fiber->theta !=-999) {

    //Need to go ahead and do Egain and doppler corrections
    //Get fiber angle
    float si23theta = 0;
    float si23phi = 0;

    si23theta =  Janus->Fiber->theta;
    si23phi = Janus->Fiber->phi;

    //Get some necessary values
    double si23mass = Gobbi->S800->Pid->getMass(14,23);
    //cout << Z << " " << A << " " << mass << endl;
    float si23ekin = s800->track->energy;
    float si23pc = sqrt(pow(si23ekin+si23mass,2) - pow(si23mass,2));
    float si23rigidity = si23pc/1000.*3.3356/14;

    //Egain time
    float si23thickness = 112.66344/cos(si23theta);// mg/cm^2   103.2 is 1 mm of BC400, S800 gives 1.2 mm (123.84)
    si23ekin = losses_fiber->getEin(si23ekin,si23thickness,14,23);

    si23thickness = Gobbi->S800->TargetThickness/2./cos(si23theta);//*1.80; 
      
    si23ekin = losses_target->getEin(si23ekin,si23thickness,14,23);

    float si23momentum = sqrt(pow(si23ekin+si23mass,2) - pow(si23mass,2));

    float si23energyTot = si23mass + si23ekin;

    float si23velocity = si23momentum/si23energyTot;

    //Get si23 unit vector
    xH = sin(si23theta)*cos(si23phi);
    yH = sin(si23theta)*sin(si23phi);
    zH = cos(si23theta);

    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Time gates
      //if (Ceasar->added[i].time < 600 || Ceasar->added[i].time > 650) continue;

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,si23velocity);

     // cout << doppE << endl;

	   Histo_sort->enAddback_Si23->Fill(doppE * 1000.); //E stored as MeV
	   if (Ceasar->Nadded == 1) Histo_sort->enAddback_Si23_mult1->Fill(doppE * 1000.); //E stored as MeV
     Histo_sort->enAddback_Si23vsT->Fill(doppE * 1000.,Ceasar->added[i].time);
    }
  }

 //For Al22
  if (S800_results.Zresidue==13 && S800_results.Aresidue==22 && Janus->Fiber->theta !=-999) {

    //Need to go ahead and do Egain and doppler corrections
    //Al22 does not always come with protons so it won't matter.
    //Get fiber angle
    float al22theta = 0;
    float al22phi = 0;

    al22theta =  Janus->Fiber->theta;
    al22phi = Janus->Fiber->phi;

    //Get some necessary values
    double al22mass = Gobbi->S800->Pid->getMass(13,22);
    //cout << Z << " " << A << " " << mass << endl;
    float al22ekin = s800->track->energy;
    float al22pc = sqrt(pow(al22ekin+al22mass,2) - pow(al22mass,2));
    float al22rigidity = al22pc/1000.*3.3356/13;

    //Egain time
    float al22thickness = 112.66344/cos(al22theta);// mg/cm^2   103.2 is 1 mm of BC400, S800 gives 1.2 mm (123.84)
    al22ekin = losses_fiber->getEin(al22ekin,al22thickness,13,22);

    al22thickness = Gobbi->S800->TargetThickness/2./cos(al22theta);//*1.80; 
      
    al22ekin = losses_target->getEin(al22ekin,al22thickness,13,22);

    float al22momentum = sqrt(pow(al22ekin+al22mass,2) - pow(al22mass,2));

    float al22energyTot = al22mass + al22ekin;

    float al22velocity = al22momentum/al22energyTot;

    //Get al22 unit vector
    xH = sin(al22theta)*cos(al22phi);
    yH = sin(al22theta)*sin(al22phi);
    zH = cos(al22theta);

    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Time gates
      //if (Ceasar->added[i].time < 600 || Ceasar->added[i].time > 650) continue;

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,al22velocity);

     // cout << doppE << endl;

	   Histo_sort->enAddback_Al22->Fill(doppE * 1000.); //E stored as MeV
	   if (Ceasar->Nadded == 1) Histo_sort->enAddback_Al22_mult1->Fill(doppE * 1000.); //E stored as MeV
     Histo_sort->enAddback_Al22vsT->Fill(doppE * 1000.,Ceasar->added[i].time);
	   if (Ceasar->Nadded == 1) Histo_sort->enAddback_Al22vsT_mult1->Fill(doppE * 1000.,Ceasar->added[i].time);

    //Look for events outside of the time gates one would use. "anti time gated"
    if (Ceasar->added[i].time < 600 && Ceasar->added[i].time > 670) {
      Histo_sort->enAddback_Al22_antiTG->Fill(doppE * 1000.); //E stored as MeV
	    if (Ceasar->Nadded == 1) Histo_sort->enAddback_Al22_mult1_antiTG->Fill(doppE * 1000.); //E stored as MeV
    }

    }

    //Very low E gamma. Also do select
    for(int i = 0; i < Ceasar->Nselect; i++) {

      //Time gates
      //if (Ceasar->added[i].time < 600 || Ceasar->added[i].time > 650) continue;

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->select[i].theta)*cos(Ceasar->select[i].phi);
      double yg = sin(Ceasar->select[i].theta)*sin(Ceasar->select[i].phi);
      double zg = cos(Ceasar->select[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->select[i].energy,thetarel,al22velocity);

     // cout << doppE << endl;

	   Histo_sort->enSelect_Al22->Fill(doppE); //E stored as keV for select
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
  int Mult = 0;  
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
      Mult++;
      solnZ++;

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
  for (int tele=0;tele<4;tele++)
  {
    for(int i=0;i<Gobbi->Telescope[tele]->Nsolution;i++)
    {
      if(Gobbi->Telescope[tele]->Solution[i].iZ >0)
      {
        Correl.load(&Gobbi->Telescope[tele]->Solution[i]);
        Mult++;
        solnZ++;
        //cout << "here" << endl;

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

  //Fill proton KE for different beams
  if(Correl.proton.mult > 0)
  {
    for (int i=0;i<Correl.proton.mult;i++)
    {
      int beamZ = S800_results.Zbeam;
      int beamA = S800_results.Abeam;
      if (beamZ == 14 && beamA == 23) Histo_sort->protonKE_Si23Beam->Fill(Correl.proton.Sol[i]->Ekin);
      if (beamZ == 13 && beamA == 22) Histo_sort->protonKE_Al22Beam->Fill(Correl.proton.Sol[i]->Ekin);
      if (beamZ == 12 && beamA == 21) Histo_sort->protonKE_Mg21Beam->Fill(Correl.proton.Sol[i]->Ekin);
    }
  }

    
  if(solnZ > 0)
  {

    //Loop for plotting correlation combinations
    //Sends fragment data to corrcomb.cpp and extracts channel value
    //Plots every combinations up to 4 types of light fragments
    int const lightmult_max = 4;
    int lightmult = 0;
    float lightcombs[lightmult_max];
    float combX = 0;
    float combY = 0;
    for (int i=0;i<Correl.Nparticles;i++) {

      if (Correl.particle[i]->mult < 1) continue;

      if (Correl.particle[i]->Z > 5) {
        //Only accept one heavy right now, take first
        if (combX != 0) continue;
        combX = Corrcomb->getX(S800_results.Zbeam,Correl.particle[i]->Z,Correl.particle[i]->A);
      }

      if (Correl.particle[i]->Z <= 5) {
        combY = Corrcomb->getY(S800_results.Zbeam,Correl.particle[i]->Z,Correl.particle[i]->A,Correl.particle[i]->mult);

        lightcombs[lightmult] = combY;
        lightmult++;
      }

    }


    for (int i=0;i<lightmult;i++) {
      if (combX == 0 || lightcombs[i] == 0) break;
      if (S800_results.Zbeam == 14) Histo_sort->SiBeam_CorrComb->Fill(combX,lightcombs[i]);
      if (S800_results.Zbeam == 13) Histo_sort->AlBeam_CorrComb->Fill(combX,lightcombs[i]);
    }


    //Calculates doppler if heavy arrived with protons and gammas
    if (Correl.proton.mult > 0 && Ceasar->Nadded > 0)
    {
      xH = 0;
      yH = 0;
      zH = 0;
      EH = 0;
      //Search through to find the heavy particle, should only be 1
      for (int i=0;i<Correl.Nparticles;i++)
      {
        if (Correl.particle[i]->mult > 0 && Correl.particle[i]->Z > 5) //If heavy fragment
        {
          //cout << Correl.particle[i]->mult << " " << Correl.particle[i]->Z << endl;
          //Calculate unit vector for heavy fragment
          xH = sin(Correl.particle[i]->Sol[0]->theta)*cos(Correl.particle[i]->Sol[0]->phi);
          yH = sin(Correl.particle[i]->Sol[0]->theta)*sin(Correl.particle[i]->Sol[0]->phi);
          zH = cos(Correl.particle[i]->Sol[0]->theta);
          //Calculate EH later, only if it has gamma, don't want to do a bunch of sqrts
          break;
        }
      }
    }

    /* //only works if loaded to correl
    // Si24 and Al22 mostly do not come with correlations, must be treated separately prior to loadtos800

    //Calculates doppler if Si24 or Al22 arrives with gammas 
    //This allows me to circumvent the light part. coincidence while not wasting time on other nuclei
    //Calculates doppler if heavy arrived with protons and gammas
    if (Correl.proton.mult == 0 && Ceasar->Nadded > 0 && (Correl.Si24.mult > 0 || Correl.Al22.mult > 0))
    {
      xH = 0;
      yH = 0;
      zH = 0;
      EH = 0;
      //Search through to find the heavy particle, should only be 1
      for (int i=0;i<Correl.Nparticles;i++)
      {
        if (Correl.particle[i]->mult > 0 && Correl.particle[i]->Z > 5) //If heavy fragment
        {
          //cout << Correl.particle[i]->mult << " " << Correl.particle[i]->Z << endl;
          //Calculate unit vector for heavy fragment
          xH = sin(Correl.particle[i]->Sol[0]->theta)*cos(Correl.particle[i]->Sol[0]->phi);
          yH = sin(Correl.particle[i]->Sol[0]->theta)*sin(Correl.particle[i]->Sol[0]->phi);
          zH = cos(Correl.particle[i]->Sol[0]->theta);
          //Calculate EH later, only if it has gamma, don't want to do a bunch of sqrts
          break;
        }
      }
    }

    //If heavy fragment is Si-24, calculate doppler and plot
    for (int i=0;i<Correl.Nparticles;i++) {
      if (Correl.particle[i]->mult > 0 && Correl.particle[i]->Z == 14 && Correl.particle[i]->A == 24) { //If Si-24
        for(int j = 0; j < Ceasar->Nadded; j++) {
          //Calculate unit vector for gamma
          double xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
          double yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
          double zg = cos(Ceasar->added[j].theta);

          //Euclidean vector lengths
          double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

          //Dot product
          double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

          //Find angle between LF and gamma
          float thetarel = acos(Ldotg/(EH*Eg));     

          //Find doppler corrected energy 
          double doppE;
          doppE = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.particle[i]->Sol[0]->velocity);

	        Histo_sort->enAddback_Si24->Fill(doppE * 1000.); //E stored as MeV

        }
    }


    //If heavy fragment is Al-22, calculate doppler and plot
      if (Correl.particle[i]->mult > 0 && Correl.particle[i]->Z == 13 && Correl.particle[i]->A == 22) {

        //For addback
        for(int j = 0; j < Ceasar->Nadded; j++) {
          //Calculate unit vector for gamma
          double xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
          double yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
          double zg = cos(Ceasar->added[j].theta);

          //Euclidean vector lengths
          double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

          //Dot product
          double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

          //Find angle between LF and gamma
          float thetarel = acos(Ldotg/(EH*Eg));     

          //Find doppler corrected energy 
          double doppE;
          doppE = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.particle[i]->Sol[0]->velocity);

	        Histo_sort->enAddback_Al22->Fill(doppE * 1000.); //E stored as MeV
          Histo_sort->enAddback_Al22vsCh->Fill(doppE * 1000.,Ceasar->added[j].id);

        }

        //For select
        for(int j = 0; j < Ceasar->Nselect; j++) {
          //Calculate unit vector for gamma
          double xg = sin(Ceasar->select[j].theta)*cos(Ceasar->select[j].phi);
          double yg = sin(Ceasar->select[j].theta)*sin(Ceasar->select[j].phi);
          double zg = cos(Ceasar->select[j].theta);

          //Euclidean vector lengths
          double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

          //Dot product
          double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

          //Find angle between LF and gamma
          float thetarel = acos(Ldotg/(EH*Eg));     

          //Find doppler corrected energy 
          double doppE;
          doppE = Doppler->correct(Ceasar->select[j].energy,thetarel,Correl.particle[i]->Sol[0]->velocity);

	        Histo_sort->enSelect_Al22->Fill(doppE * 1000.); //E stored as MeV

        }
      }
    }*/

    corr_1H();
    corr_12N();
    corr_13N();
    corr_14O();
    corr_15O();
    corr_15F();
    corr_16F();
    corr_17F();
    corr_19Ne();
    corr_18Ne();
    corr_17Ne();
    corr_16Ne();
    corr_19Na();
    corr_20Na();
    corr_21Na();
    corr_19Mg();
    corr_20Mg();
    corr_21Mg();
    corr_22Mg();
    corr_20Al();
    corr_21Al();
    corr_22Al();
    corr_23Al();
    corr_22Si();
    corr_23Si();
    corr_24Si();
    corr_23P();
    corr_24P();

    //Alphas
    corr_a15O();

  }

  //TODO make sure to save janus info where we need it before clearing, clear() should really happen at the end
  Janus->janusevts.clear();
	Janus->Histo->clear();
  Janus->nevts++;
}


//TODO much of this should remain the same
//If I make a 5th "Telescope" for the S800, I can load heavy frag solution into Telescope[4]->Solution[0]
//The S800 frag will also get x-y from Janus now. Janus is not included in Gobbi
//Frag eloss will need to have a variable target thickness and 1 mm of BC400
bool det::LoadS800toSolution()
{
  //************
  //Code originally written to do mass*m0 because Ca36 returned masses in amu
  //New version of PID returns masses in MeV

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

    if (Z == 9 && A == 17) {
      //cout << " Janus s800 theta: " << thetaf << " " << s800->track->theta << endl;
      //cout << " Janus s800 phi: " << phif << " " << s800->track->phi << endl;
      Histo_sort->ThetavsPhi_F17->Fill(s800->track->thetadeg,s800->track->phideg);
      Histo_sort->ThetavsPhi_F17_fiber->Fill(thetaf*180./acos(-1.),phif*180./acos(-1.));
    }

    Gobbi->S800->Solution[Nsol].theta_s800 = s800->track->theta; //get S800 values for fine angle comparison
    Gobbi->S800->Solution[Nsol].phi_s800 = s800->track->phi;

  }
  else
  {
    cout << "No janus, tossed" << endl;
    return false; //throw out events that use S800
    //thetaf = s800->track->theta;
    //phif = s800->track->phi;
    //Gobbi->S800->Solution[Nsol].theta = thetaf;
    //Gobbi->S800->Solution[Nsol].phi = phif;
  }
 

  double mass = Gobbi->S800->Pid->getMass(Z,A);
  //cout << Z << " " << A << " " << mass << endl;
  float ekin = s800->track->energy;
  float pc = sqrt(pow(ekin+mass,2) - pow(mass,2));
  float rigidity = pc/1000.*3.3356/Z;
  Gobbi->S800->Solution[Nsol].rigidity = rigidity;
  //cout << "Fragment track energy " << ekin << " ";
  //cout << "residue ekin/A " << ekin/A << endl;
  //cout << "rigidity " << rigidity << endl;
  //cout << "theta " << thetaf << " phi " << phif << endl;
  //cout << "theta " << Janus->Fiber->thetadeg << " phi " << Janus->Fiber->phideg << endl;
  //cout << "track values " << s800->track->theta << " " << s800->track->phi << endl;

  //Fill residue energy and rigidity

  if (Z == 12 && A == 20)
  {

  }

  //cout << "Det E " << ekin << " cos(thetaf) " << cos(thetaf) << " ";

  float thickness = 112.66344/cos(thetaf);// mg/cm^2   103.2 is 1 mm of BC400, S800 gives 1.2 mm (123.84)
  ekin = losses_fiber->getEin(ekin,thickness,Z,A);
  Histo_sort->energySi23_Si23_Fib->Fill(ekin); //Only use these for beam calibration runs
  Histo_sort->energyMg21_Mg21_Fib->Fill(ekin);
  //cout << "Fragment track energy " << ekin << " ";
  //changed to TargetThickness/2 -ND 8/17/21
  //cout << "Before Fib " << ekin << " ";
  thickness = Gobbi->S800->TargetThickness/2./cos(thetaf);//*1.80; 
  
  ekin = losses_target->getEin(ekin,thickness,Z,A);
  Histo_sort->energySi23_Si23_FibTarg->Fill(ekin); //Only use these for beam calibration runs
  Histo_sort->energyMg21_Mg21_FibTarg->Fill(ekin);
  //cout << "before target " << ekin << " Thickness " <<  endl;
  Gobbi->S800->Solution[Nsol].Ekin = ekin;

  Gobbi->S800->Solution[Nsol].mass = mass; //pid already returns total mass in MeV, don't need to convert

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

/////////////////////////////////////////////////
//Just writing these to trees
/////////////////////////////////////////////////
void det::corr_12N()
{
  //12N -> 11C + p
  if(Correl.C11.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q12N = mass_12N - (mass_11C + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.C11.mask[0]=1;
    Correl.makeArray(1);

    float Erel_12N = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_12N - Q12N;

    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p11C->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.C11.Sol[0]->velocity);

      
      p11C->Egamma[i] = doppE;
      p11C->Tgamma[i] = Ceasar->added[i].time;
      p11C->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    p11C->itele1 = Correl.proton.Sol[0]->itele;
    p11C->id1 = Correl.proton.Sol[0]->iCsI;
    p11C->ifront1 = Correl.proton.Sol[0]->ifront;
    p11C->iback1 = Correl.proton.Sol[0]->iback;
    p11C->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p11C->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p11C->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p11C->et1 = Correl.proton.Sol[0]->energyTot;
    p11C->time1 = Correl.proton.Sol[0]->CsITime;
    p11C->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p11C->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p11C->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p11C->theta1 = Correl.proton.Sol[0]->theta;
    p11C->phi1 = Correl.proton.Sol[0]->phi;

    p11C->M2[0] = Correl.C11.Sol[0]->Mvect[0];
    p11C->M2[1] = Correl.C11.Sol[0]->Mvect[1];
    p11C->M2[2] = Correl.C11.Sol[0]->Mvect[2];
    p11C->et2 = Correl.C11.Sol[0]->energyTot;
    p11C->energy_p2 = Correl.C11.Sol[0]->Ekin;
    p11C->theta2 = Correl.C11.Sol[0]->theta;
    p11C->phi2 = Correl.C11.Sol[0]->phi;

    p11C->Erel = Erel_12N;
    p11C->Ex = Ex;
    p11C->Vcm = Correl.velocityCM;
    p11C->thetaCM = thetaCM;
    p11C->cos_thetaH = Correl.cos_thetaH;

    p11C->runnum = runnum;
    p11C->beamZ = Correl.C11.Sol[0]->beamZ;

    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    p11C->T->Fill();

  }
}

void det::corr_13N()
{
  //13N -> 12C + p
  if(Correl.C12.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q13N = mass_13N - (mass_12C + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.C12.mask[0]=1;
    Correl.makeArray(1);

    float Erel_13N = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_13N - Q13N;

    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p12C->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.C12.Sol[0]->velocity);


      p12C->Egamma[i] = doppE;
      p12C->Tgamma[i] = Ceasar->added[i].time;
      p12C->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    p12C->itele1 = Correl.proton.Sol[0]->itele;
    p12C->id1 = Correl.proton.Sol[0]->iCsI;
    p12C->ifront1 = Correl.proton.Sol[0]->ifront;
    p12C->iback1 = Correl.proton.Sol[0]->iback;
    p12C->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p12C->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p12C->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p12C->et1 = Correl.proton.Sol[0]->energyTot;
    p12C->time1 = Correl.proton.Sol[0]->CsITime;
    p12C->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p12C->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p12C->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p12C->theta1 = Correl.proton.Sol[0]->theta;
    p12C->phi1 = Correl.proton.Sol[0]->phi;

    p12C->M2[0] = Correl.C12.Sol[0]->Mvect[0];
    p12C->M2[1] = Correl.C12.Sol[0]->Mvect[1];
    p12C->M2[2] = Correl.C12.Sol[0]->Mvect[2];
    p12C->et2 = Correl.C12.Sol[0]->energyTot;
    p12C->energy_p2 = Correl.C12.Sol[0]->Ekin;
    p12C->theta2 = Correl.C12.Sol[0]->theta;
    p12C->phi2 = Correl.C12.Sol[0]->phi;

    p12C->Erel = Erel_13N;
    p12C->Ex = Ex;
    p12C->Vcm = Correl.velocityCM;
    p12C->thetaCM = thetaCM;
    p12C->cos_thetaH = Correl.cos_thetaH;

    p12C->runnum = runnum;
    p12C->beamZ = Correl.C12.Sol[0]->beamZ;

    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    p12C->T->Fill();

  }
}

void det::corr_14O()
{
  //14O -> 13N + p
  if(Correl.N13.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q14O = mass_14O - (mass_13N + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.N13.mask[0]=1;
    Correl.makeArray(1);

    float Erel_14O = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_14O - Q14O;

    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p13N->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.N13.Sol[0]->velocity);


      p13N->Egamma[i] = doppE;
      p13N->Tgamma[i] = Ceasar->added[i].time;
      p13N->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    p13N->itele1 = Correl.proton.Sol[0]->itele;
    p13N->id1 = Correl.proton.Sol[0]->iCsI;
    p13N->ifront1 = Correl.proton.Sol[0]->ifront;
    p13N->iback1 = Correl.proton.Sol[0]->iback;
    p13N->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p13N->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p13N->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p13N->et1 = Correl.proton.Sol[0]->energyTot;
    p13N->time1 = Correl.proton.Sol[0]->CsITime;
    p13N->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p13N->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p13N->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p13N->theta1 = Correl.proton.Sol[0]->theta;
    p13N->phi1 = Correl.proton.Sol[0]->phi;

    p13N->M2[0] = Correl.N13.Sol[0]->Mvect[0];
    p13N->M2[1] = Correl.N13.Sol[0]->Mvect[1];
    p13N->M2[2] = Correl.N13.Sol[0]->Mvect[2];
    p13N->et2 = Correl.N13.Sol[0]->energyTot;
    p13N->energy_p2 = Correl.N13.Sol[0]->Ekin;
    p13N->theta2 = Correl.N13.Sol[0]->theta;
    p13N->phi2 = Correl.N13.Sol[0]->phi;

    p13N->Erel = Erel_14O;
    p13N->Ex = Ex;
    p13N->Vcm = Correl.velocityCM;
    p13N->thetaCM = thetaCM;
    p13N->cos_thetaH = Correl.cos_thetaH;

    p13N->beamZ = Correl.N13.Sol[0]->beamZ;

    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    p13N->runnum = runnum;
    p13N->T->Fill();

  }
}

void det::corr_15O()
{
  //15O -> 14N + p
  if(Correl.N14.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q15O = mass_15O - (mass_14N + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.N14.mask[0]=1;
    Correl.makeArray(1);

    float Erel_15O = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_15O - Q15O;

    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p14N->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.N14.Sol[0]->velocity);


      p14N->Egamma[i] = doppE;
      p14N->Tgamma[i] = Ceasar->added[i].time;
      p14N->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    p14N->itele1 = Correl.proton.Sol[0]->itele;
    p14N->id1 = Correl.proton.Sol[0]->iCsI;
    p14N->ifront1 = Correl.proton.Sol[0]->ifront;
    p14N->iback1 = Correl.proton.Sol[0]->iback;
    p14N->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p14N->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p14N->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p14N->et1 = Correl.proton.Sol[0]->energyTot;
    p14N->time1 = Correl.proton.Sol[0]->CsITime;
    p14N->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p14N->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p14N->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p14N->theta1 = Correl.proton.Sol[0]->theta;
    p14N->phi1 = Correl.proton.Sol[0]->phi;

    p14N->M2[0] = Correl.N14.Sol[0]->Mvect[0];
    p14N->M2[1] = Correl.N14.Sol[0]->Mvect[1];
    p14N->M2[2] = Correl.N14.Sol[0]->Mvect[2];
    p14N->et2 = Correl.N14.Sol[0]->energyTot;
    p14N->energy_p2 = Correl.N14.Sol[0]->Ekin;
    p14N->theta2 = Correl.N14.Sol[0]->theta;
    p14N->phi2 = Correl.N14.Sol[0]->phi;

    p14N->Erel = Erel_15O;
    p14N->Ex = Ex;
    p14N->Vcm = Correl.velocityCM;
    p14N->thetaCM = thetaCM;
    p14N->cos_thetaH = Correl.cos_thetaH;

    p14N->beamZ = Correl.N14.Sol[0]->beamZ;
    p14N->runnum = runnum;

    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    p14N->T->Fill();

  }
}

void det::corr_15F()
{
  //15F -> 14O + p
  if(Correl.O14.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q15F = mass_15F - (mass_14O + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.O14.mask[0]=1;
    Correl.makeArray(1);

    float Erel_15F = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_15F - Q15F;

    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p14O->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.O14.Sol[0]->velocity);


      p14O->Egamma[i] = doppE;
      p14O->Tgamma[i] = Ceasar->added[i].time;
      p14O->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    p14O->itele1 = Correl.proton.Sol[0]->itele;
    p14O->id1 = Correl.proton.Sol[0]->iCsI;
    p14O->ifront1 = Correl.proton.Sol[0]->ifront;
    p14O->iback1 = Correl.proton.Sol[0]->iback;
    p14O->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p14O->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p14O->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p14O->et1 = Correl.proton.Sol[0]->energyTot;
    p14O->time1 = Correl.proton.Sol[0]->CsITime;
    p14O->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p14O->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p14O->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p14O->theta1 = Correl.proton.Sol[0]->theta;
    p14O->phi1 = Correl.proton.Sol[0]->phi;

    p14O->M2[0] = Correl.O14.Sol[0]->Mvect[0];
    p14O->M2[1] = Correl.O14.Sol[0]->Mvect[1];
    p14O->M2[2] = Correl.O14.Sol[0]->Mvect[2];
    p14O->et2 = Correl.O14.Sol[0]->energyTot;
    p14O->energy_p2 = Correl.O14.Sol[0]->Ekin;
    p14O->theta2 = Correl.O14.Sol[0]->theta;
    p14O->phi2 = Correl.O14.Sol[0]->phi;

    p14O->Erel = Erel_15F;
    p14O->Ex = Ex;
    p14O->Vcm = Correl.velocityCM;
    p14O->thetaCM = thetaCM;
    p14O->cos_thetaH = Correl.cos_thetaH;

    p14O->beamZ = Correl.O14.Sol[0]->beamZ;
    p14O->runnum = runnum;
    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    p14O->T->Fill();

  }
}


void det::corr_16Ne()
{
  //16Ne -> 14O + 2p
  if(Correl.O14.mult == 1 && Correl.proton.mult == 2)
  {
    float const Q16Ne = mass_16Ne - (mass_14O + 2.*mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.O14.mask[0]=1;
    Correl.makeArray(1);

    float Erel_16Ne = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_16Ne - Q16Ne;

    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    pp14O->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.O14.Sol[0]->velocity);


      pp14O->Egamma[i] = doppE;
      pp14O->Tgamma[i] = Ceasar->added[i].time;
      pp14O->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    pp14O->itele1 = Correl.proton.Sol[0]->itele;
    pp14O->id1 = Correl.proton.Sol[0]->iCsI;
    pp14O->ifront1 = Correl.proton.Sol[0]->ifront;
    pp14O->iback1 = Correl.proton.Sol[0]->iback;
    pp14O->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    pp14O->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    pp14O->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    pp14O->et1 = Correl.proton.Sol[0]->energyTot;
    pp14O->time1 = Correl.proton.Sol[0]->CsITime;
    pp14O->energy_p1 = Correl.proton.Sol[0]->Ekin;
    pp14O->denergy1_R = Correl.proton.Sol[0]->denergyR;
    pp14O->energy_p1_R = Correl.proton.Sol[0]->energyR;
    pp14O->theta1 = Correl.proton.Sol[0]->theta;
    pp14O->phi1 = Correl.proton.Sol[0]->phi;


    pp14O->itele2 = Correl.proton.Sol[1]->itele;
    pp14O->id2 = Correl.proton.Sol[1]->iCsI;
    pp14O->ifront2 = Correl.proton.Sol[1]->ifront;
    pp14O->iback2 = Correl.proton.Sol[1]->iback;
    pp14O->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    pp14O->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    pp14O->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    pp14O->et2 = Correl.proton.Sol[1]->energyTot;
    pp14O->time2 = Correl.proton.Sol[1]->CsITime;
    pp14O->energy_p2 = Correl.proton.Sol[1]->Ekin;
    pp14O->denergy2_R = Correl.proton.Sol[1]->denergyR;
    pp14O->energy_p2_R = Correl.proton.Sol[1]->energyR;
    pp14O->theta2 = Correl.proton.Sol[1]->theta;
    pp14O->phi2 = Correl.proton.Sol[1]->phi;


    pp14O->M3[0] = Correl.O14.Sol[0]->Mvect[0];
    pp14O->M3[1] = Correl.O14.Sol[0]->Mvect[1];
    pp14O->M3[2] = Correl.O14.Sol[0]->Mvect[2];
    pp14O->et3 = Correl.O14.Sol[0]->energyTot;
    pp14O->energy_p3 = Correl.O14.Sol[0]->Ekin;
    pp14O->theta3 = Correl.O14.Sol[0]->theta;
    pp14O->phi3 = Correl.O14.Sol[0]->phi;

    pp14O->Erel = Erel_16Ne;
    pp14O->Ex = Ex;
    pp14O->Vcm = Correl.velocityCM;
    pp14O->thetaCM = thetaCM;
    pp14O->cos_thetaH = Correl.cos_thetaH;

    pp14O->beamZ = Correl.O14.Sol[0]->beamZ;
    pp14O->runnum = runnum;

    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    pp14O->T->Fill();

  }
}


void det::corr_17Ne()
{
  //17Ne -> 15O + 2p
  if(Correl.O15.mult == 1 && Correl.proton.mult == 2)
  {
    float const Q17Ne = mass_17Ne - (mass_15O + 2.*mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.O15.mask[0]=1;
    Correl.makeArray(1);

    float Erel_17Ne = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_17Ne - Q17Ne;

    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    pp15O->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.O15.Sol[0]->velocity);


      pp15O->Egamma[i] = doppE;
      pp15O->Tgamma[i] = Ceasar->added[i].time;
      pp15O->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    pp15O->itele1 = Correl.proton.Sol[0]->itele;
    pp15O->id1 = Correl.proton.Sol[0]->iCsI;
    pp15O->ifront1 = Correl.proton.Sol[0]->ifront;
    pp15O->iback1 = Correl.proton.Sol[0]->iback;
    pp15O->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    pp15O->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    pp15O->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    pp15O->et1 = Correl.proton.Sol[0]->energyTot;
    pp15O->time1 = Correl.proton.Sol[0]->CsITime;
    pp15O->energy_p1 = Correl.proton.Sol[0]->Ekin;
    pp15O->denergy1_R = Correl.proton.Sol[0]->denergyR;
    pp15O->energy_p1_R = Correl.proton.Sol[0]->energyR;
    pp15O->theta1 = Correl.proton.Sol[0]->theta;
    pp15O->phi1 = Correl.proton.Sol[0]->phi;

    pp15O->itele2 = Correl.proton.Sol[1]->itele;
    pp15O->id2 = Correl.proton.Sol[1]->iCsI;
    pp15O->ifront2 = Correl.proton.Sol[1]->ifront;
    pp15O->iback2 = Correl.proton.Sol[1]->iback;
    pp15O->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    pp15O->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    pp15O->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    pp15O->et2 = Correl.proton.Sol[1]->energyTot;
    pp15O->time2 = Correl.proton.Sol[1]->CsITime;
    pp15O->energy_p2 = Correl.proton.Sol[1]->Ekin;
    pp15O->denergy2_R = Correl.proton.Sol[1]->denergyR;
    pp15O->energy_p2_R = Correl.proton.Sol[1]->energyR;
    pp15O->theta2 = Correl.proton.Sol[1]->theta;
    pp15O->phi2 = Correl.proton.Sol[1]->phi;

    pp15O->M3[0] = Correl.O15.Sol[0]->Mvect[0];
    pp15O->M3[1] = Correl.O15.Sol[0]->Mvect[1];
    pp15O->M3[2] = Correl.O15.Sol[0]->Mvect[2];
    pp15O->et3 = Correl.O15.Sol[0]->energyTot;
    pp15O->energy_p3 = Correl.O15.Sol[0]->Ekin;
    pp15O->theta3 = Correl.O15.Sol[0]->theta;
    pp15O->phi3 = Correl.O15.Sol[0]->phi;

    pp15O->Erel = Erel_17Ne;
    pp15O->Ex = Ex;
    pp15O->Vcm = Correl.velocityCM;
    pp15O->thetaCM = thetaCM;
    pp15O->cos_thetaH = Correl.cos_thetaH;

    pp15O->beamZ = Correl.O15.Sol[0]->beamZ;
    pp15O->runnum = runnum;
    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    pp15O->T->Fill();

  }
}


/////////////////////////////////////////////////
//Hists and writes to trees
/////////////////////////////////////////////////

void det::corr_16F()
{
  //16F -> 15O + p
  if(Correl.O15.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q16F = mass_16F - (mass_15O + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.O15.mask[0]=1;
    Correl.makeArray(1);

    float Erel_16F = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_16F - Q16F;


    Histo_sort->Erel_16F_p15O->Fill(Erel_16F);
    Histo_sort->Ex_16F_p15O->Fill(Ex);

    Histo_sort->ThetaCM_16F_p15O->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_16F_p15O->Fill(Correl.velocityCM);
    Histo_sort->Erel_p15O_costhetaH->Fill(Erel_16F,Correl.cos_thetaH);
    Histo_sort->p15O_VCMvsErel->Fill(Correl.velocityCM,Erel_16F);

    if (fabs(Correl.cos_thetaH) <= 0.1)
    {
      Histo_sort->Erel_16F_p15O_trans->Fill(Erel_16F);
      Histo_sort->Ex_16F_p15O_trans->Fill(Ex);
    }
    if (fabs(Correl.cos_thetaH) >= 0.8)
    {
      Histo_sort->Erel_16F_p15O_long->Fill(Erel_16F);
      Histo_sort->Ex_16F_p15O_long->Fill(Ex);
    }
    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p15O->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.O15.Sol[0]->velocity);

	    Histo_sort->p15O_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p15O_gammasADDvsErel->Fill(Erel_16F,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.O15.Sol[0]->velocity);

        Histo_sort->p15O_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p15O->Egamma[i] = doppE;
      p15O->Tgamma[i] = Ceasar->added[i].time;
      p15O->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    p15O->itele1 = Correl.proton.Sol[0]->itele;
    p15O->id1 = Correl.proton.Sol[0]->iCsI;
    p15O->ifront1 = Correl.proton.Sol[0]->ifront;
    p15O->iback1 = Correl.proton.Sol[0]->iback;
    p15O->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p15O->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p15O->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p15O->et1 = Correl.proton.Sol[0]->energyTot;
    p15O->time1 = Correl.proton.Sol[0]->CsITime;
    p15O->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p15O->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p15O->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p15O->theta1 = Correl.proton.Sol[0]->theta;
    p15O->phi1 = Correl.proton.Sol[0]->phi;

    p15O->M2[0] = Correl.O15.Sol[0]->Mvect[0];
    p15O->M2[1] = Correl.O15.Sol[0]->Mvect[1];
    p15O->M2[2] = Correl.O15.Sol[0]->Mvect[2];
    p15O->et2 = Correl.O15.Sol[0]->energyTot;
    p15O->energy_p2 = Correl.O15.Sol[0]->Ekin;
    p15O->theta2 = Correl.O15.Sol[0]->theta;
    p15O->phi2 = Correl.O15.Sol[0]->phi;

    p15O->Erel = Erel_16F;
    p15O->Ex = Ex;
    p15O->Vcm = Correl.velocityCM;
    p15O->thetaCM = thetaCM;
    p15O->cos_thetaH = Correl.cos_thetaH;

    p15O->beamZ = Correl.O15.Sol[0]->beamZ;
    p15O->runnum = runnum;
    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    p15O->T->Fill();

  }
}


void det::corr_17F()
{
  //17F -> 16O + p
  if(Correl.O16.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q17F = mass_17F - (mass_16O + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.O16.mask[0]=1;
    Correl.makeArray(1);

    float Erel_17F = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_17F - Q17F;


    Histo_sort->Erel_17F_p16O->Fill(Erel_17F);
    Histo_sort->Ex_17F_p16O->Fill(Ex);

    Histo_sort->ThetaCM_17F_p16O->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_17F_p16O->Fill(Correl.velocityCM);
    Histo_sort->Erel_p16O_costhetaH->Fill(Erel_17F,Correl.cos_thetaH);
    Histo_sort->p16O_VCMvsErel->Fill(Correl.velocityCM,Erel_17F);

    if (fabs(Correl.cos_thetaH) <= 0.1)
    {
      Histo_sort->Erel_17F_p16O_trans->Fill(Erel_17F);
      Histo_sort->Ex_17F_p16O_trans->Fill(Ex);
    }
    if (fabs(Correl.cos_thetaH) >= 0.8)
    {
      Histo_sort->Erel_17F_p16O_long->Fill(Erel_17F);
      Histo_sort->Ex_17F_p16O_long->Fill(Ex);
    }
    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p16O->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.O16.Sol[0]->velocity);

	    Histo_sort->p16O_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p16O_gammasADDvsErel->Fill(Erel_17F,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.O16.Sol[0]->velocity);

        Histo_sort->p16O_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p16O->Egamma[i] = doppE;
      p16O->Tgamma[i] = Ceasar->added[i].time;
      p16O->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    p16O->itele1 = Correl.proton.Sol[0]->itele;
    p16O->id1 = Correl.proton.Sol[0]->iCsI;
    p16O->ifront1 = Correl.proton.Sol[0]->ifront;
    p16O->iback1 = Correl.proton.Sol[0]->iback;
    p16O->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p16O->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p16O->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p16O->et1 = Correl.proton.Sol[0]->energyTot;
    p16O->time1 = Correl.proton.Sol[0]->CsITime;
    p16O->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p16O->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p16O->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p16O->theta1 = Correl.proton.Sol[0]->theta;
    p16O->phi1 = Correl.proton.Sol[0]->phi;

    p16O->M2[0] = Correl.O16.Sol[0]->Mvect[0];
    p16O->M2[1] = Correl.O16.Sol[0]->Mvect[1];
    p16O->M2[2] = Correl.O16.Sol[0]->Mvect[2];
    p16O->et2 = Correl.O16.Sol[0]->energyTot;
    p16O->energy_p2 = Correl.O16.Sol[0]->Ekin;
    p16O->theta2 = Correl.O16.Sol[0]->theta;
    p16O->phi2 = Correl.O16.Sol[0]->phi;

    p16O->Erel = Erel_17F;
    p16O->Ex = Ex;
    p16O->Vcm = Correl.velocityCM;
    p16O->thetaCM = thetaCM;
    p16O->cos_thetaH = Correl.cos_thetaH;

    p16O->beamZ = Correl.O16.Sol[0]->beamZ;
    p16O->runnum = runnum;
    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    p16O->T->Fill();

  }
}


void det::corr_18Ne()
{

  //18Ne -> 16O + 2p
  if(Correl.O16.mult == 1 && Correl.proton.mult == 2)
  {
    float const Q18Ne = mass_18Ne - (mass_16O + 2.*mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.O16.mask[0]=1;
    Correl.makeArray(1);

    float Erel_18Ne = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_18Ne - Q18Ne;

    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    pp16O->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.O16.Sol[0]->velocity);


      pp16O->Egamma[i] = doppE;
      pp16O->Tgamma[i] = Ceasar->added[i].time;
      pp16O->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    pp16O->itele1 = Correl.proton.Sol[0]->itele;
    pp16O->id1 = Correl.proton.Sol[0]->iCsI;
    pp16O->ifront1 = Correl.proton.Sol[0]->ifront;
    pp16O->iback1 = Correl.proton.Sol[0]->iback;
    pp16O->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    pp16O->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    pp16O->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    pp16O->et1 = Correl.proton.Sol[0]->energyTot;
    pp16O->time1 = Correl.proton.Sol[0]->CsITime;
    pp16O->energy_p1 = Correl.proton.Sol[0]->Ekin;
    pp16O->denergy1_R = Correl.proton.Sol[0]->denergyR;
    pp16O->energy_p1_R = Correl.proton.Sol[0]->energyR;
    pp16O->theta1 = Correl.proton.Sol[0]->theta;
    pp16O->phi1 = Correl.proton.Sol[0]->phi;

    pp16O->itele2 = Correl.proton.Sol[1]->itele;
    pp16O->id2 = Correl.proton.Sol[1]->iCsI;
    pp16O->ifront2 = Correl.proton.Sol[1]->ifront;
    pp16O->iback2 = Correl.proton.Sol[1]->iback;
    pp16O->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    pp16O->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    pp16O->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    pp16O->et2 = Correl.proton.Sol[1]->energyTot;
    pp16O->time2 = Correl.proton.Sol[1]->CsITime;
    pp16O->energy_p2 = Correl.proton.Sol[1]->Ekin;
    pp16O->denergy2_R = Correl.proton.Sol[1]->denergyR;
    pp16O->energy_p2_R = Correl.proton.Sol[1]->energyR;
    pp16O->theta2 = Correl.proton.Sol[1]->theta;
    pp16O->phi2 = Correl.proton.Sol[1]->phi;

    pp16O->M3[0] = Correl.O16.Sol[0]->Mvect[0];
    pp16O->M3[1] = Correl.O16.Sol[0]->Mvect[1];
    pp16O->M3[2] = Correl.O16.Sol[0]->Mvect[2];
    pp16O->et3 = Correl.O16.Sol[0]->energyTot;
    pp16O->energy_p3 = Correl.O16.Sol[0]->Ekin;
    pp16O->theta3 = Correl.O16.Sol[0]->theta;
    pp16O->phi3 = Correl.O16.Sol[0]->phi;

    pp16O->Erel = Erel_18Ne;
    pp16O->Ex = Ex;
    pp16O->Vcm = Correl.velocityCM;
    pp16O->thetaCM = thetaCM;
    pp16O->cos_thetaH = Correl.cos_thetaH;

    pp16O->beamZ = Correl.O16.Sol[0]->beamZ;
    pp16O->runnum = runnum;
    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    pp16O->T->Fill();

  }


  //18Ne -> 17F + p
  if(Correl.F17.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q18Ne = mass_18Ne - (mass_17F + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.F17.mask[0]=1;
    Correl.makeArray(1);

    float Erel_18Ne = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_18Ne - Q18Ne;


    Histo_sort->Erel_18Ne_p17F->Fill(Erel_18Ne);
    Histo_sort->Ex_18Ne_p17F->Fill(Ex);

    Histo_sort->ThetaCM_18Ne_p17F->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_18Ne_p17F->Fill(Correl.velocityCM);
    Histo_sort->Erel_p17F_costhetaH->Fill(Erel_18Ne,Correl.cos_thetaH);
    Histo_sort->p17F_VCMvsErel->Fill(Correl.velocityCM,Erel_18Ne);

    if (fabs(Correl.cos_thetaH) <= 0.1)
    {
      Histo_sort->Erel_18Ne_p17F_trans->Fill(Erel_18Ne);
      Histo_sort->Ex_18Ne_p17F_trans->Fill(Ex);
    }
    if (fabs(Correl.cos_thetaH) >= 0.8)
    {
      Histo_sort->Erel_18Ne_p17F_long->Fill(Erel_18Ne);
      Histo_sort->Ex_18Ne_p17F_long->Fill(Ex);
    }
    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p17F->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.F17.Sol[0]->velocity);

	    Histo_sort->p17F_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p17F_gammasADDvsErel->Fill(Erel_18Ne,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.F17.Sol[0]->velocity);

        Histo_sort->p17F_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p17F->Egamma[i] = doppE;
      p17F->Tgamma[i] = Ceasar->added[i].time;
      p17F->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    p17F->itele1 = Correl.proton.Sol[0]->itele;
    p17F->id1 = Correl.proton.Sol[0]->iCsI;
    p17F->ifront1 = Correl.proton.Sol[0]->ifront;
    p17F->iback1 = Correl.proton.Sol[0]->iback;
    p17F->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p17F->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p17F->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p17F->et1 = Correl.proton.Sol[0]->energyTot;
    p17F->time1 = Correl.proton.Sol[0]->CsITime;
    p17F->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p17F->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p17F->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p17F->theta1 = Correl.proton.Sol[0]->theta;
    p17F->phi1 = Correl.proton.Sol[0]->phi;

    p17F->M2[0] = Correl.F17.Sol[0]->Mvect[0];
    p17F->M2[1] = Correl.F17.Sol[0]->Mvect[1];
    p17F->M2[2] = Correl.F17.Sol[0]->Mvect[2];
    p17F->et2 = Correl.F17.Sol[0]->energyTot;
    p17F->energy_p2 = Correl.F17.Sol[0]->Ekin;
    p17F->theta2 = Correl.F17.Sol[0]->theta;
    p17F->phi2 = Correl.F17.Sol[0]->phi;
    p17F->theta2 = Correl.F17.Sol[0]->theta;
    p17F->phi2 = Correl.F17.Sol[0]->phi;

    p17F->Erel = Erel_18Ne;
    p17F->Ex = Ex;
    p17F->Vcm = Correl.velocityCM;
    p17F->thetaCM = thetaCM;
    p17F->cos_thetaH = Correl.cos_thetaH;

    p17F->runnum = runnum;
    p17F->beamZ = Correl.F17.Sol[0]->beamZ;

    p17F->theta_s800 = Correl.F17.Sol[0]->theta_s800;
    p17F->phi_s800 = Correl.F17.Sol[0]->phi_s800;

    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    p17F->T->Fill();

  }
}


void det::corr_19Ne()
{
  //19Ne -> 18F + p
  if(Correl.F18.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q19Ne = mass_19Ne - (mass_18F + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.F18.mask[0]=1;
    Correl.makeArray(1);

    float Erel_19Ne = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_19Ne - Q19Ne;


    Histo_sort->Erel_19Ne_p18F->Fill(Erel_19Ne);
    Histo_sort->Ex_19Ne_p18F->Fill(Ex);

    Histo_sort->ThetaCM_19Ne_p18F->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_19Ne_p18F->Fill(Correl.velocityCM);
    Histo_sort->Erel_p18F_costhetaH->Fill(Erel_19Ne,Correl.cos_thetaH);
    Histo_sort->p18F_VCMvsErel->Fill(Correl.velocityCM,Erel_19Ne);

    if (fabs(Correl.cos_thetaH) <= 0.1)
    {
      Histo_sort->Erel_19Ne_p18F_trans->Fill(Erel_19Ne);
      Histo_sort->Ex_19Ne_p18F_trans->Fill(Ex);
    }
    if (fabs(Correl.cos_thetaH) >= 0.8)
    {
      Histo_sort->Erel_19Ne_p18F_long->Fill(Erel_19Ne);
      Histo_sort->Ex_19Ne_p18F_long->Fill(Ex);
    }
    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p18F->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.F18.Sol[0]->velocity);

	    Histo_sort->p18F_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p18F_gammasADDvsErel->Fill(Erel_19Ne,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.F18.Sol[0]->velocity);

        Histo_sort->p18F_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p18F->Egamma[i] = doppE;
      p18F->Tgamma[i] = Ceasar->added[i].time;
      p18F->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    p18F->itele1 = Correl.proton.Sol[0]->itele;
    p18F->id1 = Correl.proton.Sol[0]->iCsI;
    p18F->ifront1 = Correl.proton.Sol[0]->ifront;
    p18F->iback1 = Correl.proton.Sol[0]->iback;
    p18F->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p18F->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p18F->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p18F->et1 = Correl.proton.Sol[0]->energyTot;
    p18F->time1 = Correl.proton.Sol[0]->CsITime;
    p18F->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p18F->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p18F->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p18F->theta1 = Correl.proton.Sol[0]->theta;
    p18F->phi1 = Correl.proton.Sol[0]->phi;

    p18F->M2[0] = Correl.F18.Sol[0]->Mvect[0];
    p18F->M2[1] = Correl.F18.Sol[0]->Mvect[1];
    p18F->M2[2] = Correl.F18.Sol[0]->Mvect[2];
    p18F->et2 = Correl.F18.Sol[0]->energyTot;
    p18F->energy_p2 = Correl.F18.Sol[0]->Ekin;
    p18F->theta2 = Correl.F18.Sol[0]->theta;
    p18F->phi2 = Correl.F18.Sol[0]->phi;

    p18F->Erel = Erel_19Ne;
    p18F->Ex = Ex;
    p18F->Vcm = Correl.velocityCM;
    p18F->thetaCM = thetaCM;
    p18F->cos_thetaH = Correl.cos_thetaH;

    p18F->beamZ = Correl.F18.Sol[0]->beamZ;
    p18F->runnum = runnum;
    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    p18F->T->Fill();

  }
}

void det::corr_19Na()
{

  //19Na -> 17F + 2p
  if(Correl.F17.mult == 1 && Correl.proton.mult == 2)
  {
    float const Q19Na = mass_19Na - (mass_17F + 2.*mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.F17.mask[0]=1;
    Correl.makeArray(1);

    float Erel_19Na = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_19Na - Q19Na;

    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    pp17F->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.F17.Sol[0]->velocity);


      pp17F->Egamma[i] = doppE;
      pp17F->Tgamma[i] = Ceasar->added[i].time;
      pp17F->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    pp17F->itele1 = Correl.proton.Sol[0]->itele;
    pp17F->id1 = Correl.proton.Sol[0]->iCsI;
    pp17F->ifront1 = Correl.proton.Sol[0]->ifront;
    pp17F->iback1 = Correl.proton.Sol[0]->iback;
    pp17F->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    pp17F->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    pp17F->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    pp17F->et1 = Correl.proton.Sol[0]->energyTot;
    pp17F->time1 = Correl.proton.Sol[0]->CsITime;
    pp17F->energy_p1 = Correl.proton.Sol[0]->Ekin;
    pp17F->denergy1_R = Correl.proton.Sol[0]->denergyR;
    pp17F->energy_p1_R = Correl.proton.Sol[0]->energyR;
    pp17F->theta1 = Correl.proton.Sol[0]->theta;
    pp17F->phi1 = Correl.proton.Sol[0]->phi;


    pp17F->itele2 = Correl.proton.Sol[1]->itele;
    pp17F->id2 = Correl.proton.Sol[1]->iCsI;
    pp17F->ifront2 = Correl.proton.Sol[1]->ifront;
    pp17F->iback2 = Correl.proton.Sol[1]->iback;
    pp17F->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    pp17F->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    pp17F->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    pp17F->et2 = Correl.proton.Sol[1]->energyTot;
    pp17F->time2 = Correl.proton.Sol[1]->CsITime;
    pp17F->energy_p2 = Correl.proton.Sol[1]->Ekin;
    pp17F->denergy2_R = Correl.proton.Sol[1]->denergyR;
    pp17F->energy_p2_R = Correl.proton.Sol[1]->energyR;
    pp17F->theta2 = Correl.proton.Sol[1]->theta;
    pp17F->phi2 = Correl.proton.Sol[1]->phi;

    pp17F->M3[0] = Correl.F17.Sol[0]->Mvect[0];
    pp17F->M3[1] = Correl.F17.Sol[0]->Mvect[1];
    pp17F->M3[2] = Correl.F17.Sol[0]->Mvect[2];
    pp17F->et3 = Correl.F17.Sol[0]->energyTot;
    pp17F->energy_p3 = Correl.F17.Sol[0]->Ekin;
    pp17F->theta3 = Correl.F17.Sol[0]->theta;
    pp17F->phi3 = Correl.F17.Sol[0]->phi;

    pp17F->Erel = Erel_19Na;
    pp17F->Ex = Ex;
    pp17F->Vcm = Correl.velocityCM;
    pp17F->thetaCM = thetaCM;
    pp17F->cos_thetaH = Correl.cos_thetaH;

    pp17F->beamZ = Correl.F17.Sol[0]->beamZ;
    pp17F->runnum = runnum;
    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    pp17F->T->Fill();

  }

  //19Na -> 18Ne + p
  if(Correl.Ne18.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q19Na = mass_19Na - (mass_18Ne + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Ne18.mask[0]=1;
    Correl.makeArray(1);

    float Erel_19Na = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_19Na - Q19Na;


    Histo_sort->Erel_19Na_p18Ne->Fill(Erel_19Na);
    Histo_sort->Ex_19Na_p18Ne->Fill(Ex);

    Histo_sort->ThetaCM_19Na_p18Ne->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_19Na_p18Ne->Fill(Correl.velocityCM);
    Histo_sort->Erel_p18Ne_costhetaH->Fill(Erel_19Na,Correl.cos_thetaH);
    Histo_sort->p18Ne_VCMvsErel->Fill(Correl.velocityCM,Erel_19Na);

    if (fabs(Correl.cos_thetaH) <= 0.1)
    {
      Histo_sort->Erel_19Na_p18Ne_trans->Fill(Erel_19Na);
      Histo_sort->Ex_19Na_p18Ne_trans->Fill(Ex);
    }
    if (fabs(Correl.cos_thetaH) >= 0.8)
    {
      Histo_sort->Erel_19Na_p18Ne_long->Fill(Erel_19Na);
      Histo_sort->Ex_19Na_p18Ne_long->Fill(Ex);
    }
    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p18Ne->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Ne18.Sol[0]->velocity);

	    Histo_sort->p18Ne_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p18Ne_gammasADDvsErel->Fill(Erel_19Na,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Ne18.Sol[0]->velocity);

        Histo_sort->p18Ne_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->p18Ne_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }

      p18Ne->Egamma[i] = doppE;
      p18Ne->Tgamma[i] = Ceasar->added[i].time;
      p18Ne->Chgamma[i] = Ceasar->added[i].id;
    }

    //Run dependent histograms. Set 1 had gobbi farther away than set 1a
    if (runnum < 215) Histo_sort->Erel_19Na_p18Ne_set1->Fill(Erel_19Na);
    else Histo_sort->Erel_19Na_p18Ne_set1a->Fill(Erel_19Na);

    //Set tree variables
    p18Ne->itele1 = Correl.proton.Sol[0]->itele;
    p18Ne->id1 = Correl.proton.Sol[0]->iCsI;
    p18Ne->ifront1 = Correl.proton.Sol[0]->ifront;
    p18Ne->iback1 = Correl.proton.Sol[0]->iback;
    p18Ne->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p18Ne->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p18Ne->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p18Ne->et1 = Correl.proton.Sol[0]->energyTot;
    p18Ne->time1 = Correl.proton.Sol[0]->CsITime;
    p18Ne->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p18Ne->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p18Ne->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p18Ne->theta1 = Correl.proton.Sol[0]->theta;
    p18Ne->phi1 = Correl.proton.Sol[0]->phi;

    p18Ne->M2[0] = Correl.Ne18.Sol[0]->Mvect[0];
    p18Ne->M2[1] = Correl.Ne18.Sol[0]->Mvect[1];
    p18Ne->M2[2] = Correl.Ne18.Sol[0]->Mvect[2];
    p18Ne->et2 = Correl.Ne18.Sol[0]->energyTot;
    p18Ne->energy_p2 = Correl.Ne18.Sol[0]->Ekin;
    p18Ne->theta2 = Correl.Ne18.Sol[0]->theta;
    p18Ne->phi2 = Correl.Ne18.Sol[0]->phi;

    p18Ne->Erel = Erel_19Na;
    p18Ne->Ex = Ex;
    p18Ne->Vcm = Correl.velocityCM;
    p18Ne->thetaCM = thetaCM;
    p18Ne->cos_thetaH = Correl.cos_thetaH;

    p18Ne->runnum = runnum;
    p18Ne->beamZ = Correl.Ne18.Sol[0]->beamZ;

    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    p18Ne->T->Fill();

  }
}

void det::corr_20Na()
{

  //20Na -> 18F + 2p
  if(Correl.F18.mult == 1 && Correl.proton.mult == 2)
  {
    float const Q20Na = mass_20Na - (mass_18F + 2.*mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.F18.mask[0]=1;
    Correl.makeArray(1);

    float Erel_20Na = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_20Na - Q20Na;

    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    pp18F->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.F18.Sol[0]->velocity);


      pp18F->Egamma[i] = doppE;
      pp18F->Tgamma[i] = Ceasar->added[i].time;
      pp18F->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    pp18F->itele1 = Correl.proton.Sol[0]->itele;
    pp18F->id1 = Correl.proton.Sol[0]->iCsI;
    pp18F->ifront1 = Correl.proton.Sol[0]->ifront;
    pp18F->iback1 = Correl.proton.Sol[0]->iback;
    pp18F->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    pp18F->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    pp18F->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    pp18F->et1 = Correl.proton.Sol[0]->energyTot;
    pp18F->time1 = Correl.proton.Sol[0]->CsITime;
    pp18F->energy_p1 = Correl.proton.Sol[0]->Ekin;
    pp18F->denergy1_R = Correl.proton.Sol[0]->denergyR;
    pp18F->energy_p1_R = Correl.proton.Sol[0]->energyR;
    pp18F->theta1 = Correl.proton.Sol[0]->theta;
    pp18F->phi1 = Correl.proton.Sol[0]->phi;

    pp18F->itele2 = Correl.proton.Sol[1]->itele;
    pp18F->id2 = Correl.proton.Sol[1]->iCsI;
    pp18F->ifront2 = Correl.proton.Sol[1]->ifront;
    pp18F->iback2 = Correl.proton.Sol[1]->iback;
    pp18F->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    pp18F->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    pp18F->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    pp18F->et2 = Correl.proton.Sol[1]->energyTot;
    pp18F->time2 = Correl.proton.Sol[1]->CsITime;
    pp18F->energy_p2 = Correl.proton.Sol[1]->Ekin;
    pp18F->denergy2_R = Correl.proton.Sol[1]->denergyR;
    pp18F->energy_p2_R = Correl.proton.Sol[1]->energyR;
    pp18F->theta2 = Correl.proton.Sol[1]->theta;
    pp18F->phi2 = Correl.proton.Sol[1]->phi;

    pp18F->M3[0] = Correl.F18.Sol[0]->Mvect[0];
    pp18F->M3[1] = Correl.F18.Sol[0]->Mvect[1];
    pp18F->M3[2] = Correl.F18.Sol[0]->Mvect[2];
    pp18F->et3 = Correl.F18.Sol[0]->energyTot;
    pp18F->energy_p3 = Correl.F18.Sol[0]->Ekin;
    pp18F->theta3 = Correl.F18.Sol[0]->theta;
    pp18F->phi3 = Correl.F18.Sol[0]->phi;

    pp18F->Erel = Erel_20Na;
    pp18F->Ex = Ex;
    pp18F->Vcm = Correl.velocityCM;
    pp18F->thetaCM = thetaCM;
    pp18F->cos_thetaH = Correl.cos_thetaH;

    pp18F->beamZ = Correl.F18.Sol[0]->beamZ;
    pp18F->runnum = runnum;
    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    pp18F->T->Fill();

  }


  //19Na -> 18Ne + p
  if(Correl.Ne19.mult == 1 && Correl.proton.mult == 1)
  {

    float const Q20Na = mass_20Na - (mass_19Ne + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Ne19.mask[0]=1;
    Correl.makeArray(1);

    float Erel_20Na = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_20Na - Q20Na;


    Histo_sort->Erel_20Na_p19Ne->Fill(Erel_20Na);
    Histo_sort->Ex_20Na_p19Ne->Fill(Ex);

    Histo_sort->ThetaCM_20Na_p19Ne->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_20Na_p19Ne->Fill(Correl.velocityCM);
    Histo_sort->Erel_p19Ne_costhetaH->Fill(Erel_20Na,Correl.cos_thetaH);
    Histo_sort->p19Ne_VCMvsErel->Fill(Correl.velocityCM,Erel_20Na);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_20Na_p19Ne_trans->Fill(Erel_20Na);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p19Ne->Ngamma = Ceasar->Nadded; //Num gammas for tree
  
    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Ne19.Sol[0]->velocity);

	    Histo_sort->p19Ne_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p19Ne_gammasADDvsErel->Fill(Erel_20Na,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Ne19.Sol[0]->velocity);

        Histo_sort->p19Ne_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p19Ne->Egamma[i] = doppE;
      p19Ne->Tgamma[i] = Ceasar->added[i].time;
      p19Ne->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->p19Ne_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }

    }

    //Set tree variables
    p19Ne->itele1 = Correl.proton.Sol[0]->itele;
    p19Ne->id1 = Correl.proton.Sol[0]->iCsI;
    p19Ne->ifront1 = Correl.proton.Sol[0]->ifront;
    p19Ne->iback1 = Correl.proton.Sol[0]->iback;
    p19Ne->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p19Ne->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p19Ne->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p19Ne->et1 = Correl.proton.Sol[0]->energyTot;
    p19Ne->time1 = Correl.proton.Sol[0]->CsITime;
    p19Ne->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p19Ne->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p19Ne->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p19Ne->theta1 = Correl.proton.Sol[0]->theta;
    p19Ne->phi1 = Correl.proton.Sol[0]->phi;

    p19Ne->M2[0] = Correl.Ne19.Sol[0]->Mvect[0];
    p19Ne->M2[1] = Correl.Ne19.Sol[0]->Mvect[1];
    p19Ne->M2[2] = Correl.Ne19.Sol[0]->Mvect[2];
    p19Ne->et2 = Correl.Ne19.Sol[0]->energyTot;
    p19Ne->energy_p2 = Correl.Ne19.Sol[0]->Ekin;
    p19Ne->theta2 = Correl.Ne19.Sol[0]->theta;
    p19Ne->phi2 = Correl.Ne19.Sol[0]->phi;

    p19Ne->Erel = Erel_20Na;
    p19Ne->Ex = Ex;
    p19Ne->Vcm = Correl.velocityCM;
    p19Ne->thetaCM = thetaCM;
    p19Ne->cos_thetaH = Correl.cos_thetaH;

    p19Ne->beamZ = Correl.Ne19.Sol[0]->beamZ;
    p19Ne->runnum = runnum;
    p19Ne->T->Fill();

  }
}

void det::corr_21Na()
{
  //21Na -> 20Ne + p
  if(Correl.Ne20.mult == 1 && Correl.proton.mult == 1)
  {

    float const Q21Na = mass_21Na - (mass_20Ne + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Ne20.mask[0]=1;
    Correl.makeArray(1);

    float Erel_21Na = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_21Na - Q21Na;


    Histo_sort->Erel_21Na_p20Ne->Fill(Erel_21Na);
    Histo_sort->Ex_21Na_p20Ne->Fill(Ex);

    Histo_sort->ThetaCM_21Na_p20Ne->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_21Na_p20Ne->Fill(Correl.velocityCM);
    Histo_sort->Erel_p20Ne_costhetaH->Fill(Erel_21Na,Correl.cos_thetaH);
    Histo_sort->p20Ne_VCMvsErel->Fill(Correl.velocityCM,Erel_21Na);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_21Na_p20Ne_trans->Fill(Erel_21Na);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p20Ne->Ngamma = Ceasar->Nadded; //Num gammas for tree
  
    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Ne20.Sol[0]->velocity);

	    Histo_sort->p20Ne_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p20Ne_gammasADDvsErel->Fill(Erel_21Na,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Ne20.Sol[0]->velocity);

        Histo_sort->p20Ne_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p20Ne->Egamma[i] = doppE;
      p20Ne->Tgamma[i] = Ceasar->added[i].time;
      p20Ne->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->p20Ne_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }

    }

    //Set tree variables
    p20Ne->itele1 = Correl.proton.Sol[0]->itele;
    p20Ne->id1 = Correl.proton.Sol[0]->iCsI;
    p20Ne->ifront1 = Correl.proton.Sol[0]->ifront;
    p20Ne->iback1 = Correl.proton.Sol[0]->iback;
    p20Ne->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p20Ne->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p20Ne->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p20Ne->et1 = Correl.proton.Sol[0]->energyTot;
    p20Ne->time1 = Correl.proton.Sol[0]->CsITime;
    p20Ne->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p20Ne->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p20Ne->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p20Ne->theta1 = Correl.proton.Sol[0]->theta;
    p20Ne->phi1 = Correl.proton.Sol[0]->phi;

    p20Ne->M2[0] = Correl.Ne20.Sol[0]->Mvect[0];
    p20Ne->M2[1] = Correl.Ne20.Sol[0]->Mvect[1];
    p20Ne->M2[2] = Correl.Ne20.Sol[0]->Mvect[2];
    p20Ne->et2 = Correl.Ne20.Sol[0]->energyTot;
    p20Ne->energy_p2 = Correl.Ne20.Sol[0]->Ekin;
    p20Ne->theta2 = Correl.Ne20.Sol[0]->theta;
    p20Ne->phi2 = Correl.Ne20.Sol[0]->phi;

    p20Ne->Erel = Erel_21Na;
    p20Ne->Ex = Ex;
    p20Ne->Vcm = Correl.velocityCM;
    p20Ne->thetaCM = thetaCM;
    p20Ne->cos_thetaH = Correl.cos_thetaH;

    p20Ne->beamZ = Correl.Ne20.Sol[0]->beamZ;
    p20Ne->runnum = runnum;
    p20Ne->T->Fill();

  }
}


void det::corr_19Mg()
{
  //19Mg -> 17Ne + 2p
  if(Correl.Ne17.mult == 1 && Correl.proton.mult == 2)
  {
    float const Q19Mg = mass_19Mg - (mass_17Ne + 2*mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.Ne17.mask[0]=1;
    Correl.makeArray(1);

    float Erel_19Mg = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_19Mg - Q19Mg;


    Histo_sort->Erel_19Mg_2p17Ne->Fill(Erel_19Mg);
    Histo_sort->Ex_19Mg_2p17Ne->Fill(Ex);

    Histo_sort->ThetaCM_19Mg_2p17Ne->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_19Mg_2p17Ne->Fill(Correl.velocityCM);
    Histo_sort->Erel_2p17Ne_costhetaH->Fill(Erel_19Mg,Correl.cos_thetaH);
    Histo_sort->pp17Ne_VCMvsErel->Fill(Correl.velocityCM,Erel_19Mg);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_19Mg_2p17Ne_trans->Fill(Erel_19Mg);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    pp17Ne->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Ne17.Sol[0]->velocity);

	    Histo_sort->pp17Ne_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->pp17Ne_gammasADDvsErel->Fill(Erel_19Mg,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Ne17.Sol[0]->velocity);

        Histo_sort->pp17Ne_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      pp17Ne->Egamma[i] = doppE;
      pp17Ne->Tgamma[i] = Ceasar->added[i].time;
      pp17Ne->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->pp17Ne_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //Set tree variables
    pp17Ne->itele1 = Correl.proton.Sol[0]->itele;
    pp17Ne->id1 = Correl.proton.Sol[0]->iCsI;
    pp17Ne->ifront1 = Correl.proton.Sol[0]->ifront;
    pp17Ne->iback1 = Correl.proton.Sol[0]->iback;
    pp17Ne->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    pp17Ne->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    pp17Ne->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    pp17Ne->et1 = Correl.proton.Sol[0]->energyTot;
    pp17Ne->time1 = Correl.proton.Sol[0]->CsITime;
    pp17Ne->energy_p1 = Correl.proton.Sol[0]->Ekin;
    pp17Ne->denergy1_R = Correl.proton.Sol[0]->denergyR;
    pp17Ne->energy_p1_R = Correl.proton.Sol[0]->energyR;
    pp17Ne->theta1 = Correl.proton.Sol[0]->theta;
    pp17Ne->phi1 = Correl.proton.Sol[0]->phi;

    pp17Ne->itele2 = Correl.proton.Sol[1]->itele;
    pp17Ne->id2 = Correl.proton.Sol[1]->iCsI;
    pp17Ne->ifront2 = Correl.proton.Sol[1]->ifront;
    pp17Ne->iback2 = Correl.proton.Sol[1]->iback;
    pp17Ne->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    pp17Ne->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    pp17Ne->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    pp17Ne->et2 = Correl.proton.Sol[1]->energyTot;
    pp17Ne->time2 = Correl.proton.Sol[1]->CsITime;
    pp17Ne->energy_p2 = Correl.proton.Sol[1]->Ekin;
    pp17Ne->denergy2_R = Correl.proton.Sol[1]->denergyR;
    pp17Ne->energy_p2_R = Correl.proton.Sol[1]->energyR;
    pp17Ne->theta2 = Correl.proton.Sol[1]->theta;
    pp17Ne->phi2 = Correl.proton.Sol[1]->phi;

    pp17Ne->M3[0] = Correl.Ne17.Sol[0]->Mvect[0];
    pp17Ne->M3[1] = Correl.Ne17.Sol[0]->Mvect[1];
    pp17Ne->M3[2] = Correl.Ne17.Sol[0]->Mvect[2];
    pp17Ne->et3 = Correl.Ne17.Sol[0]->energyTot;
    pp17Ne->energy_p3 = Correl.Ne17.Sol[0]->Ekin;
    pp17Ne->theta3 = Correl.Ne17.Sol[0]->theta;
    pp17Ne->phi3 = Correl.Ne17.Sol[0]->phi;

    pp17Ne->Erel = Erel_19Mg;
    pp17Ne->Ex = Ex;
    pp17Ne->Vcm = Correl.velocityCM;
    pp17Ne->thetaCM = thetaCM;
    pp17Ne->cos_thetaH = Correl.cos_thetaH;

    pp17Ne->beamZ = Correl.Ne17.Sol[0]->beamZ;
    pp17Ne->runnum = runnum;
    //cout << Correl.cos_thetaH << endl;

    pp17Ne->T->Fill();
  }
}

void det::corr_20Mg()
{
  //20Mg -> 18Ne + 2p
  if(Correl.Ne18.mult == 1 && Correl.proton.mult == 2)
  {
    float const Q20Mg = mass_20Mg - (mass_18Ne + 2*mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.Ne18.mask[0]=1;
    Correl.makeArray(1);

    float Erel_20Mg = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_20Mg - Q20Mg;


    Histo_sort->Erel_20Mg_2p18Ne->Fill(Erel_20Mg);
    Histo_sort->Ex_20Mg_2p18Ne->Fill(Ex);

    Histo_sort->ThetaCM_20Mg_2p18Ne->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_20Mg_2p18Ne->Fill(Correl.velocityCM);
    Histo_sort->Erel_2p18Ne_costhetaH->Fill(Erel_20Mg,Correl.cos_thetaH);
    Histo_sort->pp18Ne_VCMvsErel->Fill(Correl.velocityCM,Erel_20Mg);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_20Mg_2p18Ne_trans->Fill(Erel_20Mg);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    pp18Ne->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Ne18.Sol[0]->velocity);

	    Histo_sort->pp18Ne_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->pp18Ne_gammasADDvsErel->Fill(Erel_20Mg,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Ne18.Sol[0]->velocity);

        Histo_sort->pp18Ne_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      pp18Ne->Egamma[i] = doppE;
      pp18Ne->Tgamma[i] = Ceasar->added[i].time;
      pp18Ne->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->pp18Ne_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //cout << Correl.cos_thetaH << endl;

    //Set tree variables
    pp18Ne->itele1 = Correl.proton.Sol[0]->itele;
    pp18Ne->id1 = Correl.proton.Sol[0]->iCsI;
    pp18Ne->ifront1 = Correl.proton.Sol[0]->ifront;
    pp18Ne->iback1 = Correl.proton.Sol[0]->iback;
    pp18Ne->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    pp18Ne->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    pp18Ne->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    pp18Ne->et1 = Correl.proton.Sol[0]->energyTot;
    pp18Ne->time1 = Correl.proton.Sol[0]->CsITime;
    pp18Ne->energy_p1 = Correl.proton.Sol[0]->Ekin;
    pp18Ne->denergy1_R = Correl.proton.Sol[0]->denergyR;
    pp18Ne->energy_p1_R = Correl.proton.Sol[0]->energyR;
    pp18Ne->theta1 = Correl.proton.Sol[0]->theta;
    pp18Ne->phi1 = Correl.proton.Sol[0]->phi;

    pp18Ne->itele2 = Correl.proton.Sol[1]->itele;
    pp18Ne->id2 = Correl.proton.Sol[1]->iCsI;
    pp18Ne->ifront2 = Correl.proton.Sol[1]->ifront;
    pp18Ne->iback2 = Correl.proton.Sol[1]->iback;
    pp18Ne->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    pp18Ne->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    pp18Ne->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    pp18Ne->et2 = Correl.proton.Sol[1]->energyTot;
    pp18Ne->time2 = Correl.proton.Sol[1]->CsITime;
    pp18Ne->energy_p2 = Correl.proton.Sol[1]->Ekin;
    pp18Ne->denergy2_R = Correl.proton.Sol[1]->denergyR;
    pp18Ne->energy_p2_R = Correl.proton.Sol[1]->energyR;
    pp18Ne->theta2 = Correl.proton.Sol[1]->theta;
    pp18Ne->phi2 = Correl.proton.Sol[1]->phi;

    pp18Ne->M3[0] = Correl.Ne18.Sol[0]->Mvect[0];
    pp18Ne->M3[1] = Correl.Ne18.Sol[0]->Mvect[1];
    pp18Ne->M3[2] = Correl.Ne18.Sol[0]->Mvect[2];
    pp18Ne->et3 = Correl.Ne18.Sol[0]->energyTot;
    pp18Ne->energy_p3 = Correl.Ne18.Sol[0]->Ekin;
    pp18Ne->theta3 = Correl.Ne18.Sol[0]->theta;
    pp18Ne->phi3 = Correl.Ne18.Sol[0]->phi;

    pp18Ne->Erel = Erel_20Mg;
    pp18Ne->Ex = Ex;
    pp18Ne->Vcm = Correl.velocityCM;
    pp18Ne->thetaCM = thetaCM;
    pp18Ne->cos_thetaH = Correl.cos_thetaH;

    pp18Ne->beamZ = Correl.Ne18.Sol[0]->beamZ;
    pp18Ne->runnum = runnum;
    pp18Ne->T->Fill();
  }
}

void det::corr_21Mg()
{
  //21Mg -> 20Na + p
  if(Correl.Na20.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q21Mg = mass_21Mg - (mass_20Na + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Na20.mask[0]=1;
    Correl.makeArray(1);

    float Erel_21Mg = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_21Mg - Q21Mg;


    Histo_sort->Erel_21Mg_p20Na->Fill(Erel_21Mg);
    Histo_sort->Ex_21Mg_p20Na->Fill(Ex);

    Histo_sort->ThetaCM_21Mg_p20Na->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_21Mg_p20Na->Fill(Correl.velocityCM);
    Histo_sort->Erel_p20Na_costhetaH->Fill(Erel_21Mg,Correl.cos_thetaH);
    Histo_sort->p20Na_VCMvsErel->Fill(Correl.velocityCM,Erel_21Mg);

    if (fabs(Correl.cos_thetaH) <= 0.1)
    {
      Histo_sort->Erel_21Mg_p20Na_trans->Fill(Erel_21Mg);
      Histo_sort->Ex_21Mg_p20Na_trans->Fill(Ex);
    }
    if (fabs(Correl.cos_thetaH) >= 0.8)
    {
      Histo_sort->Erel_21Mg_p20Na_long->Fill(Erel_21Mg);
      Histo_sort->Ex_21Mg_p20Na_long->Fill(Ex);
    }

    Histo_sort->Mg21_pKE->Fill(Correl.proton.Sol[0]->Ekin);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p20Na->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Na20.Sol[0]->velocity);

	    Histo_sort->p20Na_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p20Na_gammasADDvsErel->Fill(Erel_21Mg,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Na20.Sol[0]->velocity);

        Histo_sort->p20Na_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p20Na->Egamma[i] = doppE;
      p20Na->Tgamma[i] = Ceasar->added[i].time;
      p20Na->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->p20Na_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }
  
    //Run dependent histograms. Set 1 had gobbi farther away than set 1a
    if (runnum < 215) Histo_sort->Erel_21Mg_p20Na_set1->Fill(Erel_21Mg);
    else Histo_sort->Erel_21Mg_p20Na_set1a->Fill(Erel_21Mg);

    //cout << Correl.cos_thetaH << endl;

    //Set tree variables
    p20Na->itele1 = Correl.proton.Sol[0]->itele;
    p20Na->id1 = Correl.proton.Sol[0]->iCsI;
    p20Na->ifront1 = Correl.proton.Sol[0]->ifront;
    p20Na->iback1 = Correl.proton.Sol[0]->iback;
    p20Na->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p20Na->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p20Na->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p20Na->et1 = Correl.proton.Sol[0]->energyTot;
    p20Na->time1 = Correl.proton.Sol[0]->CsITime;
    p20Na->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p20Na->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p20Na->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p20Na->theta1 = Correl.proton.Sol[0]->theta;
    p20Na->phi1 = Correl.proton.Sol[0]->phi;

    p20Na->M2[0] = Correl.Na20.Sol[0]->Mvect[0];
    p20Na->M2[1] = Correl.Na20.Sol[0]->Mvect[1];
    p20Na->M2[2] = Correl.Na20.Sol[0]->Mvect[2];
    p20Na->et2 = Correl.Na20.Sol[0]->energyTot;
    p20Na->energy_p2 = Correl.Na20.Sol[0]->Ekin;
    p20Na->theta2 = Correl.Na20.Sol[0]->theta;
    p20Na->phi2 = Correl.Na20.Sol[0]->phi;

    p20Na->Erel = Erel_21Mg;
    p20Na->Ex = Ex;
    p20Na->Vcm = Correl.velocityCM;
    p20Na->thetaCM = thetaCM;
    p20Na->cos_thetaH = Correl.cos_thetaH;

    p20Na->beamZ = Correl.Na20.Sol[0]->beamZ;
    p20Na->runnum = runnum;
    p20Na->T->Fill();
  }

  //21Mg -> 19Ne + 2p
  if(Correl.Ne19.mult == 1 && Correl.proton.mult == 2)
  {
    float const Q21Mg = mass_21Mg - (mass_19Ne + 2*mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.Ne19.mask[0]=1;
    Correl.makeArray(1);

    float Erel_21Mg = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_21Mg - Q21Mg;


    Histo_sort->Erel_21Mg_2p19Ne->Fill(Erel_21Mg);
    Histo_sort->Ex_21Mg_2p19Ne->Fill(Ex);

    Histo_sort->ThetaCM_21Mg_2p19Ne->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_21Mg_2p19Ne->Fill(Correl.velocityCM);
    Histo_sort->Erel_2p19Ne_costhetaH->Fill(Erel_21Mg,Correl.cos_thetaH);
    Histo_sort->pp19Ne_VCMvsErel->Fill(Correl.velocityCM,Erel_21Mg);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_21Mg_2p19Ne_trans->Fill(Erel_21Mg);

    Histo_sort->Mg21_2p_pKE->Fill(Correl.proton.Sol[0]->Ekin);
    Histo_sort->Mg21_2p_pKE->Fill(Correl.proton.Sol[1]->Ekin);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    pp19Ne->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Ne19.Sol[0]->velocity);

	    Histo_sort->pp19Ne_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->pp19Ne_gammasADDvsErel->Fill(Erel_21Mg,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Ne19.Sol[0]->velocity);

        Histo_sort->pp19Ne_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      pp19Ne->Egamma[i] = doppE;
      pp19Ne->Tgamma[i] = Ceasar->added[i].time;
      pp19Ne->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->pp19Ne_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //cout << Correl.cos_thetaH << endl;

    //Set tree variables
    pp19Ne->itele1 = Correl.proton.Sol[0]->itele;
    pp19Ne->id1 = Correl.proton.Sol[0]->iCsI;
    pp19Ne->ifront1 = Correl.proton.Sol[0]->ifront;
    pp19Ne->iback1 = Correl.proton.Sol[0]->iback;
    pp19Ne->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    pp19Ne->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    pp19Ne->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    pp19Ne->et1 = Correl.proton.Sol[0]->energyTot;
    pp19Ne->time1 = Correl.proton.Sol[0]->CsITime;
    pp19Ne->energy_p1 = Correl.proton.Sol[0]->Ekin;
    pp19Ne->denergy1_R = Correl.proton.Sol[0]->denergyR;
    pp19Ne->energy_p1_R = Correl.proton.Sol[0]->energyR;
    pp19Ne->theta1 = Correl.proton.Sol[0]->theta;
    pp19Ne->phi1 = Correl.proton.Sol[0]->phi;

    pp19Ne->itele2 = Correl.proton.Sol[1]->itele;
    pp19Ne->id2 = Correl.proton.Sol[1]->iCsI;
    pp19Ne->ifront2 = Correl.proton.Sol[1]->ifront;
    pp19Ne->iback2 = Correl.proton.Sol[1]->iback;
    pp19Ne->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    pp19Ne->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    pp19Ne->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    pp19Ne->et2 = Correl.proton.Sol[1]->energyTot;
    pp19Ne->time2 = Correl.proton.Sol[1]->CsITime;
    pp19Ne->energy_p2 = Correl.proton.Sol[1]->Ekin;
    pp19Ne->denergy2_R = Correl.proton.Sol[1]->denergyR;
    pp19Ne->energy_p2_R = Correl.proton.Sol[1]->energyR;
    pp19Ne->theta2 = Correl.proton.Sol[1]->theta;
    pp19Ne->phi2 = Correl.proton.Sol[1]->phi;

    pp19Ne->M3[0] = Correl.Ne19.Sol[0]->Mvect[0];
    pp19Ne->M3[1] = Correl.Ne19.Sol[0]->Mvect[1];
    pp19Ne->M3[2] = Correl.Ne19.Sol[0]->Mvect[2];
    pp19Ne->et3 = Correl.Ne19.Sol[0]->energyTot;
    pp19Ne->energy_p3 = Correl.Ne19.Sol[0]->Ekin;
    pp19Ne->theta3 = Correl.Ne19.Sol[0]->theta;
    pp19Ne->phi3 = Correl.Ne19.Sol[0]->phi;

    pp19Ne->Erel = Erel_21Mg;
    pp19Ne->Ex = Ex;
    pp19Ne->Vcm = Correl.velocityCM;
    pp19Ne->thetaCM = thetaCM;
    pp19Ne->cos_thetaH = Correl.cos_thetaH;

    pp19Ne->beamZ = Correl.Ne19.Sol[0]->beamZ;
    pp19Ne->runnum = runnum;
    pp19Ne->T->Fill();
  }
}

void det::corr_22Mg()
{
  //22Mg -> 21Na + p
  if(Correl.Na21.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q22Mg = mass_22Mg - (mass_21Na + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Na21.mask[0]=1;
    Correl.makeArray(1);

    float Erel_22Mg = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_22Mg - Q22Mg;


    Histo_sort->Erel_22Mg_p21Na->Fill(Erel_22Mg);
    Histo_sort->Ex_22Mg_p21Na->Fill(Ex);

    Histo_sort->ThetaCM_22Mg_p21Na->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_22Mg_p21Na->Fill(Correl.velocityCM);
    Histo_sort->Erel_p21Na_costhetaH->Fill(Erel_22Mg,Correl.cos_thetaH);
    Histo_sort->p21Na_VCMvsErel->Fill(Correl.velocityCM,Erel_22Mg);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_22Mg_p21Na_trans->Fill(Erel_22Mg);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p21Na->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Na21.Sol[0]->velocity);

	    Histo_sort->p21Na_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p21Na_gammasADDvsErel->Fill(Erel_22Mg,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Na21.Sol[0]->velocity);

        Histo_sort->p21Na_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p21Na->Egamma[i] = doppE;
      p21Na->Tgamma[i] = Ceasar->added[i].time;
      p21Na->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->p21Na_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //cout << Correl.cos_thetaH << endl;

    //Set tree variables
    p21Na->itele1 = Correl.proton.Sol[0]->itele;
    p21Na->id1 = Correl.proton.Sol[0]->iCsI;
    p21Na->ifront1 = Correl.proton.Sol[0]->ifront;
    p21Na->iback1 = Correl.proton.Sol[0]->iback;
    p21Na->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p21Na->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p21Na->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p21Na->et1 = Correl.proton.Sol[0]->energyTot;
    p21Na->time1 = Correl.proton.Sol[0]->CsITime;
    p21Na->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p21Na->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p21Na->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p21Na->theta1 = Correl.proton.Sol[0]->theta;
    p21Na->phi1 = Correl.proton.Sol[0]->phi;

    p21Na->M2[0] = Correl.Na21.Sol[0]->Mvect[0];
    p21Na->M2[1] = Correl.Na21.Sol[0]->Mvect[1];
    p21Na->M2[2] = Correl.Na21.Sol[0]->Mvect[2];
    p21Na->et2 = Correl.Na21.Sol[0]->energyTot;
    p21Na->energy_p2 = Correl.Na21.Sol[0]->Ekin;
    p21Na->theta2 = Correl.Na21.Sol[0]->theta;
    p21Na->phi2 = Correl.Na21.Sol[0]->phi;

    p21Na->Erel = Erel_22Mg;
    p21Na->Ex = Ex;
    p21Na->Vcm = Correl.velocityCM;
    p21Na->thetaCM = thetaCM;
    p21Na->cos_thetaH = Correl.cos_thetaH;

    p21Na->beamZ = Correl.Na21.Sol[0]->beamZ;
    p21Na->runnum = runnum;
    p21Na->T->Fill();
  }
}

void det::corr_20Al()
{
  //2Al -> 17Ne + 3p
  if(Correl.Ne17.mult == 1 && Correl.proton.mult == 3)
  {
    //float const Q21Al = mass_21Al - (mass_20Mg + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.proton.mask[2]=1;
    Correl.Ne17.mask[0]=1;
    Correl.makeArray(1);

    float Erel_20Al = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    //float Ex = Erel_21Al - Q21Al;


    Histo_sort->Erel_20Al_3p17Ne->Fill(Erel_20Al);
    //Histo_sort->Ex_20Al_3p17Ne->Fill(Ex);

    Histo_sort->ThetaCM_20Al_3p17Ne->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_20Al_3p17Ne->Fill(Correl.velocityCM);
    Histo_sort->Erel_3p17Ne_costhetaH->Fill(Erel_20Al,Correl.cos_thetaH);
    Histo_sort->ppp17Ne_VCMvsErel->Fill(Correl.velocityCM,Erel_20Al);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_20Al_3p17Ne_trans->Fill(Erel_20Al);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    ppp17Ne->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Ne17.Sol[0]->velocity);

	    Histo_sort->ppp17Ne_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->ppp17Ne_gammasADDvsErel->Fill(Erel_20Al,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Ne17.Sol[0]->velocity);

        Histo_sort->ppp17Ne_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      ppp17Ne->Egamma[i] = doppE;
      ppp17Ne->Tgamma[i] = Ceasar->added[i].time;
      ppp17Ne->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->ppp17Ne_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //cout << Correl.cos_thetaH << endl;

    //Set tree variables
    ppp17Ne->itele1 = Correl.proton.Sol[0]->itele;
    ppp17Ne->id1 = Correl.proton.Sol[0]->iCsI;
    ppp17Ne->ifront1 = Correl.proton.Sol[0]->ifront;
    ppp17Ne->iback1 = Correl.proton.Sol[0]->iback;
    ppp17Ne->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    ppp17Ne->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    ppp17Ne->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    ppp17Ne->et1 = Correl.proton.Sol[0]->energyTot;
    ppp17Ne->time1 = Correl.proton.Sol[0]->CsITime;
    ppp17Ne->energy_p1 = Correl.proton.Sol[0]->Ekin;
    ppp17Ne->denergy1_R = Correl.proton.Sol[0]->denergyR;
    ppp17Ne->energy_p1_R = Correl.proton.Sol[0]->energyR;
    ppp17Ne->theta1 = Correl.proton.Sol[0]->theta;
    ppp17Ne->phi1 = Correl.proton.Sol[0]->phi;

    ppp17Ne->itele2 = Correl.proton.Sol[1]->itele;
    ppp17Ne->id2 = Correl.proton.Sol[1]->iCsI;
    ppp17Ne->ifront2 = Correl.proton.Sol[1]->ifront;
    ppp17Ne->iback2 = Correl.proton.Sol[1]->iback;
    ppp17Ne->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    ppp17Ne->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    ppp17Ne->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    ppp17Ne->et2 = Correl.proton.Sol[1]->energyTot;
    ppp17Ne->time2 = Correl.proton.Sol[1]->CsITime;
    ppp17Ne->energy_p2 = Correl.proton.Sol[1]->Ekin;
    ppp17Ne->denergy2_R = Correl.proton.Sol[1]->denergyR;
    ppp17Ne->energy_p2_R = Correl.proton.Sol[1]->energyR;
    ppp17Ne->theta2 = Correl.proton.Sol[1]->theta;
    ppp17Ne->phi2 = Correl.proton.Sol[1]->phi;

    ppp17Ne->itele3 = Correl.proton.Sol[2]->itele;
    ppp17Ne->id3 = Correl.proton.Sol[2]->iCsI;
    ppp17Ne->ifront3 = Correl.proton.Sol[2]->ifront;
    ppp17Ne->iback3 = Correl.proton.Sol[2]->iback;
    ppp17Ne->M3[0] = Correl.proton.Sol[2]->Mvect[0];
    ppp17Ne->M3[1] = Correl.proton.Sol[2]->Mvect[1];
    ppp17Ne->M3[2] = Correl.proton.Sol[2]->Mvect[2];
    ppp17Ne->et3 = Correl.proton.Sol[2]->energyTot;
    ppp17Ne->et3 = Correl.proton.Sol[2]->energyTot;
    ppp17Ne->energy_p3 = Correl.proton.Sol[2]->Ekin;

    ppp17Ne->M4[0] = Correl.Ne17.Sol[0]->Mvect[0];
    ppp17Ne->M4[1] = Correl.Ne17.Sol[0]->Mvect[1];
    ppp17Ne->M4[2] = Correl.Ne17.Sol[0]->Mvect[2];
    ppp17Ne->et4 = Correl.Ne17.Sol[0]->energyTot;
    ppp17Ne->energy_p4 = Correl.Ne17.Sol[0]->Ekin;
    ppp17Ne->theta4 = Correl.Ne17.Sol[0]->theta;
    ppp17Ne->phi4 = Correl.Ne17.Sol[0]->phi;


    ppp17Ne->Erel = Erel_20Al;
    //ppp17Ne->Ex = Ex; //No Ex
    ppp17Ne->Vcm = Correl.velocityCM;
    ppp17Ne->thetaCM = thetaCM;
    ppp17Ne->cos_thetaH = Correl.cos_thetaH;

    ppp17Ne->beamZ = Correl.Ne17.Sol[0]->beamZ;
    ppp17Ne->runnum = runnum;
    ppp17Ne->T->Fill();
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

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_21Al_p20Mg_trans->Fill(Erel_21Al);

    Histo_sort->Al21_pKE->Fill(Correl.proton.Sol[0]->Ekin);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p20Mg->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Mg20.Sol[0]->velocity);

      //cout << doppE << endl;
	    Histo_sort->p20Mg_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p20Mg_gammasADDvsErel->Fill(Erel_21Al,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Mg20.Sol[0]->velocity);

        Histo_sort->p20Mg_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p20Mg->Egamma[i] = doppE;
      p20Mg->Tgamma[i] = Ceasar->added[i].time;
      p20Mg->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->p20Mg_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //cout << Correl.cos_thetaH << endl;

    //Set tree variables
    p20Mg->itele1 = Correl.proton.Sol[0]->itele;
    p20Mg->id1 = Correl.proton.Sol[0]->iCsI;
    p20Mg->ifront1 = Correl.proton.Sol[0]->ifront;
    p20Mg->iback1 = Correl.proton.Sol[0]->iback;
    p20Mg->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p20Mg->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p20Mg->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p20Mg->et1 = Correl.proton.Sol[0]->energyTot;
    p20Mg->time1 = Correl.proton.Sol[0]->CsITime;
    p20Mg->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p20Mg->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p20Mg->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p20Mg->theta1 = Correl.proton.Sol[0]->theta;
    p20Mg->phi1 = Correl.proton.Sol[0]->phi;

    p20Mg->M2[0] = Correl.Mg20.Sol[0]->Mvect[0];
    p20Mg->M2[1] = Correl.Mg20.Sol[0]->Mvect[1];
    p20Mg->M2[2] = Correl.Mg20.Sol[0]->Mvect[2];
    p20Mg->et2 = Correl.Mg20.Sol[0]->energyTot;
    p20Mg->energy_p2 = Correl.Mg20.Sol[0]->Ekin;
    p20Mg->theta2 = Correl.Mg20.Sol[0]->theta;
    p20Mg->phi2 = Correl.Mg20.Sol[0]->phi;

    p20Mg->Erel = Erel_21Al;
    p20Mg->Ex = Ex;
    p20Mg->Vcm = Correl.velocityCM;
    p20Mg->thetaCM = thetaCM;
    p20Mg->cos_thetaH = Correl.cos_thetaH;

    p20Mg->beamZ = Correl.Mg20.Sol[0]->beamZ;
    p20Mg->runnum = runnum;
    p20Mg->T->Fill();
  }

  //21Al -> 18Ne + 3p
  if(Correl.Ne18.mult == 1 && Correl.proton.mult == 3)
  {
    //float const Q21Al = mass_21Al - (mass_20Mg + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.proton.mask[2]=1;
    Correl.Ne18.mask[0]=1;
    Correl.makeArray(1);

    float Erel_21Al = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    //float Ex = Erel_21Al - Q21Al;


    Histo_sort->Erel_21Al_3p18Ne->Fill(Erel_21Al);
    //Histo_sort->Ex_20Al_3p17Ne->Fill(Ex);

    Histo_sort->ThetaCM_21Al_3p18Ne->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_21Al_3p18Ne->Fill(Correl.velocityCM);
    Histo_sort->Erel_3p18Ne_costhetaH->Fill(Erel_21Al,Correl.cos_thetaH);
    Histo_sort->ppp18Ne_VCMvsErel->Fill(Correl.velocityCM,Erel_21Al);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_21Al_3p18Ne_trans->Fill(Erel_21Al);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    ppp18Ne->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Ne18.Sol[0]->velocity);

	    Histo_sort->ppp18Ne_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->ppp18Ne_gammasADDvsErel->Fill(Erel_21Al,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Ne18.Sol[0]->velocity);

        Histo_sort->ppp18Ne_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      ppp18Ne->Egamma[i] = doppE;
      ppp18Ne->Tgamma[i] = Ceasar->added[i].time;
      ppp18Ne->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->ppp18Ne_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //cout << Correl.cos_thetaH << endl;

    //Set tree variables
    ppp18Ne->itele1 = Correl.proton.Sol[0]->itele;
    ppp18Ne->id1 = Correl.proton.Sol[0]->iCsI;
    ppp18Ne->ifront1 = Correl.proton.Sol[0]->ifront;
    ppp18Ne->iback1 = Correl.proton.Sol[0]->iback;
    ppp18Ne->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    ppp18Ne->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    ppp18Ne->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    ppp18Ne->et1 = Correl.proton.Sol[0]->energyTot;
    ppp18Ne->time1 = Correl.proton.Sol[0]->CsITime;
    ppp18Ne->energy_p1 = Correl.proton.Sol[0]->Ekin;
    ppp18Ne->denergy1_R = Correl.proton.Sol[0]->denergyR;
    ppp18Ne->energy_p1_R = Correl.proton.Sol[0]->energyR;
    ppp18Ne->theta1 = Correl.proton.Sol[0]->theta;
    ppp18Ne->phi1 = Correl.proton.Sol[0]->phi;

    ppp18Ne->itele2 = Correl.proton.Sol[1]->itele;
    ppp18Ne->id2 = Correl.proton.Sol[1]->iCsI;
    ppp18Ne->ifront2 = Correl.proton.Sol[1]->ifront;
    ppp18Ne->iback2 = Correl.proton.Sol[1]->iback;
    ppp18Ne->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    ppp18Ne->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    ppp18Ne->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    ppp18Ne->et2 = Correl.proton.Sol[1]->energyTot;
    ppp18Ne->time2 = Correl.proton.Sol[1]->CsITime;
    ppp18Ne->energy_p2 = Correl.proton.Sol[1]->Ekin;
    ppp18Ne->denergy2_R = Correl.proton.Sol[1]->denergyR;
    ppp18Ne->energy_p2_R = Correl.proton.Sol[1]->energyR;
    ppp18Ne->theta2 = Correl.proton.Sol[1]->theta;
    ppp18Ne->phi2 = Correl.proton.Sol[1]->phi;

    ppp18Ne->itele3 = Correl.proton.Sol[2]->itele;
    ppp18Ne->id3 = Correl.proton.Sol[2]->iCsI;
    ppp18Ne->ifront3 = Correl.proton.Sol[2]->ifront;
    ppp18Ne->iback3 = Correl.proton.Sol[2]->iback;
    ppp18Ne->M3[0] = Correl.proton.Sol[2]->Mvect[0];
    ppp18Ne->M3[1] = Correl.proton.Sol[2]->Mvect[1];
    ppp18Ne->M3[2] = Correl.proton.Sol[2]->Mvect[2];
    ppp18Ne->et3 = Correl.proton.Sol[2]->energyTot;
    ppp18Ne->time3 = Correl.proton.Sol[2]->CsITime;
    ppp18Ne->energy_p3 = Correl.proton.Sol[2]->Ekin;
    ppp18Ne->denergy3_R = Correl.proton.Sol[2]->denergyR;
    ppp18Ne->energy_p3_R = Correl.proton.Sol[2]->energyR;
    ppp18Ne->theta3 = Correl.proton.Sol[2]->theta;
    ppp18Ne->phi3 = Correl.proton.Sol[2]->phi;

    ppp18Ne->M4[0] = Correl.Ne18.Sol[0]->Mvect[0];
    ppp18Ne->M4[1] = Correl.Ne18.Sol[0]->Mvect[1];
    ppp18Ne->M4[2] = Correl.Ne18.Sol[0]->Mvect[2];
    ppp18Ne->et4 = Correl.Ne18.Sol[0]->energyTot;
    ppp18Ne->energy_p4 = Correl.Ne18.Sol[0]->Ekin;
    ppp18Ne->theta4 = Correl.Ne18.Sol[0]->theta;
    ppp18Ne->phi4 = Correl.Ne18.Sol[0]->phi;

    ppp18Ne->Erel = Erel_21Al;
    //ppp18Ne->Ex = Ex; //No Ex
    ppp18Ne->Vcm = Correl.velocityCM;
    ppp18Ne->thetaCM = thetaCM;
    ppp18Ne->cos_thetaH = Correl.cos_thetaH;

    ppp18Ne->beamZ = Correl.Ne18.Sol[0]->beamZ;
    ppp18Ne->runnum = runnum;
    ppp18Ne->T->Fill();
  }
}

void det::corr_22Al()
{
  //22Al -> 21Mg + p
  if(Correl.Mg21.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q22Al = mass_22Al - (mass_21Mg + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Mg21.mask[0]=1;
    Correl.makeArray(1);

    float Erel_22Al = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_22Al - Q22Al;


    Histo_sort->Erel_22Al_p21Mg->Fill(Erel_22Al);
    Histo_sort->Ex_22Al_p21Mg->Fill(Ex);

    Histo_sort->ThetaCM_22Al_p21Mg->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_22Al_p21Mg->Fill(Correl.velocityCM);
    Histo_sort->Erel_p21Mg_costhetaH->Fill(Erel_22Al,Correl.cos_thetaH);
    Histo_sort->p21Mg_VCMvsErel->Fill(Correl.velocityCM,Erel_22Al);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_22Al_p21Mg_trans->Fill(Erel_22Al);

    Histo_sort->Al22_pKE->Fill(Correl.proton.Sol[0]->Ekin);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p21Mg->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Currently skip events greater than 10
      if (Ceasar->Nadded > 15) {
        p21Mg->Ngamma = 0;
        break;
      }

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Mg21.Sol[0]->velocity);

	    Histo_sort->p21Mg_gammasADD_nodopp->Fill(Ceasar->added[i].energy * 1000.); //E stored as MeV
	    Histo_sort->p21Mg_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p21Mg_gammasADDvsErel->Fill(Erel_22Al,doppE * 1000.);
      if (Erel_22Al >= 1 && Erel_22Al <= 1.5) Histo_sort->p21Mg_gammasADD_peak2->Fill(doppE * 1000.);
      if (Erel_22Al >= 1.7 && Erel_22Al <= 2.4) Histo_sort->p21Mg_gammasADD_peak3->Fill(doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Mg21.Sol[0]->velocity);

        Histo_sort->p21Mg_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p21Mg->Egamma[i] = doppE;
      p21Mg->Tgamma[i] = Ceasar->added[i].time;
      p21Mg->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->p21Mg_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
 	      if (Erel_22Al >= 1 && Erel_22Al <= 1.5) Histo_sort->p21Mg_gammasADD_peak2_tgate->Fill(doppE * 1000.); //E stored as MeV
 	      if (Erel_22Al >= 1.7 && Erel_22Al <= 2.4) Histo_sort->p21Mg_gammasADD_peak3_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

      p21Mg->Ngamma_Select = Ceasar->Nselect; //Num gammas for tree

      //For non addback gammas
      for(int i = 0; i < Ceasar->Nselect; i++) {

      //Currently skip events greater than 10
      if (Ceasar->Nselect > 15) {
        p21Mg->Ngamma_Select = 0;
        break;
      }

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->select[i].theta)*cos(Ceasar->select[i].phi);
      double yg = sin(Ceasar->select[i].theta)*sin(Ceasar->select[i].phi);
      double zg = cos(Ceasar->select[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->select[i].energy,thetarel,Correl.Mg21.Sol[0]->velocity);
      //cout << "here " << i << " " << Ceasar->Nselect << endl;
      //cout << doppE << " " << Ceasar->select[i].time << " " << Ceasar->select[i].id << endl;
      p21Mg->Egamma_Select[i] = doppE;
      p21Mg->Tgamma_Select[i] = Ceasar->select[i].time;
      p21Mg->Chgamma_Select[i] = Ceasar->select[i].id;

    }

    //Set tree variables
    p21Mg->itele1 = Correl.proton.Sol[0]->itele;
    p21Mg->id1 = Correl.proton.Sol[0]->iCsI;
    p21Mg->ifront1 = Correl.proton.Sol[0]->ifront;
    p21Mg->iback1 = Correl.proton.Sol[0]->iback;
    p21Mg->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p21Mg->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p21Mg->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p21Mg->et1 = Correl.proton.Sol[0]->energyTot;
    p21Mg->time1 = Correl.proton.Sol[0]->CsITime;
    p21Mg->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p21Mg->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p21Mg->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p21Mg->theta1 = Correl.proton.Sol[0]->theta;
    p21Mg->phi1 = Correl.proton.Sol[0]->phi;


    p21Mg->M2[0] = Correl.Mg21.Sol[0]->Mvect[0];
    p21Mg->M2[1] = Correl.Mg21.Sol[0]->Mvect[1];
    p21Mg->M2[2] = Correl.Mg21.Sol[0]->Mvect[2];
    p21Mg->et2 = Correl.Mg21.Sol[0]->energyTot;
    p21Mg->energy_p2 = Correl.Mg21.Sol[0]->Ekin;
    p21Mg->theta2 = Correl.Mg21.Sol[0]->theta;
    p21Mg->phi2 = Correl.Mg21.Sol[0]->phi;


    p21Mg->Erel = Erel_22Al;
    p21Mg->Ex = Ex;
    p21Mg->Vcm = Correl.velocityCM;
    p21Mg->thetaCM = thetaCM;
    p21Mg->cos_thetaH = Correl.cos_thetaH;

    p21Mg->beamZ = Correl.Mg21.Sol[0]->beamZ;
    p21Mg->runnum = runnum;
    p21Mg->T->Fill();

  }
}

void det::corr_23Al()
{
  //23Al -> 22Mg + p
  if(Correl.Mg22.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q23Al = mass_23Al - (mass_22Mg + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Mg22.mask[0]=1;
    Correl.makeArray(1);

    float Erel_23Al = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_23Al - Q23Al;


    Histo_sort->Erel_23Al_p22Mg->Fill(Erel_23Al);
    Histo_sort->Ex_23Al_p22Mg->Fill(Ex);

    Histo_sort->ThetaCM_23Al_p22Mg->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_23Al_p22Mg->Fill(Correl.velocityCM);
    Histo_sort->Erel_p22Mg_costhetaH->Fill(Erel_23Al,Correl.cos_thetaH);
    Histo_sort->p22Mg_VCMvsErel->Fill(Correl.velocityCM,Erel_23Al);

    if (fabs(Correl.cos_thetaH) <= 0.1)
    {
      Histo_sort->Erel_23Al_p22Mg_trans->Fill(Erel_23Al);
      Histo_sort->Ex_23Al_p22Mg_trans->Fill(Erel_23Al);
    }
    if (fabs(Correl.cos_thetaH) >= 0.8)
    {
      Histo_sort->Erel_23Al_p22Mg_long->Fill(Erel_23Al);
      Histo_sort->Ex_23Al_p22Mg_long->Fill(Erel_23Al);
    }

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p22Mg->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Mg22.Sol[0]->velocity);

	    Histo_sort->p22Mg_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p22Mg_gammasADDvsErel->Fill(Erel_23Al,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Mg22.Sol[0]->velocity);

        Histo_sort->p22Mg_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p22Mg->Egamma[i] = doppE;
      p22Mg->Tgamma[i] = Ceasar->added[i].time;
      p22Mg->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->p22Mg_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //Run dependent histograms. Set 1 had gobbi farther away than set 1a
    if (runnum < 215) Histo_sort->Erel_23Al_p22Mg_set1->Fill(Erel_23Al);
    else Histo_sort->Erel_23Al_p22Mg_set1a->Fill(Erel_23Al);

    //Set tree variables
    p22Mg->itele1 = Correl.proton.Sol[0]->itele;
    p22Mg->id1 = Correl.proton.Sol[0]->iCsI;
    p22Mg->ifront1 = Correl.proton.Sol[0]->ifront;
    p22Mg->iback1 = Correl.proton.Sol[0]->iback;
    p22Mg->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p22Mg->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p22Mg->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p22Mg->et1 = Correl.proton.Sol[0]->energyTot;
    p22Mg->time1 = Correl.proton.Sol[0]->CsITime;
    p22Mg->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p22Mg->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p22Mg->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p22Mg->theta1 = Correl.proton.Sol[0]->theta;
    p22Mg->phi1 = Correl.proton.Sol[0]->phi;

    p22Mg->M2[0] = Correl.Mg22.Sol[0]->Mvect[0];
    p22Mg->M2[1] = Correl.Mg22.Sol[0]->Mvect[1];
    p22Mg->M2[2] = Correl.Mg22.Sol[0]->Mvect[2];
    p22Mg->et2 = Correl.Mg22.Sol[0]->energyTot;
    p22Mg->energy_p2 = Correl.Mg22.Sol[0]->Ekin;
    p22Mg->theta2 = Correl.Mg22.Sol[0]->theta;
    p22Mg->phi2 = Correl.Mg22.Sol[0]->phi;

    p22Mg->Erel = Erel_23Al;
    p22Mg->Ex = Ex;
    p22Mg->Vcm = Correl.velocityCM;
    p22Mg->thetaCM = thetaCM;
    p22Mg->cos_thetaH = Correl.cos_thetaH;

    p22Mg->beamZ = Correl.Mg22.Sol[0]->beamZ;
    p22Mg->runnum = runnum;
    p22Mg->T->Fill();

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

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_22Si_2p20Mg_trans->Fill(Erel_22Si);

    Histo_sort->Si22_pKE->Fill(Correl.proton.Sol[0]->Ekin);
    Histo_sort->Si22_pKE->Fill(Correl.proton.Sol[1]->Ekin);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    pp20Mg->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Mg20.Sol[0]->velocity);

	    Histo_sort->pp20Mg_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->pp20Mg_gammasADDvsErel->Fill(Erel_22Si,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Mg20.Sol[0]->velocity);

        Histo_sort->pp20Mg_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      pp20Mg->Egamma[i] = doppE;
      pp20Mg->Tgamma[i] = Ceasar->added[i].time;
      pp20Mg->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->pp20Mg_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //Peek at telescope 0 csi 0
    //if (Correl.proton.Sol[0]->itele == 0 && Correl.proton.Sol[0]->iCsI == 0) cout << Erel_22Si << endl;
    //if (Correl.proton.Sol[1]->itele == 0 && Correl.proton.Sol[1]->iCsI == 0) cout << Erel_22Si << endl;

    //Set tree variables
    pp20Mg->itele1 = Correl.proton.Sol[0]->itele;
    pp20Mg->id1 = Correl.proton.Sol[0]->iCsI;
    pp20Mg->ifront1 = Correl.proton.Sol[0]->ifront;
    pp20Mg->iback1 = Correl.proton.Sol[0]->iback;
    pp20Mg->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    pp20Mg->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    pp20Mg->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    pp20Mg->et1 = Correl.proton.Sol[0]->energyTot;
    pp20Mg->time1 = Correl.proton.Sol[0]->CsITime;
    pp20Mg->energy_p1 = Correl.proton.Sol[0]->Ekin;
    pp20Mg->denergy1_R = Correl.proton.Sol[0]->denergyR;
    pp20Mg->energy_p1_R = Correl.proton.Sol[0]->energyR;
    pp20Mg->theta1 = Correl.proton.Sol[0]->theta;
    pp20Mg->phi1 = Correl.proton.Sol[0]->phi;

    pp20Mg->itele2 = Correl.proton.Sol[1]->itele;
    pp20Mg->id2 = Correl.proton.Sol[1]->iCsI;
    pp20Mg->ifront2 = Correl.proton.Sol[1]->ifront;
    pp20Mg->iback2 = Correl.proton.Sol[1]->iback;
    pp20Mg->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    pp20Mg->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    pp20Mg->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    pp20Mg->et2 = Correl.proton.Sol[1]->energyTot;
    pp20Mg->time2 = Correl.proton.Sol[1]->CsITime;
    pp20Mg->energy_p2 = Correl.proton.Sol[1]->Ekin;
    pp20Mg->denergy2_R = Correl.proton.Sol[1]->denergyR;
    pp20Mg->energy_p2_R = Correl.proton.Sol[1]->energyR;
    pp20Mg->theta2 = Correl.proton.Sol[1]->theta;
    pp20Mg->phi2 = Correl.proton.Sol[1]->phi;

    pp20Mg->M3[0] = Correl.Mg20.Sol[0]->Mvect[0];
    pp20Mg->M3[1] = Correl.Mg20.Sol[0]->Mvect[1];
    pp20Mg->M3[2] = Correl.Mg20.Sol[0]->Mvect[2];
    pp20Mg->et3 = Correl.Mg20.Sol[0]->energyTot;
    pp20Mg->energy_p3 = Correl.Mg20.Sol[0]->Ekin;
    pp20Mg->theta3 = Correl.Mg20.Sol[0]->theta;
    pp20Mg->phi3 = Correl.Mg20.Sol[0]->phi;

    pp20Mg->Erel = Erel_22Si;
    //pp20Mg->Ex = Ex; //No Ex
    pp20Mg->Vcm = Correl.velocityCM;
    pp20Mg->thetaCM = thetaCM;
    pp20Mg->cos_thetaH = Correl.cos_thetaH;

    pp20Mg->beamZ = Correl.Mg20.Sol[0]->beamZ;
    pp20Mg->runnum = runnum;
    pp20Mg->T->Fill();

    //Jacobi, do after so it doesn't mess it up tree
    if (Erel_22Si >= 1.5 && Erel_22Si <= 2.5)
    {
      Correl.getJacobi();

      //Mask out fragments to get the relative energy
      Correl.Mg20.mask[0] = 0;
      Correl.makeArray(1);
      float Erel_pp = Correl.findErel()/Erel_22Si;

      Histo_sort->Si22_JacobiT_xy_s->Fill(Erel_pp,Correl.cosThetaT);
      Histo_sort->Si22_JacobiT_xy_s->Fill(Erel_pp,Correl.cosThetaT);

      Correl.Mg20.mask[1] = 1;
      Correl.proton.mask[0] = 1;
      Correl.proton.mask[1] = 0;
      Correl.makeArray(1);
      float Erel_pp1 = Correl.findErel()/Erel_22Si;

      Histo_sort->Si22_JacobiY_xy_s->Fill(Erel_pp1,Correl.cosThetaY[0]);

      Correl.proton.mask[0] = 0;
      Correl.proton.mask[1] = 1;
      Correl.makeArray(1);
      float Erel_pp2 = Correl.findErel()/Erel_22Si;

      Histo_sort->Si22_JacobiY_xy_s->Fill(Erel_pp2,Correl.cosThetaY[1]);

    }

  }


    //22Si -> 18Ne + 4p
  if(Correl.Ne18.mult == 1 && Correl.proton.mult == 4)
  {
    //TODO Si22 and P23 masses are not known
    //     Use calculations?
    //float const Q22Si = mass_22Si - (mass_20Mg + 2.*mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.proton.mask[2]=1;
    Correl.proton.mask[3]=1;
    Correl.Ne18.mask[0]=1;
    Correl.makeArray(1);

    float Erel_22Si = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    //float Ex = Erel_22Si - Q22Si;

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    pppp18Ne->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Ne18.Sol[0]->velocity);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Ne18.Sol[0]->velocity);

      }

      pppp18Ne->Egamma[i] = doppE;
      pppp18Ne->Tgamma[i] = Ceasar->added[i].time;
      pppp18Ne->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    pppp18Ne->itele1 = Correl.proton.Sol[0]->itele;
    pppp18Ne->id1 = Correl.proton.Sol[0]->iCsI;
    pppp18Ne->ifront1 = Correl.proton.Sol[0]->ifront;
    pppp18Ne->iback1 = Correl.proton.Sol[0]->iback;
    pppp18Ne->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    pppp18Ne->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    pppp18Ne->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    pppp18Ne->et1 = Correl.proton.Sol[0]->energyTot;
    pppp18Ne->time1 = Correl.proton.Sol[0]->CsITime;
    pppp18Ne->energy_p1 = Correl.proton.Sol[0]->Ekin;
    pppp18Ne->denergy1_R = Correl.proton.Sol[0]->denergyR;
    pppp18Ne->energy_p1_R = Correl.proton.Sol[0]->energyR;
    pppp18Ne->theta1 = Correl.proton.Sol[0]->theta;
    pppp18Ne->phi1 = Correl.proton.Sol[0]->phi;

    pppp18Ne->itele2 = Correl.proton.Sol[1]->itele;
    pppp18Ne->id2 = Correl.proton.Sol[1]->iCsI;
    pppp18Ne->ifront2 = Correl.proton.Sol[1]->ifront;
    pppp18Ne->iback2 = Correl.proton.Sol[1]->iback;
    pppp18Ne->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    pppp18Ne->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    pppp18Ne->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    pppp18Ne->et2 = Correl.proton.Sol[1]->energyTot;
    pppp18Ne->time2 = Correl.proton.Sol[1]->CsITime;
    pppp18Ne->energy_p2 = Correl.proton.Sol[1]->Ekin;
    pppp18Ne->denergy2_R = Correl.proton.Sol[1]->denergyR;
    pppp18Ne->energy_p2_R = Correl.proton.Sol[1]->energyR;
    pppp18Ne->theta2 = Correl.proton.Sol[1]->theta;
    pppp18Ne->phi2 = Correl.proton.Sol[1]->phi;

    pppp18Ne->itele3 = Correl.proton.Sol[2]->itele;
    pppp18Ne->id3 = Correl.proton.Sol[2]->iCsI;
    pppp18Ne->ifront3 = Correl.proton.Sol[2]->ifront;
    pppp18Ne->iback3 = Correl.proton.Sol[2]->iback;
    pppp18Ne->M3[0] = Correl.proton.Sol[2]->Mvect[0];
    pppp18Ne->M3[1] = Correl.proton.Sol[2]->Mvect[1];
    pppp18Ne->M3[2] = Correl.proton.Sol[2]->Mvect[2];
    pppp18Ne->et3 = Correl.proton.Sol[2]->energyTot;
    pppp18Ne->time3 = Correl.proton.Sol[2]->CsITime;
    pppp18Ne->energy_p3 = Correl.proton.Sol[2]->Ekin;
    pppp18Ne->denergy3_R = Correl.proton.Sol[2]->denergyR;
    pppp18Ne->energy_p3_R = Correl.proton.Sol[2]->energyR;
    pppp18Ne->theta3 = Correl.proton.Sol[2]->theta;
    pppp18Ne->phi3 = Correl.proton.Sol[2]->phi;

    pppp18Ne->itele4 = Correl.proton.Sol[3]->itele;
    pppp18Ne->id4 = Correl.proton.Sol[3]->iCsI;
    pppp18Ne->ifront4 = Correl.proton.Sol[3]->ifront;
    pppp18Ne->iback4 = Correl.proton.Sol[3]->iback;
    pppp18Ne->M4[0] = Correl.proton.Sol[3]->Mvect[0];
    pppp18Ne->M4[1] = Correl.proton.Sol[3]->Mvect[1];
    pppp18Ne->M4[2] = Correl.proton.Sol[3]->Mvect[2];
    pppp18Ne->et4 = Correl.proton.Sol[3]->energyTot;
    pppp18Ne->time4 = Correl.proton.Sol[3]->CsITime;
    pppp18Ne->energy_p4 = Correl.proton.Sol[3]->Ekin;
    pppp18Ne->denergy4_R = Correl.proton.Sol[3]->denergyR;
    pppp18Ne->energy_p4_R = Correl.proton.Sol[3]->energyR;
    pppp18Ne->theta4 = Correl.proton.Sol[3]->theta;
    pppp18Ne->phi4 = Correl.proton.Sol[3]->phi;

    pppp18Ne->M5[0] = Correl.Ne18.Sol[0]->Mvect[0];
    pppp18Ne->M5[1] = Correl.Ne18.Sol[0]->Mvect[1];
    pppp18Ne->M5[2] = Correl.Ne18.Sol[0]->Mvect[2];
    pppp18Ne->et5 = Correl.Ne18.Sol[0]->energyTot;
    pppp18Ne->energy_p5 = Correl.Ne18.Sol[0]->Ekin;
    pppp18Ne->theta5 = Correl.Ne18.Sol[0]->theta;
    pppp18Ne->phi5 = Correl.Ne18.Sol[0]->phi;

    pppp18Ne->Erel = Erel_22Si;
    //pp20Mg->Ex = Ex; //No Ex
    pppp18Ne->Vcm = Correl.velocityCM;
    pppp18Ne->thetaCM = thetaCM;
    pppp18Ne->cos_thetaH = Correl.cos_thetaH;

    pppp18Ne->beamZ = Correl.Ne18.Sol[0]->beamZ;
    pppp18Ne->runnum = runnum;
    pppp18Ne->T->Fill();


  }

}

void det::corr_23Si()
{
  //23Si -> 22Al + p
  if(Correl.Al22.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q23Si = mass_23Si - (mass_22Al + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Al22.mask[0]=1;
    Correl.makeArray(1);

    float Erel_23Si = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_23Si - Q23Si;


    Histo_sort->Erel_23Si_p22Al->Fill(Erel_23Si);
    Histo_sort->Ex_23Si_p22Al->Fill(Ex);

    Histo_sort->ThetaCM_23Si_p22Al->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_23Si_p22Al->Fill(Correl.velocityCM);
    Histo_sort->Erel_p22Al_costhetaH->Fill(Erel_23Si,Correl.cos_thetaH);
    Histo_sort->p22Al_VCMvsErel->Fill(Correl.velocityCM,Erel_23Si);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_23Si_p22Al_trans->Fill(Erel_23Si);

    Histo_sort->Si23_pKE->Fill(Correl.proton.Sol[0]->Ekin);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p22Al->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Al22.Sol[0]->velocity);

	    Histo_sort->p22Al_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p22Al_gammasADDvsErel->Fill(Erel_23Si,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Al22.Sol[0]->velocity);

        Histo_sort->p22Al_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p22Al->Egamma[i] = doppE;
      p22Al->Tgamma[i] = Ceasar->added[i].time;
      p22Al->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->p22Al_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }

    }

        //Set tree variables
      p22Al->itele1 = Correl.proton.Sol[0]->itele;
      p22Al->id1 = Correl.proton.Sol[0]->iCsI;
      p22Al->ifront1 = Correl.proton.Sol[0]->ifront;
      p22Al->iback1 = Correl.proton.Sol[0]->iback;
      p22Al->M1[0] = Correl.proton.Sol[0]->Mvect[0];
      p22Al->M1[1] = Correl.proton.Sol[0]->Mvect[1];
      p22Al->M1[2] = Correl.proton.Sol[0]->Mvect[2];
      p22Al->et1 = Correl.proton.Sol[0]->energyTot;
      p22Al->time1 = Correl.proton.Sol[0]->CsITime;
      p22Al->energy_p1 = Correl.proton.Sol[0]->Ekin;
      p22Al->denergy1_R = Correl.proton.Sol[0]->denergyR;
      p22Al->energy_p1_R = Correl.proton.Sol[0]->energyR;
      p22Al->theta1 = Correl.proton.Sol[0]->theta;
      p22Al->phi1 = Correl.proton.Sol[0]->phi;

      p22Al->M2[0] = Correl.Al22.Sol[0]->Mvect[0];
      p22Al->M2[1] = Correl.Al22.Sol[0]->Mvect[1];
      p22Al->M2[2] = Correl.Al22.Sol[0]->Mvect[2];
      p22Al->et2 = Correl.Al22.Sol[0]->energyTot;
      p22Al->energy_p2 = Correl.Al22.Sol[0]->Ekin;
      p22Al->theta2 = Correl.Al22.Sol[0]->theta;
      p22Al->phi2 = Correl.Al22.Sol[0]->phi;

      p22Al->Erel = Erel_23Si;
      p22Al->Ex = Ex;
      p22Al->Vcm = Correl.velocityCM;
      p22Al->thetaCM = thetaCM;
      p22Al->cos_thetaH = Correl.cos_thetaH;

      p22Al->beamZ = Correl.Al22.Sol[0]->beamZ;
      p22Al->runnum = runnum;
      p22Al->T->Fill();

  }

  //23Si -> 21Mg + 2p
  if(Correl.Mg21.mult == 1 && Correl.proton.mult == 2)
  {
    float const Q23Si = mass_23Si - (mass_21Mg + 2.*mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.Mg21.mask[0]=1;
    Correl.makeArray(1);

    float Erel_23Si = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_23Si - Q23Si;


    Histo_sort->Erel_23Si_2p21Mg->Fill(Erel_23Si);
    Histo_sort->Ex_23Si_2p21Mg->Fill(Ex);

    Histo_sort->ThetaCM_23Si_2p21Mg->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_23Si_2p21Mg->Fill(Correl.velocityCM);
    Histo_sort->Erel_2p21Mg_costhetaH->Fill(Erel_23Si,Correl.cos_thetaH);
    Histo_sort->pp21Mg_VCMvsErel->Fill(Correl.velocityCM,Erel_23Si);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_23Si_2p21Mg_trans->Fill(Erel_23Si);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    pp21Mg->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Mg21.Sol[0]->velocity);

	    Histo_sort->pp21Mg_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->pp21Mg_gammasADDvsErel->Fill(Erel_23Si,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Mg21.Sol[0]->velocity);

        Histo_sort->pp21Mg_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      pp21Mg->Egamma[i] = doppE;
      pp21Mg->Tgamma[i] = Ceasar->added[i].time;
      pp21Mg->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->pp21Mg_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //Set tree variables
    pp21Mg->itele1 = Correl.proton.Sol[0]->itele;
    pp21Mg->id1 = Correl.proton.Sol[0]->iCsI;
    pp21Mg->ifront1 = Correl.proton.Sol[0]->ifront;
    pp21Mg->iback1 = Correl.proton.Sol[0]->iback;
    pp21Mg->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    pp21Mg->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    pp21Mg->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    pp21Mg->et1 = Correl.proton.Sol[0]->energyTot;
    pp21Mg->time1 = Correl.proton.Sol[0]->CsITime;
    pp21Mg->energy_p1 = Correl.proton.Sol[0]->Ekin;
    pp21Mg->denergy1_R = Correl.proton.Sol[0]->denergyR;
    pp21Mg->energy_p1_R = Correl.proton.Sol[0]->energyR;
    pp21Mg->theta1 = Correl.proton.Sol[0]->theta;
    pp21Mg->phi1 = Correl.proton.Sol[0]->phi;

    pp21Mg->itele2 = Correl.proton.Sol[1]->itele;
    pp21Mg->id2 = Correl.proton.Sol[1]->iCsI;
    pp21Mg->ifront2 = Correl.proton.Sol[1]->ifront;
    pp21Mg->iback2 = Correl.proton.Sol[1]->iback;
    pp21Mg->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    pp21Mg->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    pp21Mg->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    pp21Mg->et2 = Correl.proton.Sol[1]->energyTot;
    pp21Mg->time2 = Correl.proton.Sol[1]->CsITime;
    pp21Mg->energy_p2 = Correl.proton.Sol[1]->Ekin;
    pp21Mg->denergy2_R = Correl.proton.Sol[1]->denergyR;
    pp21Mg->energy_p2_R = Correl.proton.Sol[1]->energyR;
    pp21Mg->theta2 = Correl.proton.Sol[1]->theta;
    pp21Mg->phi2 = Correl.proton.Sol[1]->phi;

    pp21Mg->M3[0] = Correl.Mg21.Sol[0]->Mvect[0];
    pp21Mg->M3[1] = Correl.Mg21.Sol[0]->Mvect[1];
    pp21Mg->M3[2] = Correl.Mg21.Sol[0]->Mvect[2];
    pp21Mg->et3 = Correl.Mg21.Sol[0]->energyTot;
    pp21Mg->energy_p3 = Correl.Mg21.Sol[0]->Ekin;
    pp21Mg->theta3 = Correl.Mg21.Sol[0]->theta;
    pp21Mg->phi3 = Correl.Mg21.Sol[0]->phi;

    pp21Mg->Erel = Erel_23Si;
    pp21Mg->Ex = Ex;
    pp21Mg->Vcm = Correl.velocityCM;
    pp21Mg->thetaCM = thetaCM;
    pp21Mg->cos_thetaH = Correl.cos_thetaH;

    pp21Mg->beamZ = Correl.Mg21.Sol[0]->beamZ;
    pp21Mg->runnum = runnum;
    pp21Mg->T->Fill();

  }
}

void det::corr_24Si()
{
  //24Si -> 23Al + p
  if(Correl.Al23.mult == 1 && Correl.proton.mult == 1)
  {
    float const Q24Si = mass_24Si - (mass_23Al + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Al23.mask[0]=1;
    Correl.makeArray(1);

    float Erel_24Si = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_24Si - Q24Si;


    Histo_sort->Erel_24Si_p23Al->Fill(Erel_24Si);
    Histo_sort->Ex_24Si_p23Al->Fill(Ex);

    Histo_sort->ThetaCM_24Si_p23Al->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_24Si_p23Al->Fill(Correl.velocityCM);
    Histo_sort->Erel_p23Al_costhetaH->Fill(Erel_24Si,Correl.cos_thetaH);
    Histo_sort->p23Al_VCMvsErel->Fill(Correl.velocityCM,Erel_24Si);

    if (fabs(Correl.cos_thetaH) <= 0.1)
    {
      Histo_sort->Erel_24Si_p23Al_trans->Fill(Erel_24Si);
      //Histo_sort->Ex_24Si_p23Al_trans->Fill(Erel_24Si);
    }
    //if (fabs(Correl.cos_thetaH) >= 0.8)
    //{
      //Histo_sort->Erel_24Si_p23Al_long->Fill(Erel_24Si);
      //Histo_sort->Ex_24Si_p23Al_long->Fill(Erel_24Si);
    //}

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p23Al->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Al23.Sol[0]->velocity);

	    Histo_sort->p23Al_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p23Al_gammasADDvsErel->Fill(Erel_24Si,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Al23.Sol[0]->velocity);

        Histo_sort->p23Al_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p23Al->Egamma[i] = doppE;
      p23Al->Tgamma[i] = Ceasar->added[i].time;
      p23Al->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->p23Al_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //Set tree variables
    p23Al->itele1 = Correl.proton.Sol[0]->itele;
    p23Al->id1 = Correl.proton.Sol[0]->iCsI;
    p23Al->ifront1 = Correl.proton.Sol[0]->ifront;
    p23Al->iback1 = Correl.proton.Sol[0]->iback;
    p23Al->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p23Al->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p23Al->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p23Al->et1 = Correl.proton.Sol[0]->energyTot;
    p23Al->time1 = Correl.proton.Sol[0]->CsITime;
    p23Al->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p23Al->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p23Al->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p23Al->theta1 = Correl.proton.Sol[0]->theta;
    p23Al->phi1 = Correl.proton.Sol[0]->phi;

    p23Al->M2[0] = Correl.Al23.Sol[0]->Mvect[0];
    p23Al->M2[1] = Correl.Al23.Sol[0]->Mvect[1];
    p23Al->M2[2] = Correl.Al23.Sol[0]->Mvect[2];
    p23Al->et2 = Correl.Al23.Sol[0]->energyTot;
    p23Al->energy_p2 = Correl.Al23.Sol[0]->Ekin;
    p23Al->theta2 = Correl.Al23.Sol[0]->theta;
    p23Al->phi2 = Correl.Al23.Sol[0]->phi;

    p23Al->Erel = Erel_24Si;
    p23Al->Ex = Ex;
    p23Al->Vcm = Correl.velocityCM;
    p23Al->thetaCM = thetaCM;
    p23Al->cos_thetaH = Correl.cos_thetaH;

    p23Al->beamZ = Correl.Al23.Sol[0]->beamZ;
    p23Al->runnum = runnum;
    p23Al->T->Fill();
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

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_23P_3p20Mg_trans->Fill(Erel_23P);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    ppp20Mg->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Mg20.Sol[0]->velocity);

	    Histo_sort->ppp20Mg_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->ppp20Mg_gammasADDvsErel->Fill(Erel_23P,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Mg20.Sol[0]->velocity);

        Histo_sort->ppp20Mg_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      ppp20Mg->Egamma[i] = doppE;
      ppp20Mg->Tgamma[i] = Ceasar->added[i].time;
      ppp20Mg->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->ppp20Mg_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //Set tree variables
    ppp20Mg->itele1 = Correl.proton.Sol[0]->itele;
    ppp20Mg->id1 = Correl.proton.Sol[0]->iCsI;
    ppp20Mg->ifront1 = Correl.proton.Sol[0]->ifront;
    ppp20Mg->iback1 = Correl.proton.Sol[0]->iback;
    ppp20Mg->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    ppp20Mg->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    ppp20Mg->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    ppp20Mg->et1 = Correl.proton.Sol[0]->energyTot;
    ppp20Mg->time1 = Correl.proton.Sol[0]->CsITime;
    ppp20Mg->energy_p1 = Correl.proton.Sol[0]->Ekin;
    ppp20Mg->denergy1_R = Correl.proton.Sol[0]->denergyR;
    ppp20Mg->energy_p1_R = Correl.proton.Sol[0]->energyR;
    ppp20Mg->theta1 = Correl.proton.Sol[0]->theta;
    ppp20Mg->phi1 = Correl.proton.Sol[0]->phi;

    ppp20Mg->itele2 = Correl.proton.Sol[1]->itele;
    ppp20Mg->id2 = Correl.proton.Sol[1]->iCsI;
    ppp20Mg->ifront2 = Correl.proton.Sol[1]->ifront;
    ppp20Mg->iback2 = Correl.proton.Sol[1]->iback;
    ppp20Mg->M2[0] = Correl.proton.Sol[1]->Mvect[0];
    ppp20Mg->M2[1] = Correl.proton.Sol[1]->Mvect[1];
    ppp20Mg->M2[2] = Correl.proton.Sol[1]->Mvect[2];
    ppp20Mg->et2 = Correl.proton.Sol[1]->energyTot;
    ppp20Mg->time2 = Correl.proton.Sol[1]->CsITime;
    ppp20Mg->energy_p2 = Correl.proton.Sol[1]->Ekin;
    ppp20Mg->denergy2_R = Correl.proton.Sol[1]->denergyR;
    ppp20Mg->energy_p2_R = Correl.proton.Sol[1]->energyR;
    ppp20Mg->theta2 = Correl.proton.Sol[1]->theta;
    ppp20Mg->phi2 = Correl.proton.Sol[1]->phi;

    ppp20Mg->itele3 = Correl.proton.Sol[2]->itele;
    ppp20Mg->id3 = Correl.proton.Sol[2]->iCsI;
    ppp20Mg->ifront3 = Correl.proton.Sol[2]->ifront;
    ppp20Mg->iback3 = Correl.proton.Sol[2]->iback;
    ppp20Mg->M3[0] = Correl.proton.Sol[2]->Mvect[0];
    ppp20Mg->M3[1] = Correl.proton.Sol[2]->Mvect[1];
    ppp20Mg->M3[2] = Correl.proton.Sol[2]->Mvect[2];
    ppp20Mg->et3 = Correl.proton.Sol[2]->energyTot;
    ppp20Mg->time3 = Correl.proton.Sol[2]->CsITime;
    ppp20Mg->energy_p3 = Correl.proton.Sol[2]->Ekin;
    ppp20Mg->denergy3_R = Correl.proton.Sol[2]->denergyR;
    ppp20Mg->energy_p3_R = Correl.proton.Sol[2]->energyR;
    ppp20Mg->theta3 = Correl.proton.Sol[2]->theta;
    ppp20Mg->phi3 = Correl.proton.Sol[2]->phi;

    ppp20Mg->M4[0] = Correl.Mg20.Sol[0]->Mvect[0];
    ppp20Mg->M4[1] = Correl.Mg20.Sol[0]->Mvect[1];
    ppp20Mg->M4[2] = Correl.Mg20.Sol[0]->Mvect[2];
    ppp20Mg->et4 = Correl.Mg20.Sol[0]->energyTot;
    ppp20Mg->energy_p4 = Correl.Mg20.Sol[0]->Ekin;
    ppp20Mg->theta4 = Correl.Mg20.Sol[0]->theta;
    ppp20Mg->phi4 = Correl.Mg20.Sol[0]->phi;

    ppp20Mg->Erel = Erel_23P;
    //ppp20Mg->Ex = Ex; //No Ex
    ppp20Mg->Vcm = Correl.velocityCM;
    ppp20Mg->thetaCM = thetaCM;
    ppp20Mg->cos_thetaH = Correl.cos_thetaH;

    ppp20Mg->beamZ = Correl.Mg20.Sol[0]->beamZ;
    ppp20Mg->runnum = runnum;
    ppp20Mg->T->Fill();
  }

  //23P -> 22Si + p
  if(Correl.Si22.mult == 1 && Correl.proton.mult == 1)
  {
    //float const Q23P = mass_23P - (mass_22Si + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Si22.mask[0]=1;
    Correl.makeArray(1);
    corr_23P_counter++;
    //cout << "corr count " << corr_23P_counter << endl;
    float Erel_23P = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    //float Ex = Erel_22Si - Q22Si;


    Histo_sort->Erel_23P_p22Si->Fill(Erel_23P);
    //Histo_sort->Ex_23P_p22Si->Fill(Ex);

    Histo_sort->ThetaCM_23P_p22Si->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_23P_p22Si->Fill(Correl.velocityCM);
    Histo_sort->Erel_p22Si_costhetaH->Fill(Erel_23P,Correl.cos_thetaH);
    Histo_sort->p22Si_VCMvsErel->Fill(Correl.velocityCM,Erel_23P);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_23P_p22Si_trans->Fill(Erel_23P);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p22Si->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Si22.Sol[0]->velocity);

	    Histo_sort->p22Si_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p22Si_gammasADDvsErel->Fill(Erel_23P,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Si22.Sol[0]->velocity);

        Histo_sort->p22Si_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p22Si->Egamma[i] = doppE;
      p22Si->Tgamma[i] = Ceasar->added[i].time;
      p22Si->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->p22Si_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //Set tree variables
    p22Si->itele1 = Correl.proton.Sol[0]->itele;
    p22Si->id1 = Correl.proton.Sol[0]->iCsI;
    p22Si->ifront1 = Correl.proton.Sol[0]->ifront;
    p22Si->iback1 = Correl.proton.Sol[0]->iback;
    p22Si->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p22Si->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p22Si->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p22Si->et1 = Correl.proton.Sol[0]->energyTot;
    p22Si->time1 = Correl.proton.Sol[0]->CsITime;
    p22Si->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p22Si->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p22Si->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p22Si->theta1 = Correl.proton.Sol[0]->theta;
    p22Si->phi1 = Correl.proton.Sol[0]->phi;

    p22Si->M2[0] = Correl.Si22.Sol[0]->Mvect[0];
    p22Si->M2[1] = Correl.Si22.Sol[0]->Mvect[1];
    p22Si->M2[2] = Correl.Si22.Sol[0]->Mvect[2];
    p22Si->et2 = Correl.Si22.Sol[0]->energyTot;
    p22Si->energy_p2 = Correl.Si22.Sol[0]->Ekin;
    p22Si->theta2 = Correl.Si22.Sol[0]->theta;
    p22Si->phi2 = Correl.Si22.Sol[0]->phi;

    p22Si->Erel = Erel_23P;
    //p22Si->Ex = Ex; //No Ex
    p22Si->Vcm = Correl.velocityCM;
    p22Si->thetaCM = thetaCM;
    p22Si->cos_thetaH = Correl.cos_thetaH;

    p22Si->beamZ = Correl.Si22.Sol[0]->beamZ;
    p22Si->runnum = runnum;
    p22Si->T->Fill();

  }
}

void det::corr_24P()
{
  //24P -> 23Si + p
  if(Correl.Si23.mult == 1 && Correl.proton.mult == 1)
  {
    //float const Q23P = mass_23P - (mass_22Si + mass_p);
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.Si23.mask[0]=1;
    Correl.makeArray(1);
    //cout << "corr count " << corr_23P_counter << endl;
    float Erel_24P = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    //float Ex = Erel_22Si - Q22Si;


    Histo_sort->Erel_24P_p23Si->Fill(Erel_24P);
    //Histo_sort->Ex_23P_p22Si->Fill(Ex);

    Histo_sort->ThetaCM_24P_p23Si->Fill(thetaCM*180./acos(-1));
    Histo_sort->VCM_24P_p23Si->Fill(Correl.velocityCM);
    Histo_sort->Erel_p23Si_costhetaH->Fill(Erel_24P,Correl.cos_thetaH);
    Histo_sort->p23Si_VCMvsErel->Fill(Correl.velocityCM,Erel_24P);

    if (fabs(Correl.cos_thetaH) <= 0.1) Histo_sort->Erel_24P_p23Si_trans->Fill(Erel_24P);

    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    p23Si->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {

      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.Si23.Sol[0]->velocity);

	    Histo_sort->p23Si_gammasADD->Fill(doppE * 1000.); //E stored as MeV
	    Histo_sort->p23Si_gammasADDvsErel->Fill(Erel_24P,doppE * 1000.);

      //Make 2D gamma correlation plot
      for (int j=i+1;j<Ceasar->Nadded;j++) {
        xg = sin(Ceasar->added[j].theta)*cos(Ceasar->added[j].phi);
        yg = sin(Ceasar->added[j].theta)*sin(Ceasar->added[j].phi);
        zg = cos(Ceasar->added[j].theta);
        Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));
        Ldotg = (xH*xg) + (yH*yg) + (zH*zg);
        thetarel = acos(Ldotg/(EH*Eg));
        double doppE_2 = Doppler->correct(Ceasar->added[j].energy,thetarel,Correl.Si23.Sol[0]->velocity);

        Histo_sort->p23Si_gammasADDvsgammasADD->Fill(doppE * 1000.,doppE_2 * 1000.);
      }

      p23Si->Egamma[i] = doppE;
      p23Si->Tgamma[i] = Ceasar->added[i].time;
      p23Si->Chgamma[i] = Ceasar->added[i].id;

      //Tgate
      if (Ceasar->added[i].itime >= 2000 && Ceasar->added[i].itime <= 3700)
      {
 	      Histo_sort->p23Si_gammasADD_tgate->Fill(doppE * 1000.); //E stored as MeV
      }
    }

    //Set tree variables
    p23Si->itele1 = Correl.proton.Sol[0]->itele;
    p23Si->id1 = Correl.proton.Sol[0]->iCsI;
    p23Si->ifront1 = Correl.proton.Sol[0]->ifront;
    p23Si->iback1 = Correl.proton.Sol[0]->iback;
    p23Si->M1[0] = Correl.proton.Sol[0]->Mvect[0];
    p23Si->M1[1] = Correl.proton.Sol[0]->Mvect[1];
    p23Si->M1[2] = Correl.proton.Sol[0]->Mvect[2];
    p23Si->et1 = Correl.proton.Sol[0]->energyTot;
    p23Si->time1 = Correl.proton.Sol[0]->CsITime;
    p23Si->energy_p1 = Correl.proton.Sol[0]->Ekin;
    p23Si->denergy1_R = Correl.proton.Sol[0]->denergyR;
    p23Si->energy_p1_R = Correl.proton.Sol[0]->energyR;
    p23Si->theta1 = Correl.proton.Sol[0]->theta;
    p23Si->phi1 = Correl.proton.Sol[0]->phi;

    p23Si->M2[0] = Correl.Si23.Sol[0]->Mvect[0];
    p23Si->M2[1] = Correl.Si23.Sol[0]->Mvect[1];
    p23Si->M2[2] = Correl.Si23.Sol[0]->Mvect[2];
    p23Si->et2 = Correl.Si23.Sol[0]->energyTot;
    p23Si->energy_p2 = Correl.Si23.Sol[0]->Ekin;
    p23Si->theta2 = Correl.Si23.Sol[0]->theta;
    p23Si->phi2 = Correl.Si23.Sol[0]->phi;

    p23Si->Erel = Erel_24P;
    //p22Si->Ex = Ex; //No Ex
    p23Si->Vcm = Correl.velocityCM;
    p23Si->thetaCM = thetaCM;
    p23Si->cos_thetaH = Correl.cos_thetaH;

    p23Si->beamZ = Correl.Si23.Sol[0]->beamZ;
    p23Si->runnum = runnum;
    p23Si->T->Fill();

  }
}

//Alphas
void det::corr_a15O()
{
  //12N -> 11C + p
  if(Correl.O15.mult == 1 && Correl.alpha.mult == 1)
  {
    float const Q19Ne = mass_19Ne - (mass_15O + mass_alpha);
    Correl.zeroMask();
    Correl.alpha.mask[0]=1;
    Correl.O15.mask[0]=1;
    Correl.makeArray(1);

    float Erel_19Ne = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_19Ne - Q19Ne;

    
    //Euclidean vector length for heavy
    double EH = 0;
    if (Ceasar->Nadded > 0) EH = sqrt((xH*xH) + (yH*yH) + (zH*zH));

    a15O->Ngamma = Ceasar->Nadded; //Num gammas for tree

    for(int i = 0; i < Ceasar->Nadded; i++) {
      if (Ceasar->Nadded > 100) cout << "Huge caesar " << Ceasar->Nadded << endl;
      //Calculate unit vector for gamma
      double xg = sin(Ceasar->added[i].theta)*cos(Ceasar->added[i].phi);
      double yg = sin(Ceasar->added[i].theta)*sin(Ceasar->added[i].phi);
      double zg = cos(Ceasar->added[i].theta);

      //Euclidean vector lengths
      double Eg = sqrt((xg*xg) + (yg*yg) + (zg*zg));

      //Dot product
      double Ldotg = (xH*xg) + (yH*yg) + (zH*zg);

      //Find angle between LF and gamma
      float thetarel = acos(Ldotg/(EH*Eg));     

      //Find doppler corrected energy 
      double doppE;
      doppE = Doppler->correct(Ceasar->added[i].energy,thetarel,Correl.O15.Sol[0]->velocity);


      a15O->Egamma[i] = doppE;
      a15O->Tgamma[i] = Ceasar->added[i].time;
      a15O->Chgamma[i] = Ceasar->added[i].id;
    }

    //Set tree variables
    a15O->itele1 = Correl.alpha.Sol[0]->itele;
    a15O->id1 = Correl.alpha.Sol[0]->iCsI;
    a15O->ifront1 = Correl.alpha.Sol[0]->ifront;
    a15O->iback1 = Correl.alpha.Sol[0]->iback;
    a15O->M1[0] = Correl.alpha.Sol[0]->Mvect[0];
    a15O->M1[1] = Correl.alpha.Sol[0]->Mvect[1];
    a15O->M1[2] = Correl.alpha.Sol[0]->Mvect[2];
    a15O->et1 = Correl.alpha.Sol[0]->energyTot;
    a15O->time1 = Correl.alpha.Sol[0]->CsITime;
    a15O->energy_p1 = Correl.alpha.Sol[0]->Ekin;
    a15O->denergy1_R = Correl.alpha.Sol[0]->denergyR;
    a15O->energy_p1_R = Correl.alpha.Sol[0]->energyR;
    a15O->theta1 = Correl.alpha.Sol[0]->theta;
    a15O->phi1 = Correl.alpha.Sol[0]->phi;

    a15O->M2[0] = Correl.O15.Sol[0]->Mvect[0];
    a15O->M2[1] = Correl.O15.Sol[0]->Mvect[1];
    a15O->M2[2] = Correl.O15.Sol[0]->Mvect[2];
    a15O->et2 = Correl.O15.Sol[0]->energyTot;
    a15O->energy_p2 = Correl.O15.Sol[0]->Ekin;
    a15O->theta2 = Correl.O15.Sol[0]->theta;
    a15O->phi2 = Correl.O15.Sol[0]->phi;

    a15O->Erel = Erel_19Ne;
    a15O->Ex = Ex;
    a15O->Vcm = Correl.velocityCM;
    a15O->thetaCM = thetaCM;
    a15O->cos_thetaH = Correl.cos_thetaH;

    a15O->runnum = runnum;
    a15O->beamZ = Correl.O15.Sol[0]->beamZ;

    //cout << Correl.cos_thetaH << endl;

    //Fill tree
    a15O->T->Fill();

  }
}
