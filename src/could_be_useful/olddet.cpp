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
  Hira = new hira(ran,Histo_sort);
  //last value is the s800 setting, 0=Ca36, 1=K35
  s800 = new S800(ran,Histo_sort,setting); 
  Ceasar = new ceasar(ran,Histo_sort,Histo_read, 13);
  Doppler = new doppler(0.3323); //beta for Ca-36

  losses_fiber = new CLosses(20,"_fiber.loss",true);
  losses_target = new CLosses(20,"_Be.loss",true);

  Egamma0 = 0.;
  thetarel0 = 0.;
  thetagamma = 0.;
  thetares = 0.;
  res_Vel0 = 0.;
  detID = -1;
  Ring = -1;
  Loc = -1;
  Edop = 0.;

  Histo_sort->tree->Branch("Egamma0", &Egamma0);
  Histo_sort->tree->Branch("thetarel0", &thetarel0);
  Histo_sort->tree->Branch("thetagamma", &thetagamma);
  Histo_sort->tree->Branch("thetares", &thetares);
  Histo_sort->tree->Branch("res_Vel0", &res_Vel0);
  Histo_sort->tree->Branch("detID", &detID);
  Histo_sort->tree->Branch("Ring", &Ring);
  Histo_sort->tree->Branch("Loc", &Loc);
  Histo_sort->tree->Branch("Edop", &Edop);


  //string name_a8B("Wood_a8B.root");
  //Wood_a8B = new wood(2,false,&name_a8B);
}

//************************************************
/**
 * destructor
 */
det::~det()
{

  delete losses_fiber;
  delete losses_target;
  
  delete Hira;
  delete ran;
  delete s800;
  delete Ceasar;

  cout << "You made it!" << endl;
}

void det::Reset()
{
  Hira->reset();
  s800->Reset();
}

//*************************************************************
/**
 * unpacks a physics event from the data stream
 \param point0 - pointer to location in data stream
*/
bool det::unpack(unsigned short *point,int runno,int sourceID)
{

  unsigned short  words = *point++;
  // S Gillespie 2020/11/2
  // Defining a new value to store tdc position to pass to CEASAR
  unsigned short *tdcpoint;

  bool stat = false;

  //cout << "sourceID = " << sourceID << endl;
  if(sourceID == SiID)
  {
    point += 6; //This skips the SIS3820 information which is already in the timestamp
    stat = Hira->unpack(point, tdcpoint, runno); // SG 2020/11/02

    if(!stat)
    {
      //  cout << "Didn't read hira right" <<endl;
      return stat;
    }
    stat = Ceasar->unpack(point, tdcpoint, runno); // SG 2020/11/02
  }
  else if (sourceID == S800ID)
  {
    stat = s800->unpack(point,runno);
    //if(stat)
      //NS800++;
  }
  else
  {
    cout << "found unexpected sourceID = " << endl;
    return stat;
  }
  
  return stat;
}
//*********************************

void det::analyze(int event, int run)
{
  Correl.reset();

  bool foundresidue = false;

  s800_results S800_results;
  S800_results.Reset();
  // S800_results.trig_coin=false;
  // S800_results.trig_singles=false;
  // S800_results.trig_s800_singles=false;


  
  if (s800->Trig.registr & 1) S800_results.trig_s800_singles =true;
  if (s800->Trig.registr & 2) S800_results.trig_coin = true;
  if (s800->Trig.registr & 16) S800_results.trig_singles = true;
  

  //ND adding in a counter for these events
  //if (s800->Trig.registr & 1) N_s800_singles++;
  //if (s800->Trig.registr & 2) N_coin++;
  //if (s800->Trig.registr & 16) N_singles++;

  S800_results.Zbeam = -1;
  S800_results.Abeam = -1;
  S800_results.Zresidue = -1;
  S800_results.Aresidue = -1;

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
  release = Hira->RingCounter->analysis(S800_results);  // SG 2020/10/27

  
  if(s800->BeamID <0) return; //Need a good beam id
  
  S800_results.Zresidue = s800->residue_pid[s800->S800Setting][s800->BeamID]->Z;
  S800_results.Aresidue = s800->residue_pid[s800->S800Setting][s800->BeamID]->A;


  //ND, investigative counting

  if (S800_results.Zbeam == 20 && S800_results.Abeam==37)
  {
    if (S800_results.Zresidue==19 && S800_results.Aresidue==35)
    {

      if (release == 0)
        release0++;
      if (release == 3)
        release3++;
      if (release == 4)
        release4++;
      if (release == 5)
        release5++;
      if (release == 6)
        release6++;
    }
  }
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
  //fiber position maps for only K35 fragments
  if (S800_results.Zresidue==19 && S800_results.Aresidue==35)
  {
    //NK35residue++;t
    if(Hira->XY_mon->has_data)
    {
      //NK35_withfiber++;
      Histo_sort->VertHitMap_35Kgate->Fill(Hira->XY_mon->vert->pmtx,Hira->XY_mon->vert->pmty);
      Histo_sort->HorzHitMap_35Kgate->Fill(Hira->XY_mon->horz->pmtx,Hira->XY_mon->horz->pmty);
    }
  }

  if (S800_results.trig_coin) Histo_sort->S800_Csi_time->Fill(Hira->T_CsiTrig/10.);
  if (S800_results.trig_singles) Histo_sort->singles_trig_time->Fill(Hira->T_CsiTrig/10.);

  //  cout << "before fiber" << endl;
  //demand the fibre array gave information of both x and y position
  //if(!Hira->XY_mon->has_data)return;
  //cout << "after fiber" << endl;

  //cout << S800_results.Zbeam << " " << S800_results.Abeam << endl;
 
  if (S800_results.trig_coin && Hira->RingCounter->proton_present)
    Histo_sort->S800_Csi_time_with_proton->Fill(Hira->T_CsiTrig/10.);


  if (Hira->RingCounter->multAlpha == 1)
  {
    if((int)s800->mTDC.objCorrected.size()>0)
    {
      if((int)s800->mTDC.xfp.size()>0)
      {
        Histo_sort->ObjvsXFPwithAlpha1->Fill(s800->mTDC.objCorrected.at(0),
                                             s800->mTDC.xfp.at(0));
      }
    }
  }


  if (Hira->RingCounter->multProton == 1)
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
  }

  //relative angles for multProton == 1
  if (Hira->RingCounter->multProton == 1 && Hira->XY_mon->has_data &&
      S800_results.trig_coin)
  {
    float xp,yp,zp;
    bool good = false;
    for (int i=0;i<Hira->RingCounter->Nsolution;i++)
    {
      if (Hira->RingCounter->Solution[i].ipid == 1)
      {
        zp = cos(Hira->RingCounter->Solution[i].theta);
        xp = sin(Hira->RingCounter->Solution[i].theta)*
                  cos(Hira->RingCounter->Solution[i].phi);
        yp = sin (Hira->RingCounter->Solution[i].theta)*
                  sin(Hira->RingCounter->Solution[i].phi);
        good = true;
        exit;
      }
    }
    if (good)
    {
      float xr,yr,zr;
      zr = cos(Hira->XY_mon->theta);
      xr = sin(Hira->XY_mon->theta)*cos(Hira->XY_mon->phi);
      yr = sin(Hira->XY_mon->theta)*sin(Hira->XY_mon->phi);

      float dot = xp*xr + yp*yr + zp*zr;
      float deltaTheta = acos(dot)*180./acos(-1);

      //if (S800_results.Zbeam == 16 && S800_results.Abeam == 29
      //    && S800_results.Zresidue == 15 && S800_results.Aresidue == 27)
      //{
      //  Histo_sort->Si28_deltaTheta->Fill(deltaTheta);
      //}
    }
  }


  if(S800_results.Zresidue == 19 && S800_results.Aresidue == 35)
  {
    Histo_sort->CornerMult->Fill(Hira->XY_mon->vert->Mult + Hira->XY_mon->horz->Mult);
    Histo_sort->CornerMult_Vert->Fill(Hira->XY_mon->vert->Mult);
    Histo_sort->CornerMult_Hori->Fill(Hira->XY_mon->horz->Mult);
    Histo_sort->CornerMult_VvsH->Fill(Hira->XY_mon->horz->Mult,Hira->XY_mon->vert->Mult);
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
  
  //  if(S800_results.Zbeam ==18 && S800_results.Abeam ==31)
  if(S800_results.Zbeam >0 && S800_results.Abeam >0)
  {
    if(S800_results.Zresidue==19 && S800_results.Aresidue ==35)
    {
      Histo_sort->ThetavsPhi_fiber->Fill(Hira->XY_mon->thetadeg,Hira->XY_mon->phideg);
      Histo_sort->ThetaFibervsThetaS800->Fill(Hira->XY_mon->thetadeg,s800->track->thetadeg);
      Histo_sort->PhiFibervsPhiS800->Fill(Hira->XY_mon->phideg,s800->track->phideg);
    }
    //cout << "looking for residue" << endl;
    //Lets add the S800 to the Solution Class
    if(S800_results.Zresidue >0 && S800_results.Aresidue >0)
    {
      foundresidue = true;
      foundresidue = LoadS800toSolution();
      //cout << "found residue" << endl;
      // if(S800_results.Zresidue==16 && S800_results.Aresidue==28)
      //   {
      //     cout << "mult proton = " << Hira->RingCounter->multProton << endl;
      //   }
    }
  }
  
  if(foundresidue)
  {
    int Nsol = Hira->RingCounter->Nsolution;

    for (int n=0;n<Hira->RingCounter->Csi.Nstore;n++)
    {
      //cout << "n " << n << endl;
      //cout << " energy " << Hira->RingCounter->Csi.Order[n].energyR << "    time " <<  Hira->RingCounter->Csi.Order[n].time << endl;
      Histo_sort->ET_csi_res->Fill(Hira->RingCounter->Csi.Order[n].energyR,Hira->RingCounter->Csi.Order[n].time/10.);
    }


    float res_Vel = Hira->RingCounter->Solution[Nsol-1].velocity;
    float phiR = Hira->RingCounter->Solution[Nsol-1].phi;//s800->track->phi;//Hira->RingCounter->Solution[Nsol-1].phi;
    float thetaR = Hira->RingCounter->Solution[Nsol-1].theta;//s800->track->theta;//Hira->RingCounter->Solution[Nsol-1].theta;
    float xR = sin(thetaR)*cos(phiR);
    float yR = sin(thetaR)*sin(phiR);
    float zR = cos(thetaR);

    //track velocity distributions of all heavy fragments
    //velocities should be similar to beam velocity
    if(S800_results.Zbeam == 20 && S800_results.Abeam == 37)
    {
      if(S800_results.Zresidue == 20 && S800_results.Aresidue == 38)
        Histo_sort->Vlab_Ca38->Fill(res_Vel);
      if(S800_results.Zresidue == 20 && S800_results.Aresidue == 37)
        Histo_sort->Vlab_Ca37->Fill(res_Vel);      

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
      if(S800_results.Zresidue == 18 && S800_results.Aresidue == 32)
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
    }

    //This section has gamma spectra with no add back
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

        //gamma spec for K35 s800 setting
        //K35        
        if(S800_results.Zresidue == 20 && S800_results.Aresidue == 37)
          if(dopp > 0)
          {                   
            Histo_read->TEC_Ca37->Fill(dopp);
            Histo_read->dopp_timing_Ca37->Fill(dopp, Ceasar->added[i].time);
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

            Histo_read->TEC_Ca36_simpletheta->Fill(dopp_simple);
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
            Histo_read->gamma_K36_K36beam->Fill(dopp);
        }
      }

    }

  }

  //////// end of gamma code

 
  //cout << "found residue = " << foundresidue << " " << S800_results.Zbeam << endl;

  int Mult = 0;
  
  
  for(int i=0;i<Hira->RingCounter->Nsolution;i++)
  {
    if(Hira->RingCounter->Solution[i].iZ >0)
    {
      Correl.load(&Hira->RingCounter->Solution[i]);
      Mult++;
      solnZ++;
    }
    else
    {
      solnnoZ++;
    }
  }
  
  
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
  }

  if(foundresidue)
  {
    corr_36Ca();
    corr_37Ca();
    corr_37Sc();




    //2p emitters
    corr_38Ti2p();
    corr_37Sc2p();
    corr_35Ca2p();
    corr_36Ca2p();

    corr_34K();
    corr_35K();
    corr_36K();

    corr_33Ar();
    corr_34Ar();
    corr_35Ar();

    corr_32Cl();
    corr_33Cl();
    corr_34Cl();

    corr_31S();
    corr_32S();
    corr_33S();

    corr_29P();


  }
}


bool det::LoadS800toSolution()
{

  int Nsol = Hira->RingCounter->Nsolution;
  Hira->RingCounter->Nsolution++;

  int Z = s800->residue_pid[s800->S800Setting][s800->BeamID]->Z;
  int A = s800->residue_pid[s800->S800Setting][s800->BeamID]->A;

  Hira->RingCounter->Solution[Nsol].iZ = Z;
  Hira->RingCounter->Solution[Nsol].iA = A;

  //cout << "theta fiber = " << Hira->XY_mon->theta << " theta S800 = " << s800->track->theta << " ratio = " << Hira->XY_mon->theta/s800->track->theta << endl;

  //using the theta phi from the fiber
  float thetaf = 0;
  float phif = 0;

  //cout << Hira->XY_mon->theta << endl;
  if(Hira->XY_mon->theta !=-999)
  {
    thetaf =  Hira->XY_mon->theta;
    Hira->RingCounter->Solution[Nsol].theta = thetaf;
    phif = Hira->XY_mon->phi;
    Hira->RingCounter->Solution[Nsol].phi = phif;
  }
  else
  {
    thetaf = s800->track->theta;
    
    phif = s800->track->phi;
    return false;
  }
  
  double mass = Hira->RingCounter->getMass(Z,A);

  float ekin = s800->track->energy;
  float pc = sqrt(pow(ekin+mass*m0,2) - pow(mass*m0,2));
  float rigidity = pc/1000.*3.3356/Z;

  if (Z == 19 && A == 35)
  {
    Histo_sort->rigidityK35->Fill(rigidity);
    Histo_read->Vlab_HF_p35K_before->Fill(pc/(ekin + mass*m0));
  }
  if (Z == 19 && A == 36)
  {
    Histo_read->Vlab_HF_p36K_before->Fill(pc/(ekin + mass*m0));
  }

  float thickness = 51.1/cos(thetaf);//*1.30; // NB: 0.5mm -> 51.1 mg/cm2
  ekin = losses_fiber->getEin(ekin,thickness,Z,A);
  //changed to TargetThickness/2 -ND 8/17/21
  thickness = Hira->RingCounter->TargetThickness/2/cos(thetaf);//*1.80;

  ekin = losses_target->getEin(ekin,thickness,Z,A);
  Hira->RingCounter->Solution[Nsol].Ekin = ekin;
  
  Hira->RingCounter->Solution[Nsol].mass = mass*m0; //mass*(amu->MeV)

  float momentum = Hira->RingCounter->Solution[Nsol].Kinematics.
    getMomentum(Hira->RingCounter->Solution[Nsol].Ekin,Hira->RingCounter->Solution[Nsol].mass);

  Hira->RingCounter->Solution[Nsol].momentum = momentum;

  Hira->RingCounter->Solution[Nsol].Mvect[0] = momentum*sin(thetaf)*cos(phif);
  Hira->RingCounter->Solution[Nsol].Mvect[1] = momentum*sin(thetaf)*sin(phif);
  Hira->RingCounter->Solution[Nsol].Mvect[2] = momentum*cos(thetaf);



  Hira->RingCounter->Solution[Nsol].energyTot = Hira->RingCounter->Solution[Nsol].mass
                                            + Hira->RingCounter->Solution[Nsol].Ekin;

  Hira->RingCounter->Solution[Nsol].velocity = momentum/Hira->RingCounter->Solution[Nsol].energyTot;

  return true;
}


/////////////////////////////////////////////////
//Correlation functions
/////////////////////////////////////////////////


void det::corr_36Ca()
{
  // p + 35K
  if(Correl.proton.mult == 1 && Correl.K35.mult == 1)
  {
 
    float const Q36Ca = -2.567; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.K35.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_36Ca = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_36Ca - Q36Ca;
    //angle between beam axis and the momentum vector between the center
    // of mass and the heavy fragment
    float mag2 = 0;
    for (int j=0;j<3;j++) mag2 += pow(Correl.frag[1]->MomCM[j],2);
    float cosbeamCMtoHF = Correl.frag[1]->MomCM[2]/sqrt(mag2);
    //cout << "cosbeamCMtoHF " << cosbeamCMtoHF << endl;

    Histo_read->Erel_36Ca_p35K->Fill(Erel_36Ca);
    Histo_read->Erel_36Ca_p35K_lowres->Fill(Erel_36Ca);
    Histo_read->Ex_36Ca_p35K->Fill(Ex);
    Histo_read->Ex_36Ca_p35K_lowres->Fill(Ex);    
    Histo_read->ThetaCM_36Ca_p35K->Fill(thetaCM*180./acos(-1));
    Histo_read->VCM_36Ca_p35K->Fill(Correl.velocityCM);
    Histo_read->Vlab_LF_p35K->Fill(Correl.proton.Sol[0]->velocity);
    Histo_read->Vlab_HF_p35K->Fill(Correl.K35.Sol[0]->velocity);
    Histo_read->Vlab_LF_p35K_before->Fill(Hira->RingCounter->velocity_before);
    Histo_read->cosbeamCMtoHF_Ex_p35K->Fill(Ex, cosbeamCMtoHF);
    Histo_sort->CRDC1X_K35->Fill(s800->CRDC[0].x_gravity);

    //printing out s800 variables
    //cout << "x_gravity " << s800->CRDC[1].x_gravity << " " << s800->CRDC[0].x_gravity << endl;
 
    /*
    //printing out particle info
    cout << "in correl velocity prot " << Correl.proton.Sol[0]->velocity << endl;
    cout << "in correl momentum prot " << Correl.proton.Sol[0]->momentum << endl;
    cout << "in correl energyTot prot " << Correl.proton.Sol[0]->energyTot << endl;
    cout << "in correl vel calc prot " << Correl.proton.Sol[0]->momentum/Correl.proton.Sol[0]->energyTot << endl;

    cout << "in correl velocity K35 " << Correl.K35.Sol[0]->velocity << endl;
    cout << "in correl momentum K35 " << Correl.K35.Sol[0]->momentum << endl;
    cout << "in correl energyTot K35 " << Correl.K35.Sol[0]->energyTot << endl;
    cout << "in correl vel calc K35 " << Correl.K35.Sol[0]->momentum/Correl.K35.Sol[0]->energyTot << endl;

    */  

  }
}

void det::corr_37Ca()
{
  // p + 36K
  if(Correl.proton.mult == 1 && Correl.K36.mult == 1)
  {
    float const Q37Ca = -3.008; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.K36.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_37Ca = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_37Ca - Q37Ca;

    //angle between beam axis and the momentum vector between the center
    // of mass and the heavy fragment
    float mag2 = 0;
    for (int j=0;j<3;j++) mag2 += pow(Correl.frag[1]->MomCM[j],2);
    float cosbeamCMtoHF = Correl.frag[1]->MomCM[2]/sqrt(mag2);
    //cout << "cosbeamCMtoHF " << cosbeamCMtoHF << endl;

    Histo_read->Erel_37Ca_p36K->Fill(Erel_37Ca);
    Histo_read->Erel_37Ca_p36K_lowres->Fill(Erel_37Ca);
    Histo_read->Ex_37Ca_p36K->Fill(Ex);
    Histo_read->Ex_37Ca_p36K_lowres->Fill(Ex);    
    Histo_read->ThetaCM_37Ca_p36K->Fill(thetaCM*180./acos(-1));
    Histo_read->VCM_37Ca_p36K->Fill(Correl.velocityCM);
    Histo_read->Vlab_LF_p36K->Fill(Correl.proton.Sol[0]->velocity);
    Histo_read->Vlab_HF_p36K->Fill(Correl.K36.Sol[0]->velocity);
    Histo_read->Vlab_LF_p36K_before->Fill(Hira->RingCounter->velocity_before);

    Histo_read->cosbeamCMtoHF_Ex_p36K->Fill(Ex, cosbeamCMtoHF);

    Histo_sort->CRDC1X_K36->Fill(s800->CRDC[0].x_gravity);

/*
    //cout << "x_gravity " << s800->CRDC[1].x_gravity << " " << s800->CRDC[0].x_gravity << endl;

    //printing out particle info
    cout << "in correl velocity prot " << Correl.proton.Sol[0]->velocity << endl;
    cout << "in correl momentum prot " << Correl.proton.Sol[0]->momentum << endl;
    cout << "in correl energyTot prot " << Correl.proton.Sol[0]->energyTot << endl;
    cout << "in correl vel calc prot " << Correl.proton.Sol[0]->momentum/Correl.proton.Sol[0]->energyTot << endl;

    cout << "in correl velocity K36 " << Correl.K36.Sol[0]->velocity << endl;
    cout << "in correl momentum K36 " << Correl.K36.Sol[0]->momentum << endl;
    cout << "in correl energyTot K36 " << Correl.K36.Sol[0]->energyTot << endl;
    cout << "in correl vel calc K36 " << Correl.K36.Sol[0]->momentum/Correl.K36.Sol[0]->energyTot << endl;
*/

  }
}

void det::corr_37Sc()
{
  // p + 36Ca
  if(Correl.proton.mult == 1 && Correl.Ca36.mult == 1)
  {
    float const Q37Sc = 2.682; //Ame2016, theory estimate for 37Sc mass
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Ca36.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_37Sc = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_37Sc - Q37Sc;

    Histo_read->Erel_37Sc_p36Ca->Fill(Erel_37Sc);
    Histo_read->Ex_37Sc_p36Ca->Fill(Ex);
    Histo_read->ThetaCM_37Sc_p36Ca->Fill(thetaCM*180./acos(-1));
    Histo_read->VCM_37Sc_p36Ca->Fill(Correl.velocityCM);
  }
}


////////////////////////////
//K isotopes
void det::corr_34K()
{
  // p + 33Ar
  if(Correl.proton.mult == 1 && Correl.Ar33.mult == 1)
  {
    float const Q34K = 0.875; //Ame2016, theory estimates for 34K
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Ar33.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q34K;

    Histo_read->Erel_34K_p33Ar->Fill(Erel);
    Histo_read->Ex_34K_p33Ar->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}
void det::corr_35K()
{
  // p + 34Ar
  if(Correl.proton.mult == 1 && Correl.Ar34.mult == 1)
  {
    float const Q35K = -83.57; //Ame2016
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Ar34.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q35K;

    Histo_read->Erel_35K_p34Ar->Fill(Erel);
    Histo_read->Ex_35K_p34Ar->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}
void det::corr_36K()
{
  // p + 35Ar
  if(Correl.proton.mult == 1 && Correl.Ar35.mult == 1)
  {
    float const Q36K = -1.658; //Ame2016
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Ar35.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q36K;

    Histo_read->Erel_36K_p35Ar->Fill(Erel);
    Histo_read->Ex_36K_p35Ar->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}
//////////////////////////////////
//Ar isotopes
void det::corr_33Ar()
{
  // p + 32Cl
  if(Correl.proton.mult == 1 && Correl.Cl32.mult == 1)
  {
    float const Q33Ar = -3.338; //Ame2016
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Cl32.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q33Ar;

    Histo_read->Erel_33Ar_p32Cl->Fill(Erel);
    Histo_read->Ex_33Ar_p32Cl->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}
void det::corr_34Ar()
{
  // p + 33Cl
  if(Correl.proton.mult == 1 && Correl.Cl33.mult == 1)
  {
    float const Q34Ar = -4.664; //Ame2016
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Cl33.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q34Ar;

    Histo_read->Erel_34Ar_p33Cl->Fill(Erel);
    Histo_read->Ex_34Ar_p33Cl->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}
void det::corr_35Ar()
{
  // p + 34Cl
  if(Correl.proton.mult == 1 && Correl.Cl34.mult == 1)
  {
    float const Q35Ar = -5.896; //Ame2016
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Cl34.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q35Ar;

    Histo_read->Erel_35Ar_p34Cl->Fill(Erel);
    Histo_read->Ex_35Ar_p34Cl->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}

/////////////////////////////////////
//Cl isotopes
void det::corr_32Cl()
{
  // p + 31S
  if(Correl.proton.mult == 1 && Correl.S31.mult == 1)
  {
    float const Q32Cl = -1.581; //Ame2016
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.S31.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q32Cl;

    Histo_read->Erel_32Cl_p31S->Fill(Erel);
    Histo_read->Ex_32Cl_p31S->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}
void det::corr_33Cl()
{
  // p + 32S
  if(Correl.proton.mult == 1 && Correl.S32.mult == 1)
  {
    float const Q = -2.276; //Ame2016
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.S32.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q;

    Histo_read->Erel_33Cl_p32S->Fill(Erel);
    Histo_read->Ex_33Cl_p32S->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}
void det::corr_34Cl()
{
  // p + 33S
  if(Correl.proton.mult == 1 && Correl.S33.mult == 1)
  {
    float const Q = -5.143; //Ame2016
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.S33.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q;

    Histo_read->Erel_34Cl_p33S->Fill(Erel);
    Histo_read->Ex_34Cl_p33S->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}

/////////////////////////////////
//S isotopes
void det::corr_31S()
{
  // p + 30P
  if(Correl.proton.mult == 1 && Correl.P30.mult == 1)
  {
    float const Q = 0; //Ame2016
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.P30.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q;

    Histo_read->Erel_31S_p30P->Fill(Erel);
    Histo_read->Ex_31S_p30P->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}
void det::corr_32S()
{
  // p + 31P
  if(Correl.proton.mult == 1 && Correl.P31.mult == 1)
  {
    float const Q = 0; //Ame2016
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.P31.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q;

    Histo_read->Erel_32S_p31P->Fill(Erel);
    Histo_read->Ex_32S_p31P->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}
void det::corr_33S()
{
  // p + 32P
  if(Correl.proton.mult == 1 && Correl.P32.mult == 1)
  {
    float const Q = 0; //Ame2016
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.P31.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q;

    Histo_read->Erel_33S_p32P->Fill(Erel);
    Histo_read->Ex_33S_p32P->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}

////////////////////////////////////
//P isotopes
void det::corr_29P()
{
  // p + 28Si
  if(Correl.proton.mult == 1 && Correl.Si28.mult == 1)
  {
    float const Q = 0; //Ame2016
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Si28.mask[0]=1;   
    Correl.makeArray(1);
    float Erel = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel - Q;

    Histo_read->Erel_29P_p28Si->Fill(Erel);
    Histo_read->Ex_29P_p28Si->Fill(Ex);
    //Histo_read->ThetaCM_34K_p33Ar->Fill(thetaCM*180./acos(-1));
    //Histo_read->VCM_34K_p33Ar->Fill(Correl.velocityCM);
  }
}


////////////////////////////////
//2p emitters
void det::corr_38Ti2p()
{
  //2p + 36Ca
  if(Correl.proton.mult ==2 && Correl.Ca36.mult ==1)
  {
 
    float const Q38Ti2p = 2.74; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.Ca36.mask[0]=1;   
    Correl.makeArray(1);

    float Erel_38Ti = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_38Ti - Q38Ti2p;
    //      cout << "Erel = " << Erel_30Ar << endl;
    Histo_read->Erel_38Ti_2p36Ca->Fill(Erel_38Ti);
    Histo_read->Ex_38Ti_2p36Ca->Fill(Ex);
    Histo_read->ThetaCM_38Ti_2p36Ca->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_38Ti_2p36Ca->Fill(vCM);


    /*
    //Get relative energy between the core and protons
    Correl.proton.mask[0]=0;
    if (S800_results.Zbeam == 20 && S800_results.Abeam==37)
    {
      Correl.makeArray(1);
      float ExY1 = Correl.findErel()/Erel_30Ar;

      Correl.proton.mask[0]=1;
      Correl.proton.mask[1]=0;
      Correl.makeArray(1);
      float ExY2 = Correl.findErel()/Erel_30Ar;

      //Get the relative energy betweeen the protons
      Correl.proton.mask[1]=1;
      Correl.S29.mask[0]=0;   
      Correl.makeArray(1);
      float ExT  = Correl.findErel()/Erel_30Ar;

      Correl.S29.mask[0] =1;
      Correl.makeArray(1);
      Correl.findErel();
      Correl.getJacobi();

      //Gate on Ground state
      if(Erel_30Ar <1.0)
      {
        Histo_read->Jacobi_T_30Ar_gs->Fill(ExT,Correl.cosThetaT);
        Histo_read->Jacobi_T_30Ar_gs->Fill(ExT,-Correl.cosThetaT);   
        Histo_read->Jacobi_Y_30Ar_gs->Fill(ExY1,Correl.cosThetaY[0]);
        Histo_read->Jacobi_Y_30Ar_gs->Fill(ExY2,Correl.cosThetaY[1]);
      } 
    }
    */
  }
}
void det::corr_37Sc2p()
{
  //2p + 35K
  if(Correl.proton.mult ==2 && Correl.K35.mult ==1)
  {
 
    float const Q37Sc2p = -4.606; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.K35.mask[0]=1;   
    Correl.makeArray(1);

    float Erel_37Sc = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_37Sc - Q37Sc2p;
    //      cout << "Erel = " << Erel_30Ar << endl;
    Histo_read->Erel_37Sc_2p35K->Fill(Erel_37Sc);
    Histo_read->Ex_37Sc_2p35K->Fill(Ex);
    Histo_read->ThetaCM_37Sc_2p35K->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_37Sc_2p35K->Fill(vCM);
  }
}
void det::corr_35Ca2p()
{
  //2p + 33Ar
  if(Correl.proton.mult ==2 && Correl.Ar33.mult ==1)
  {
 
    float const Q35Ca2p = -0.406; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.Ar33.mask[0]=1;   
    Correl.makeArray(1);

    float Erel_35Ca = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_35Ca - Q35Ca2p;
    //      cout << "Erel = " << Erel_30Ar << endl;
    Histo_read->Erel_35Ca_2p33Ar->Fill(Erel_35Ca);
    Histo_read->Ex_35Ca_2p33Ar->Fill(Ex);
    Histo_read->ThetaCM_35Ca_2p33Ar->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_35Ca_2p33Ar->Fill(vCM);
  }
}

void det::corr_36Ca2p()
{
  //2p + 34Ar
  if(Correl.proton.mult ==2 && Correl.Ar34.mult ==1)
  {
 
    float const Q36Ca2p = -2.651; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.Ar34.mask[0]=1;   
    Correl.makeArray(1);

    float Erel_36Ca = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_36Ca - Q36Ca2p;
    //      cout << "Erel = " << Erel_30Ar << endl;
    Histo_read->Erel_36Ca_2p34Ar->Fill(Erel_36Ca);
    Histo_read->Ex_36Ca_2p34Ar->Fill(Ex);
    Histo_read->ThetaCM_36Ca_2p34Ar->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_36Ca_2p34Ar->Fill(vCM);
  }
}



/*

void det::corr_6Li()
{
  //d + alpha
  if(Correl.alpha.mult ==1 && Correl.H2.mult ==1)
  {
 
    float const Q6Li = 1.4738;
    Correl.zeroMask();
    Correl.alpha.mask[0]=1;
    Correl.H2.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_6Li = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_6Li + Q6Li;
    Histo_read->Erel_6Li->Fill(Erel_6Li);
    Histo_read->Ex_6Li_da->Fill(Ex);

    float vCM = Correl.velocityCM;
    Histo_read->vel_6Li->Fill(vCM);

    Histo_read->Costheta_6Li->Fill(thetaCM);
  }
}

void det::corr_30Ar()
{
  //2p + 28S
  if(Correl.proton.mult ==2 && Correl.S28.mult ==1)
  {
 
    float const Q30Ar = -2.28; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.S28.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_30Ar = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_30Ar + Q30Ar;
    //      cout << "Erel = " << Erel_30Ar << endl;
    Histo_read->Erel_30Ar_2p28S->Fill(Erel_30Ar);
    Histo_read->Ex_30Ar_2p28S->Fill(Ex);
    Histo_read->ThetaCM_30Ar_2p28S->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_30Ar_2p28S->Fill(vCM);


      
    //Get relative energy between the core and protons
    Correl.proton.mask[0]=0;if (S800_results.Zbeam == 20 && S800_results.Abeam==37)
  {
    Correl.makeArray(1);
    float ExY1 = Correl.findErel()/Erel_30Ar;

    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=0;
    Correl.makeArray(1);
    float ExY2 = Correl.findErel()/Erel_30Ar;

    //Get the relative energy betweeen the protons
    Correl.proton.mask[1]=1;
    Correl.S29.mask[0]=0;   
    Correl.makeArray(1);
    float ExT  = Correl.findErel()/Erel_30Ar;

    Correl.S29.mask[0] =1;
    Correl.makeArray(1);
    Correl.findErel();
    Correl.getJacobi();

    //Gate on Ground state
    if(Erel_30Ar <1.0)
    {
      Histo_read->Jacobi_T_30Ar_gs->Fill(ExT,Correl.cosThetaT);
      Histo_read->Jacobi_T_30Ar_gs->Fill(ExT,-Correl.cosThetaT);   
      Histo_read->Jacobi_Y_30Ar_gs->Fill(ExY1,Correl.cosThetaY[0]);
      Histo_read->Jacobi_Y_30Ar_gs->Fill(ExY2,Correl.cosThetaY[1]);
    } 
  }
}

void det::corr_28S()
{  float ekin = s800->track->energy;
  float thickness = 51.1/cos(thetaf); // NB: 0.5mm -> 51.1 mg/cm2

  ekin = losses_fiber->getEin(ekin,thickness,Z,A);
  //changed to TargetThickness/2 -ND 8/17/21
  thickness = Hira->RingCounter->TargetThickness/2/cos(thetaf);


  ekin = losses_target->getEin(ekin,thickness,Z,A);

  Hira->RingCounter->Solution[Nsol].Ekin = ekin;

  
  double mass = Hira->RingCounter->getMass(Z,A);
  Hira->RingCounter->Solution[Nsol].mass = mass*m0; //mass*(amu->MeV)
  cout << "mass " << Hira->RingCounter->Solution[Nsol].mass << endl;

  float momentum = Hira->RingCounter->Solution[Nsol].Kinematics.
    getMomentum(Hira->RingCounter->Solution[Nsol].Ekin,Hira->RingCounter->Solution[Nsol].mass);

  Hira->RingCounter->Solution[Nsol].momentum = momentum;

  
  Hira->RingCounter->Solution[Nsol].Mvect[0] = momentum*sin(thetaf)*cos(phif);
  Hira->RingCounter->Solution[Nsol].Mvect[1] = momentum*sin(thetaf)*sin(phif);
  Hira->RingCounter->Solution[Nsol].Mvect[2] = momentum*cos(thetaf);



  Hira->RingCounter->Solution[Nsol].energyTot = Hira->RingCounter->Solution[Nsol].mass + 
  Hira->RingCounter->Solution[Nsol].Ekin;
  Hira->RingCounter->Solution[Nsol].velocity = momentum/Hira->RingCounter->Solution[Nsol].energyTot;
  // p + 27P
  if(Correl.proton.mult ==1 && Correl.P27.mult ==1)
  {
 
    float const Q28S = -2.49; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.P27.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_28S = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_28S + Q28S;
    Histo_read->Erel_28S_1p27P->Fill(Erel_28S);
    Histo_read->Ex_28S_1p27P->Fill(Ex);
    Histo_read->ThetaCM_28S_1p27P->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_28S_1p27P->Fill(vCM);

  }

}

void det::corr_29S()
{
  // p + 28P
  if(Correl.proton.mult ==1 && Correl.P28.mult ==1)
  {
 
    float const Q29S = -3.30; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.P28.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_29S = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_29S + Q29S;
    Histo_read->Erel_29S_1p28P->Fill(Erel_29S);
    Histo_read->Ex_29S_1p28P->Fill(Ex);
    Histo_read->ThetaCM_29S_1p28P->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_29S_1p28P->Fill(vCM);

  }

}
// void det::corr_29S()
// {
//   // p + 28P
//   if(Correl.proton.mult ==1 && Correl.P28.mult ==1)
//     {
 
//       float const Q29S = -3.30; //Taken from Ame2016 
//       Correl.zeroMask();
//       Correl.proton.mask[0]=1;
//       //      Correl.proton.mask[1]=1;
//       Correl.P28.mask[0]=1;   
//       Correl.makeArray(1);
//       float Erel_29S = Correl.findErel();
//       float thetaCM = Correl.thetaCM;
//       float Ex = Erel_29S + Q29S;
//       Histo_read->Erel_29S_1p28P->Fill(Erel_29S);
//       Histo_read->Ex_29S_1p28P->Fill(Ex);
//       Histo_read->ThetaCM_29S_1p28P->Fill(thetaCM*180./acos(-1));
//       float vCM = Correl.velocityCM;
//       Histo_read->VCM_29S_1p28P->Fill(vCM);

//     }

// }

void det::corr_28P() ////
{
  // p + 27Si
  if(Correl.proton.mult ==1 && Correl.Si27.mult ==1)
  {
    float const Q28P = -2.052; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Si27.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_28P = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_28P + Q28P;
    Histo_read->Erel_28P_1p27Si->Fill(Erel_28P);
    Histo_read->Ex_28P_1p27Si->Fill(Ex);
    Histo_read->ThetaCM_28P_1p27Si->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_28P_1p27Si->Fill(vCM);

  }

}

void det::corr_29Cl()
{
  // p + 28S
  if(Correl.proton.mult ==1 && Correl.S28.mult ==1)
  {
 
    float const Q29Cl = 1.8; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.S28.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_29Cl = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_29Cl + Q29Cl;
    Histo_read->Erel_29Cl_1p28S->Fill(Erel_29Cl);
    Histo_read->Ex_29Cl_1p28S->Fill(Ex);
    Histo_read->ThetaCM_29Cl_1p28S->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_29Cl_1p28S->Fill(vCM);

  }

}


void det::corr_30Cl()
{
  //p + 29S
  if(Correl.proton.mult ==1 && Correl.S29.mult ==1)
  {
 
    float const Q30Cl = -0.31; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.S29.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_30Cl = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_30Cl + Q30Cl;
    // cout << "Erel = " << Erel_30Cl << endl;
    Histo_read->Erel_30Cl_p29S->Fill(Erel_30Cl);
    Histo_read->Ex_30Cl_p29S->Fill(Ex);
    Histo_read->ThetaCM_30Cl_p29S->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_30Cl_p29S->Fill(vCM);


      
      
  }
}


void det::corr_31Cl()
{
  //p + 30S
  //  cout << Correl.proton.mult << "   " << Correl.S30.mult << endl;
  if(Correl.proton.mult ==1 && Correl.S30.mult ==1)
  {
    float const Q31Cl = -0.264; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.S30.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_31Cl = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_31Cl + Q31Cl;
    Histo_read->Erel_31Cl_1p30S->Fill(Erel_31Cl);
    Histo_read->Ex_31Cl_1p30S->Fill(Ex);
    Histo_read->ThetaCM_31Cl_1p30S->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_31Cl_1p30S->Fill(vCM);
      
  }
}

void det::corr_30S()
{
  // p + 29P
  if(Correl.proton.mult ==1 && Correl.P29.mult ==1)
  {
 
    float const Q30S = -4.395; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.P29.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_30S = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_30S + Q30S;
    Histo_read->Erel_30S_1p29P->Fill(Erel_30S);
    Histo_read->Ex_30S_1p29P->Fill(Ex);
    Histo_read->ThetaCM_30S_1p29P->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_30S_1p29P->Fill(vCM);

  }

}

void det::corr_27P() ////
{
  // p + 26Si
  if(Correl.proton.mult ==1 && Correl.Si26.mult ==1)
  {
    float const Q27P = -0.87; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Si26.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_27P = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_27P + Q27P;
    Histo_read->Erel_27P_1p26Si->Fill(Erel_27P);
    Histo_read->Ex_27P_1p26Si->Fill(Ex);
    Histo_read->ThetaCM_27P_1p26Si->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_27P_1p26Si->Fill(vCM);

  }

}

void det::corr_29P() 
{
  // p + 28Si
  if(Correl.proton.mult ==1 && Correl.Si28.mult ==1)
  {
    float const Q29P = -2.749; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Si28.mask[0]=1;   << "in correl velocity1 " << Correl.frag[0].velocity << endl;
    Correl.makeArray(1);
    float Erel_29P = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_29P + Q29P;
    Histo_read->Erel_29P_1p28Si->Fill(Erel_29P);
    Histo_read->Ex_29P_1p28Si->Fill(Ex);
    Histo_read->ThetaCM_29P_1p28Si->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_29P_1p28Si->Fill(vCM);

  }

}

void det::corr_29P_2p() 
{
  //2p + 27Al
  if(Correl.proton.mult ==2 && Correl.Al27.mult ==1)
  {
    float const Q29P_2p = -14.334; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.Al27.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_29P_2p = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_29P_2p + Q29P_2p;
    Histo_read->Erel_29P_2p27Al->Fill(Erel_29P_2p);
    Histo_read->Ex_29P_2p27Al->Fill(Ex);
    Histo_read->ThetaCM_29P_2p27Al->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_29P_2p27Al->Fill(vCM);

  } 

}

void det::corr_27P_2p() 
{
  //2p + 25Al
  if(Correl.proton.mult ==2 && Correl.Al25.mult ==1)
  {
    float const Q27P_2p = -6.383; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    Correl.proton.mask[1]=1;
    Correl.Al25.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_27P_2p = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_27P_2p + Q27P_2p;
    Histo_read->Erel_27P_2p25Al->Fill(Erel_27P_2p);
    Histo_read->Ex_27P_2p25Al->Fill(Ex);
    Histo_read->ThetaCM_27P_2p25Al->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_27P_2p25Al->Fill(vCM);

  } 

}

void det::corr_28Si()
{
  // p + 27Al
  if(Correl.proton.mult ==1 && Correl.Al27.mult ==1)
  {
 
    float const Q28Si = -11.584; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Al27.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_28Si = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_28Si + Q28Si;
    Histo_read->Erel_28Si_1p27Al->Fill(Erel_28Si);
    Histo_read->Ex_28Si_1p27Al->Fill(Ex);
    Histo_read->ThetaCM_28Si_1p27Al->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_28Si_1p27Al->Fill(vCM);

  }

}

void det::corr_26Si()
{
  // p + 25Al
  if(Correl.proton.mult ==1 && Correl.Al25.mult ==1)
  {
 
    float const Q26Si = -5.514; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Al25.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_26Si = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_26Si + Q26Si;
    Histo_read->Erel_26Si_1p25Al->Fill(Erel_26Si);
    Histo_read->Ex_26Si_1p25Al->Fill(Ex);
    Histo_read->ThetaCM_26Si_1p25Al->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_26Si_1p25Al->Fill(vCM);

  }

}


void det::corr_27Si()
{
  // p + 26Al
  if(Correl.proton.mult ==1 && Correl.Al26.mult ==1)
  {
 
    float const Q27Si = -7.463; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Al26.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_27Si = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_27Si + Q27Si;
    Histo_read->Erel_27Si_1p26Al->Fill(Erel_27Si);
    Histo_read->Ex_27Si_1p26Al->Fill(Ex);
    Histo_read->ThetaCM_27Si_1p26Al->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_27Si_1p26Al->Fill(vCM);

  }

}

void det::corr_25Si()
{
  // p + 24Al
  if(Correl.proton.mult ==1 && Correl.Al24.mult ==1)
  {
 
    float const Q25Si = -3.414; //Taken from Ame2016 
    Correl.zeroMask();
    Correl.proton.mask[0]=1;
    //      Correl.proton.mask[1]=1;
    Correl.Al24.mask[0]=1;   
    Correl.makeArray(1);
    float Erel_25Si = Correl.findErel();
    float thetaCM = Correl.thetaCM;
    float Ex = Erel_25Si + Q25Si;
    Histo_read->Erel_25Si_1p24Al->Fill(Erel_25Si);
    Histo_read->Ex_25Si_1p24Al->Fill(Ex);
    Histo_read->ThetaCM_25Si_1p24Al->Fill(thetaCM*180./acos(-1));
    float vCM = Correl.velocityCM;
    Histo_read->VCM_25Si_1p24Al->Fill(vCM);

  }

}
*/
