#include "pixels.h"
//#include "loss.h"
#include "TFile.h"
#include "TH1I.h"
#include "TTree.h"
#include "kin.h"
#include "calibrate.h"
#include "TH2I.h"
#include "constants.h"
//#include "calCorrection.h"
//#include "fineCal.h"

//from p+21Na 
double suppress(double E)
{
   const double c0 = .71969;
   const double c1 = 1.067;
   return ((1.+exp(c0/c1))/(1.+exp(-(E-c0)/c1))-1.)*exp(-c0/c1);
}


int main()
{
  const int nn = 200;
  save savep[nn];
  save save19[nn];
  int nwrite = 0;
  int nsave = 0.;

  bool B17[nn];
  bool B15[nn];


  kin Kin;

  TFile file("pp16O_save.root");
  TTree * tree = (TTree*)file.Get("T");
  
  int N = tree->GetEntries();
  cout << "Number of ENtries = " << N << endl;


  vel ee1;
  vel ee2;
  vel ee3;
  vel CM;

  float MM1;
  float MM2;
  float M1[3];
  float M2[3];

  int id2;
  int id1;
  int ifront2;
  int iback2;
  int ifront1;
  int iback1;

  int itele1;

  float et1;
  float et2;
  float denergy1;
  float energy_p1;


  float denergy2;
  float energy_p2;

  float timeCsI1;
  float timeCsI2;

  int Ngamma;
  float Egamma[10];
  float Tgamma[10];
  int Chgamma[10];

  float Erel;
  float Ex;
  float Vcm;
  float thetaCM;
  float cos_thetaH;

  int beamZ;  



  tree->SetBranchAddress("M1",M1);
  tree->SetBranchAddress("M2",M2);
  tree->SetBranchAddress("et1",&et1);
  tree->SetBranchAddress("et2",&et2);
  tree->SetBranchAddress("id2",&id2);
  tree->SetBranchAddress("id1",&id1);
  tree->SetBranchAddress("id1",&id1);
  tree->SetBranchAddress("itele1",&itele1);
  tree->SetBranchAddress("ifront2",&ifront2);
  tree->SetBranchAddress("iback2",&iback2);
  tree->SetBranchAddress("ifront1",&ifront1);
  tree->SetBranchAddress("iback1",&iback1);


  //tree->SetBranchAddress("denergy1",&denergy1);
  tree->SetBranchAddress("energy_p1",&energy_p1);


  //tree->SetBranchAddress("denergy2",&denergy2);
  tree->SetBranchAddress("energy_p2",&energy_p2);


  tree->SetBranchAddress("Ngamma",&Ngamma);
  tree->SetBranchAddress("Egamma",Egamma);
  tree->SetBranchAddress("Tgamma",Tgamma);
  tree->SetBranchAddress("Chgamma",Chgamma);

  tree->SetBranchAddress("Erel",&Erel);
  tree->SetBranchAddress("Ex",&Ex);
  tree->SetBranchAddress("Vcm",&Vcm);
  tree->SetBranchAddress("thetaCM",&thetaCM);
  tree->SetBranchAddress("cos_thetaH",&cos_thetaH);

  tree->SetBranchAddress("beamZ",&beamZ);


  tree->SetBranchAddress("time1",&timeCsI1);
  //tree->SetBranchAddress("timeCsI2",&timeCsI2);


  //name ="Carbon.loss";
  //CLoss loss_Core(name,Acore);
  //name="Hydrogen.loss";
  //CLoss loss_proton(name,1.);

  TFile fileOut("sort.root","RECREATE");

  TH1I * Erel_hist = new TH1I("Erel_hist","",250,0,10);
  TH1I * Erel_trans = new TH1I("Erel_trans","",250,0,10);
  TH1I * Erel_long = new TH1I("Erel_long","",250,0,10);
  TH1I * Erel_Sibeam = new TH1I("Erel_Sibeam","",250,0,10);
  TH1I * Erel_Albeam = new TH1I("Erel_Albeam","",250,0,10);
  TH1I * Erel_Mgbeam = new TH1I("Erel_Mgbeam","",250,0,10);
  TH1I * Ex_hist = new TH1I("Ex_hist","",250,0,10);
  TH1I * Ex_trans = new TH1I("Ex_trans","",250,0,10);
  TH1I * Ex_trans_0_4 = new TH1I("Ex_trans_0_4","",250,0,10);
  TH1I * Ex_long = new TH1I("Ex_long","",250,0,10);
  TH1I * Ex_narrowtime = new TH1I("Ex_narrowtime","",250,0,10);
  TH1I * Ex_0_5costhetaH = new TH1I("Ex_0_5costhetaH","",250,0,10);
  //TH1I * Erel_19Na_p18Ne_set1 = new TH1I("Erel_19Na_p18Ne_set1","",500,-5,15);
  //TH1I * Erel_19Na_p18Ne_set1a = new TH1I("Erel_19Na_p18Ne_set1a","",500,-5,15);
  TH1I * ThetaCM = new TH1I("ThetaCM","",200,0,10);
  TH1I * VCM = new TH1I("VCM","",400,7,15);
  TH2I * Erel_costhetaH = new TH2I("Erel_costhetaH","",800,0,8,50,-1,1);
  Erel_costhetaH->SetMinimum(1.0);
  TH2I * VCMvsErel = new TH2I("VCMvsErel","",400,7,15,250,0,10);
  VCMvsErel->SetMinimum(1.0);
  TH1I * gammasADD = new TH1I("gammasADD","",2048, 0, 8192);
  //TH1I * p18Ne_gammasADD_tgate = new TH1I("p18Ne_gammasADD_tgate","",2048, 0, 8192);
  TH2I * gammasADDvsErel = new TH2I("gammasADDvsErel","",250,0,10,2048/8, 0, 2048);
  gammasADDvsErel->SetMinimum(1.0);
  TH2I * gammasADDvsEx = new TH2I("gammasADDvsEx","",250,0,10,2048/8, 0, 2048);
  gammasADDvsEx->SetMinimum(1.0);
  TH2I * gammasADDvsgammasADD = new TH2I("gammasADDvsgammasADD","",2048/4, 0, 8192,2048/4, 0, 8192);
  gammasADDvsgammasADD->SetMinimum(1.0);

  //Gamma-Gamma gates
  TH2I * gammasADDvsErel_greater2 = new TH2I("gammasADDvsErel_greater2","",250,0,10,2048, 0, 8192);
  gammasADDvsErel_greater2->SetMinimum(1.0);
  TH2I * gammasADDvsErel_with200 = new TH2I("gammasADDvsErel_with200","",250,0,10,2048, 0, 8192);
  gammasADDvsErel_with200->SetMinimum(1.0);
  TH2I * gammasADDvsErel_with600 = new TH2I("gammasADDvsErel_with600","",250,0,10,2048, 0, 8192);
  gammasADDvsErel_with600->SetMinimum(1.0);

  TH2I * gammasADDvsgammasADD_with200 = new TH2I("gammasADDvsgammasADD_with200","",2048/4, 0, 8192,2048/4, 0, 8192);
  gammasADDvsgammasADD_with200->SetMinimum(1.0);
  TH2I * gammasADDvsgammasADD_with600 = new TH2I("gammasADDvsgammasADD_with600","",2048/4, 0, 8192,2048/4, 0, 8192);
  gammasADDvsgammasADD_with600->SetMinimum(1.0);

  //Erel gated on gammas
  TH1I * Erel_with500 = new TH1I("Erel_with500","",250,0,10);
  TH1I * Erel_no500 = new TH1I("Erel_no500","",250,0,10);

  //Ex gated on gammas
  TH1I * Ex_with500 = new TH1I("Ex_with500","",250,0,10);
  TH1I * Ex_no500 = new TH1I("Ex_no500","",250,0,10);

  TH1I * Ex_nogamma= new TH1I("Ex_nogamma","",250,0,10);
  TH1I * Ex_nogamma_0_5costhetaH = new TH1I("Ex_nogamma_0_5costhetaH","",250,0,10);

  //2D gamma t vs channel
  TH2I * gammaTvsChannel = new TH2I("gammaTvsChannel","",192,0,192,5000,-5000,5000);
  TH2I * gammaTvsE = new TH2I("gammaTvsE","",2048/4, 0, 8192,5000,-5000,5000);

  //2D csi time vs channel
  TH2I * CsITvsChannel = new TH2I("CsITvsChannel","",16,-0.5,15.5,2000,-2000,2000);
  //vs Ex
  TH2I * CsITvsEx = new TH2I("CsITvsEx","",250,0,10,2000,-2000,2000);

  file.cd();
  for (int i=0;i<N;i++)
  {
    if (i%3000 == 0)cout << i << endl;
    
    tree->GetEntry(i); //Pull data from tree
    cout << timeCsI1 << endl;

    //if (timeCsI1 < 960 || timeCsI1 > 1030) continue;

    Erel_hist->Fill(Erel);
    Ex_hist->Fill(Ex);
    ThetaCM->Fill(thetaCM*180./pi);
    VCM->Fill(Vcm);
    VCMvsErel->Fill(Vcm,Erel);
    Erel_costhetaH->Fill(Erel,cos_thetaH);

    //Fill time spectra 
    CsITvsChannel->Fill(itele1*4 + id1,timeCsI1);
    CsITvsEx->Fill(Ex,timeCsI1);

    //CsI time gate between 970 and 1030
    if (timeCsI1 >= 970 && timeCsI1 <= 1030) Ex_narrowtime->Fill(Ex);

    //ThetaH gates
    if (fabs(cos_thetaH) <= 0.1)
    {
      Erel_trans->Fill(Erel);
      Ex_trans->Fill(Ex);
    }
    if (fabs(cos_thetaH) >= 0.9)
    {
      Erel_long->Fill(Erel);
      Ex_long->Fill(Ex);
    }

    //ThetaH gates
    if (fabs(cos_thetaH) <= 0.4)
    {
      Ex_trans_0_4->Fill(Ex);
    }

    //Gamma coincidence gates
    int with500 = 0;
    int no500 = 1;

    for (int j=0;j<Ngamma;j++) {
      if (Tgamma[j] < 675 || Tgamma[j] > 715) continue;
      if (Egamma[j] >= 0.443 && Egamma[j] <= 0.522) with500 = 1;
      if (Egamma[j] >= 0.443 && Egamma[j] <= 0.522) no500 = 0;
    } 

    //Erel gated on gammas
    if (with500 == 1) Erel_with500->Fill(Erel);
    if (no500 == 1) Erel_no500->Fill(Erel);

    //Ex gated on gammas
    if (with500 == 1) Ex_with500->Fill(Ex);
    if (with500 == 0) Ex_no500->Fill(Ex);

    if (fabs(cos_thetaH) <= 0.5) Ex_0_5costhetaH->Fill(Ex);

    if (Ngamma == 0) Ex_nogamma->Fill(Ex);
    if (Ngamma == 0 && fabs(cos_thetaH) <= 0.5) Ex_nogamma_0_5costhetaH->Fill(Ex);

    //Gammas - already doppler shifted
    for (int j=0;j<Ngamma;j++)
    {
      if (Ngamma >= 30) {
        cout << "Ngamma >= 30, skip " << endl;
        break;
      }

      //Gamma time gate
      //if (Tgamma[j] < 675 || Tgamma[j] > 715) continue;

      gammasADD->Fill(Egamma[j]*1000.);
      gammasADDvsErel->Fill(Erel,Egamma[j]*1000.);
      gammasADDvsEx->Fill(Ex,Egamma[j]*1000.);

      gammaTvsChannel->Fill(Chgamma[j],Tgamma[j]);
      gammaTvsE->Fill(Egamma[j]*1000.,Tgamma[j]);
      
      cout << Tgamma[j] << endl;

      if (Ngamma > 2) gammasADDvsErel_greater2->Fill(Erel,Egamma[j]*1000.);

      for (int k=j+1;k<Ngamma;k++)
      {
        gammasADDvsgammasADD->Fill(Egamma[j]*1000.,Egamma[k]*1000.);

      }
    }




/*
//et1 = FineCal.cal_p(energy_p1,denergy1,ifront1,iback1,id1,193.,M1);
//double theta1 = FineCal.theta;
//double phi1 =FineCal.phi;
//et2 = FineCal.cal_O14(energy_p2,denergy2,ifront2,iback2,id2,193.,M2);
//double theta2 = FineCal.theta;
//double phi2 =FineCal.phi;

  double theta1;
  double phi1;
  
  double theta2;
  double phi2;

/*
      double theta1 = Pix.getAngle(id1/4,ifront1,iback1);
      double phi1 = Pix.phi;
      double E = energy_p1;
      E += denergy1;
        double thick = 193./2./cos(theta1);
        E = loss_p.getEin(E,thick);

        et1 = E + Mass_p;
        double  pc = sqrt(pow(et1,2) - pow(Mass_p,2));

        M1[0] = pc*sin(theta1)*cos(phi1);
        M1[1] = pc*sin(theta1)*sin(phi1);
        M1[2] = pc*cos(theta1);



      double theta2 = Pix.getAngle(id2/4,ifront2,iback2);
      double phi2 = Pix.phi;
       E = energy_p2;
       E = cal.getEnergy(0,id2,E);
       E += 2.8;
       E += denergy2;
        thick = 193./2./cos(theta2);
        E = loss_O.getEin(E,thick);

        et2 = E + Mass_14O;
        pc = sqrt(pow(et2,2) - pow(Mass_14O,2));

        M2[0] = pc*sin(theta2)*cos(phi2);
        M2[1] = pc*sin(theta2)*sin(phi2);
        M2[2] = pc*cos(theta2);


      //hist_Ek_14O->Fill(et2-Mass_14O);

      
      //CM = Kin.findCM(M1,M2,et1,et2);
      //hist_Vcm->Fill(CM.vv);
      //ee1 = Kin.trans(CM,M1,et1);
      //ee2 = Kin.trans(CM,M2,et2);
      //float Erel = ee1.vv + ee2.vv - Mass_14O - Mass_p;





      bool O15beam = false;
      bool Ne17beam = false;
      if (CM.vv > 9.39)
      {
      Ne17beam = true;
      }
     if (CM.vv < 9.27)
      {
       O15beam = true;
      }



     bool Bgamma = false;
     for (int j=0;j<Ngamma;j++)
       {
	 if (Egamma[j] > 0.01)
	   {
            hist_Egamma->Fill(Egamma[j]);
            map_Erel_Egamma->Fill(Erel,Egamma[j]);
            if(O15beam)hist_Egamma_15O->Fill(Egamma[j]);
            if(Ne17beam)hist_Egamma_17Ne->Fill(Egamma[j]);
            if (Egamma[j] > 4.16 && Egamma[j] < 5.6) Bgamma = true;
	   }
       }
*/


//if (CM.vv < 9.36) continue;    //sort only 17Ne beam contribution
/*
      hist_Vcm->Fill(CM.vv);
      hist_Erel->Fill(Erel);


/if (Erel > 2.5 && Erel < 3. && CM.vv > 9.36)
{
 hist_Vcm_peak->Fill(CM.vv);
 hist_thetaCM_peak->Fill(acos(CM.v[2]/CM.vv)*180./acos(-1.));
}
      double cosTheta = ee2.v[2]/sqrt(pow(ee2.v[0],2)+pow(ee2.v[1],2)+pow(ee2.v[2],2));

      map_Erel_cos->Fill(Erel,cosTheta);
      map_Vcm_cos->Fill(CM.vv,cosTheta);
      if (fabs(cosTheta) < .5) hist_Erel_trans->Fill(Erel);
      if (fabs(cosTheta) < .2) hist_Erel_trans_narrow->Fill(Erel);
      if (cosTheta > .7) hist_Erel_for->Fill(Erel);
      if (cosTheta < -.4 && cosTheta > -.7) hist_Erel_mid_back->Fill(Erel);
      if (cosTheta < -.7) hist_Erel_back->Fill(Erel);

      bool BB17 = false;
      bool BB15 = false;
      if (Bgamma) hist_Erel_gamma->Fill(Erel);



      if (CM.vv > 9.39)
      {

      hist_Erel_17Ne->Fill(Erel);
      if (Bgamma) hist_Erel_17Ne_gamma->Fill(Erel);
      if (fabs(cosTheta) < .5) hist_Erel_trans_17Ne->Fill(Erel);
      if (fabs(cosTheta) < .2) hist_Erel_trans_narrow_17Ne->Fill(Erel);
      BB17 = true;
      map_Erel_cos_17Ne->Fill(Erel,cosTheta);
      }
     else if (CM.vv < 9.27)
      {

       hist_Erel_15O->Fill(Erel);
      if (Bgamma) hist_Erel_15O_gamma->Fill(Erel);
      if (fabs(cosTheta) < .5) hist_Erel_trans_15O->Fill(Erel);
      if (fabs(cosTheta) < .2) hist_Erel_trans_narrow_15O->Fill(Erel);
      BB15 = true;
      map_Erel_cos_15O->Fill(Erel,cosTheta);
      }




   map_proton->Fill(theta1*cos(phi1),theta1*sin(phi1));    
   map_core->Fill(theta2*cos(phi2),theta2*sin(phi2));    
      save saveP(id1,et1,M1);
      save saveO(id2,et2,M2);
*/
      /*for (int j=0;j<nsave;j++)
        {

           if(savep[j].id != saveP.id && savep[j].id != saveO.id && BB17 && B17[j])
             {
              CM = Kin.findCM(savep[j],saveP,saveO);
              ee1 = Kin.trans(CM,savep[j]);
              ee2 = Kin.trans(CM,saveO);
              ee3 = Kin.trans(CM,saveP);
              double Erel = ee1.vv + ee2.vv + ee3.vv - 2.*Mass_p - Mass_14O;

              CM = Kin.findCM(savep[j],saveO);
              ee1 = Kin.trans(CM,savep[j]);
              ee2 = Kin.trans(CM,saveO);
              double Erel_p_O = ee1.vv + ee2.vv;
              double supp = suppress(Erel_p_O);
               hist_Erel_16Ne_mix_weight->Fill(Erel,supp);
               double cosTheta = ee3.v[2]/sqrt(pow(ee3.v[0],2)+pow(ee3.v[1],2)+pow(ee3.v[2],2));
               if(fabs(cosTheta) < .2) hist_Erel_16Ne_mix_weight_trans_narrow->Fill(Erel,supp);
             }

           if(saveP.id != savep[j].id && saveP.id != save18[j].id && BB17 && B17[j])
             {
              CM = Kin.findCM(savep[j],saveP,save18[j]);
              ee1 = Kin.trans(CM,savep[j]);
              ee2 = Kin.trans(CM,save18[j]);
              ee3 = Kin.trans(CM,saveP);
              double Erel = ee1.vv + ee2.vv + ee3.vv - 2.*Mass_p - Mass_14O;

              CM = Kin.findCM(saveP,save18[j]);
              ee1 = Kin.trans(CM,savep[j]);
              ee2 = Kin.trans(CM,saveO);
              double Erel_p_O = ee1.vv + ee2.vv;
              double supp = suppress(Erel_p_O);
              hist_Erel_16Ne_mix_weight->Fill(Erel,supp);
               double cosTheta = ee3.v[2]/sqrt(pow(ee3.v[0],2)+pow(ee3.v[1],2)+pow(ee3.v[2],2));
               if(fabs(cosTheta) < .2) hist_Erel_16Ne_mix_weight_trans_narrow->Fill(Erel,supp);

             }


           if (savep[j].id != saveO.id)
             {
              CM = Kin.findCM(savep[j],saveO);
              ee1 = Kin.trans(CM,savep[j]);
              ee2 = Kin.trans(CM,saveO);
              double Erel = ee1.vv + ee2.vv - Mass_p - Mass_14O;
              double supp = suppress(Erel);
              //cout << Erel << " " << supp << endl;
              hist_Erel_mix->Fill(Erel,supp);
               if (BB17 && B17[j])hist_Erel_17Ne_mix->Fill(Erel,supp);
               if (BB15 && B15[j])hist_Erel_15O_mix->Fill(Erel,supp);
              double cosTheta = ee2.v[2]/sqrt(pow(ee2.v[0],2)+pow(ee2.v[1],2)+pow(ee2.v[2],2));
              if (fabs(cosTheta) < .2)
                {

                 hist_Erel_trans_narrow_mix->Fill(Erel,supp);
                 if (BB17 && B17[j])hist_Erel_trans_narrow_17Ne_mix->Fill(Erel,supp);
                 if (BB15 && B15[j])hist_Erel_trans_narrow_15O_mix->Fill(Erel,supp);
                }
             }

           if (saveP.id != save18[j].id)
             {
              CM = Kin.findCM(saveP,save18[j]);
              ee1 = Kin.trans(CM,saveP);
              ee2 = Kin.trans(CM,save18[j]);
              double Erel = ee1.vv + ee2.vv - Mass_p - Mass_14O;
              double supp = suppress(Erel);
              //cout << Erel << " " << supp << endl;
              hist_Erel_mix->Fill(Erel,supp);
                 if (BB17 && B17[j])hist_Erel_17Ne_mix->Fill(Erel,supp);
                 if (BB15 && B15[j])hist_Erel_15O_mix->Fill(Erel,supp);


              double cosTheta = ee2.v[2]/sqrt(pow(ee2.v[0],2)+pow(ee2.v[1],2)+pow(ee2.v[2],2));
              if (fabs(cosTheta) < .2)
                {

                hist_Erel_trans_narrow_mix->Fill(Erel,supp);
                 if (BB17 && B17[j])hist_Erel_trans_narrow_17Ne_mix->Fill(Erel,supp);
                 if (BB15 && B15[j])hist_Erel_trans_narrow_15O_mix->Fill(Erel,supp);
             }
           }
        }*/

      //B17[nwrite] = BB17;
      //B15[nwrite] = BB15;
      //savep[nwrite] = saveP;
      //save18[nwrite] = saveO;
      //nwrite++;
      //if(nwrite == nn) nwrite = 0;
      //nsave++;
      //if(nsave == nn) nsave = nn -1;
    }

  fileOut.Write();
}
