#include "correl.h"
#include <cmath>

using namespace std;
void correl::load(solution* particle)
{

  if (particle->iZ == 1 && particle->iA == 1)
    {
      proton[multProton] = particle;
      if (multProton < 7) multProton++;
      multPro[particle->itele]++;
      //proton[multProton]->A = 1;

    }
  else if (particle->iZ == 1 && particle->iA == 2)
    {
      H2[multH2] = particle;
      if (multH2 < 5) multH2++;
      itele = particle->itele;
      //H2[multH2]->A = 2;
      multh2[itele]++;
    }
  else if (particle->iZ == 1 && particle->iA == 3)
    {

      H3[multH3] = particle;
      if (multH3 < 5) multH3++;
      itele = particle->itele;
      //H3[multH3]->A = 3;
      multh3[itele]++;
    }
 else if (particle->iZ == 2 && particle->iA == 3)
    {
      He3[mult3He] = particle;
      if (mult3He < 5) mult3He++;
      itele = particle->itele;
      //He3[mult3He]->A = 3;
      mult3he[itele]++;
    }
  else if (particle->iZ == 2 && particle->iA ==4)
    {
      alpha[multAlpha] = particle;
      if (multAlpha < 5) multAlpha++;
      itele = particle->itele;
      //alpha[multAlpha]->A = 4;
      multAlp[itele]++;
    }
  else if (particle->iZ == 3 && particle->iA == 6)
    {
      Li6[mult6Li] = particle;
      if (mult6Li < 5) mult6Li++;
      itele = particle->itele;
     
      mult6li[itele]++;
    }
  else if (particle->iZ == 4 && particle->iA == 7)
    {
      Be7[mult7Be] = particle;
      if (mult7Be < 5) mult7Be++;
      itele = particle->itele;

      mult7be[itele]++;
    }
  else if (particle->iZ == 5 && particle->iA == 8)
    {
      B8[mult8B] = particle;
      if (mult8B < 5) mult8B++;
      itele = particle->itele;
     
      mult8b[itele]++;
    }
  else if (particle->iZ == 5 && particle->iA == 10)
    {
      B10[mult10B] = particle;
      if (mult10B < 5) mult10B++;
      itele = particle->itele;
     
      mult10b[itele]++;
    }
  else if (particle->iZ == 6 && particle->iA == 9)
    {
      C9[mult9C] = particle;
      if (mult9C < 5) mult9C++;
      itele = particle->itele;
     
      mult9c[itele]++;
    }
    else if (particle->iZ == 99 && particle->iA == 99)
    {
      C9Fake[mult9CFake] = particle;
      if (mult9CFake < 5) mult9CFake++;
      itele = particle->itele;
     
      mult9cFake[itele]++;
    }
    else if (particle->iZ == 89 && particle->iA == 99)
    {
      C9Fake_Channeling[mult9CFake_Channeling] = particle;
      if (mult9CFake_Channeling < 5) mult9CFake_Channeling++;
      itele = particle->itele;
     
      mult9cFake_Channeling[itele]++;
    }
    else if (particle->iZ == 79 && particle->iA == 99)
    {
      C9Fake_CsIreaction[mult9CFake_CsIreaction] = particle;
      if (mult9CFake_CsIreaction < 5) mult9CFake_CsIreaction++;
      itele = particle->itele;
     
      mult9cFake_CsIreaction[itele]++;
    }
    else if (particle->iZ == 69 && particle->iA == 99)
    {
      C9Fake_CsIreaction_C11[mult9CFake_CsIreaction_C11] = particle;
      if (mult9CFake_CsIreaction_C11 < 5) mult9CFake_CsIreaction_C11++;
      itele = particle->itele;
      
      mult9cFake_CsIreaction_C11[itele]++;
    }
    else if (particle->iZ == 59 && particle->iA == 99)
    {
      C10Fake_CsIreaction[mult10CFake_CsIreaction] = particle;
      if (mult10CFake_CsIreaction < 5) mult10CFake_CsIreaction++;
      itele = particle->itele;
      
      mult10cFake_CsIreaction[itele]++;
    }
  else if (particle->iZ == 6 && particle->iA == 10)
    {
      C10[mult10C] = particle;
      if (mult10C < 5) mult10C++;
      itele = particle->itele;
     
      mult10c[itele]++;
    }
    else if (particle->iZ == 6 && particle->iA == 11)
    {
      C11[mult11C] = particle;
      if (mult11C < 5) mult11C++;
      itele = particle->itele;
     
      mult11c[itele]++;
    }
    else if (particle->iZ == 6 && particle->iA == 12)
    {
      C12[mult12C] = particle;
      if (mult12C < 5) mult12C++;
      itele = particle->itele;
     
      mult12c[itele]++;
    }
  else if (particle->iZ == 8 && particle->iA == 13)
    {
      O13[mult13O] = particle;
      if (mult13O < 5) mult13O++;
      itele = particle->itele;
     
      mult13o[itele]++;
    }
}
//**********************************************************
void correl::setMask()
{
  for (int i=0;i<5;i++)
    {
      maskProton[i] = 1;
      maskAlpha[i] = 1;
      mask2H[i] = 1;
      mask3He[i] = 1;
      mask6Li[i] = 1;
      mask3H[i] = 1;
    }
}
//************************************************************
void correl::zeroMask()
{
  for (int i=0;i<5;i++)
    {
      maskProton[i] = 0;
      maskAlpha[i] = 0;
      mask3He[i] = 0;
      mask2H[i] = 0;
      mask3H[i] = 0;
      mask6Li[i] = 0;
      mask7Be[i] = 0;
      mask10B[i] = 0;
      mask8B[i] = 0;
      mask9C[i] = 0;
      mask9CFake[i] = 0;
      mask9CFake_Channeling[i] = 0;
      mask9CFake_CsIreaction[i] = 0;
      mask9CFake_CsIreaction_C11[i] = 0;
      mask10CFake_CsIreaction[i] = 0;
      mask10C[i] = 0;
      mask11C[i] = 0;
      mask13O[i] = 0;
    }
}
//**********************************************************
void correl::zeroH2()
{
  for (int i=0;i<5;i++) mask2H[i] = 0;
}
//**********************************************************
void correl::zero3He()
{
  for (int i=0;i<5;i++) mask3He[i] = 0;
}
//*********************************************************
void correl::zeroProton()
{
  for (int i=0;i<5;i++) maskProton[i] = 0;
}
//************************************************
void correl::zeroAlpha()
{
  for (int i=0;i<5;i++) maskAlpha[i] = 0;
}
//*******************************************************
  /**
   * make an array of fragments for later analyisis
   */
void correl::makeArray(bool flagMask/*=0*/)
{

 N = 0;

 for (int j=0;j<multProton;j++) 
    {
    if (!flagMask || maskProton[j])
      {
	frag[N] = proton[j];
        N++;
      }
    }
  for (int j=0;j<multAlpha;j++)
    {
    if (!flagMask || maskAlpha[j]) 
      {
	frag[N] = alpha[j];
        N++;
      }
    }
  for (int j=0;j<multH2;j++)
    {
    if (!flagMask || mask2H[j]) 
      {
	frag[N] = H2[j];
	N++;
      }
    }
 for (int j=0;j<multH3;j++)
    {
    if (!flagMask || mask3H[j]) 
      {
	frag[N] = H3[j];
	N++;
      }
    }
 for (int j=0;j<mult3He;j++)
    {
    if (!flagMask || mask3He[j]) 
      {
	frag[N] = He3[j];
	N++;
      }
    } 
 for (int j=0;j<mult6Li;j++)
    {
    if (!flagMask || mask6Li[j]) 
      {
	frag[N] = Li6[j];
	N++;
      }
    } 
 for (int j=0;j<mult7Be;j++)
    {
    if (!flagMask || mask7Be[j]) 
      {
	frag[N] = Be7[j];
	N++;
      }
    } 
  for (int j=0;j<mult8B;j++)
    {
    if (!flagMask || mask8B[j]) 
      {
	frag[N] = B8[j];
	N++;
      }
    } 
  for (int j=0;j<mult10B;j++)
    {
    if (!flagMask || mask10B[j]) 
      {
	frag[N] = B10[j];
	N++;
      }
    } 
 for (int j=0;j<mult9C;j++)
    {
    if (!flagMask || mask9C[j]) 
      {
	frag[N] = C9[j];
	N++;
      }
    }
  for (int j=0;j<mult9CFake;j++)
    {
    if (!flagMask || mask9CFake[j]) 
      {
	frag[N] = C9Fake[j];
	N++;
      }
    } 
  for (int j=0;j<mult9CFake_Channeling;j++)
    {
    if (!flagMask || mask9CFake_Channeling[j]) 
      {
	frag[N] = C9Fake_Channeling[j];
	N++;
      }
    } 
  for (int j=0;j<mult9CFake_CsIreaction;j++)
    {
    if (!flagMask || mask9CFake_CsIreaction[j]) 
      {
	frag[N] = C9Fake_CsIreaction[j];
	N++;
      }
    } 
  for (int j=0;j<mult9CFake_CsIreaction_C11;j++)
    {
    if (!flagMask || mask9CFake_CsIreaction_C11[j]) 
      {
	frag[N] = C9Fake_CsIreaction_C11[j];
	N++;
      }
    } 
  for (int j=0;j<mult10CFake_CsIreaction;j++)
    {
    if (!flagMask || mask10CFake_CsIreaction[j]) 
      {
	frag[N] = C10Fake_CsIreaction[j];
	N++;
      }
    } 
 for (int j=0;j<mult10C;j++)
    {
    if (!flagMask || mask10C[j]) 
      {
	frag[N] = C10[j];
	N++;
      }
    }
  for (int j=0;j<mult11C;j++)
    {
    if (!flagMask || mask11C[j]) 
      {
	frag[N] = C11[j];
	N++;
      }
    }  

  for (int j=0;j<mult12C;j++)
    {
    if (!flagMask || mask12C[j]) 
      {
	frag[N] = C12[j];
	N++;
      }
    }  
 for (int j=0;j<mult13O;j++)
    {
    if (!flagMask || mask13O[j]) 
      {
	frag[N] = O13[j];
	N++;
      }
    } 
}

//***************************************
void correl::Reset()
{
  multProton = 0;
  multAlpha = 0;
  multH2 = 0;
  mult3He = 0;
  mult6Li = 0;
  mult7Be = 0;
  multH3 = 0;
  mult9C =0;
  mult9CFake =0;
  mult9CFake_Channeling =0;
  mult9CFake_CsIreaction =0;
  mult9CFake_CsIreaction_C11 =0;
  mult10CFake_CsIreaction =0;
  mult8B = 0;
  mult10B = 0;
  mult10C=0;
  mult11C=0;
  mult12C=0;
  mult13O =0;
  multGamma = 0;
  for (int itele=0;itele<14;itele++) 
    {
      multAlp[itele] = 0;
      multPro[itele] = 0;
      multh2[itele] = 0;
      mult3he[itele] = 0;
      mult6li[itele] = 0;
      mult7be[itele] = 0;
      multh3[itele] = 0;
      mult8b[itele] =0;
      mult10b[itele] =0;
      mult9c[itele]=0;
      mult12c[itele]=0;
      mult9cFake[itele]=0;
      mult10c[itele]=0;
      mult13o[itele]=0;
      mult_Gamma[itele] = 0;
    }
}

//**************************************************************
  /**
   * Finds the total kinetic energy of the fragments
   * in there center-of-mass frame.
   */
float correl::findErel()
{
  
  //first find total momentum
  for (int i=0;i<3;i++) Mtot[i] = 0.;
  float energyTot = 0.;   // total energy for relativity, total mass for newton



  for (int i=0;i<N;i++) 
    {
      energyTot += frag[i]->energyTot;
      for (int j=0;j<3;j++)
        Mtot[j] += frag[i]->Mvect[j];
    }

  momentumCM = 0.;
  for (int j=0;j<3;j++) momentumCM += pow(Mtot[j],2);
  momentumCM = sqrt(momentumCM);
  
  //transform to average
  
  //float velC[3]={0.,0.,10.8699};// 7Be beam
  //float velC[3]={0.,0.,10.7685};// 9C
  //float velC[3]={0.,0.,10.8381};// 9C
  float velC[3]={0.,0.,10.9687}; //13O beam

  Kinematics.transformMomentum(Mtot,velC,
        energyTot,momC);
  float mmc = 0.;
  for (int i=0;i<3;i++) mmc += pow(momC[i],2);
  mmc = sqrt(mmc);
  cosThetaC = momC[2]/mmc;
  PperpC = sqrt(pow(momC[0],2)+pow(momC[1],2));
  PparaC = momC[2];
  PtotC = sqrt(pow(momC[0],2)+pow(momC[1],2)+pow(momC[2],2));

  //velocity of cemter of mass
  velocityCM = momentumCM*Kinematics.c/energyTot;
  

  float velCM[3]={0.};
  for (int j=0;j<3;j++) velCM[j] = velocityCM/momentumCM*Mtot[j];
  thetaCM = acos(velCM[2]/velocityCM);
  phiCM = atan2(velCM[1],velCM[0]);

  float totalKE = 0.;
  for (int i=0;i<N;i++)
    {
      float eKinNew = Kinematics.transformMomentum(frag[i]->Mvect,velCM,frag[i]->energyTot,frag[i]->MomCM);
      frag[i]->energyCM = eKinNew - Kinematics.scale*frag[i]->mass;
      totalKE += eKinNew - Kinematics.scale*frag[i]->mass;
      check_ke = eKinNew - Kinematics.scale*frag[i]->mass;
      check_mass = Kinematics.scale*frag[i]->mass;
    }


  if (N == 3)
    {
      float dot = 0.;
      float mm = 0.;
      for (int j=0;j<3;j++)
	{
         dot += frag[2]->MomCM[j]*momC[j];
         mm += pow(frag[2]->MomCM[j],2);
	}
      mm = sqrt(mm);
      cosAlphaQ = dot/mm/PtotC;
    }


  return totalKE;
}
//***********************************************************
void correl::getJacobi()
{

  for (int i=0;i<3;i++)
    {
      frag[i]->momentumCM = 0.;
      for (int k=0;k<3;k++) frag[i]->momentumCM +=
	pow(frag[i]->MomCM[k],2);
      frag[i]->momentumCM = sqrt(frag[i]->momentumCM);
    }


  //alpha is the  third fragment
  //first JacobiT
  float dot = 0.;
  float pp[3] = {0.};
  float PP = 0.;
  for (int k=0;k<3;k++)
    {
      pp[k] = frag[0]->MomCM[k] - frag[1]->MomCM[k];
      PP += pow(pp[k],2);
      dot += pp[k]*frag[2]->MomCM[k];
    }
  PP = sqrt(PP);
  cosThetaT = dot/PP/frag[2]->momentumCM;


  dot = 0;
  double PP1 = 0;
  double pp1[3]={0.};
  for (int k=0;k<3;k++)
    {
      pp1[k] = frag[0]->MomCM[k]/frag[0]->mass 
            - frag[2]->MomCM[k]/frag[2]->mass;
      PP1 += pow(pp1[k],2);
      dot += pp1[k]*frag[1]->MomCM[k];
    }  
  PP1 = sqrt(PP1);
  cosThetaY[0] = -dot/PP1/frag[1]->momentumCM;


  dot = 0;
  double PP2 = 0;
  double pp2[3]={0.};
  for (int k=0;k<3;k++)
    {
      pp2[k] = frag[1]->MomCM[k]/frag[1]->mass 
            - frag[2]->MomCM[k]/frag[2]->mass;
      PP2 += pow(pp2[k],2);
      dot += pp2[k]*frag[0]->MomCM[k];
    }  
  PP2 = sqrt(PP2);
  cosThetaY[1] = -dot/PP2/frag[0]->momentumCM;


  cosThetaV = (pp1[0]*pp2[0] + pp1[1]*pp2[1] + pp1[2]*pp2[2])/PP1/PP2;
}
//*****************************************************
void correl::printP()
{

  for (int i=0;i<5;i++)
    {
      frag[i]->momentumCM = 0.;
      for (int k=0;k<3;k++) frag[i]->momentumCM +=
	pow(frag[i]->MomCM[k],2);
      frag[i]->momentumCM = sqrt(frag[i]->momentumCM);
    }

  cout << frag[0]->MomCM[0] << " " << frag[0]->MomCM[1] << " " << 
          frag[0]->MomCM[2] <<endl;

  cout << frag[1]->MomCM[0] << " " << frag[1]->MomCM[1] << " " << 
          frag[1]->MomCM[2] <<endl;


  cout << frag[2]->MomCM[0] << " " << frag[2]->MomCM[1] << " " << 
          frag[2]->MomCM[2] <<endl;

  cout << frag[3]->MomCM[0] << " " << frag[3]->MomCM[1] << " " << 
          frag[3]->MomCM[2] <<endl;


  cout << frag[4]->MomCM[0] << " " << frag[4]->MomCM[1] << " " << 
    frag[4]->MomCM[2] << " " << frag[4]->momentumCM << " " <<
    Kinematics.getKE(frag[4]->momentumCM,frag[4]->mass) <<  endl;




  cout << endl;

    
}
//****************************************************************
float correl::getAlphaMom()
{
     frag[4]->momentumCM = 0.;
      for (int k=0;k<3;k++) frag[4]->momentumCM +=
	pow(frag[4]->MomCM[k],2);
      frag[4]->momentumCM = sqrt(frag[4]->momentumCM);
      return frag[4]->momentumCM;
}
//********************************************************************
void correl::rotate()
{
  float phi = atan2(frag[4]->MomCM[1],frag[4]->MomCM[0]);

  //rotate in x-y plane
  for (int j=0;j<5;j++)
    {
      frag[j]->MomRot[0] = frag[j]->MomCM[0]*cos(phi) + frag[j]->MomCM[1]*sin(phi);
      frag[j]->MomRot[1] = frag[j]->MomCM[1]*cos(phi) - frag[j]->MomCM[0]*sin(phi);
      frag[j]->MomRot[2] = frag[j]->MomCM[2];
    }

  float theta = acos(frag[4]->MomCM[2]/frag[4]->momentumCM);


  //rotate in x-y plane
  for (int j=0;j<5;j++)
    {
      frag[j]->MomRot2[2] = frag[j]->MomRot[2]*cos(theta) + frag[j]->MomRot[0]*sin(theta);
      frag[j]->MomRot2[0] = frag[j]->MomRot[0]*cos(theta) - frag[j]->MomRot[2]*sin(theta);
      frag[j]->MomRot2[1] = frag[j]->MomRot[1];
    }
}
