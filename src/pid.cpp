#include "pid.h"

/**
 * Constructor reads in banana gates
 */
pid::pid(string file, bool S800)
{
  string name;

  if (S800 == false) name = "zline/Gobbi/"+file +".zline";
  if (S800 == true) name = "zline/"+file +".zline";
  cout << "opening zline file: " << name << endl;

  ifstream ifile(name.c_str());
  if (!ifile.is_open()) 
  {
    cout << "could not open zline file " << name << endl;
    nlines = 0;
    return;
  }
  else ifile >> nlines;
  par = new ZApar*[nlines];
  for (int i=0;i<nlines;i++)
  {
    par[i] = new ZApar(ifile);
  }
  ifile.close();
  ifile.clear();
  Z=-1;
  A=-1;
}
//************************
  /**
   * destructor
   */
pid::~pid()
{
  for (int i=0;i<nlines;i++) delete par[i];
  delete [] par;
}
/**
 * returns true if particle is in a banana gate.
 * returns false otherwise. The parameters Z and A
 * are loaded with the detected particle's values.
\param x energy of particle
\param y energy loss of particle
 */
//*********************************
bool pid::getPID(float x, float y)
{
  if (nlines == 0) return false;
  Z = 0;
  A = 0;
  for (int i=0;i<nlines;i++)
  {
    if (par[i]->inBanana(x,y))
    {
      Z = par[i]->Z;
      A = par[i]->A;
      mass = getMass(Z,A);
      return true;
    }
  }
  return false;
}

//Returns whether the F and B energy match inside the energy zgate
bool pid::getEGate(float x, float y)
{
  if (nlines == 0) return false;
  
  for (int i=0;i<nlines;i++)
  {
    if (par[i]->inBanana(x,y)) return true;
  }

  return false;
}

float pid::getMass(int iZ,int iA)
{
  //All masses are hard coded in. If adding new isotope, just copy/paste
  if (iZ == 1)
  {
    if (iA == 1)  return Mass_p;
    else if (iA == 2) return Mass_d;
    else if (iA == 3) return Mass_t;
    else if (iA == 4) return 4.*931.; //This is purely for testing an anomalous zline
    else
    {
      cout << "No mass info for Z = "<< iZ << " A =" << iA << endl;
      abort(); 
    }
  }
  else if (iZ == 2)
  {
    if (iA == 3) return Mass_3He;
    else if (iA == 4) return Mass_alpha;
    else if (iA == 6) return Mass_6He;
    else if (iA == 8) return Mass_8He;
    else 
    {
      cout << "No mass info for Z = "<< iZ << " A =" << iA << endl;
      abort(); 
    }
  }
  else if (iZ == 3)
  {
    if (iA == 4) return 4. *931.; //Purely for testing anomalous zline
    if (iA == 6) return Mass_6Li;
    else if(iA == 7) return Mass_7Li;
    else if (iA == 8) return Mass_8Li;
    else if (iA == 9) return Mass_9Li;
    else 
    {
      cout << "No mass info for Z = "<< iZ << " A =" << iA << endl;
      abort(); 
    }
  }
  else if (iZ == 4)
  {
    if (iA == 7) return Mass_7Be;
    else if (iA == 8) return Mass_8Be;
    else if (iA == 9) return Mass_9Be;
    else if (iA == 10) return Mass_10Be;
    else if (iA == 11) return Mass_11Be;
    else 
    {
      cout << "No mass info for Z = "<< iZ << " A =" << iA << endl;
      abort(); 
    }
  }
  else if (iZ==5)
  {
    if (iA == 8) return Mass_8B;
    else if (iA == 10) return Mass_10B;
    else if (iA == 11) return Mass_11B;
    else 
    {
      cout << "No mass info for Z = "<< iZ << " A =" << iA << endl;
      abort(); 
    }
  }
  else if (iZ == 6)
  {
    if (iA == 9) return Mass_9C;
    else if (iA == 10) return Mass_10C;
    else if (iA == 11) return Mass_11C;
    else if (iA == 12) return Mass_12C;
    else if (iA == 13) return Mass_13C;
    else if (iA == 14) return Mass_14C;
    else if (iA == 15) return Mass_15C;
    else 
    {
      cout << "No mass info for Z = "<< iZ << " A =" << iA << endl;
      abort(); 
    }
  }
  else if (iZ == 7)
  {
    if (iA == 11) return Mass_11N;
    else if (iA == 12) return Mass_12N;
    else if (iA == 13) return Mass_13N;
    else if (iA == 14) return Mass_14N;
    else if (iA == 15) return Mass_15N;
    else if (iA == 16) return Mass_16N;
    else if (iA == 17) return Mass_17N;
    else 
    {
      cout << "No mass info for Z = "<< iZ << " A =" << iA << endl;
      abort(); 
    }
   }
  else if (iZ == 8)
  {
    if (iA == 13) return Mass_13O;
    else if(iA == 14) return Mass_14O;
    else if (iA == 15) return Mass_15O;
    else if (iA == 16) return Mass_16O;
    else if (iA == 17) return Mass_17O;
    else if (iA == 18) return Mass_18O;
    else 
    {
      cout << "No mass info for Z = "<< iZ << " A =" << iA << endl;
      abort(); 
    }
  }
  else if (iZ == 9)
  {
    if (iA == 14) return Mass_14F;
    else if (iA == 15) return Mass_15F;
    else if (iA == 16) return Mass_16F;
    else if (iA == 17) return Mass_17F;
    else if (iA == 18) return Mass_18F;
    else if (iA == 19) return Mass_19F;
    else 
    {
      cout << "No mass info for Z = "<< iZ << " A =" << iA << endl;
      abort(); 
    }
  }
  else
  {
    cout << "No mass info for Z = "<< iZ << " A =" << iA << endl;
    abort();
    return 0;
  }
}

//Corrects each element for its pulse-height defect
double pid::PHDCorrect(int iZ,double energy)
{
  if (iZ == 8)
  {
    double a = 0.5696272;
    double b = -1.69875;
    double PHD = pow(10.,b)*pow(energy,a);
    return energy + PHD;
  }

  if (iZ == 7)
  {
    double a = 0.569293;
    double b = -1.953214;
    double PHD = pow(10.,b)*pow(energy,a);
    return energy + PHD;
  }

  if (iZ == 6)
  {
    double a = 0.5690028;
    double b = -2.7675;
    double PHD = pow(10.,b)*pow(energy,a);
    return energy + PHD;
  }

  return energy;
}

