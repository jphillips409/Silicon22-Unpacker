#include "calibrate.h"

/**
 * Constructor
  \param Nstrip0 is number of strips or csi in a telescope
  \param name is string contain the file of coefficients
  \param order is order of polynomial

 */
calibrate::calibrate(int Ntele0, int Nstrip0, string name, int order0, bool weave, bool bback/*=false*/)
{
  Nstrip = Nstrip0;
  Ntele = Ntele0;
  order = order0;

  //cout << Nstrip << " " << name << endl;

  Coeff = new coeff*[Ntele];
  for (int i=0;i<Ntele;i++)
  {
    Coeff[i] = new coeff [Nstrip];
  }

  ifstream file(name.c_str());
  if (file.is_open() != 1)
  {
    cout << " could not open calibration file " << name << endl;
    abort();
  }

  string title;
  //  getline(file,title);
  // cout << title << endl;
  //getline(file,title);
  //cout << "compare " << name.c_str() << "  " << name.compare("Cal/WVCal.dat") << endl;
  int itele,istrip;
  int board,chan;
  double slope, intercept, a2,a3,qa,qb,qc,qd,qp;
  if ((name.find("CsI") == string::npos || name.find("Time") != string::npos) && name.compare("Cal/WVCal.dat") == 0)
  {
    for(;;)
    {
      file >>  itele >> istrip >> slope >> intercept;
      //cout << itele << " " << istrip << " " << slope << " " << intercept <<endl;

      if (itele >= Ntele || itele < 0 || itele != itele)
      {
        cout <<"problem in calibration.cpp itele >= Ntele, itele = " << itele << " Ntele = " << Ntele << endl; 
        cout <<" file " << name << " " << itele << " " << istrip << " " << slope << " " << intercept << endl;
        abort();
      }
      if(istrip >= Nstrip || istrip < 0 || istrip != istrip)
      {
        cout <<"problem in calibration.cpp istrip >= Nstrip, istrip = " << itele << " Nstrip = " << Nstrip << endl;
        cout <<" file " << name << " " << itele << " " << istrip << " " << slope << " " << intercept << endl; 

        abort();
      }

      if (weave && itele < 4)
      {
        board = istrip/16+1;
        //need to undo old chip# assignment and then redo it to be position correct
        if (board%2 == 0)
        {
          chan = (istrip - 16)*2;
        }
        else
        {
          chan = (istrip)*2+1;
        }
      }
      else if(weave && itele==4 && bback)  //back side of W detor needs it own version o weaving
      {
        if (istrip < 8) chan = 7 - istrip;
        else chan = istrip; 
      }
      else
      {
        chan = istrip;
      }

      //cout << "Board# " << board << " new chip# " << chan << endl;

      if (order >=2) file >> a2;
      else a2 = 0.;
      if (order == 3) file >> a3;
      else a3 = 0.;
      if (file.eof()) break;
      if (file.bad()) break;

      if(chan >= Nstrip || chan < 0 || chan != chan)
      {
        cout <<"problem in calibration.cpp chan >= Nstrip, chan = " << chan << 
             " Nstrip = " << Nstrip << endl; 
        cout <<" file " << name << endl;
        abort();
      }

      Coeff[itele][chan].slope = slope;
      Coeff[itele][chan].intercept = intercept;
      Coeff[itele][chan].a2 = a2;
      Coeff[itele][chan].a3 = a3;
      Coeff[itele][chan].flag = "Linear";
    
    }
  }
  else
  {
    for(;;)
    {

      file >>  itele >> istrip; 
      if (istrip >= 0 && name.compare("Cal/WVCal.dat") != 0)
      {

        file >> slope >> intercept;
      //cout << itele << " " << istrip << " " << slope << " " << intercept <<endl;

        if (itele >= Ntele || itele < 0 || itele != itele)
        {
          cout <<"problem in calibration.cpp itele >= Ntele, itele = " << itele << " Ntele = " << Ntele << endl; 
          cout <<" file " << name << " " << itele << " " << istrip << " " << slope << " " << intercept << endl;
          abort();
        }
        if(istrip >= Nstrip || istrip < 0 || istrip != istrip)
        {
          cout <<"problem in calibration.cpp istrip >= Nstrip, istrip = " << itele << " Nstrip = " << Nstrip << endl;
          cout <<" file " << name << " " << itele << " " << istrip << " " << slope << " " << intercept << endl; 

          abort();
        }

        if (weave && itele < 4)
        {
          board = istrip/16+1;
          //need to undo old chip# assignment and then redo it to be position correct
          if (board%2 == 0)
          {
            chan = (istrip - 16)*2;
          }
          else
          {
            chan = (istrip)*2+1;
          }
        }
        else if(weave && itele==4 && bback)  //back side of W detor needs it own version o weaving
        {
          if (istrip < 8) chan = 7 - istrip;
          else chan = istrip; 
        }
        else
        {
          chan = istrip;
        }

        //cout << "Board# " << board << " new chip# " << chan << endl;

        if (order >=2) file >> a2;
        else a2 = 0.;
        if (order == 3) file >> a3;
        else a3 = 0.;
        if (file.eof()) break;
        if (file.bad()) break;

        if(chan >= Nstrip || chan < 0 || chan != chan)
        {
          cout <<"problem in calibration.cpp chan >= Nstrip, chan = " << chan << 
             " Nstrip = " << Nstrip << endl; 
          cout <<" file " << name << endl;
          abort();
        }

        Coeff[itele][chan].slope = slope;
        Coeff[itele][chan].intercept = intercept;
        Coeff[itele][chan].a2 = a2;
        Coeff[itele][chan].a3 = a3;
        Coeff[itele][chan].flag = "Linear";
      }
      else
      {
        //cout << "here " << name.c_str() << endl;
        file >> qa >> qb >> qc >> qp;

        if (itele >= Ntele || itele < 0 || itele != itele)
        {
          cout <<"problem in calibration.cpp itele >= Ntele, itele = " << itele << " Ntele = " << Ntele << endl; 
          cout <<" file " << name << " " << itele << " " << istrip << " " << slope << " " << intercept << endl;
          abort();
        }
        if(istrip >= Nstrip || istrip < 0 || istrip != istrip)
        {
          cout <<"problem in calibration.cpp istrip >= Nstrip, istrip = " << itele << " Nstrip = " << Nstrip << endl;
          cout <<" file " << name << " " << itele << " " << istrip << " " << slope << " " << intercept << endl; 

          abort();
        }

        if (weave && itele < 4)
        {
          board = istrip/16+1;
          //need to undo old chip# assignment and then redo it to be position correct
          if (board%2 == 0)
          {
            chan = (istrip - 16)*2;
          }
          else
          {
            chan = (istrip)*2+1;
          }
        }
        else if(weave && itele==4 && bback)  //back side of W detor needs it own version o weaving
        {
          if (istrip < 8) chan = 7 - istrip;
          else chan = istrip; 
        }
        else
        {
          chan = istrip;
        }


        if (file.eof()) break;
        if (file.bad()) break;

        if(chan >= Nstrip || chan < 0 || chan != chan)
        {
          cout <<"problem in calibration.cpp chan >= Nstrip, chan = " << chan << 
             " Nstrip = " << Nstrip << endl; 
          cout <<" file " << name << endl;
          abort();
        }
        //cout << qa << " " << qb << " " << qc << " " << qp << endl;
        Coeff[itele][chan].qa = qa;
        Coeff[itele][chan].qb = qb;
        Coeff[itele][chan].qc = qc;
        Coeff[itele][chan].qp = qp;
        Coeff[itele][chan].flag = "Quench";
        //Coeff[itele][chan].flag = "pol3";
      }
     
    }
  } 

  file.close();
  file.clear();  

}
//*****************************************************
  /**
   * destructor
   */
calibrate::~calibrate()
{
  for (int i=0;i<Ntele;i++)
  {
    delete [] Coeff[i];
  }

  delete [] Coeff;
}
//*****************************************
  /**
   * returns the calibrated energy
\param istrip - number of the strip or detector
\param channel - raw channels from the ADC, etc
  */
float calibrate::getEnergy(int itele,int istrip,float channel)
{

  if(itele < 0. || itele >= Ntele || itele != itele) 
   {
    cout << "problem in getEnergy, itele = "<< itele << "  Ntele= "  << Ntele << endl;
   return 0.;
   }
  if(istrip < 0. || istrip >= Nstrip || istrip != istrip) 
   {
    cout << "problem in getEnergy, istrip = "<< istrip << " Nstrip= " << Nstrip << endl;
   return 0.;
   }
  
  if (Coeff[itele][istrip].flag == "Linear")
  {
    float fact = channel*Coeff[itele][istrip].slope + Coeff[itele][istrip].intercept;
    if (order == 1) return fact;

    fact += pow(channel,2)*Coeff[itele][istrip].a2;
    if (order == 2) return fact;
    if (order == 3)return pow(channel,3)*Coeff[itele][istrip].a3 + fact;
    else abort();
  }
  else
  {
    float fact;
    if (Coeff[itele][istrip].flag == "Quench") fact = Coeff[itele][istrip].qa*(channel-Coeff[itele][istrip].qp)+Coeff[itele][istrip].qb*log(1+Coeff[itele][istrip].qc*(channel-Coeff[itele][istrip].qp));

    if (Coeff[itele][istrip].flag == "pol3") fact = Coeff[itele][istrip].qa + Coeff[itele][istrip].qb*channel + Coeff[itele][istrip].qc*pow(channel,2) + Coeff[itele][istrip].qp*pow(channel,3);
    return fact;
  }
}

float calibrate::getTime(int itele,int istrip,float channel)
{
  if(itele < 0. || itele >= Ntele || itele != itele) 
   {
    cout << "problem in getTime, itele = "<< itele << "  Ntele= "  << Ntele << endl;
   return 0.;
   }
  if(istrip < 0. || istrip >= Nstrip || istrip != istrip) 
   {
    cout << "problem in getTime, istrip = "<< istrip << " Nstrip= " << Nstrip << endl;
   return 0.;
   }
  return channel + Coeff[itele][istrip].intercept;
}

float calibrate::reverseCal(int itele, int istrip, float energy)
{
  float fact = (energy - Coeff[itele][istrip].intercept)/Coeff[itele][istrip].slope;
  return fact;
}

