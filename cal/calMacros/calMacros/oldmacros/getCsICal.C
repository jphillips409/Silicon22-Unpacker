void getCsICal()
{
  
  ifstream Lfile("CsIdats/H_AL_Centroids.dat");
  ifstream Hfile("CsIdats/H_Be_Centroids.dat");
  ifstream Efile("CsIdats/CsI_H_Energy.txt");
  if(!Lfile || !Hfile || !Efile)
    {
      cout << "Couldn't find the data" << endl;
      return;
    }
  
  float posL[20] = {0.};
  float posH[20] = {0.};
  float EL[20] = {0.};
  float EH[20] = {0.};
  float Pslope[20] = {0.};
  float Pinter[20] = {0.};
  //float energy[2] = {52.0077,71.4806}; //Protons
  ofstream ofile("output.cal");
  int idum = 0;
  int idum2 = 0;
  float fdum =0.;
  float fdum2 =0.;

  char cdum[200];
  Efile.getline(cdum,200);
  
  for(;;)
    {
      if(Efile.eof())break;
      Efile >> idum >> fdum >> fdum2;
      EL[idum] = fdum2;
      EH[idum] = fdum;

    }
  Efile.close();
  for(;;)
    {
      if(Lfile.eof())break;
      Lfile >>idum >> fdum;
      posL[idum] = fdum;
    }
  for(int j =0;j<20;j++)
    {
      if(Hfile.eof())break;
      Hfile >>idum >> fdum;
      posH[idum] = fdum;
    }
  
  float slope = 0;
  float inter = 0;
  float N = 0;
  float aveS = 0;
  float aveI = 0;
  for(int ii = 0; ii<20;ii++)
    {

      slope = (EH[ii]-EL[ii])/(posH[ii]-posL[ii]);
      inter = EH[ii] - slope*posH[ii];
      //slope = (energy[1] - energy[0])/(posH[ii]-posL[ii]); //proton
      //inter = energy[0] - posL[ii]*slope; //proton
      N +=1;
      aveS += slope;
      aveI += inter;
      // cout << slope << " " << inter << endl;
      ofile << 0 << " "  << ii << " " << slope;
      ofile << " " << inter << endl;
      
    }
  
  cout << aveS/N << " " << aveI/N << endl;
  Lfile.close();
  Hfile.close();
  ofile.close();
  return;
}
