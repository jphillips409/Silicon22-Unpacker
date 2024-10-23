void dum()
{

  ifstream infile("Pies.cal");
  ofstream ofile("a.dat");
  float slope[128]={0};
  float inter[128]={0};

  int dum,dum2;

  for(int i=0;i<128;i++)
    {
      infile >> dum2 >> dum >> slope[i] >> inter[i];
    }

  float a = 8.632;
  float b = 5.478;
  float bprime = 6.595;

  float slope2[128] = {0};
  float inter2[128] = {0};

  for(int i=0;i<128;i++)
    {
      if(slope[i] ==-1)
	{
	  slope2[i] = -1;
	  inter2[i] = 0;
	}
      else
	{
	  slope2[i] = slope[i]*(a-bprime)/(a-b);
	  inter2[i] = bprime + slope2[i]*((inter[i]-b)/(slope[i]));
	}
      ofile <<"0 " << i << " " << slope2[i] << " " << inter2[i] << endl;
    }


  
  return;
}
