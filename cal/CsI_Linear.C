{
  ifstream fin1;
  fin1.open("CsI_roughcal_peaks.dat");

  ofstream fout;
  fout.open("CsI_roughcal.dat");

  double y1, y2, x1, x2;

  y1 = 0; //pedestal value for x1
  y2 = 100; //MeV, take max position x2

  for (int j=0;j<4;j++)
  {
    for (int i=0;i<4;i++)
    {
      fin1 >> x1 >> x2;

      double slope = (y2-y1)/(x2-x1);
      double intercept = -x1*slope;

      cout << " " << j << " " << i << " " << slope << " " << intercept << endl;
      fout << j << " " << i << " " << slope << " " << intercept << endl;
    }
} 
}
