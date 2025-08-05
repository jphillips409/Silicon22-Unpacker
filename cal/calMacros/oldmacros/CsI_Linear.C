{
  ifstream file("../CsI_roughcal_peaks.dat");
  ofstream fout;
  fout.open("../CsI_roughcal.cal");

  double y1, y2, x1, x2;

  y1 = 0; //pedestal value for x1
  y2 = 100; //MeV, take max position x2

  for (int j=0;j<4;j++)
  {
    for (int i=0;i<4;i++)
    {
      file >> x1 >> x2;

      double slope = (y2-y1)/(x2-x1);
      double intercept = -x1*slope;

      cout << " " << j << " " << i << " " << slope << " " << intercept << endl;
      fout << slope << " " << intercept << endl;
    }
} 
}
