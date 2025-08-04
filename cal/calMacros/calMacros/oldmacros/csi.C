{
  ifstream file("csiPunch.dat");

  double y1, y2;
  int ii;
  double x1 = 0.;
  double x2 = 115.8;
  for (int i=0;i<56;i++)
    {
      file >> ii >> y1 >> y2;

      double slope = (x2-x1)/(y2-y1);
      double intercept = -y1*slope;

      cout << 0 << " " << ii << " " << slope << " " << intercept << endl;
    }
}
