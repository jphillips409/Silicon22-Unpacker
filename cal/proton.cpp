#include <sstream>
#include <fstream>

int main()
{
  ifstream fileIn("proton.dat");
  ofstream fileOut("proton2.dat");
  int id;
  float low,high;
  float Elow = 69.81;
  float Ehigh = 78.67;
 
  
  for (int i=0;i<20;i++)
    {
      fileIn >> id >> low >> high;
      float slope = (78.67-69.81)/(high-low);
      float intercept = low*slope - Elow;
      fileOut << 0 << " " << i << " " << slope << " " << intercept << endl;
    }
  
}
