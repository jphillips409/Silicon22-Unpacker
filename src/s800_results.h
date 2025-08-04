#ifndef s800_resultsH
#define s800_resultsH

#include <cmath>

//Struct needed to pass s800 results in for time gates
struct s800_results
{
	int run;
  bool trig_coin;
  bool trig_singles;
  bool trig_s800_singles;
	double obj;
	double XFP;
	double ICsum;
  int Zbeam;
  int Abeam;
  int Zresidue;
  int Aresidue;
  double tstamp;
  void Reset()
  {
		run = 0;
    trig_coin = false;
    trig_singles = false;
    trig_s800_singles = false;
		obj = NAN;
		XFP = NAN;
		ICsum = NAN;
    Zbeam =-1;
    Abeam=-1;
    Zresidue=-1;
    Aresidue=-1;
    tstamp = -1;
  }
}; 

#endif
