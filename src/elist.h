#ifndef elist_
#define elist_
#include <string>

using namespace std;

struct order
{
  float energy;
  float energyR;  // high gain channels
  float energylowR;  //low channels
  float energylow;
  float energyMax;
	float qdcflag; //only used for CsI
	float qdc; //only used for CsI
  int strip;
  int neighbours; //I smell australian here
  float time;
  int CsIFlag; //flag for a CsI match
  int SiFlag; //flag for a Si-Si match
};

int const nnn=60;

/**
 * !\brief Energy ordered list
 *
 * This class creates an energy ordered list of the strips
 * read out from a strip detector, keeping track of the strip 
 * numbers that fired.
 */

class elist
{
 public:

  int Nstore = 0; //number stored in list
  order Order[nnn];
  int mult;
  
  elist();
  void Add(int,float,float, int,int,float,float,int);
  void Add(int, float, int, int,float,int);
  void Remove(int);
  void Copy(int, int);
  int  Reduce(const char*);
  void reset();
  void Neighbours(int);
  void Threshold(float);
  float threshold0;
};
#endif
