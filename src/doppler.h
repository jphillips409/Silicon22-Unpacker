#ifndef doppler_
#define doppler_
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
using namespace std;


/**
 *!\brief performed dopler correction to gamma energies
 */

class doppler
{
 public:
  doppler(float VelSource);
  doppler();
  float correct(float energy, float theta);
  float correct(float energy, float theta,float velocity);
  float correct(float energy, int id);
  float correct(float energy, int id, float velocity);
  float velSource; //!< velocity of gamma source in units of c
  float gamma; //!< 1/sqrt(1-pow(v/c,2))
  float cosTheta[10]; //!< array of cosine of angle angle to each detector
  int Nnai; //!< number of gamma detectors
  void read();

};

#endif
