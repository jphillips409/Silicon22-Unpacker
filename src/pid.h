#ifndef pid_
#define pid_

#include <string>
#include <fstream>
#include <iostream>
#include "ZApar.h"
#include "constants.h"

using namespace std;
/**
 * !\brief detemine PID (Particle identification from E-DE map
 *
 * stores banana gate for each particle type
 */


class pid
{
 public:
  pid(string file); 
  ~pid();
  ZApar ** par;  //!< individual banana gates
  int nlines=0;    //!< number of banana gated stored 
  bool getPID(float x, float y);
  bool getEGate(float x, float y);
  float getMass(int iZ,int iA);
  double PHDCorrect(int iZ, double energy);
  int Z; //!< Z of particle in gate
  int A; //!< A of particle in gate
  float mass; //!< mass of particle in amu
};


#endif
