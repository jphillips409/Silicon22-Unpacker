#ifndef corrcomb_
#define corrcomb_
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "math.h"
#include "TSpline.h"
using namespace std;

/**
 * storage of correlation combinations
 */

class corrcomb
{
 public:
  corrcomb(int setting);
  ~corrcomb();

  float getX(int beam, int Z, int A);
  float getY(int beam, int Z, int A, int mult);

  int s800_set;


};
#endif
