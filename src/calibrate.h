#ifndef calibrate_
#define calibrate_
#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "math.h"
#include "TSpline.h"
using namespace std;

/**
 * storge of calibration coefficients
 */
struct coeff
{
  float slope; //!< slope for calibration
  float intercept; //!< intercept for calibration
  float a2;  //!< quadratic coeff if needed
  float a3; //!< cubic coeff if needed

  string flag; // flag for CsI vs silicon calibrations, CsI uses different equation

  //If using the CsI quenching equation
  double qa;
  double qb;
  double qc;
  double qd;
  double qp;
};

class calibrate
{
 public:
  calibrate(int Ntele,int Nstrip,string file,int order,bool weave, bool bback=false);
  ~calibrate();
  float getEnergy(int itele,int istrip,float channel);
  float getTime(int itele,int istrip,float channel);
  float reverseCal(int itele, int istrip, float energy);
  int order;
  int Nstrip;  //!< number of strips
  int Ntele;   //!<number of telescopes
  coeff ** Coeff;  //!< array with calibration coefficients for each strip

};
#endif
