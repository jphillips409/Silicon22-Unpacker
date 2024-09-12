#include "doppler.h"

/**
 * constructor
\param velSource0 - velocity of gamma source in units of c
*/
doppler::doppler(float velSource0)
{
  velSource = velSource0;
  gamma = 1./sqrt(1.-pow(velSource,2));
  // read();
}
//*****************************************************
  /**
   * alternate constructor
   */
doppler::doppler()
{
  read();
}
  //**************************************************
    /**
     * reds in the detector angles
     */
void doppler::read()
{
  ifstream f("cal/Nai.ang");
  if (!f.is_open())
  {
    cout << "NaI angles not found " << endl;
    abort();
  }
  int k;
  Nnai = 0;
  float thetaDeg;
  for (int i=0;i<10;i++)
  {
    f >> k >> thetaDeg;
    if (f.eof()) break;
    if (f.bad()) break;
    if (k >= 10) break;
    cosTheta[k] = cos(thetaDeg*0.01745);
    Nnai++;
    //cout << k << " " << Nnai << endl;
  }
}
//*****************************************************
  /**
   * returns the gamma energy in the frame of the source
   * specified in the constructor
\param energy of detected gamma ray
\param theta is angle in radians of detected gamma ray
  */
float doppler::correct(float energy, float theta)
{
  return gamma*energy*(1.-velSource*cos(theta));
}
//****************************************************
  /**
   * returns the gamma energy in the frame specified by the 
   * velocity.
\param energy is the energy of the detected gamma ray
\param theta is the angle of teh detected gamma ray
\param velocity is the velocity of the source in units of c
  */
float doppler::correct(float energy, float theta, float velocity)
{
  float gamma_v = 1./sqrt(1.-pow(velocity,2));
 
   return gamma_v*energy*(1.-velocity*cos(theta));
}
//*****************************************************
  /**
   * returns the gamma energy in the source frame spectified
   * in the constructor
   \param energy is the energy of the detected gamma ray
   \param id is the number of the detector
  */
float doppler::correct(float energy, int id)
{
  return gamma*energy*(1.-velSource*cosTheta[id]);
}
//****************************************************
  /**
   * returns the gamma energy in the frame specified by the 
   * velocity.
\param energy is the energy of the detected gamma ray
\param id is the number of the gamma-ray detector
\param velocity is the velocity of the source in units of c
  */
float doppler::correct(float energy, int id, float velocity)
{
  float gamma_v = 1./sqrt(1.-pow(velocity,2));

  return gamma_v*energy*(1.-velocity*cosTheta[id]);
}
