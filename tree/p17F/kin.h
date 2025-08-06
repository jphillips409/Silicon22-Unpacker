#if !defined(kin_)
#define kin_ 1
#include "save.h"
#include <cmath>

struct vel
{
  double vv;
  double v[3];
};


struct dvel
{
  double vv;
  double v[3];
};

class kin
{
 public:
  kin(){};
  vel findCM(save,save);
  vel findCM(save,save,save);
  vel findCM(save,save,save,save);
  vel findCM(save,save,save,save,save);
  vel findCM(float[3],float[3],float,float);
  vel findCM(float[3],float[3],float[3],float,float,float);
  vel findCM(float[3],float[3],float[3],float[3],float);
  vel findCM(float[3],float[3],float[3],float[3],float[3],float);
  vel findCM(float[3],float[3],float[3],float[3],float[3],float[3],float);
  vel trans(vel,float[3],float);
  vel trans(vel,save);
  vel transV(vel,vel);
  float getMinv(save,save);

  dvel findCM(double[3],double[3],double,double);
  dvel findCM(double[3],double[3],double[3],double,double,double);
  dvel findCM(double[3],double[3],double[3],double[3],double);
  dvel findCM(double[3],double[3],double[3],double[3],double[3],double);
  dvel findCM(double[3],double[3],double[3],double[3],double[3],double[3],double);
  dvel trans(dvel,double[3],double);
  dvel transV(dvel,dvel);


};
#endif
