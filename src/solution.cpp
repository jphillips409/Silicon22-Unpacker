#include "solution.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


using namespace std;

void solution::reset()
{
  //variables filled matching events in a telescope
  energy = -1;
  energyR = -1;
  energylow = -1;
  energylowR = -1;
  benergy = -1;
  benergyR = -1;
  benergylow = -1;
  benergylowR = -1;
  denergy = -1;
  qdc = -1;
  ifront = -1;
  iback = -1;
  ide = -1;
  iCsI = -1;
  itele = -1;
  CsITime = -1;
  timediff = -100000.;
  fTime = -1;
  isSiCsI = false;
  
  //variables filled after getPID() from telescope.cpp
  ipid = 0;
  iZ = 0;
  iA = 0;
  mass = 0;
  
  //variables filled after position() and calcEloss() from telescope.cpp
  Xpos = -1;
  Ypos = -1;
  theta = -1;
  phi = -1;
  theta_s800 = -1;
  phi_s800 = -1;
  energyTot = -1;
  Ekin = -1;
  velocity = -1;
  rigidity = -1; //Only used for S800
  beamZ = -1; //Only S800
  beamA = -1; //Only S800
}

void solution::SetTargetDistance(double dist0)
{
  distTarget = (float)dist0;
}

float solution::angle()
{
  float XYZ2 = pow(Xpos , 2) + pow(Ypos, 2) + pow(distTarget,2);
  theta = acos(distTarget/sqrt(XYZ2));
  phi = atan2(Ypos , Xpos);

  return theta;
}
//********************************************************
void solution::getMomentum()
{
  momentum = Kinematics.getMomentum(Ekin,mass);
  Mvect[0] = momentum*sin(theta)*cos(phi);
  Mvect[1] = momentum*sin(theta)*sin(phi);
  Mvect[2] = momentum*cos(theta);



  //scale = 1 einstein, 0 for newton
  energyTot = Ekin*Kinematics.scale + mass;

  velocity = momentum/energyTot;
}
