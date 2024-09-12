#include "tele.h"

tele::tele(float * r_center[3], float * r_front[3], float * r_back[3], float active_x0, float active_y0)
{
  center[0] = r_center[0];
  center[1] = r_center[1];
  center[2] = r_center[2];

  diagonal_mag = sqrt(pow(r_front[0],2) + pow(r_front[1],2) + pow(r_front[2],2));
  r_back_mag = sqrt(pow(r_back[0],2) + pow(r_back[1],2) + pow(r_back[2],2))

  diagonal[0] = r_front[0]/diagonal_mag;
  diagonal[1] = r_front[1]/diagonal_mag;
  diagonal[2] = r_front[2]/diagonal_mag;

  unit_back[0] = r_back[0]/r_back_mag;
  unit_back[1] = r_back[1]/r_back_mag;
  unit_back[2] = r_back[2]/r_back_mag;

  active_x = active_x0;
  active_y = active_y0;

  normal[0] = diagonal[1]*unit_back[2] - diagonal[2]*unit_back[1];
  normal[1] = diagonal[2]*unit_back[0] - diagonal[0]*unit_back[2];
  normal[2] = diagonal[0]*unit_back[1] - diagonal[1]*unit_back[0];

  unit_front[0] = normal[1]*unit_back[2] - normal[2]*unit_back[1];
  unit_front[1] = normal[2]*unit_back[0] - normal[0]*unit_back[2];
  unit_front[2] = normal[0]*unit_back[1] - normal[1]*unit_back[0];


}
