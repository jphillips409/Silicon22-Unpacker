#include "corrcomb.h"

/**
 * Constructor reads in s800 setting, only 1 is valid currently
 */
corrcomb::corrcomb(int setting)
{

  s800_set = setting;

}
//************************
  /**
   * destructor
   */
corrcomb::~corrcomb()
{

}
/**
 * returns 0 if s800_set != 1, 0 won't be plotted
 * returns x and y values for 2D plots of correlation combinations
 */
//*********************************
float corrcomb::getX(int beam, int Z, int A)
{

  if (s800_set != 1) return 0;

  //Returns x value for different residues
  if (beam == 14) {

    //Currently 12 options for Si-23 beam

    if (Z == 10) {
      if (A == 17) return 0.5;    
      if (A == 18) return 1.5; 
      if (A == 19) return 2.5; 
      else return 0;
    }

    if (Z == 11) {
      if (A == 20) return 3.5; 
      if (A == 21) return 4.5;
      else return 0; 
    }

    if (Z == 12) {
      if (A == 20) return 5.5; 
      if (A == 21) return 6.5;
      if (A == 22) return 7.5;
      else return 0; 
    }

    if (Z == 13) {
      if (A == 22) return 8.5;
      if (A == 23) return 9.5;
      else return 0; 
    }

    if (Z == 14) {
      if (A == 23) return 10.5;
      if (A == 24) return 11.5;
      else return 0; 
    }

  }

  if (beam == 13) {

    //Currently 12 options for Si-23 beam

    if (Z == 13) {
      if (A == 22) return 0.5;
      else return 0; 
    }

    if (Z == 12) {
      if (A == 21) return 1.5;
      if (A == 20) return 2.5;
      else return 0; 
    }

    if (Z == 11) {
      if (A == 21) return 3.5;
      if (A == 20) return 4.5;
      else return 0; 
    }

    if (Z == 10) {
      if (A == 19) return 5.5;
      if (A == 18) return 6.5;
      if (A == 17) return 7.5;
      else return 0; 
    }

  }
  
  return 0;

}


float corrcomb::getY(int beam, int Z, int A, int mult)
{

  if (s800_set != 1) return 0;

  //Returns y value for different light fragment combinations
  if (beam == 14 || beam == 13) {

    //Currently 7 options

    if (Z == 1) {
      if (A == 1 && mult == 1) return 0.5;    
      if (A == 1 && mult == 2) return 1.5;    
      if (A == 1 && mult == 3) return 2.5;    
      if (A == 1 && mult == 4) return 3.5;    
      if (A == 2 && mult == 1) return 4.5;    
      if (A == 3 && mult == 1) return 5.5;    
      else return 0;
    }

    if (Z == 2) {
      if (A == 4 && mult == 1) return 6.5;       
      else return 0;
    }

  }
  
  return 0;

}
