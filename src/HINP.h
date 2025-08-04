#ifndef HINP_H
#define HINP_H
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "constants.h"

/**
 * !\brief handles the readout from Washu HIMP chips
 * This class deals with the read out from HIMP chips for silicon
 */

class HINP
{
 public:
  
  HINP();
  
  unsigned short NWords;
  //int32_t NWords;
  unsigned short NWords2;
  int NstripsRead;
  
  //max length of the packet. increase if event rate is high or coincidence is high
  static const int maxlen = 500;
  
  unsigned short board[maxlen];
  unsigned short chan[maxlen];
  
  unsigned short high[maxlen];
  unsigned short low[maxlen];
  unsigned short time[maxlen];

  bool unpackSi_HINP4(unsigned short *&point);
};
#endif
