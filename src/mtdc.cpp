#include "mtdc.h"
#include <iostream>
#include <iomanip>
using namespace std;

//MTDC-32 resolution values
double mr[10] = {0.0, 0.0, 3.9, 7.8, 15.6, 31.3, 62.5, 125.0, 250.0, 500.0};
/***********************************************************************/
/* Unpacker for Mesytec modules ****************************************/
/* Written primarily for MTDC-32 ***************************************/
/***********************************************************************/
unsigned short * mtdc::read(unsigned short * point) {

  // read header
  unsigned short header2 = *point++;
  unsigned short header1 = *point++;
  unsigned short hsig = (header1 & 0xC000) >> 14;
  //cout << hex << header1 << " " << header2 << " " << hsig << dec << endl;
  if(hsig != 1) {
    cout << "Bad header Word" << endl;
    return point;
  }
  moduleid = header1 & 0xFF;
  tdc_res = (header2 & 0xF000 ) >> 12;
  res = mr[tdc_res];
  nwords = header2 & 0xFFF;
  number = nwords - 1;
  // loop over each channel
  int j=0;
  tstampxtend = 0;
  for (int i = 0; i < nwords-1; i++) {
    unsigned short data2 = *point++;
    unsigned short data1 = *point++;
    if( ((data1 & 0xFF00) >> 8) != 4) {
      number--;
      continue; //Dummy words
    }
    if( ((data1 & 0xFF) >> 4) == 8) {
      tstampxtend = data2 & 0xFFFF;
      channel[j] = 34;
    } else {
      unsigned short Ichan = data1 & 0x3F;
      unsigned short value = data2 & 0xFFFF;
      if(number < 100) channel[j] = Ichan;
      if(number < 100) data[j] = value;
    }
    j++;
  }
  //read end of block
  header2 = *point++;
  header1 = *point++;
  tstamplow = header2;
  tstamphigh = header1 & 0x3FFF;
  
  //cout << header2 << " " << header1 << endl;
  return point;
}
