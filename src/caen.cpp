#include "caen.h"
#include <iostream>
#include <iomanip>
using namespace std;



//*********************************************************************
unsigned short * caen::read(unsigned short * point)
{

  // read header
  unsigned short header2 = *point++;
  unsigned short header1 = *point++;
  //cout << "header 2 " << hex << header2 << dec << endl;
 // cout << "header 1 " << hex << header1 << dec << endl;

  //extract information
  unsigned short bit =  (header1>>9) & 0x1;
  //if (bit == 0) cout << "header not found" << endl;
  geographic = header1>>11;
  //cout << "geographic " << hex << geographic << dec << "  " << geographic << endl;

  crate = (header1 >> 3) & 0x7F;
  number = header2 >> 8;

  //cout << "crate " << hex << crate << dec << "  " << crate << endl;
  //cout << "number " << hex << header1 << dec << "  " << number << endl;

  // loop over each channel
  int j=0;
  for (int i=0;i<number;i++)
  {
    unsigned short data2 = *point++;
    unsigned short data1 = *point++;
    bit = (data1>>10) & 0x1;
    if (bit)
	  {
      //this is a trailer block
      //cout << "missing caen data" << endl;
  	  point -= 2;
  	  number = j;
  	  break;
	  }
    unsigned short Ichan = data1 & 0x1F;
    unsigned short value = data2 & 0xFFF;
    channel[j] = Ichan;
    underflow[j] = (data2>>13)&0x1;
    overflow[j] = (data2>>12)&0x1;
    //if (overflow[j]) cout << "overflow "<< value << endl;
    //if (underflow[j]) cout << "underflow " << value << endl;

    data[j] = value;
    j++;

    //cout << Ichan << " " << value << endl;
  }

  //read end of block
  header2 = *point++;
  header1 = *point++;
  //cout << "trailer " << hex << header1 << " " << header2 << endl;
  bit = (header1>>10) & 0x1;
  if (bit == 0) 
  {
   cout << "caen End-of-Block not found for slot "  
		     << geographic << endl;
   //abort();
  }
  //cout << " ADC point " << *point << " point 2 " << *(point+1) << endl;
  return point;
}

