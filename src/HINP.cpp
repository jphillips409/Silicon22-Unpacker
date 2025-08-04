#include "HINP.h"

using namespace std;

//unpacking the XLM with ADC on the CHIP BOARDS (HINP 4)
//constructor sets multiplicities to 0
HINP::HINP()
{
  NstripsRead = 0;
}

/* Example in the buffer
13 1ff3 10 0 8 0 0 0 0 28 412 ee 38b1 6b 30f d5 3815 0 a2 d

13          words=19 after this
1ff3        marker, unique to the motherboard
10 0        Nwords=16 (don't know the use of this)
8         NstripsRead, must be divisible by 4
0 0 0 0       *point+5, some junk that needs to be skipped, could be time stamp
28 412 ee 38b1      id,Ehigh,Elow,time
6b 30f d5 3815      id,Ehigh,Elow,time

0 a2 d    just some junk included in the buffer?
*/

bool HINP::unpackSi_HINP4(unsigned short *&point)
{

/*
  unsigned short *endpos = point;
  //point++;
  //point++;
  unsigned short words = *point++;
  endpos += words;

  unsigned short marker;

  marker = *point++;
 
 
  if (marker != 0x1ff3)
  { 
    cout << "Did not read the proper XLM marker. Was " << hex << marker << " expected 0x1ff3" << dec <<endl;
    return false;
  }
  
  NstripsRead = 0;
  NWords = *point++;
  NWords2 = *point++;
  
  //cout << NWords << endl;
  point += NWords;
  return true;
 */

  //point++;
  //point++;
  unsigned short words = *point++;
  //endpos += words;

  unsigned short marker;

  marker = *point++;
 
  //cout << hex << marker << endl;
  if (marker != 0x1ff3)
  { 
    cout << "Did not read the proper XLM marker. Was " << hex << marker << " expected 0x1ff3" << dec <<endl;
    return false;
  }

  unsigned short *endpos = point;  
  NstripsRead = 0;
  NWords = *point;
  //NWords2 = *point++;
  
  //cout << "HINP NWords " << NWords << endl;
  //point += NWords;
  //cout << "point " << *point << " point 2 " << hex << *(point+1) << dec << endl;
  //return true;
  

  //cout << hex << marker << " "  << NWords << dec << endl;
  //unsigned short * endHINP = point;

  if (NWords >= 2048)
  {
    point+=10;
    return false;
  }
  
  endpos += NWords+2;

  point++;
  point++;

  NstripsRead = *point;



  if(NstripsRead %4 !=0)
  {
    cout << "NstripsRead %4 !=0" << endl;
    point+=8;
    return false;
  }
  //cout << "Nstripsread " << hex << NstripsRead << dec << endl;
  NstripsRead /= 4;

  if (NstripsRead > 512) 
  {
    cout << "NstripsRead %4 !=0" << endl;
    point +=8;
    return false; // bad buffer
    //return true; //need to return to this KB - 10/26/20
  }

  if (NstripsRead > maxlen)
  {
    cout << "please increase the value of maxlen" << endl;
    cout << "must be at least " << NstripsRead << endl;
    return false;
  }

  if (Verb)
  {
    cout << "HINP words: " << words << endl;  
    cout << "  Nwords: " << NWords << endl;  
    cout << "  NstripsRead: " << NstripsRead << endl;  
  }

  point += 5;
  //cout << "point after skip " << hex << *point << dec << endl;
  
  //cout << "NStrips " << NstripsRead <<  endl;
  for (int istrip = 0;istrip < NstripsRead;istrip++)
  {
    unsigned short id = *point++;
    unsigned short chipNum = (id&0x1FE0)>>5;
    unsigned short chanNum = id & 0x1F;
    unsigned short ienergy = *point++;
    unsigned short ilowenergy = *point++;
    unsigned short itime =  *point++;
    unsigned short boardNum;
    //cout << hex << id << " " << chipNum << " " << chanNum << " " << ienergy << " " << ilowenergy << dec << endl;

    //cout << "id " << id << endl;
    //cout << "chipNum " <<chipNum << endl;
    //cout << "chanNum " << chanNum << endl;
    //cout << "ienergy " << ienergy << endl;
    //cout << "ilowenergy " << ilowenergy << endl;
    //cout << "itime " << itime << endl;

    //Gobbi telescopes
    //each board has 2 chips, convert so that it is per Board
    if (chipNum%2 == 0)
    { //for even chipnum (second half of chip)
      chanNum = 2*chanNum;//nick
      //chanNum = 16+chanNum;
      //chipNum /= 2;
      boardNum = chipNum/2;
      //cout << "b Num " << boardNum << endl;
    }
    else
    { //for odd chipnum (first half of chip)
      chanNum = 2*chanNum+1; //nick
      //chanNum = chanNum;
      boardNum = chipNum/2 + 1;
      //cout << "b Num " << boardNum << endl;
    }
    //cout << "  boardNum " << boardNum << endl;
    //cout << "  chanNum " << chanNum << endl;

    if (chanNum > 31)
    {
      cout << "chanNum too big" << endl;
      return false;
    }
    if (boardNum > 8)
    {
      cout << "boardNum too big" << endl;
      return false;
    }

    //save information to arrays
    board[istrip] = boardNum;
    chan[istrip] = chanNum;
    high[istrip] = ienergy;
    low[istrip] = ilowenergy;
    time[istrip] = itime;

  } //end loop over strips
    
  point = endpos;

  return true;
}
