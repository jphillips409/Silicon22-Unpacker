//creates energy order (descending) lists of strips

#include "elist.h"
#include <algorithm>
#include <iostream>
using namespace std;

elist::elist()
{
 Nstore = 0;
}

//***********************************************************************
//places a new strip energy in an order list from max to min energy
//updated for high/low gain (HINP4) chips
void elist::Add(int StripNum, float energy, float energylow, int energyRlow, int energyR, float time, float qdc, int qdcflag)
{

  //first find place in list
  int i = 0;
  for (;;)
  {
    if (i == Nstore) break;
    if(energy >0)
    {
      if (energy > Order[i].energy) break;
    }
    else
    {
      if(energyR > Order[i].energyR) break;
    }
    i++;
  }
  if (i == nnn) return; // not enougth room in list 

  //new list length
  int N = min(nnn,Nstore+1);

  // move those in energy below the new value down the list 
  for (int j=N-1;j>i;j--) Order[j] = Order[j-1];

  //add present energy to list
  Order[i].energy = energy;
  Order[i].energyR = energyR;
  Order[i].energylow = energylow;
  Order[i].energylowR = energyRlow;
  Order[i].strip = StripNum;
  Order[i].time = time;
	Order[i].qdc = qdc;

  // increase list length
  Nstore = N;
  mult = N;
}

void elist::Remove(int entry)
{
  if (entry >= Nstore)
  {
    cout << "entry does not exist: " << entry << endl;
    return;
  }
  if (entry == Nstore -1)
  {
    Nstore--;
    return;
  }
  for (int i=entry;i<Nstore;i++)
  {
    Order[i] = Order[i+1];
  }
  Nstore--;
  return;
}

void elist::Add(int StripNum, float energy, int rawenergy, int time, float qdc, int qdcflag)
{
  Order[Nstore].energy = energy;
  Order[Nstore].strip = StripNum;
  Order[Nstore].time = time;
  Order[Nstore].energylow = 0;
  Order[Nstore].energyR = rawenergy;
  Order[Nstore].qdcflag = qdcflag;
  Order[Nstore].qdc = qdc;
  Nstore++;
}

//Copies an entry's qualities onto another index
void elist::Copy(int icopy, int iplace)
{
  Order[iplace].energy = Order[icopy].energy;
  Order[iplace].energyR = Order[icopy].energyR;
  Order[iplace].energylow = Order[icopy].energylow;
  Order[iplace].energylowR = Order[icopy].energylowR;
  Order[iplace].strip = Order[icopy].strip;
  Order[iplace].time = Order[icopy].time;
  Order[iplace].qdcflag = Order[icopy].qdcflag;
  Order[iplace].qdc = Order[icopy].qdc;
}

//**********************************************************************
//remove silicon strip cross talk from list 
int elist::Reduce(const char*face)
{

  if (Nstore <= 1) return Nstore;
  for (int ii = Nstore-1;ii > 0;ii--) // ii location in list of interest
  {
    int xtalk = 0;
    if (Nstore <= 1) return Nstore;
    if (Nstore > nnn) 
    {
      cout << "problem in Reduce" << endl;
      return 0;
    }
    if (Nstore < 0)
    {
      cout << "problem in Reduce()" << endl;
      return 0;
    }
    for (int i=ii-1;i>=0;i--)
    {
      //look for neighboring strips
      if (abs(Order[ii].strip-Order[i].strip) == 1)
      {
        if (*face == 'F')
        {
          if (Order[ii].energy < Order[i].energy*0.044 && Order[i].energy > 10.)
          {
            xtalk = 1; // signal xtalk
            break;
          }
        }
        else if (*face == 'B')
        {
          if (Order[ii].energy < Order[i].energy*0.30 && Order[i].energy > 10.)
          {
            xtalk = 1; // signal xtalk
            break;
          }
        }
      }
    }
    if (xtalk == 0) continue;
    //remove from list
    if (ii+1 < Nstore) 
    {
      for (int j=ii;j<Nstore;j++) {Order[j] = Order[j+1];}
    }

    Nstore--;
  }

  return Nstore;
}
//*********************************************************************
void elist::reset()
{
  for(int i =0;i<Nstore;i++)
  {
    Order[i].energy = 0;
    Order[i].energyR = 0;
    Order[i].energylow = 0;
    Order[i].energylowR = 0;
    Order[i].strip = 0;
    Order[i].time = 0;
    Order[i].qdcflag = 0;
    Order[i].qdc = -1;
    Order[i].neighbours=0;
    Order[i].energyMax=0;
    Order[i].CsIFlag=0;
    Order[i].SiFlag=0;
  }
  Nstore = 0;
  mult = 0;
  
}
//*********************************************************************
  //looks for cross talks events
void elist::Neighbours(int id)
{
  if (Nstore <1) return; // nothing to look at 

/*
  for (int i = 0; i<Nstore; i++)
  {
    cout << "  strip " << Order[i].strip << " E " << Order[i].energy << " Eraw " << Order[i].energyR << " time " << Order[i].time << endl;
  }
*/

  if (Nstore ==1) // no neighbors to look for
  {
    Order[0].energyMax = Order[0].energy;
    Order[0].neighbours = 0;
    return ;
  }
  int i=-1;
  for (;;)
  {
    i++;
    if (i >= Nstore) break;
    Order[i].energyMax = Order[i].energy;
    Order[i].neighbours = 0;

    int j = i;
    for (;;)
    {
      j++;
      if (j >= Nstore) break;
      if (abs(Order[i].strip - Order[j].strip) == 1) //neighboring strips
      {

        //cout << "!!!!!  neighbour addback " << endl;
        Order[i].energy += Order[j].energy; // add energy from adjacent strip
        Order[i].energylow += Order[j].energylow; //adds low gain energy as well
        Order[i].neighbours++;

        //remove this strip from the list
        if (j != Nstore-1)
        {
          //shift all Order objects back one
          for (int k=j+1;k<Nstore;k++) {Order[k-1] = Order[k];}
        }
        Nstore--;
        j--;
      }
    }
  }
}

//cut threshold
void elist::Threshold(float threshold)
{
  for(int i=0;i<Nstore;i++)
  {
    threshold0 = threshold;
    if(Order[i].energy<threshold0)
    {
      Remove(i);
      i--;
    }
  }
}
