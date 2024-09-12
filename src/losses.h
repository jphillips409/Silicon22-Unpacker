#ifndef losses_
#define losses_
#include "loss2.h"

class CLosses
{
 private:
   CLoss2 ** loss;
   int Zmax;
 public:
   CLosses(int,string);
   CLosses(int,string,bool);
   ~CLosses();
   float getEin(float,float,int,int);
   float getEout(float,float,int,int);

   
   
};

#endif
