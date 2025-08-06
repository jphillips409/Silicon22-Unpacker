#include "save.h"

save::save(short id0, float et0, float M0[])
{
  id = id0;
  et = et0;
  for (int k=0;k<3;k++) M[k] = M0[k];

}
save::save(){};
