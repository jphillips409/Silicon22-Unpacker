#if !defined(save_)
#define save_ 1

class save{
public:
float M[3];
float weight;
float et;
short  id;
bool elast;
save(short,float,float[3]);
save();
inline save operator=(save a)
   {
      id = a.id;
      weight = a.weight;
      et = a.et;
      for (int k=0;k<3;k++) M[k] = a.M[k];
      return a;
   }


};




#endif
