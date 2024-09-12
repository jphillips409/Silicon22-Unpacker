#ifndef targthick_
#define targthick_
#include <string>

class targthick
{
 protected:

  static  targthick* fInstance; //!< instance member to make tis a singleton
 targthick();


 public:
 
  static targthick* instance();
  ~targthick();

  void SetThick(int numb);
  double TargetThickness;
  std::string TargType;

};


#endif
