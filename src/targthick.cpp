#include "targthick.h"

targthick* targthick::fInstance = 0;

targthick* targthick::instance() 
{
    if (fInstance == 0) {
        fInstance = new targthick;
    }
    return fInstance;
}

/**
 * COnstructor for CRandom
 */
targthick::targthick()
{


}

void targthick::SetThick(int numb)
{
   
   if (numb >= 1155 && numb <= 1158)
   {
     TargType = std::string("_Polypropylene.loss");
     TargetThickness = 3.04;
   }
   else
   {
     TargType = std::string("_Be.loss");
     if (numb < 1050) TargetThickness = 9.472;
     else if (numb > 1073) TargetThickness = 9.472;
     else TargetThickness = 4.625;
   }
}

targthick::~targthick()
{

}

