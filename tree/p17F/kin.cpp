#include "kin.h"
#include "constants.h"

vel kin::findCM(float M1[3], float M2[3], float et1, float et2)
{
      vel CM;
      double Mcm[3];
      double MMcm = 0.;
      for (int j=0;j<3;j++)
	{
	  Mcm[j] = M1[j]+M2[j];
	  MMcm += pow(Mcm[j],2);
	}
      MMcm = sqrt(MMcm);

      float et = et1+et2;
      CM.vv = MMcm*c/et;


      for (int j=0;j<3;j++)
	{
	  CM.v[j] = Mcm[j]/MMcm*CM.vv;
	}
      return CM;
}

//*********************************
vel kin::findCM(save save1,save save2)
{
  return findCM(save1.M,save2.M,save1.et,save2.et);
}
//*********************************
vel kin::findCM(save save1,save save2,save save3)
{
  return findCM(save1.M,save2.M,save3.M,save1.et,save2.et,save3.et);
}
//********************************************************************************


vel kin::findCM(save save1,save save2,save save3, save save4)
{
  return findCM(save1.M,save2.M,save3.M,save4.M,save1.et+save2.et+save3.et+save4.et);
}
//**********************************************************************
vel kin::findCM(save save1,save save2,save save3, save save4, save save5)
{
  return findCM(save1.M,save2.M,save3.M,save4.M,save5.M,save1.et+save2.et+save3.et+save4.et+save5.et);
}
//**********************************************************************
vel kin::findCM(float M1[3], float M2[3], float M3[3],float et1, float et2, float et3)
{
      vel CM;
      double Mcm[3];
      double MMcm = 0.;
      for (int j=0;j<3;j++)
	{
	  Mcm[j] = M1[j]+M2[j]+M3[j];
	  MMcm += pow(Mcm[j],2);
	}
      MMcm = sqrt(MMcm);

      double et = et1+et2+et3;
      CM.vv = MMcm*c/et;


      for (int j=0;j<3;j++)
	{
	  CM.v[j] = Mcm[j]/MMcm*CM.vv;
	}
      return CM;
}


vel kin::findCM(float M1[3], float M2[3], float M3[3],float M4[3],float et)
{
      vel CM;
      double Mcm[3];
      double MMcm = 0.;
      for (int j=0;j<3;j++)
	{
	  Mcm[j] = M1[j]+M2[j]+M3[j]+M4[j];
	  MMcm += pow(Mcm[j],2);
	}
      MMcm = sqrt(MMcm);


      CM.vv = MMcm*c/et;


      for (int j=0;j<3;j++)
	{
	  CM.v[j] = Mcm[j]/MMcm*CM.vv;
	}
      return CM;
}


vel kin::findCM(float M1[3], float M2[3], float M3[3],float M4[3],float M5[3],float et)
{
      vel CM;
      double Mcm[3];
      double MMcm = 0.;
      for (int j=0;j<3;j++)
	{
	  Mcm[j] = M1[j]+M2[j]+M3[j]+M4[j]+M5[j];
	  MMcm += pow(Mcm[j],2);
	}
      MMcm = sqrt(MMcm);


      CM.vv = MMcm*c/et;


      for (int j=0;j<3;j++)
	{
	  CM.v[j] = Mcm[j]/MMcm*CM.vv;
	}
      return CM;
}



vel kin::findCM(float M1[3], float M2[3], float M3[3],float M4[3],float M5[3],float M6[3],float et)
{
      vel CM;
      double Mcm[3];
      double MMcm = 0.;
      for (int j=0;j<3;j++)
	{
	  Mcm[j] = M1[j]+M2[j]+M3[j]+M4[j]+M5[j]+M6[j];
	  MMcm += pow(Mcm[j],2);
	}
      MMcm = sqrt(MMcm);


      CM.vv = MMcm*c/et;


      for (int j=0;j<3;j++)
	{
	  CM.v[j] = Mcm[j]/MMcm*CM.vv;
	}
      return CM;
}


vel kin::trans(vel CM,save save1)
{
  return trans(CM,save1.M,save1.et);
}


vel kin::trans(vel CM,float M[3], float et)
{
      //transform alpha1 to com frame
      double dot = 0.;
      for (int j=0;j<3;j++)
	{
	  dot += M[j]*CM.v[j];
	}



      double para[3];
      double perp[3];
      double Para = dot/CM.vv;
      for (int j=0;j<3;j++)
	{
	  para[j] = Para/CM.vv*CM.v[j];
          perp[j] = M[j] - para[j];

	}


      double gamma = 1./sqrt(1.-pow(CM.vv/c,2));
      double paraNew = (Para - et*CM.vv/c)*gamma;

     
      vel ee;  //stores energy

      for (int j= 0;j<3;j++)
	{
	  ee.v[j] = perp[j] + paraNew/Para*para[j]; 
	}
     
      ee.vv = gamma*(et - CM.vv*Para/c);
      //ee.vv -=  m0*(double)A;

      return ee;


}

vel kin::transV(vel CM,vel start)
{
      //transform alpha1 to com frame
      double dot = 0.;
      for (int j=0;j<3;j++)
	{
	  dot += start.v[j]*CM.v[j];
	}



      double para[3];
      double perp[3];
      double Para = dot/CM.vv;
      double Perp = 0.;
      for (int j=0;j<3;j++)
	{
	  para[j] = Para/CM.vv*CM.v[j];
          perp[j] = start.v[j] - para[j];
          Perp += pow(perp[j],2);
	}
      Perp = sqrt(Perp);

      double gamma = 1./sqrt(1.-pow(CM.vv/c,2));

      vel stop;
      stop.vv = 0.;
      for (int j= 0;j<3;j++)
	{

          perp[j] = perp[j]/(1.-Para*CM.vv/pow(c,2))/gamma;

	  para[j] = (para[j]-CM.v[j])/(1-Para*CM.vv/pow(c,2));

	  stop.v[j] = perp[j] + para[j];
          stop.vv += pow(stop.v[j],2);
	}
     
      stop.vv = sqrt(stop.vv);
      return stop;


}

dvel kin::findCM(double M1[3], double M2[3], double et1, double et2)
{
      dvel CM;
      double Mcm[3];
      double MMcm = 0.;
      for (int j=0;j<3;j++)
	{
	  Mcm[j] = M1[j]+M2[j];
	  MMcm += pow(Mcm[j],2);
	}
      MMcm = sqrt(MMcm);

      double et = et1+et2;
      CM.vv = MMcm*c/et;


      for (int j=0;j<3;j++)
	{
	  CM.v[j] = Mcm[j]/MMcm*CM.vv;
	}
      return CM;
}



dvel kin::findCM(double M1[3], double M2[3], double M3[3],double et1, double et2, double et3)
{
      dvel CM;
      double Mcm[3];
      double MMcm = 0.;
      for (int j=0;j<3;j++)
	{
	  Mcm[j] = M1[j]+M2[j]+M3[j];
	  MMcm += pow(Mcm[j],2);
	}
      MMcm = sqrt(MMcm);

      double et = et1+et2+et3;
      CM.vv = MMcm*c/et;


      for (int j=0;j<3;j++)
	{
	  CM.v[j] = Mcm[j]/MMcm*CM.vv;
	}
      return CM;
}


dvel kin::findCM(double M1[3], double M2[3], double M3[3],double M4[3],double et)
{
      dvel CM;
      double Mcm[3];
      double MMcm = 0.;
      for (int j=0;j<3;j++)
	{
	  Mcm[j] = M1[j]+M2[j]+M3[j]+M4[j];
	  MMcm += pow(Mcm[j],2);
	}
      MMcm = sqrt(MMcm);


      CM.vv = MMcm*c/et;


      for (int j=0;j<3;j++)
	{
	  CM.v[j] = Mcm[j]/MMcm*CM.vv;
	}
      return CM;
}


dvel kin::findCM(double M1[3], double M2[3], double M3[3],double M4[3],double M5[3],double et)
{
      dvel CM;
      double Mcm[3];
      double MMcm = 0.;
      for (int j=0;j<3;j++)
	{
	  Mcm[j] = M1[j]+M2[j]+M3[j]+M4[j]+M5[j];
	  MMcm += pow(Mcm[j],2);
	}
      MMcm = sqrt(MMcm);


      CM.vv = MMcm*c/et;


      for (int j=0;j<3;j++)
	{
	  CM.v[j] = Mcm[j]/MMcm*CM.vv;
	}
      return CM;
}


dvel kin::findCM(double M1[3], double M2[3], double M3[3],double M4[3],double M5[3],double M6[3],double et)
{
      dvel CM;
      double Mcm[3];
      double MMcm = 0.;
      for (int j=0;j<3;j++)
	{
	  Mcm[j] = M1[j]+M2[j]+M3[j]+M4[j]+M5[j]+M6[j];
	  MMcm += pow(Mcm[j],2);
	}
      MMcm = sqrt(MMcm);


      CM.vv = MMcm*c/et;


      for (int j=0;j<3;j++)
	{
	  CM.v[j] = Mcm[j]/MMcm*CM.vv;
	}
      return CM;
}



dvel kin::trans(dvel CM,double M[3], double et)
{
      //transform alpha1 to com frame
      double dot = 0.;
      for (int j=0;j<3;j++)
	{
	  dot += M[j]*CM.v[j];
	}



      double para[3];
      double perp[3];
      double Para = dot/CM.vv;
      for (int j=0;j<3;j++)
	{
	  para[j] = Para/CM.vv*CM.v[j];
          perp[j] = M[j] - para[j];

	}


      double gamma = 1./sqrt(1.-pow(CM.vv/c,2));
      double paraNew = (Para - et*CM.vv/c)*gamma;

     
      dvel ee;  //stores energy

      for (int j= 0;j<3;j++)
	{
	  ee.v[j] = perp[j] + paraNew/Para*para[j]; 
	}
     
      ee.vv = gamma*(et - CM.vv*Para/c);
      //ee.vv -=  931.478*(double)A;

      return ee;


}

dvel kin::transV(dvel CM,dvel start)
{
      //transform alpha1 to com frame
      double dot = 0.;
      for (int j=0;j<3;j++)
	{
	  dot += start.v[j]*CM.v[j];
	}



      double para[3];
      double perp[3];
      double Para = dot/CM.vv;
      double Perp = 0.;
      for (int j=0;j<3;j++)
	{
	  para[j] = Para/CM.vv*CM.v[j];
          perp[j] = start.v[j] - para[j];
          Perp += pow(perp[j],2);
	}
      Perp = sqrt(Perp);

      double gamma = 1./sqrt(1.-pow(CM.vv/c,2));

      dvel stop;
      stop.vv = 0.;
      for (int j= 0;j<3;j++)
	{

          perp[j] = perp[j]/(1.-Para*CM.vv/pow(c,2))/gamma;

	  para[j] = (para[j]-CM.v[j])/(1-Para*CM.vv/pow(c,2));

	  stop.v[j] = perp[j] + para[j];
          stop.vv += pow(stop.v[j],2);
	}
     
      stop.vv = sqrt(stop.vv);
      return stop;


}


float kin::getMinv(save save1, save save2)
{
    if (save1.id == save2.id) return -1.;

    vel CM = findCM(save1,save2);
    vel ee1 = trans(CM,save1);
    vel ee2 = trans(CM,save2);
    return ee1.vv + ee2.vv;
}
