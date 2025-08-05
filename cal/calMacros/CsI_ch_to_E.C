void CsI_ch_to_E()
{

  //tele num, crystal, target
  float CsIChannels[4][4][3] = {{0.}};
  // Be, Thin Al, Thick Al (with 5mm Al blocker and 1.5mm Si)
  //float CsIEnergies[3] = {103.591,85.556,60.952}; //MeV
  float CsIEnergies[4][3] = {{102.314,84.199,59.192},
                             {102.314,84.199,59.192},
                             {102.314,84.199,59.192},
                             {102.314,84.199,59.192}}; //MeV takes angle into account for each CsI

  ofstream cal_out("CsI_cal_out.dat");

  ifstream inCsIBe("CsI_p_Be.dat");
  if(!inCsIBe.is_open())
    {
      cout << "could not open Be csi channel text file" << endl;
      return;
    }
  
  ifstream inCsIThinAl("CsI_p_ThinAl.dat");
  if(!inCsIBe.is_open())
    {
      cout << "could not open Thin Al csi channel text file" << endl;
      return;
    }
  
  ifstream inCsIThickAl("CsI_p_ThickAl.dat");
  if(!inCsIBe.is_open())
    {
      cout << "could not open ThickAl csi channel text file" << endl;
      return;
    }
  

  TCanvas *mycan = (TCanvas*)gROOT->FindObject("mycan");
  if(!mycan)
    {
      mycan = new TCanvas("mycan","",1200,800);
      mycan->Divide(4,4,0.00001,0.00001);
    }

  int quad, crystal;
  float channel;
  float slope[4][4]={{0.}};
  float inter[4][4]={{0.}};
  for(int i=0;i<4;i++)
    {
      for(int j=0;j<4;j++)
	{
	  inCsIBe >> quad >> crystal >> channel;
	  CsIChannels[i][j][0] = channel;
	  inCsIThinAl >> quad >> crystal >> channel;
	  CsIChannels[i][j][1] = channel;
	  inCsIThickAl >> quad >> crystal >> channel;
	  CsIChannels[i][j][2] = channel;
	  mycan->cd(i*4+j+1);
	  //TGraph *mygraph = new TGraph(3,CsIChannels[i][j],CsIEnergies);
	  TGraph *mygraph = new TGraph(3,CsIChannels[i][j],CsIEnergies[crystal]);
	  mygraph->SetMarkerStyle(20);
	  mygraph->Draw("AP");


	  TF1 *myfit = new TF1("myfit","pol1",300,1300);
	  mygraph->Fit(myfit,"Q");
	  inter[i][j] = myfit->GetParameter(0);
	  slope[i][j] = myfit->GetParameter(1);
	    
	}
    }

    for(int i=0;i<4;i++)
    {
      for(int j=0;j<4;j++)
	{
	  cout << i << " "<< j <<  " " <<slope[i][j] << " " << inter[i][j] << endl;
	 
	}
    }

  inCsIBe.close();
  inCsIThinAl.close();
  inCsIThickAl.close();

  return;
}
