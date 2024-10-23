void Zcal()
{


  TFile *myfile = new TFile("~/root/read.root");
  if(!myfile->IsOpen())
    {
      cout << "cannot open the root file, try harder" << endl;
      return;
    }
  gROOT->cd();

  ofstream ofile("z.dat");
  if(!ofile.is_open())
    {
      cout <<"Could not open output txt file" << endl;
      return;
    }

  TCanvas * mycan = (TCanvas*)gROOT->FindObject("mycan");
  if(!mycan)
    {
      mycan = new TCanvas("mycan","mycan",1350,900);
    }
  
  mycan->SetCrosshair(1);
  int NCsI = 56;
  for(int icsi = 0;icsi<NCsI;icsi++)
    {
      TH1I *Light = (TH1I*)myfile->Get(Form("CsI/CsIGate/Light_%i",icsi));
      Light->Rebin(2);
      Light->Draw();

      Light->GetXaxis()->SetRangeUser(200,500);

      ofile <<icsi;
      TMarker * mark;

      for(int peak = 0;peak<4;peak++)
	{
      
	  mycan->Modified();
	  mycan->Update();
	  
	  mark=(TMarker*)mycan->WaitPrimitive("TMarker"); //Get the Background limits
	  int bkg_lo = mark->GetX();
	  delete mark;  
	  mark=(TMarker*)mycan->WaitPrimitive("TMarker");
	  int bkg_hi = mark->GetX();
	  delete mark;
	  mark=(TMarker*)mycan->WaitPrimitive("TMarker"); // Get the 1st peak initial guess
	  int peak1 = mark->GetX();
	  delete mark;
	  
	  
	  TF1 *myfit = new TF1("g1","gaus",bkg_lo,bkg_hi);
	  myfit->SetParameter(1,peak1);
	  Light->Fit(myfit,"R");

	  ofile << " " <<myfit->GetParameter(1);
	  
	}
      ofile << endl;
      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); // Get the 1st peak initial guess
      delete mark;      
      

    }
  ofile.close();
  
  return;
}
