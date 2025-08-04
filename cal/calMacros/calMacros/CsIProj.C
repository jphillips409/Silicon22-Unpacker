void CsIProj()
{
  TFile *file = new TFile("../root/NZ_55_New.root");
  TFile *gates = new TFile("../gates/zlines.root");
  TFile *gates2 = new TFile("../gates/zlines_new.root");

  ofstream ofile("CsI_55A_New.dat");
  
  TCanvas *mycan = (TCanvas*)gROOT->FindObjectAny("mycan");
  if(!mycan)
    {
      mycan = new TCanvas("mycan","mycan");
      mycan->Divide(1,2);
    }
  
  
  ostringstream outstring;
  string name;
  int p1= 30, p2=50; //+- fit limits up to 2 peaks. May be different.
  int const num_par = 5; //number of peaks times 2(pol1)+3(gaus).
  
  for(int ic =0;ic<56;ic++)
    {
      
      outstring.str("");
      outstring << "dEE/dEE_" << ic;
      name = outstring.str();
      
      mycan->cd(1);
      TH2I *hist = (TH2I*)file->Get(name.c_str());
      hist->Draw("col");  
      hist->GetXaxis()->SetRangeUser(200.,1800.); 
      hist->GetYaxis()->SetRangeUser(5.,50.); 
   
      
      if(ic <16 || ic > 31)
	TCutG *mycut = (TCutG*)gates->Get(Form("Zline_%i_2_4",ic));
      else 
	TCutG *mycut = (TCutG*)gates2->Get(Form("Zline_%i_2_4",ic));
      mycut->Draw();

      file->cd();
      outstring.str("");
      outstring << "CsI/CsIGate/ECsI_" << ic << "_Gate";
      name = outstring.str();
      gPad->SetLogz();

      mycan->cd(2);
      TH1I * proj = (TH1I*)file->Get(name.c_str());
      proj->Draw();
      proj->Rebin(4);
      proj->GetXaxis()->SetRangeUser(700.,1800.);

      mycan->Modified();
      mycan->Update();

      TMarker * mark;
      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); //Get the Background limits
      int bkg_lo = mark->GetX();
      delete mark;  
      mark=(TMarker*)mycan->WaitPrimitive("TMarker");
      int bkg_hi = mark->GetX();
      delete mark;
      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); // Get the 1st peak initial guess
      int peak1 = mark->GetX();
      delete mark;
      
      
      double par[num_par] = {0.};
      double out[num_par] = {0.}; 
      int peak1_lo = peak1 - p1, peak1_hi = peak1 + p1; // Peak center and limits
      
      
      TF1 *l1 = new TF1("l1", "pol1", bkg_lo, bkg_hi);
      TF1 *g1 = new TF1("g1", "gaus", peak1_lo,peak1_hi);
      
      TF1 *total = new TF1("total", "pol1(0)+gaus(2)", bkg_lo,bkg_hi);
      
      proj->Fit(l1,"R");
      proj->Fit(g1,"R+");
      
      l1->GetParameters(&par[0]);
      g1->GetParameters(&par[2]);
      
      total->SetParameters(par);
      proj->Fit(total,"R");
      total->GetParameters(out);
      
      
      ofile << ic << " " << out[3] << endl;
      
      outstring.str("");
      outstring << "55A_" << ic;
      name = outstring.str();
      total->SetName(name.c_str());
      total->Draw("same");
      mycan->Modified();
      mycan->Update();
      
      bool IsGood = 0;

      cout << "Good fit?" << endl;
      cin >> IsGood;
  

      if(IsGood)
	{
      ofile << ic << " " << out[3] << endl;
	}      
      else
	ofile << ic << " " << -1 << endl;      

   

    }
  
  
  return;
}
