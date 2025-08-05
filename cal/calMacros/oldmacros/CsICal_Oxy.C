void CsICal_Oxy()
{
  TFile *file = new TFile("../root/Oxygen55_new.root");
  TFile *file2 = new TFile("CsIfits/CsI_55O_new.root","Update");
   
  ofstream ofile("CsI_55O_New.dat");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  ostringstream outstring;
  string name;
  
  int p1= 80, p2=50; //+- fit limits up to 2 peaks. May be different.
  int const num_par = 5; //number of peaks times 2(pol1)+3(gaus).

  int ent = 0;
  
  gROOT->cd();
  TCanvas *mycan =(TCanvas*)gROOT->FindObjectAny("mycan");
  if(!mycan)
    {
      mycan = new TCanvas("mycan","mycan");
    }
  
  for(int ic=0;ic<56;ic++)
    {
      cout << ic << endl;
      outstring.str("");
      outstring << "CsI/CsISum/CsIOver_" << ic;
      name = outstring.str();
      TH1I * hist = (TH1I*)file->Get(name.c_str());
      
      hist->Draw("");
      hist->GetXaxis()->SetRangeUser(500,2000);
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
      
      hist->Fit(l1,"R");
      hist->Fit(g1,"R+");
      
      l1->GetParameters(&par[0]);
      g1->GetParameters(&par[2]);
      
      total->SetParameters(par);
      hist->Fit(total,"R");
      total->GetParameters(out);
      
      
      ofile << ic << " " << out[3] << endl;
      
      outstring.str("");
      outstring << "55O_" << ic;
      name = outstring.str();
      total->SetName(name.c_str());
      total->Draw("same");
      mycan->Modified();
      mycan->Update();
      
      //   total->Write();
      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); // Get the 1st peak initial guess
      delete mark;      
      
      
   }
  
  ofile.close();
  file2->Close();
  file->Close();
  Gates->Close();
  
  return;
}
