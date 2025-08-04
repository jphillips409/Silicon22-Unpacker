void CsICal()
{
  gROOT->Reset();
  TFile *file = new TFile("../root/NZ_75_Old.root");
  TFile *file2 = new TFile("CsIfits/CsI_75A_old.root","Update");
  TFile *Gates = new TFile("../gates/zlines.root");
  
  
  ofstream ofile("CsI_75A_Old.dat");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1);
  
  ostringstream outstring;
  string name;
  
  int p1= 30, p2=50; //+- fit limits up to 2 peaks. May be different.
  int const num_par = 5; //number of peaks times 2(pol1)+3(gaus).
  
  int ent = 0;
  
  gROOT->cd();
  TCanvas *mycan =(TCanvas*)gROOT->FindObjectAny("mycan");
  if(!mycan)
    {
      mycan = new TCanvas("mycan","mycan");
      mycan->Divide(1,2);
    }
  
  for(int ic=0;ic<1;ic++)
    {
      cout << ic << endl;
      outstring.str("");
      outstring << "dEE/dEE_" << ic;
      name = outstring.str();
      TH2I * hist = (TH2I*)file->Get(name.c_str())->Clone(Form("dee_%i",ic));
      outstring.str("");
      outstring << "Zline_" << ic << "_2_4";
      name = outstring.str();
      TCutG *mycut = (TCutG*)Gates->Get(name.c_str());
      if(!mycut)
      	{
      	  cout << "Can't find the cut for " << ic << endl;
      	  cout << "continuing search for " << ic+1 << endl;
      	  continue;
      	}
      
      mycan->cd(1);
      gPad->SetLogz();
      hist->Draw("zcol");

      outstring.str("");
      outstring << "[Zline_" << ic << "_2_4]";
      name = outstring.str();

      TH1D *projx = (TH1D*)hist->ProjectionX("projx",1,500,name.c_str());
      hist->GetXaxis()->SetRangeUser(0,2000);
      hist->GetYaxis()->SetRangeUser(0,100);
      mycut->Draw();
      mycan->cd(2);
      projx->Draw();
      projx->Rebin(2);
      projx->GetXaxis()->SetRangeUser(0,2000);

      
      
      
      mycan->Modified();
      mycan->Update();
      cout << "tele = "<< ic << endl;
      //      gPad->SetLogy();
      
      
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
      
      projx->Fit(l1,"R");
      projx->Fit(g1,"R+");
      
      l1->GetParameters(&par[0]);
      g1->GetParameters(&par[2]);
      
      total->SetParameters(par);
      projx->Fit(total,"R");
      total->GetParameters(out);
      
      
      outstring.str("");
      outstring << "55A_" << ic;
      name = outstring.str();
      total->SetName(name.c_str());
      total->Draw("same");
      mycan->Modified();
      mycan->Update();
      
       //   total->Write();
      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); // Get the 1st peak initial guess
      delete mark;      
      

    }
  
  cout << "hey" << endl;
  ofile.close();
  file2->Close();
  file->Close();
  Gates->Close();
  
  return;
}
