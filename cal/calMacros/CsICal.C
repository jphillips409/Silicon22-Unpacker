void CsICal()
{
  gROOT->Reset();
//  TFile *file = new TFile("../../Proton_Be.root");
  TFile *file = new TFile("../../Proton_Al.root");
  
  
//  ofstream ofile("CsI_p_Be.dat");
  ofstream ofile("CsI_p_Al.dat");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1);
  
  ostringstream outstring;
  string name;
  
  int p1= 10, p2=50; //+- fit limits up to 2 peaks. May be different.
  int const num_par = 3; //number of peaks times 2(pol1)+3(gaus).
  
  int ent = 0;
  
  gROOT->cd();
  TCanvas *mycan =(TCanvas*)gROOT->FindObjectAny("mycan");
  if(!mycan)
    {
      mycan = new TCanvas("mycan","mycan");
    }
  
  for(int ic=0;ic<20;ic++)
    {
      cout << ic << endl;
      outstring.str("");
      outstring << "CsI/CsIRaw/ECsI_" << ic;
      name = outstring.str();
      TH1I * hist = (TH1I*)file->Get(name.c_str())->Clone(Form("ECsIRaw_%i",ic));
      hist->Draw("");

      outstring.str("");
      outstring << "[Zline_" << ic << "_2_4]";
      name = outstring.str();

//      hist->Rebin(2);
      hist->GetXaxis()->SetRangeUser(0,400);

      
      
      
      mycan->Modified();
      mycan->Update();
      cout << "tele = "<< ic << endl;
      //      gPad->SetLogy();
      
      
      TMarker * mark;
      // mark=(TMarker*)mycan->WaitPrimitive("TMarker"); //Get the Background limits
      // int bkg_lo = mark->GetX();
      // delete mark;  
      // mark=(TMarker*)mycan->WaitPrimitive("TMarker");
      // int bkg_hi = mark->GetX();
      // delete mark;
      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); // Get the 1st peak initial guess
      int peak1 = mark->GetX();
      delete mark;
      
      
      double par[num_par] = {0.};
      double out[num_par] = {0.}; 
      int peak1_lo = peak1 - p1, peak1_hi = peak1 + p1; // Peak center and limits
      
      
      // TF1 *l1 = new TF1("l1", "pol1", bkg_lo, bkg_hi);
      // TF1 *g1 = new TF1("g1", "gaus", peak1_lo,peak1_hi);
      
      // TF1 *total = new TF1("total", "pol1(0)+gaus(2)", bkg_lo,bkg_hi);
      TF1 *total = new TF1("total", "gaus(0)", peak1_lo,peak1_hi);
      
      // hist->Fit(l1,"R");
      // hist->Fit(g1,"R+");
      
      // l1->GetParameters(&par[0]);
      // g1->GetParameters(&par[2]);
      
      total->SetParameters(1,peak1);
      hist->Fit(total,"R+");
      total->GetParameters(out);
      
      
       //   total->Write();
      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); // Get the 1st peak initial guess
      delete mark;      
      
      ofile << ic << " " << out[1] << endl;
      cout << ic << " " << out[1] << endl;
    }
  
  cout << "hey" << endl;
  ofile.close();
  file->Close();
  
  return;
}
