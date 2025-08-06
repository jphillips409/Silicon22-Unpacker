void gammafit()
{
  //For fitting the gamma spectrum
  gROOT->Reset();
  TStyle * Sty = new TStyle("MyStyle","MyStyle");
  Sty->SetOptTitle(0);
  Sty->SetOptStat(0);
  Sty->SetLineWidth(4);
  Sty->SetPalette(55);
  Sty->SetCanvasColor(10);
  Sty->SetCanvasBorderMode(0);
  Sty->SetFrameLineWidth(5);
  Sty->SetFrameFillColor(10);
  Sty->SetPadColor(10);
  Sty->SetPadTickX(1);
  Sty->SetPadTickY(1);
  Sty->SetPadBottomMargin(.2);
  Sty->SetPadTopMargin(.05);
  Sty->SetPadLeftMargin(.15);
  Sty->SetPadRightMargin(.03);
  Sty->SetHistLineWidth(3);
  Sty->SetFuncWidth(3);
  Sty->SetFuncColor(kRed);
  Sty->SetLineWidth(3);
  Sty->SetLabelSize(0.06,"xyz");
  Sty->SetLabelOffset(0.01,"y");
  Sty->SetLabelOffset(0.02,"x");
  Sty->SetLabelColor(kBlack,"xyz");
  Sty->SetTitleSize(0.07,"xyz");
  Sty->SetTitleOffset(1.1,"y");
  Sty->SetTitleOffset(1.2,"x");
  Sty->SetTitleFillColor(10);
  Sty->SetTitleTextColor(kBlack);
  Sty->SetTickLength(.05,"xz");
  Sty->SetTickLength(.025,"y");
  Sty->SetNdivisions(410,"xyz");
  Sty->SetEndErrorSize(0);
  gROOT->SetStyle("MyStyle");
  gROOT->ForceStyle();


  TFile *f;
  f = new TFile("sort.root");

  gROOT->cd();
  TCanvas *mycan = (TCanvas*)gROOT->FindObject("mycan"); //Make canvas for fitting inv hists
  if(!mycan)
    {
      mycan = new TCanvas("mycan","",800,600);
    }

  TH1I * gamma = (TH1I*) f->Get("SortCode/p20Mg_gammasADD")->Clone(Form("gamma"));

  mycan->cd();

  gamma->Rebin(12);

  gamma->GetXaxis()->SetRangeUser(0,2500);
  gamma->GetXaxis()->SetTitle("E_{#gamma} (keV)");
  gamma->GetXaxis()->CenterTitle();
  gamma->GetYaxis()->SetTitle("Counts/48 keV");
  gamma->GetYaxis()->CenterTitle();
  gamma->GetYaxis()->SetRangeUser(0,350);
  gamma->SetMarkerStyle(20);
  gamma->SetMarkerColor(1);
  gamma->SetLineColor(1);
  gamma->SetMarkerSize(2);
  gamma->Draw("same P");

  //Array of parameters
  double par[6];
  TF1 *pol = new TF1("pol", "pol2",200, 2400);
  TF1 *g = new TF1("g", "gaus", 1000, 2200);

  TF1 *total = new TF1("total", "pol2(0)+gaus(3)", 200, 2400);

  gamma->Fit(pol, "0R+");
  gamma->Fit(g, "0R+");

  pol->GetParameters(&par[0]);
  g->GetParameters(&par[3]);

  total->SetParameters(par);
  gamma->Fit(total, "0R+");

  total->GetParameters(&par[0]);

  pol->SetParameters(&par[0]);
  g->SetParameters(&par[3]);

  total->SetLineColor(kRed);
  pol->SetLineColor(kGreen+2);
  g->SetLineColor(kBlue);

  pol->Draw("Csame");
  g->Draw("Csame");
  total->Draw("Csame");

  TLatex tt;
  tt.SetTextSize(.07);
  tt.DrawLatex(2150, 300,"(d)");

  //Grab the errors from each point
  int nbins = gamma->GetNbinsX(); //Same bins for all
  TAxis *xaxis = gamma->GetXaxis(); //Same x-axis for all
  double binpos[nbins];
  double vals[nbins];
  double err[nbins];
  double xerr[nbins];


   for (int j=0;j<nbins;j++)
   {
     binpos[j] = xaxis->GetBinCenter(j);
     vals[j] = gamma->GetBinContent(j);
     err[j] = sqrt(vals[j]);
     xerr[j] = 0;
   }



  //Integrate functions for across Mg-20 peak
  //Mg-20 peak bounds
  double lowerb = 1450.;
  double upperb = 1700.;
  double binw = gamma->GetXaxis()->GetBinWidth(1);

  double bkg_int = pol->Integral(lowerb,upperb)/binw;
  double peak_int = g->Integral(lowerb,upperb)/binw;
  double spect_int = total->Integral(200, 2400)/binw;

  double peak_tot = bkg_int + peak_int;
  double bkg_ratio = bkg_int/peak_tot;

  cout << "Background integral for peak: " << bkg_int << endl;
  cout << "Peak integral: " << peak_int << endl;
  cout << "Total counts in peak region: " << peak_tot << endl;
  cout << "Bkg/Total: " << bkg_ratio << endl;
  cout << "Total spectrum integral " << spect_int << endl;

  //Plot errors
  TGraphErrors errors(nbins,binpos,vals,xerr,err);
  errors.SetMarkerStyle(20);
  errors.SetMarkerSize(1.2);
  errors.SetLineWidth(2);
  errors.Draw("same EP");

  mycan->Print("p20Mg_gammafit_d.png", "png");
  mycan->Print("p20Mg_gammafit_d.eps", "eps");

  f->Close();

}
