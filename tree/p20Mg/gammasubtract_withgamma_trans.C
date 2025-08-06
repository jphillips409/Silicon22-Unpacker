void gammasubtract_withgamma_trans()
{
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
  Sty->SetFuncColor(kGreen);
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
      mycan = new TCanvas("mycan","",1200,800);
      mycan->Divide(2,2,0.00001,0.00001);
    }

  TH1I * Nogamma = (TH1I*) f->Get("SortCode/Erel_21Al_p20Mg_trans")->Clone(Form("Nogamma"));
  TH1I * with = (TH1I*) f->Get("SortCode/Erel_21Al_p20Mg_withgamma_trans")->Clone(Form("Erel_21Al_p20Mg_withgamma"));

  mycan->cd(1);
  TH1I * Nogamma_cop = (TH1I*)Nogamma->DrawCopy();
  Nogamma_cop->GetXaxis()->SetRangeUser(0,7);
  Nogamma_cop->SetLineColor(kRed);
  Nogamma_cop->GetXaxis()->SetTitle("Erel (MeV)");
  Nogamma_cop->GetXaxis()->CenterTitle();
  Nogamma_cop->GetYaxis()->SetTitle("Counts/40 keV");
  Nogamma_cop->GetYaxis()->CenterTitle();
  Nogamma_cop->Draw("same hist");

  mycan->cd(2);
  TH1I * with_cop = (TH1I*)with->DrawCopy();
  with_cop->GetYaxis()->SetRangeUser(0,200);
  with_cop->GetXaxis()->SetRangeUser(0,7);
  with_cop->GetXaxis()->SetTitle("Erel (MeV)");
  with_cop->GetXaxis()->CenterTitle();
  with_cop->GetYaxis()->SetTitle("Counts/40 keV");
  with_cop->GetYaxis()->CenterTitle();
  with_cop->Draw("same hist");

  //Find the number of counts inside of first peak for both spectra
  int nbins = Nogamma->GetNbinsX(); //Same bins for all
  TAxis *xaxis = Nogamma->GetXaxis(); //Same x-axis for all

  //Integrate first peak for both spectra
  double peak1_nogamma = 0;
  double peak1_with = 0;
  
  double xmin = 0.24;
  double xmax = 0.56;

  for (int i=0;i<nbins;i++)
  {
    if (xaxis->GetBinLowEdge(i) >= xmin && xaxis->GetBinLowEdge(i) <= xmax)
    {
      peak1_nogamma += Nogamma->GetBinContent(i);
      peak1_with += with->GetBinContent(i);
    }
  }

  cout << "Peak 1: No gamma = " << peak1_nogamma << "  with = " << peak1_with << endl;

  //Scale no gamma down
  double scale_nogamma = peak1_with / peak1_nogamma;

  Nogamma->Scale(0.32);

  cout << "scale factor " << scale_nogamma << endl;

  mycan->cd(3);
  Nogamma->GetXaxis()->SetRangeUser(0,7);
  Nogamma->GetYaxis()->SetRangeUser(0,200);
  Nogamma->SetLineColor(kRed);
  Nogamma->GetXaxis()->SetTitle("Erel (MeV)");
  Nogamma->GetXaxis()->CenterTitle();
  Nogamma->GetYaxis()->SetTitle("Counts/40 keV");
  Nogamma->GetYaxis()->CenterTitle();
  Nogamma->Draw("same hist");
  with_cop->Draw("same hist");

  //Subtract nogamma from with1650
  with->Add(Nogamma,-1);
  mycan->cd(4);
  with->GetYaxis()->SetRangeUser(0,30);
  with->GetXaxis()->SetRangeUser(0,7);
  with->GetXaxis()->SetTitle("Erel (MeV)");
  with->GetXaxis()->CenterTitle();
  with->GetYaxis()->SetTitle("Counts/40 keV");
  with->GetYaxis()->CenterTitle();
  with->Draw("same hist");
  
  mycan->Print("Erel_21Al_p20Mg_gsubtract1600.png", "png");

  f->Close();

}
