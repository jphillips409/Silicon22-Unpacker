void gammasubtract_500()
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

  TH1I * Nogamma = (TH1I*) f->Get("Ex_hist")->Clone(Form("Nogamma"));
  TH1I * with = (TH1I*) f->Get("Ex_with500")->Clone(Form("with"));

  mycan->cd(1);
  TH1I * Nogamma_cop = (TH1I*)Nogamma->DrawCopy();
  Nogamma_cop->GetXaxis()->SetRangeUser(4,7);
  Nogamma_cop->SetLineColor(kRed);
  Nogamma_cop->GetXaxis()->SetTitle("Ex (MeV)");
  Nogamma_cop->GetXaxis()->CenterTitle();
  Nogamma_cop->GetYaxis()->SetTitle("Counts/40 keV");
  Nogamma_cop->GetYaxis()->CenterTitle();
  Nogamma_cop->Draw("same hist");

  mycan->cd(2);
  TH1I * with_cop = (TH1I*)with->DrawCopy();
  with_cop->GetYaxis()->SetRangeUser(0,20);
  with_cop->GetXaxis()->SetRangeUser(4,7);
  with_cop->GetXaxis()->SetTitle("Ex (MeV)");
  with_cop->GetXaxis()->CenterTitle();
  with_cop->GetYaxis()->SetTitle("Counts/40 keV");
  with_cop->GetYaxis()->CenterTitle();
  with_cop->Draw("same hist");


  //Scale no gamma down
  double scale_nogamma = 0.015;

  Nogamma->Scale(scale_nogamma);

  cout << "scale factor " << scale_nogamma << endl;

  mycan->cd(3);
  Nogamma->GetXaxis()->SetRangeUser(4,7);
  Nogamma->GetYaxis()->SetRangeUser(0,20);
  Nogamma->SetLineColor(kRed);
  Nogamma->GetXaxis()->SetTitle("Ex (MeV)");
  Nogamma->GetXaxis()->CenterTitle();
  Nogamma->GetYaxis()->SetTitle("Counts/40 keV");
  Nogamma->GetYaxis()->CenterTitle();
  Nogamma->Draw("same hist");
  with_cop->Draw("same hist");

  //Subtract nogamma from with1650
  with->Add(Nogamma,-1);
  mycan->cd(4);
  with->GetYaxis()->SetRangeUser(0,20);
  with->GetXaxis()->SetRangeUser(4,7);
  with->GetXaxis()->SetTitle("Ex (MeV)");
  with->GetXaxis()->CenterTitle();
  with->GetYaxis()->SetTitle("Counts/40 keV");
  with->GetYaxis()->CenterTitle();
  with->Draw("same hist");
  
  mycan->Print("Erel_18Ne_p17F_gsubtract500.png", "png");

  f->Close();

}
