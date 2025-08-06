void gammasubtract_1600()
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

  TH1I * Nogamma = (TH1I*) f->Get("SortCode/Erel_21Al_p20Mg")->Clone(Form("Nogamma"));
  TH1I * with = (TH1I*) f->Get("SortCode/Erel_21Al_p20Mg_with1600")->Clone(Form("with1600"));

  mycan->cd(1);
  TH1I * Nogamma_cop = (TH1I*)Nogamma->DrawCopy();
  Nogamma_cop->GetXaxis()->SetRangeUser(0,5);
  Nogamma_cop->SetLineColor(kRed);
  Nogamma_cop->GetXaxis()->SetTitle("Erel (MeV)");
  Nogamma_cop->GetXaxis()->CenterTitle();
  Nogamma_cop->GetYaxis()->SetTitle("Counts/40 keV");
  Nogamma_cop->GetYaxis()->CenterTitle();
  Nogamma_cop->Draw("same hist");

  mycan->cd(2);
  TH1I * with_cop = (TH1I*)with->DrawCopy();
  with_cop->GetYaxis()->SetRangeUser(0,40);
  with_cop->GetXaxis()->SetRangeUser(0,5);
  with_cop->GetXaxis()->SetTitle("Erel (MeV)");
  with_cop->GetXaxis()->CenterTitle();
  with_cop->GetYaxis()->SetTitle("Counts/40 keV");
  with_cop->GetYaxis()->CenterTitle();
  with_cop->Draw("same hist");

  //Find the number of counts inside of first peak for both spectra
  int nbins = Nogamma->GetNbinsX(); //Same bins for all
  TAxis *xaxis = Nogamma->GetXaxis(); //Same x-axis for all

  //Get total counts in each spectra;
  //Reduce ungated spectrum to 0.275097 of the gated and subtract
  //Find the scaling factor needed
  double ratio = 0.275097;
  double nogamma_counts = Nogamma_cop->Integral();
  double with_counts = with_cop->Integral();
  double scale_factor = (with_counts/nogamma_counts)*ratio;

  cout << "No gamma scaling factor = " << scale_factor << endl;;

  //Nogamma->Scale(0.023);
  Nogamma->Scale(scale_factor);

  cout << "scale factor " << scale_factor << endl;

  mycan->cd(3);
  Nogamma->GetXaxis()->SetRangeUser(0,5);
  Nogamma->GetYaxis()->SetRangeUser(0,40);
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
  with->GetYaxis()->SetRangeUser(0,40);
  with->GetXaxis()->SetRangeUser(0,5);
  with->GetXaxis()->SetTitle("Erel (MeV)");
  with->GetXaxis()->CenterTitle();
  with->GetYaxis()->SetTitle("Counts/40 keV");
  with->GetYaxis()->CenterTitle();
  with->Draw("same hist");
  
  mycan->Print("Erel_21Al_p20Mg_gsubtract1600.png", "png");

  f->Close();


  //Output the subtracted spectrum to a root file
  //Adds simulated peaks together at some specified ratio
  ostringstream filename;
  filename.str("");
  filename << "subtracted.root";

  TFile * file_read = new TFile(filename.str().c_str(),"RECREATE");


  //Write the subtracted spectrum
  with->GetXaxis()->SetRangeUser(-5,15);
  with->SetLineColor(kBlue+2);
  with->Sumw2(0);
  with->Write();


  file_read->Write();
  cout << "histo written" << endl;
  file_read->Close();

}
