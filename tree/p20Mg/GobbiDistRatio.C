void GobbiDistRatio()
{
  gROOT->Reset();

  //Adds up the counts for setting 1 and setting 1a and finds the ratio of set1/tot for each peak
  //Need ratio for simulation
  //Setting 1a is closer than set1
  TFile *f;
  f = new TFile("sort.root");

  gROOT->cd();
  TCanvas *mycan = (TCanvas*)gROOT->FindObject("mycan"); //Make canvas for fitting inv hists
  if(!mycan)
    {
      mycan = new TCanvas("mycan","",1800,600);
      mycan->Divide(3,1,0.00001,0.00001);
    }

  TH1I * htot = (TH1I*) f->Get("ReCalc/Erel_re")->Clone(Form("Tot"));
  TH1I * hset1 = (TH1I*) f->Get("ReCalc/Erel_re_set1")->Clone(Form("Set 1"));
  TH1I * hset1a = (TH1I*) f->Get("ReCalc/Erel_re_set1a")->Clone(Form("Set 1a"));

  //Draw hists
  mycan->cd(1);
  htot->GetXaxis()->SetRangeUser(0,5);
  htot->Draw("same");

  mycan->cd(2);
  hset1->GetXaxis()->SetRangeUser(0,5);
  hset1->Draw("same");

  mycan->cd(3);
  hset1a->GetXaxis()->SetRangeUser(0,5);
  hset1a->Draw("same");

  //Make cuts for summing
  int nbins = htot->GetNbinsX(); //Same bins for all
  int numpeaks = 4;
  TAxis *xaxis = htot->GetXaxis(); //Same x-axis for all
  double Peak_xmin[numpeaks];
  double Peak_xmax[numpeaks];
  double Peak_int[2][numpeaks]; //Set 1 and Set 1a

  //Assign peak limits
  Peak_xmin[0] = 0.7;
  Peak_xmin[1] = 1.1;
  Peak_xmin[2] = 1.7;
  Peak_xmin[3] = 3.5;

  Peak_xmax[0] = 1;
  Peak_xmax[1] = 1.6;
  Peak_xmax[2] = 2.1;
  Peak_xmax[3] = 4;

  for (int n=0;n<2;n++)
  {
    for (int i=0;i<numpeaks;i++)
    {
      Peak_int[n][i] = 0; //reset integral

      for (int j=0;j<nbins;j++)
      {
        if (xaxis->GetBinLowEdge(j) >= Peak_xmin[i] && xaxis->GetBinLowEdge(j) <= Peak_xmax[i])
        {
          cout << hset1->GetBinContent(j) << endl;
          if (n == 0) Peak_int[n][i] += hset1->GetBinContent(j);
          if (n == 1) Peak_int[n][i] += hset1a->GetBinContent(j);
        }
      }
    }
  }

  //Calc ratios
  double Peak_ratio[numpeaks]; //Do ratio of set1/tot
  for (int i=0;i<numpeaks;i++)
  {
    Peak_ratio[i] = 0;
    
    Peak_ratio[i] = Peak_int[0][i]/(Peak_int[0][i] + Peak_int[1][i]);
    cout << "Peak set1/tot for peak " << i << " " << Peak_ratio[i] << endl;
  }

  f->Close();
}
