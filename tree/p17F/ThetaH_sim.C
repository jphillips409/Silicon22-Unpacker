{
TCanvas canvas("ThetaH");

  //For comparing cos_thetaH
  TFile f("sort.root");
  TH2I * ThetaH_sort = (TH2I*) f.Get("ReCalc/Erel_costhetaH_re");
  TH2I * ThetaH_sorts1 = (TH2I*) f.Get("ReCalc/Erel_costhetaH_re_set1");
  TH2I * ThetaH_sorts1a = (TH2I*) f.Get("ReCalc/Erel_costhetaH_re_set1a");

  TFile fsim("/home/Silicon22/FittingSims/p17F/p17F_sim/rootout/sim_p17F_gradient.root");
  TH2I * ThetaH_sim = (TH2I*) fsim.Get("hist_Erel_thetaH");
  TH2I * ThetaH_sims1 = (TH2I*) fsim.Get("hist_Erel_thetaH_set1");
  TH2I * ThetaH_sims1a = (TH2I*) fsim.Get("hist_Erel_thetaH_set1a");

  TCanvas *mycan = (TCanvas*)gROOT->FindObject("mycan"); //Make canvas for fitting inv hists
  if(!mycan)
    {
      mycan = new TCanvas("mycan","",600,800);
      mycan->Divide(1,2,0.00001,0.00001);
    }
  mycan->SetTitle("Total");
  mycan->cd(1);
  ThetaH_sort->SetStats(kFALSE);
  ThetaH_sort->Draw("colz");
  
  mycan->cd(2);
  ThetaH_sim->SetStats(kFALSE);
  ThetaH_sim->Draw("colz");



  TCanvas *mycans1 = (TCanvas*)gROOT->FindObject("mycans1"); //Make canvas for fitting inv hists
  if(!mycans1)
    {
      mycans1 = new TCanvas("mycans1","",600,800);
      mycans1->Divide(1,2,0.00001,0.00001);
    }
  mycans1->SetTitle("set 1");
  mycans1->cd(1);
  ThetaH_sorts1->SetStats(kFALSE);
  ThetaH_sorts1->Draw("colz");
  
  mycans1->cd(2);
  ThetaH_sims1->SetStats(kFALSE);
  ThetaH_sims1->Draw("colz");



  TCanvas *mycans1a = (TCanvas*)gROOT->FindObject("mycans1a"); //Make canvas for fitting inv hists
  if(!mycans1a)
    {
      mycans1a = new TCanvas("mycans1a","",600,800);
      mycans1a->Divide(1,2,0.00001,0.00001);
    }
  mycans1a->SetTitle("set 1a");
  mycans1a->cd(1);
  ThetaH_sorts1a->SetStats(kFALSE);
  ThetaH_sorts1a->Draw("colz");
  
  mycans1a->cd(2);
  ThetaH_sims1a->SetStats(kFALSE);
  ThetaH_sims1a->Draw("colz");

}
