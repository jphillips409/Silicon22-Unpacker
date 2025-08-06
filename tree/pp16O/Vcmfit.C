{
TCanvas canvas("Vcm");

  //For fitting the 19Na -> p + 18Ne Vcm
  TFile f("sort.root");
  TH1I * Vcm_B = (TH1I*) f.Get("VCM");

  //Array of parameters
  double par[9];

  //par[1] = 4.16;
  //par[4] = 5.15;

  TF1 *g1 = new TF1("g1", "gaus",11.2, 11.6);
  TF1 *g2 = new TF1("g2", "gaus", 11.6, 12.1);

  TF1 *total = new TF1("total", "gaus(0)+gaus(3)", 11.2, 12.1);


 /* TH2I frame("frame","",1000,0,250,1000,0,100);
  frame.GetXaxis()->SetTitle("E");
  frame.GetYaxis()->SetTitle("dE");
  frame.GetYaxis()->CenterTitle();
  frame.GetXaxis()->CenterTitle();
  frame.Draw();*/

   Vcm_B->Draw();
   Vcm_B->Fit(g1, "0R+");
   Vcm_B->Fit(g2, "0R+");

   g1->GetParameters(&par[0]);
   g2->GetParameters(&par[3]);

   total->SetParameters(par);
   Vcm_B->Fit(total, "0R+");

   total->GetParameters(&par[0]);

   g1->SetParameters(&par[0]);
   g2->SetParameters(&par[3]);


   total->Draw("same");
   g1->Draw("same");
   g2->Draw("same");
    
   double inttot = (par[0]*fabs(par[2]) + par[3]*fabs(par[5]))*sqrt(2*M_PI);
   cout << "Total Integral " << inttot << endl;
   cout << "amp ratio 1 " << par[0]*fabs(par[2])*sqrt(2*M_PI)/inttot << endl;
   cout << "amp ratio 2 " << par[3]*fabs(par[5])*sqrt(2*M_PI)/inttot << endl;

}
