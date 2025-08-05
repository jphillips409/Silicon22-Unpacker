void ceasarCal()
{
  TFile file("../../sort.root");

  TCanvas *mycan = new TCanvas("mycan","mycan");

  ofstream ofile("ceasar_Cs137.dat");
  TFile file2("ceasar_fits_Cs137.root","recreate");


  ostringstream outstring;
  string name;

  int p1= 50, p2=50; //+- fit limits up to 2 peaks. May be different.
  
  int const num_par = 5; //number of peaks times 2(pol1)+3(gaus).

  int ent = 0;


  // for (int ic=0;ic<13;ic++)
  for(int ic=0;ic<192;ic++)
    {
   
      outstring.str("");
      if (ic < 10)      outstring << "ceasar/raw/EC00" << ic;
      else if (ic < 100) outstring << "ceasar/raw/EC0" << ic;
      else outstring << "ceasar/raw/EC" << ic;
      name = outstring.str();
      TH1I * hist = (TH1I*)file.Get(name.c_str());
      hist->GetXaxis()->SetRangeUser(000,600);
      ent = hist->GetEntries();

      cout << ent << endl;
      hist->Draw();
      mycan->Modified();
      mycan->Update();
      if(ent <50)
	{
	  ofile << ic << " " << -1 << " "  << -1 << endl;
	  cout << ic << " Only has " << ent << " Entries. " << endl;
	  continue;
	}
      else
	{
        
	  TMarker * mark;
	  mark=(TMarker*)mycan->WaitPrimitive("TMarker"); //Get the Background limits
	  float bkg_lo = mark->GetX();
	  delete mark;  
	  mark=(TMarker*)mycan->WaitPrimitive("TMarker");
	  float bkg_hi = mark->GetX();
	  delete mark;
	  mark=(TMarker*)mycan->WaitPrimitive("TMarker"); // Get the 1st peak initial guess
	  float peak1 = mark->GetX();
	  delete mark;

	  double par[num_par] = {0.};
          double out[num_par] = {0.}; 
	  
	  int peak1_lo = peak1 - p1, peak1_hi = peak1 + p1; // Peak center and limits
	  // int peak2_lo = peak2 - p2, peak2_hi = peak2 + p2; // Peak center and limits
 	
	  TF1 *l1 = new TF1("l1", "pol1", bkg_lo, bkg_hi);
	  TF1 *g1 = new TF1("g1", "gaus", peak1_lo,peak1_hi);

	  TF1 *total = new TF1("total", "pol1(0)+gaus(2)", bkg_lo,bkg_hi);

	  hist->Fit(l1,"R");
	  hist->Fit(g1,"R+");

	  l1->GetParameters(&par[0]);
	  g1->GetParameters(&par[2]);
	  cout << "Peak1 = " << peak1_lo << " " << peak1_hi << endl;
	  total->SetParLimits(3,peak1_lo,peak1_hi);
	  total->SetParameters(par);
	  hist->Fit(total,"BR");
	  total->GetParameters(out);
	  cout << total << endl;



	  ofile << ic << " " << out[3]  << endl;

	  mycan->cd();
	  outstring.str("");
	  outstring << "Cs137_" << ic;
	  name = outstring.str();
	  total->SetName(name.c_str());
	  total->Draw("same L");

 	  file2.cd();
	  total->Write();

	  mycan->Modified();
	  mycan->Update();
	  
	  mark=(TMarker*)mycan->WaitPrimitive("TMarker"); //Click to see the next spectrum


        }
      
      
   }
  
  ofile.close();
  file.Close();
  return;
}
