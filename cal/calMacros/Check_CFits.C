void Check_CFits()
{

  TFile *file = new TFile("../root/na22.root");

  TFile *file2 = new TFile("ceasar_fits_Na22.root","UPDATE");
  string name;
  string name2;


  TH1I *hist = new TH1I();
  TF1 *fit1 = new TF1();
  TF1 *fit2 = new TF1();
  ostringstream outstring;
  ostringstream outstring2;

  TH1I *hist2 = new TH1I();

  gROOT->cd();
  int ent = 0;
  const int num_par = 5;
  int p1 = 15;





  TCanvas *mycan= new TCanvas("mycan","mycan");
  // mycan->Divide(10,20,0.000001,0.00001);
  for(int ic =0; ic<196;ic++)
    {

      //  mycan->cd(ic+1);
      outstring.str("");
      if (ic < 10)  outstring << "ceasar/raw/EC00" << ic;
      else if (ic < 100) outstring << "ceasar/raw/EC0" << ic;
      else outstring << "ceasar/raw/EC" << ic;
      name = outstring.str();
      outstring2.str("");
      outstring2 << "ceasar/raw/EC" << ic << "_temp";
      name2 = outstring2.str();
      hist = (TH1I*)file->Get(name.c_str());
      hist2 = (TH1I*)hist->Clone(name2.c_str());
      hist2->GetXaxis()->SetRangeUser(100,1000);
      hist2->Draw();
      ent = hist2->GetEntries();

      if(ent < 50) 
	{
	 cout << "Only found " << ent << " entries for" << ic << "." << endl;
	 continue;
	 
	}
      mycan->Modified();
      mycan->Update();
   

      outstring.str("");
      outstring << "Na22_" << ic;
      name = outstring.str();


      fit1 = (TF1*)file2->Get(name.c_str());

      if (!fit1) continue;

      fit1->Draw("SAME");
      fit1->Print();

      mycan->Modified();
      mycan->Update();

      bool IsGood = 0;
      bool IsFixed = 0;
      cout << "Is the fit good?" << endl;
      cin >> IsGood;

      if(!IsGood)
	{
	  for(;;)
	    {
	      if(IsFixed) break;
	      TMarker * mark;
	      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); //Get the Background limits
	      int bkg_lo = mark->GetX();
	      delete mark;  
	      mark=(TMarker*)mycan->WaitPrimitive("TMarker");
	      int bkg_hi = mark->GetX();
	      delete mark;
	      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); // Get the 1st peak 
	      int peak1 = mark->GetX();
	      delete mark;
	      
	      
	      double par[num_par] = {0.};
	      double out[num_par] = {0.}; 
	      
	      int peak1_lo = peak1 - p1, peak1_hi = peak1 + p1; // Peak center and limits
	      
	      
	      TF1 *l1 = new TF1("l1", "pol1", bkg_lo, bkg_hi);
	      TF1 *g1 = new TF1("g1", "gaus", peak1_lo,peak1_hi);
	      
	      TF1 *total = new TF1("total", "pol1(0)+gaus(2)", bkg_lo,bkg_hi);
	      
	      hist2->Fit(l1,"R");
	      hist2->Fit(g1,"R+");
	      
	      l1->GetParameters(&par[0]);
	      g1->GetParameters(&par[2]);
	      
	      total->SetParameters(par);
	      hist2->Fit(total,"R");
	      total->GetParameters(out);
	      hist2->Draw();
	      total->Draw("same");
	      gPad->Modified();
	      gPad->Update();

	      cout << "Did you fix it?" << endl;
	      cin >> IsFixed;
	      
	      if(IsFixed)
		{
		  
		  cout << ic << " " << out[3] << endl;
		  
		  outstring.str("");
		  outstring << "Na22_" << ic;
		  name = outstring.str();
		  total->SetName(name.c_str());
		  file2->cd();
		  total->Write();
		  gROOT->cd();
		}
	    }
	}

    }
  //  mycan->cd(0);
  //  mycan->Update();

  cout << "hey" << endl;

  file->Close();
  file2->Close();
  return;
}
