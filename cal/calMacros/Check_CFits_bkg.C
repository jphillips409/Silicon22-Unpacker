void Check_CFits_bkg()
{

  TFile *file = new TFile("../bkg.root");

  TFile *file2 = new TFile("ceasar_fits.root_bkg");
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
      hist2->GetXaxis()->SetRangeUser(150,1000);
      hist2->Draw();
      gPad->Modified();
      gPad->Update();
      ent = hist2->GetEntries();

      if(ent < 50) 
	{
	 cout << "Only found " << ent << " entries for" << ic << "." << endl;
	 continue;
	 
	}
   

      outstring.str("");
      outstring << "K40_" << ic;
      name = outstring.str();


      fit1 = (TF1*)file2->Get(name.c_str());
      outstring.str("");
      outstring << "Th228_" << ic;
      name = outstring.str();
      fit2 =  (TF1*)file2->Get(name.c_str());

      if (!fit1 && !fit2) continue;

      fit1->Draw("SAME");
      fit1->Print();
      fit2->Draw("SAME");

      gPad->Modified();
      gPad->Update();

      TMarker * mark;
      mark=(TMarker*)mycan->WaitPrimitive("TMarker"); //Get the Background limits
      delete mark;


    }
  //  mycan->cd(0);
  //  mycan->Update();

  cout << "hey" << endl;

  file->Close();
  file2->Close();
  return;
}
