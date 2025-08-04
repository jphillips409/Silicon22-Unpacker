{  // Script: c:\root\macros\draw1_fitN.c 
  gROOT -> Reset();
  TCanvas *fits= new TCanvas("Fits", "Energy Histograms", 1);
  TH1F*spectrum;

  // FIT STUFF: Define the Fit Parameters Here:
  string spec;
  cout << "Which spectrum ?  "; cin >>spec;
  
  TString hist(spec);
  spectrum=(TH1F*)gROOT->FindObject(hist);
  
  //spectrum   =  E1;  //string spec= "E1";
  int Ecal;
  float comp;// = 1.0;  // ROOT's compression of the data  = channel-dim/hist-chan 
  int p1= 20, p2=20, p3=20; //+- fit limits up to 3 peaks. May be different.
  int num_peaks;
  int peak1=50, peak2=75, peak3=100;
   
  cout << "How many peaks  = "; cin >>num_peaks;
  cout << "Root compression= "; cin >> comp;
  cout << "Ecal=0, no cal, 1 or 2 = "; cin >> Ecal ;

  int num_par  = num_peaks*3 + 2; // Upto 3 peaks Maximum
   
  //added by Chris----------
  spectrum->Draw();
 
  mark=(TMarker*)gPad->WaitPrimitive("TMarker","Marker"); //Get the Background limits
  int bkg_lo = mark->GetX();
  delete (TMarker*)gPad->FindObject("TMarker");  
  mark=(TMarker*)gPad->WaitPrimitive("TMarker","Marker");
  int bkg_hi = mark->GetX();
  delete (TMarker*)gPad->FindObject("TMarker");
  
  mark=(TMarker*)gPad->WaitPrimitive("TMarker","Marker"); // Get the 1st peak initial guess
  peak1 = mark->GetX();
  delete (TMarker*)gPad->FindObject("TMarker");
    
  if (num_peaks == 2)  // Get the 2nd-peak initial position guess
  {	mark=(TMarker*)gPad->WaitPrimitive("TMarker","Marker");
        int peak2 = mark->GetX();
  	delete (TMarker*)gPad->FindObject("TMarker"); 
  }
  	
  if (num_peaks == 3)  // Get the 2nd- and 3rd-initial peak position guesses
  {     mark=(TMarker*)gPad->WaitPrimitive("TMarker","Marker");
  	peak2 = mark->GetX();
  	delete (TMarker*)gPad->FindObject("TMarker");  	
  	mark=(TMarker*)gPad->WaitPrimitive("TMarker","Marker");
  	peak3 = mark->GetX();
  	delete (TMarker*)gPad->FindObject("TMarker");
  }
  //end of the stuff added by Chris------

  int peak1_lo = peak1 - p1, peak1_hi = peak1 + p1; // Peak center and limits
  int peak2_lo = peak2 - p2, peak2_hi = peak2 + p2; // Peak center and limits
  int peak3_lo = peak3 - p3, peak3_hi = peak3 + p3; // Peak center and limits
 
  // Calibration coefficients
  int Ecal;
  if (Ecal == 0){
     float a= 0, b = 1;
  }
  if (Ecal == 1){
     float a = -51.29, b = 1.08174; // Analyzer
     //float a = -45.23, b = 1.1125; // cal at end
  }
  if (Ecal == 2){
     float a = -82.1001, b =1.12346 ;     //Scatterer
  }                                                                                                                                                                                                                                                                                                                                                                                                                                                          ; // For No calibration, Energy = a + b*chan_num
  //float a = -12.5367, b = 1.3268; // Gamma calibration
  //float a = -8.56, b = 1.34605; // Det 1, dgs
  // float a= 310.32, b= 4.756;  //Alphas

  // Efficiency Input Choice
  int get_eff= 0;  //If = 1 Get and use efficiencies, If = 0 Skip Efficiencies
  // Eff. coefs in the "Ge_eff.c" , Load it first by typing: .L Ge_eff.c <ret> 
  //double efit =C0 + C1*ee + C2*pow(ee,2) + C3*pow(ee,3);
  //double Efficiency =exp(efit);
  double c1=1.00158, c2=0.00081172527, c3=0.40473091/1000000;//Coeffs for FWHM

  // Define the file name (+directory) for the output of the fit data.
  // The results of each fit are appended at the end of this file. 
  // Bad fits should be removed by editing the file. 
  ofstream outfile("e:/322/2013/Compton/Peaks.fit",ios::app);
  
  //---------------------------------------------
  fits->Divide(1,1);
  fits->cd(1);
  g1 = new TF1("m1", "pol1", bkg_lo, bkg_hi);
  g2 = new TF1("m2", "gaus", peak1_lo,peak1_hi);
  if (num_peaks >= 2){g3 = new TF1("m3", "gaus", peak2_lo,peak2_hi);}
  if (num_peaks == 3){g4 = new TF1("m4", "gaus", peak3_lo,peak3_hi);}
  
  if (num_peaks == 1)
    {total = new TF1("total","pol1(0)+gaus(2)", bkg_lo,bkg_hi);}
  if (num_peaks == 2)
    {total = new TF1("total","pol1(0)+gaus(2)+gaus(5)", bkg_lo,bkg_hi);}
  if (num_peaks == 3)
    {total = new TF1("total","pol1(0)+gaus(2)+gaus(5)+gaus(8)", bkg_lo,bkg_hi);}
  
  Double_t par[num_par], out[num_par], Area[num_peaks]; 
  
  spectrum->Fit(g1, "R");
  spectrum->Fit(g2, "R+");
  if (num_peaks >= 2){spectrum->Fit(g3, "R+");}
  if (num_peaks == 3){spectrum->Fit(g4, "R+");}
  g1->GetParameters(&par[0]);
  g2->GetParameters(&par[2]);
  if (num_peaks >= 2){g3->GetParameters(&par[5]);}
  if (num_peaks == 3){g4->GetParameters(&par[8]);}
  total->SetParameters(par);
  spectrum->Fit(total,"R");
  total.GetParameters(out);
  
  // Calculate the areas: A = p2*p4*sqrt(2 Pi), etc.
  if (num_peaks >= 1){Area[0]=fabs(out[2]*out[4]*2.5066/comp);}
  if (num_peaks >= 2){Area[1]=fabs(out[5]*out[7]*2.5066/comp);}
  if (num_peaks == 3){Area[2]=fabs(out[8]*out[10]*2.5066/comp);}
  
  if (num_peaks >= 1){float en1 = a + b*out[3];}
  if (num_peaks >= 2){float en2 = a + b*out[6];}
  if (num_peaks == 3){float en3 = a + b*out[9];}

  if (num_peaks >= 1){float res1=fabs(out[4]*2.3548*b);}
  if (num_peaks >= 2){float res2=fabs(out[7]*2.3548*b);}
  if (num_peaks == 3){float res3=fabs(out[10]*2.3548*b);}
 
  if (get_eff == 1) // Get Efficiencies
    {
      if (num_peaks >= 1){double Eff1=Efficiency(en1/1000);}
      if (num_peaks >= 2){double Eff2=Efficiency(en2/1000);}
      if (num_peaks == 3){double Eff3=Efficiency(en3/1000);}

      if (num_peaks >= 1){double fwhm1=c1 + c2*en1 + c3*en1*en1;}
      if (num_peaks >= 2){double fwhm2=c1 + c2*en2 + c3*en2*en2;}
      if (num_peaks == 3){double fwhm3=c1 + c2*en3 + c3*en3*en3;}
    } 
  cout    <<"-----------------------" << endl;
  if (get_eff == 1)  // Write also efficiencies
    {
      cout <<spec<<",  Chan1="<<out[3]<<", E_1(keV)="<<en1<<", Net_Area1=" << Area[0];
      cout <<", FWHM(keV)="<<res1<<endl;
      cout <<" Efficiency="<<Eff1<<", Intensity= "<<Area[0]/Eff1<<", Expect FWHM="<<fwhm1<<endl;
      cout  << "  " << endl;
      outfile <<"-----------------------" << endl;
      outfile <<spec<<",   Chan     En(keV)  Counts     FWHM    Eff    Intensity  Exp_FWHM " << endl;
      outfile <<spec<<",  "<<out[3]<<", " <<en1<<", "<<Area[0]<<", " <<res1<<", " <<Eff1<<", ";
      outfile <<Area[0]/Eff1<<", "<<fwhm1<<endl;
      
      if (num_peaks >= 2)
	{
	  cout <<spec<<",  Chan2="<<out[6]<<", E_2(keV)="<<en2<<", Net_Area2=" <<Area[1];
	  cout <<", FWHM(keV)="<<res2<<endl;
	  cout <<" Efficiency="<<Eff2<<", Intensity= "<<Area[1]/Eff2<<", Expect FWHM="<<fwhm2<<endl;
	  cout << "  " << endl;
	  outfile <<spec<<",  "<<out[6]<<", "<<en2<<", "<<Area[1]<<", " <<res2<<", " <<Eff2<<", ";
	  outfile <<Area[1]/Eff2<<", "<<fwhm2<<endl; 	  
	  //outfile.close();
	}
      
      if (num_peaks == 3)
	{
	  cout <<spec<<",  Chan3="<< out[9]<<", E_3(keV)="<<en3<<", Net_Area3=";
	  cout <<Area[2]<<", FWHM(keV)="<<res3<<endl;
	  cout <<" Efficiency="<<Eff3<<", Intensity= "<<Area[2]/Eff3<<", Expect FWHM="<<fwhm3<<endl;
	  cout <<"-----------------------" << endl;
	  outfile <<spec<<",  "<<out[9]<<", "<<en3<<", "<<Area[2]<<", " <<res3<<", " <<Eff3<<", ";
	  outfile <<Area[2]/Eff3<<", "<<fwhm2<<endl;	  
	}
      outfile.close();     
    }
  if (get_eff == 0) // No efficiencies  output just counts
    {      
      cout <<spec<<",  Chan1="<<out[3]<<", E_1(keV)="<<en1 <<", Net_Area1=" <<Area[0];
      cout <<", FWHM(keV)="<<res1<<endl;     
      outfile <<"-----------------------" << endl; 
      outfile <<spec<<",   Chan    En (keV)   Counts      FWHM " << endl;
      outfile <<spec<<",  "<<out[3]<<", "<<en1<<", "<<  Area[0]<<", " <<res1<<endl;        
      
      if (num_peaks >= 2)
	{
	  cout <<spec<<",  Chan2="<<out[6]<<", E_2(keV)="<<en2 <<", Net_Area2=" <<Area[1];
	  cout <<", FWHM(keV)="<<res2<<endl;
	  //cout << "  " << endl;
	  outfile <<spec<<",  "<<out[6]<<", "<<en2<<", "<<  Area[1]<<", " <<res2<<endl;	  
	}
      
      if(num_peaks == 3) 
	{
	  cout <<spec<<",  Chan3="<<out[9]<<", E_3(keV)="<<en3 <<", Net_Area3=" <<Area[2];
	  cout <<", FWHM(keV)="<<res3<<endl;
	  //cout << "  " << endl;     
	  cout <<"-----------------------" << endl;
	  outfile <<spec<<",  "<<out[9]<<", "<<en3<<", "<<Area[2]<<", " <<res3<<endl; 
	}
      outfile.close();	
    }
  
}
