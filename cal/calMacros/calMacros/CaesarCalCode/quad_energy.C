Int_t num_cores = 192; //number of cores in TIGRESS detector
Int_t num_known_sources = 7; //number of sources that can be used for calibration

//names of usable sources
string source_names[] = {
"eu152", "ba133", "co56", "co60", "y88", "ambe", "Online"
};

//number of peaks used for each source's calibration
Int_t num_peaks[] = {
  3, 4, 4, 2, 2, 1, 1
};

//adds second array to the end of the first one
//note: first array must have enough space for second
void move_into_array(Double_t to[], Double_t from[], int start, int len) {
  for (int i = start; i < start + len; i++) {
    to[i] = from[i - start];
  }//for
}//move_into_array

//loads the histograms from root file
void load_histograms(const char histogram_filepath[], TH1F *hist[], Int_t source_count, Int_t num_cores) {
  TFile *infile = new TFile(histogram_filepath);
  for (int i = 0; i < num_cores; i++) {
    hist[i] = (TH1F*)infile->Get(Form("hist%i_%i",source_count,i));
  }
}

//fits an equation on the centroids vs energies for each core of the TIGRESS detector
void fit_equation(TList *list, Double_t energy[], Double_t energy_er[], Double_t** centroids, Double_t** centroids_er, Int_t num_peaks_used, Double_t thirds[num_cores], Double_t gains[num_cores], Double_t offsets[num_cores]) {

  ofstream quad_fit("quad_fit_points.txt");
  for (int i = 0; i < num_peaks_used; i++) {
    quad_fit << energy[i] << " " << energy_er[i] << endl;
  }//for
  for (int i = 0; i < num_cores; i++) {
    for (int j = 0; j < num_peaks_used; j++) {
      quad_fit << centroids[i][j] << " " << centroids_er[i][j] << endl;
    }//for

    //if this centroid is null, then skip this core
    if (centroids[i][0] == -1) {
      thirds[i] = 0;
      gains[i] = 0;
      offsets[i] = -1;
      continue;
    }//if
    //creating graph of centroids vs energy
    TGraphErrors *gr = new TGraphErrors(num_peaks_used, centroids[i], energy, centroids_er[i], energy_er);
    gr->SetName(Form("Equation Fit - %i", i));
    gr->Draw("AP");

    list->Add(gr); //saving graph for future
    //fitting quadratic on data points
    TF1 *coeffFit;
    coeffFit = new TF1("coeffFit", "[0] + [1]*x + [2]*x*x");
    gr->Fit(coeffFit, "Q");
    thirds[i] = coeffFit->GetParameter(2);
    gains[i] = coeffFit->GetParameter(1);
    offsets[i] = coeffFit->GetParameter(0);

    double x[num_peaks_used];
    double resid[num_peaks_used];
    for (int j = 0; j < num_peaks_used; j++) {
      x[j] = j;
      resid[j] = energy[j] - coeffFit->Eval(centroids[i][j]);
    }//for

    TGraph *resid_g = new TGraph(num_peaks_used, x, resid);
    resid_g->SetName(Form("Fit Residuals - %i", i));
    list->Add(resid_g);
  } //for
}//fit_quadratic

void fit_lin_equation(TList *list, Double_t energy[], Double_t energy_er[], Double_t** centroids, Double_t** centroids_er, Int_t num_peaks_used, Double_t gains[num_cores], Double_t offsets[num_cores]) {
  for (int i = 0; i < num_cores; i++) {
    //if this centroid is null, then skip this core
    if (centroids[i][0] == -1 || centroids_er[i][0] == -1) {
      gains[i] = 0;
      offsets[i] = -1;
      continue;
    } //if
    //creating graph of centroids vs energy
    TGraphErrors *gr = new TGraphErrors(num_peaks_used, centroids[i], energy, centroids_er[i], energy_er);
    gr->SetName(Form("Equation Fit Linear - %i", i));
    gr->Draw("AP");

    list->Add(gr); //saving graph for future

    TF1 *lincoeffFit;
    lincoeffFit = new TF1("lincoeffFit", "[0] + [1]*x");
    lincoeffFit->SetParameters(0,1.5);
    lincoeffFit->SetParLimits(0,-500,500);
    lincoeffFit->SetParLimits(1,0.5,4);
    gr->Fit(lincoeffFit, "Q");
    gains[i] = lincoeffFit->GetParameter(1)/1000.0;
    offsets[i] = lincoeffFit->GetParameter(0)/1000.0;
  }//for
}//fit_linear

//using linear calibration, estimates the location of the centroids
//and then refines that estimate
void find_centroids(TH1F *hist[], Double_t energy[], Double_t energy_er[], Double_t lin_gains[], Double_t lin_offsets[],
					Int_t start_pos, Int_t num_peaks_used, Double_t** centroids, Double_t** centroids_er) {

  for (int i = 0; i < num_cores; i++) {
    //if this linear equation is null, skip this core
    if ((lin_gains[i] == -1 && lin_offsets[i] == -1) || (lin_gains[i] == -1 && lin_offsets[i] == 0) || (lin_gains[i] == -0 && lin_offsets[i] == -1)) {
      for (int j = 0; j < num_peaks_used; j++) {
	centroids[i][j] = -1;
	centroids_er[i][j] = -1;
      }//for
      continue;
    }//if
hist[i]->RebinX(2);

    for (int j = 0; j < num_peaks_used; j++) {
      int max_bin = hist[i]->GetXaxis()->GetXmax();
      int min_bin = hist[i]->GetXaxis()->GetXmin();
      //using linear calibration
      //guess where this peak would be in this core's gamma ray spectrum
      Double_t x_guess = (energy[j] - lin_offsets[i]) / lin_gains[i];
      Double_t init_guess = (energy[j] - lin_offsets[i]) / lin_gains[i];
      Int_t bin_guess = hist[i]->GetXaxis()->FindBin(x_guess);
      hist[i]->SetAxisRange(x_guess - (45 * 2000/1000), x_guess + (45 * 2000/1000),"X");
      //to help with the peak fitting, move initial centroid guess to bin with most counts
      Int_t bin = hist[i]->GetMaximumBin();
      x_guess = hist[i]->GetXaxis()->GetBinCenter(bin);
      hist[i]->SetAxisRange(min_bin, max_bin,"X");

double minrange = init_guess - (45 * init_guess/1000);
double maxrange = init_guess + (45 * init_guess/1000);

      hist[i]->SetAxisRange(minrange, init_guess-8, "X");

      Int_t leftbin = hist[i]->GetMaximumBin();
      double lefty = hist[i]->GetBinContent(leftbin);

      hist[i]->SetAxisRange(min_bin, max_bin,"X");
      hist[i]->SetAxisRange(init_guess+8, maxrange,"X");

      Int_t rightbin = hist[i]->GetMaximumBin();
      double righty = hist[i]->GetBinContent(rightbin);
cout << (righty / lefty) << "\t" << righty << "\t" << lefty << "\t" << rightbin << "\t" << leftbin << "\t" << bin << endl;
      if((righty / lefty) > 0.5 && lefty>5) bin = rightbin;

      x_guess = hist[i]->GetXaxis()->GetBinCenter(bin);

      Double_t y_guess = hist[i]->GetBinContent(bin);
      hist[i]->SetAxisRange(min_bin, max_bin,"X");

      hist[i]->SetAxisRange(x_guess - (45 * x_guess/1000), x_guess + (45 * x_guess/1000),"X");
      //to help with the peak fitting, move initial centroid guess to bin with most counts
      bin = hist[i]->GetMaximumBin();
      x_guess = hist[i]->GetXaxis()->GetBinCenter(bin);
      y_guess = hist[i]->GetBinContent(bin);
      hist[i]->SetAxisRange(min_bin, max_bin,"X");

      TF1 *fit = new TF1(Form("fit %i-%i", i, j),"[0]*(exp(-((x-[1])^2/(2*[2]^2)))) + [3] + [4]*x",x_guess - (50 * x_guess/1000), x_guess + (50 * x_guess/1000));
      fit->SetParameters(y_guess, x_guess, 15, 0, 0);
      fit->SetParLimits(0, y_guess*0.5, y_guess*1.5); //area
      fit->SetParLimits(1, x_guess - (50 * x_guess/1000), x_guess + (50 * x_guess/1000)); //centroid
      fit->SetParLimits(2, 3.0, 30); //sigma of gaussian distribution
      //fitting equation and saving centroids
      hist[i]->Fit(fit, "QR+");
      centroids[i][start_pos + j] = fit->GetParameter(1);
      centroids_er[i][start_pos + j] = fit->GetParError(1);
      cout << "Graph " << i << " Guess: " << init_guess << "\t" << x_guess  << " Actual: " << centroids[i][start_pos + j] << endl; // << " " << fit->GetChisquare()/fit->GetNDF() << endl;
    }//for
  }//for
}//find_centroids

void quad_energy() {

  /*******Initialization of Source Energy Peak Data********/
  Double_t** source_energy = new Double_t*[num_known_sources];
  source_energy[0] = new Double_t[num_peaks[0]]{121.77,244.6976, 344.2785}; //eu152
  source_energy[1] = new Double_t[num_peaks[1]]{276.4, 302.851, 356.013, 383.848}; //ba133
  source_energy[2] = new Double_t[num_peaks[2]]{1771.3567, 2034.7097, 2598.500, 3253.5030}; //co56
  source_energy[3] = new Double_t[num_peaks[3]]{1173.240, 1332.508}; //co60
  source_energy[4] = new Double_t[num_peaks[4]]{898.042, 1836.063}; //y88
  source_energy[5] = new Double_t[num_peaks[5]]{4438.94}; //ambe
  source_energy[6] = new Double_t[num_peaks[6]]{4438.94}; //online

  Double_t** source_energy_er = new Double_t*[num_known_sources];
  source_energy_er[0] = new Double_t[num_peaks[0]]{0.0008,0.0008, 0.0012}; //eu152
  source_energy_er[1] = new Double_t[num_peaks[1]]{0.0021, 0.0016, 0.0017, 0.0012}; //ba133
  source_energy_er[2] = new Double_t[num_peaks[2]]{0.0039, 0.0047, 0.004, 0.0044}; //co56
  source_energy_er[3] = new Double_t[num_peaks[3]]{0.003, 0.004}; //co60
  source_energy_er[4] = new Double_t[num_peaks[4]]{0.003,0.012}; //y88
  source_energy_er[5] = new Double_t[num_peaks[5]]{0.05}; //ambe
  source_energy_er[6] = new Double_t[num_peaks[6]]{0.05}; //Online
  /*************************************************************/

  Int_t num_sources = 0; //number of sources used in this calibration
  Int_t num_peaks_used = 0; //number of peaks to process in this calibration

  TList *list = new TList;
  // Loads Source list and analysistree files from sourcefile.dat
  string source[5];
  string rootfile[5];
  int iii = 0;
  string sourcefile = "sourcelist.dat";
  ifstream fp;
  fp.open(sourcefile.c_str());
  while (fp.good()) {
    fp >> source[iii] >> rootfile[iii];
    iii++;
  }
  num_sources = iii-1;
  if (num_sources == 0) {
    cout << "No sources inputted. Exiting program..." << endl;
    return;
  }//if
  for(int i = 0; i < num_sources; i++) {
    for(int j = 0; j < num_known_sources; j++) {
      if ( source[i].compare(source_names[j])==0) {
        num_peaks_used += num_peaks[j];
      }
    }
  }

  Double_t energy[num_peaks_used];
  Double_t energy_er[num_peaks_used];

  Double_t** centroids = new Double_t*[num_cores];
  Double_t** centroids_er = new Double_t*[num_cores];
  for (int i = 0; i < num_cores; i++) {
    centroids[i] = new Double_t[num_peaks_used];
    centroids_er[i] = new Double_t[num_peaks_used];
  }//for

  TH1F *quad_hist[num_sources][num_cores];
  string sources_used_files[num_sources];

  Double_t lin_gains[num_cores];
  Double_t lin_offsets[num_cores];
  ifstream coeff_fr("lin_energy_coeff.txt");

  for (int i = 0; i < num_cores; i++) {
    coeff_fr >> lin_gains[i];
    coeff_fr >> lin_offsets[i];
  }//for

  coeff_fr.close();
  Int_t num_proc_peaks = 0;
  Int_t num_proc_sources = 0;
  for (int k = 0; k < num_sources; k++) {
    for (int i = 0; i < num_known_sources; i++) {
      if ( source[k].compare(source_names[i])==0) {
 	move_into_array(energy, source_energy[i], num_proc_peaks, num_peaks[i]);
	move_into_array(energy_er, source_energy_er[i], num_proc_peaks, num_peaks[i]);
	load_histograms("CalHist.root", quad_hist[num_proc_sources], k, num_cores);
 	for (int j = 0; j < num_cores; j++) {
   	  list->Add(quad_hist[num_proc_sources][j]);
	}//for
	find_centroids(quad_hist[num_proc_sources], source_energy[i], source_energy_er[i],lin_gains, lin_offsets, num_proc_peaks, num_peaks[i],	centroids, centroids_er);
	num_proc_peaks += num_peaks[i];
	num_proc_sources++;
      }//if
    }//for
  }//for
  Double_t quad_thirds[num_cores];
  Double_t quad_gains[num_cores];
  Double_t quad_offsets[num_cores];

  Double_t linear_gains[num_cores];
  Double_t linear_offsets[num_cores];

  fit_equation(list, energy, energy_er, centroids, centroids_er, num_peaks_used, quad_thirds, quad_gains, quad_offsets);
  fit_lin_equation(list, energy, energy_er, centroids, centroids_er, num_peaks_used, linear_gains, linear_offsets);

  cout << "Writing fitted quadratic coefficients to 'quad_energy_coeff.txt'..." << endl;
  ofstream quad_out("quad_energy_coeff.txt");
  for (int i = 0; i < num_cores; i++) {
    quad_out << 0 << "\t" << i << "\t" << linear_gains[i] << "\t" << linear_offsets[i] << endl;
//    quad_out << linear_gains[i] << "\n" <<  linear_offsets[i] << endl;
  }//for

  quad_out.close();

  TFile * outfile = new TFile("quad_cal.root", "RECREATE");
  outfile->cd();
  list->Write();
  outfile->Close();
}
