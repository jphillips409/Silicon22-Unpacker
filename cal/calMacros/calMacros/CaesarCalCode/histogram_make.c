//file path to default calibration file
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <dirent.h>

Double_t num_bins = 1024; //number of bins to read in from histograms
Double_t min_bin = 0; //first bin in histogram
Double_t max_bin = 4096; //last bin in histogram
Int_t num_cores = 192;

//Loads histogram from one subrun only
void load_histograms(const char analysis_filepath[], TH1F *hist[], Int_t source_count, Int_t num_cores) {

  //opening analysis root file
  TFile *input_file = new TFile(analysis_filepath, "READ");
  if (!input_file->IsOpen()) {
    cout << "Cannot open input file!" << endl;
    return 0;
  }//if

  char hname[20];
  char hname2[20];
  for(int i = 0; i < num_cores; i++) {
    sprintf(hname2,"ceasar/raw/EC%3.3i",i);
    hist[i] = (TH1F *) input_file->Get(hname2);
    sprintf(hname,"hist%i_%i",source_count, i);
    hist[i]->SetName(hname);
    //hist[i] = new TH1F(hname, hname, num_bins, min_bin, max_bin);
  }

  cout << "Histograms created." << endl;
}
//main function
void histogram_make() {

  TList *list = new TList;

  Int_t num_sources = 0; //number of sources used in this calibration
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
  }//

  TH1F *cal_hist[num_sources][num_cores];
  for (int k = 0; k < num_sources; k++) {
    load_histograms(rootfile[k].c_str(), cal_hist[k], k, num_cores);
    for (int j = 0; j < num_cores; j++) list->Add(cal_hist[k][j]);
  }
  //writes out the histograms in the TList
  TFile *out = new TFile("CalHist.root", "RECREATE");
  out->cd();
  list->Write();
  out->Close();//closes files

}//histogram_make
