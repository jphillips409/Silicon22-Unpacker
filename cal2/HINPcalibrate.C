#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TSpectrum.h"
#include "TObject.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2.h"


void calibrate(TString detector="Front"){

  TString inputFile1;

  int polDeg = 1;

  if (detector.Contains("Front"))
  {
    inputFile1 = "peakpositions_Front_U232.txt";
  }
  if (detector.Contains("Back"))
  {
    inputFile1 = "peakpositions_Back_U232.txt";
  }
  if (detector.Contains("FrontLow"))
  {
    inputFile1 = "peakpositions_Frontlow_U232.txt";
  }
  if (detector.Contains("BackLow"))
  {
    inputFile1 = "peakpositions_Backlow_U232.txt";
  }

  cout << "Read in: inputFile = " << inputFile1 << ", detector = " << detector << ", polynomial degree = " << polDeg << endl;

  //FrontEcal.dat
  TString outFileName = "Gobbi";
  outFileName.Append(detector);
  //if (detector.Contains("Low")) outFileName.Append("Low");
  outFileName.Append("Ecal.dat");
  //outFileName.Append("EcalTEST.dat");
  
  cout << "Opening calibration results file: " << outFileName << endl;
  ofstream fout;
  fout.open(outFileName);
  
  const int numDetectorElements = 4*32;
  const int numPeaks = 5;

  Double_t alphas[numPeaks];
  Double_t alphaserr[numPeaks];

  Double_t peaks[numPeaks];
  Double_t peakserr[numPeaks];

  //*************************
  //Alpha energies for U-232 used
  //  U-232:  5.320 MeV
  //  Ra-224: 5.685 MeV
  //  Rn-220: 6.288 MeV
  //  Po-216: 6.778 MeV
  //  Po-212: 8.784 MeV
  //*************************
  
  for (Int_t i=0;i<numDetectorElements;i++)
  {
    //alphas[0] = 5340.3;
   // alphas[1] = 5685.;
    //alphas[2] = 6288.1;
   // alphas[3] = 6778.3;
    //alphas[4] = 8784.8;

    alphas[0] = 8784.8;
    alphas[1] = 6778.3;
    alphas[2] = 6288.1;
    alphas[3] = 5685.;
    alphas[4] = 5340.3;

    alphaserr[0] = 1.0;
    alphaserr[1] = 1.0;
    alphaserr[2] = 1.0;
    alphaserr[3] = 1.0;
    alphaserr[4] = 1.0;
  }

  ifstream fin1;
  fin1.open(inputFile1);

  for(Int_t i=0; i<numDetectorElements; i++)
  {
    cout << "\n\nDoing detector element " << i+1 << endl;
    Int_t ch;
    fin1 >> ch;

    for (Int_t j=0; j<numPeaks; j++) {
      fin1 >> peaks[j] >> peakserr[j];
      //peakserr[j] = 1;
    }
    
    bool zeros = false;
    for (Int_t j=0;j<numPeaks;j++) {
      if (peaks[j] == 0)
        zeros = true;
    }

    if (zeros) continue;

    cout << "read peaks: ";
    for (Int_t j=0;j<numPeaks;j++){cout << peaks[j] << " " << peakserr[j] << " " ;}
    cout << endl << "alphas peaks: ";
    for (Int_t j=0;j<numPeaks;j++){ cout << alphas[j] << " " << alphaserr[j] << " ";}

    TGraphErrors *g1err = new TGraphErrors(numPeaks,peaks,alphas,peakserr,alphaserr);
    g1err->Draw("A*");

    //TF1 *f1 = new TF1("f1","pol1",100.,16000.);
    //g1err->Fit("f1","R+");
    g1err->Fit("pol1");

    int board = i/32;
    int chip = i%32;

    double para[10];
    g1err->GetFunction("pol1")->GetParameters(para);
    fout << board << " " << chip << " " << para[1]/1000. << " " << para[0]/1000. << endl; 

  }

}
