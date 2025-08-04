/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * asicCalibrate.C
 * A ROOT script to do peak search for the updated TECSA+HINP setup.
 *
 * Author: Antti J. Saastamoinen, Cyclotron Insititute, Texas A&M University (ajsaasta@comp.tamu.edu)
 * Uses ROOT so GNU LGPL applies. See more https://root.cern/about/license/
 *
 * Based on various old scripts from who knows where.
 *
 * First version: 2021-11-04
 * This version 2021-11-11
 *
 * Usage: root[] .L asicCalibrate.C++                         <-- load this stuff
 * Possible uses:
 *        root[] calibrate("peak_positions_from_peaksearch.txt","Up","Sector","Ra226",1)    <-- calibrate whatever detector / elements with sources
 *
 *  arguments required:
 *          inputFile = file produced by asicPeakSearch.C
 *          detector = Up/Down (defaults to Up)
 *          detectorSide = Sector/Ring (defaults to Sector)
 *          source = 4peak/Ra226/pulser (defaults to 4peak)
 *          polDeg = degree of the fitted polynomial (defaults to 1, i.e. linear)
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


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

  TCanvas *cChi;
  cChi = new TCanvas("cChi","StripChiSqds");
  TCanvas *cRes;
  cRes = new TCanvas("cRes","Residuals");
  TH2F *hResiduals = new TH2F("hResiduals","Residuals",128,-0.5,127.5,4096,0,16383);

  TList *list = new TList;

  int polDeg = 1;

  TString inputFile1;
  TString inputFile2;
  TString inputFile3;

  if (detector.Contains("Front"))
  {
    inputFile1 = "peakpositions_Front_U232.txt";
  }
  if (detector.Contains("Back"))
  {
    inputFile1 = "peakpositions_Back_4peak.txt";
    inputFile2 = "peakpositions_Back_Ra226.txt";
  }

  cout << "Read in: inputFile = " << inputFile1 << ", detector = " << detector << ", polynomial degree = " << polDeg << endl;
  
  //FrontEcal.dat
  TString outFileName = "";
  outFileName.Append(detector);
  //outFileName.Append("Ecal.dat");
  outFileName.Append("EcalTEST.dat");
  
  cout << "Opening calibration results file: " << outFileName << endl;
  ofstream fout;
  fout.open(outFileName);
  
  const int numDetectorElements = 4*32;
  const int numPeaks = 5;
  
  Double_t chisqdarray[numDetectorElements], stripno[numDetectorElements], residualarray[numDetectorElements][numPeaks], peaks[numPeaks], alphastemp[numPeaks], alphas[numDetectorElements][numPeaks], peakerr[numPeaks], alphastemperr[numPeaks], alphaserr[numDetectorElements][numPeaks];
  Double_t tempcalcval[numPeaks];
  Double_t tenparray[numPeaks][numDetectorElements+1];
  Double_t tenparrayerr[numPeaks][numDetectorElements+1];
//  double alphasDL[4][24];
  // 239Pu value seems wrong in this (found this from old SAMUTAMU script)
//  alphas[0] = 5795.,alphas[1] = 5479.3,alphas[2] = 5130.4,alphas[3] = 3182.7; // lit. values
//  alphas[0] = 5719.3,alphas[1] = 5400.3,alphas[2] = 5047.4,alphas[3] = 3073.8; // after 0.5 um Al window
  
  // Lit. values as per NNDC. In case of Cm, Am and Pu energies are weighted averages of major lines as those do not resolve to separate lines with TTT.
  //  alphas[0] = 5795.,alphas[1] = 5479.3,alphas[2] = 5148.4,alphas[3] = 3182.7;
//  alphas[0] = 5719.3,alphas[1] = 5400.3,alphas[2] = 5066.4,alphas[3] = 3073.8; // after 0.5 um Al window (LISE++ phys. calc. for loss calculation)

  //TGraph *g1 = new TGraph(numPeaks,peaks,alphastemp);
  TGraph *g2 = new TGraph(numDetectorElements,stripno,chisqdarray);
  //TGraphErrors *g1err = new TGraphErrors(numPeaks,peaks,alphastemp,peakerr,alphastemperr);
  
  for (Int_t i=0; i<numDetectorElements; i++) {
   
    //alphas[i]={5.3403,5.4232,5.6854,6.059,6.2881,6.7783,8.7848};

    //*************************
    //Alpha energies for U-232 used
    //  U-232:  5.320 MeV
    //  Ra-224: 5.685 MeV
    //  Rn-220: 6.288 MeV
    //  Po-216: 6.778 MeV
    //  Po-212: 8.784 MeV
    //*************************
  
    alphas[i][0] = 5340.3;
    alphas[i][1] = 5685.;
    alphas[i][2] = 6288.1;
    alphas[i][3] = 6778.3;
    alphas[i][4] = 8784.8;


    /*alphas[i][0] = 5795.,alphas[i][1] = 5479.3, alphas[i][2] = 5148.4, alphas[i][3] = 3182.7; // Lit. values as per NNDC. In case of Cm, Am and Pu energies are weighted averages of     alphas[i][0] = 5795.,alphas[i][1] = 5479.3, alphas[i][2] = 5148.4, alphas[i][3] = 3182.7; // Lit. values as permajor lines as those do not resolve to separate lines*/
    alphaserr[i][0] = 1.0, alphaserr[i][1] = 1.0,alphaserr[i][2] = 1.0,alphaserr[i][3] = 1.0,alphaserr[i][4] = 1.0;

    //alphas[i][4] = 7686.82, alphas[i][5] = 6002.55, alphas[i][6] = 5489.48, alphas[i][7] = 5304.33, alphas[i][8] = 4784.34;
    //alphaserr[i][4] = 1.0, alphaserr[i][5] = 1.0,alphaserr[i][6] = 1.0,alphaserr[i][7] = 1.0,alphaserr[i][8] = 1.0;

    //lphas[i][9] = 34369.04*0.98;
    
    //if (detector.Contains("Delta"))
      //alphas[i][9] = 7665.0;

    //alphaserr[i][9] = 1.0;

  }


  ifstream fin1;
  fin1.open(inputFile1);

  for(Int_t i=0; i<numDetectorElements; i++){
    cout << "\n\nDoing detector element " << i+1 << endl;
    stripno[i] = i+1;
    Int_t ch;
    fin1 >> ch;

    for (Int_t j=0; j<numPeaks; j++) {
      fin1 >> peaks[j] >> peakerr[j];
      peakerr[j] = 1;
      tenparray[j][i+1] = peaks[j];
      tenparrayerr[j][i+1] = peakerr[j];
    }
    
    bool zeros = false;
    for (Int_t j=0;j<numPeaks;j++) {
      if (peaks[j] == 0)
        zeros = true;
    }

    if (zeros) continue;
    
    cout << "read peaks: ";
    for (Int_t j=0;j<numPeaks;j++){cout << peaks[j] << " " << peakerr[j] << " " << tenparray[j][i+1] << " " << tenparrayerr[j][i+1] << " " ;}

    cout << "\n";

    //Hack here to get the energies to smaller array expected by TGraph

    for (Int_t j=0;j<numPeaks;j++) {
      alphastemp[j] = alphas[i][j];
      tenparray[j][0] = alphastemp[j];
      alphastemperr[j] = alphaserr[i][j];
      tenparrayerr[j][0] = alphastemperr[j];
    }

    cout << "alphas temp peaks: ";
    for (Int_t j=0;j<numPeaks;j++){ cout << alphastemp[j] << " " << alphastemperr[j] << " ";}

    TGraph *g1 = new TGraph(numPeaks,peaks,alphastemp);
    TGraphErrors *g1err = new TGraphErrors(numPeaks,peaks,alphastemp,peakerr,alphastemperr);

    if (polDeg == 1) {
      list->Add(g1err);
     // TF1 *f1 = new TF1("f1","pol1",100.,16000.);
     // g1err->Fit("f1","R+");

      //TF1* myfit;
      //myfit = (TF1*)g1err->GetFunction("f1");

      int board = i/32;
      int chip = i%32;

      //fout << board << " " << chip << " " << myfit->GetParameter(1)/1000 << " " << myfit->GetParameter(0)/1000 << endl;

       g1err->Fit("pol1");
       double para[10];
       g1err->GetFunction("pol1")->GetParameters(para);
       fout << board << " " << chip << " " << para[0]/1000. << " " << para[1]/1000. << endl;   

      //chisqdarray[i] = myfit->GetChisquare()/myfit->GetNDF();
      //cout << "Chisqd: " << myfit->GetChisquare() << " ndf: " << myfit->GetNDF() << " \nChiSqd/NDF = " << myfit->GetChisquare()/myfit->GetNDF() << " \nProbability: " << myfit->GetProb() << endl;


    for (Int_t j=0;j<numPeaks;j++) {
      residualarray[i][j] = (( para[1]* peaks[j] + para[0] ) - alphastemp[j]);
    }
  }

  g2->SetPoint(i,peaks[i],chisqdarray[i]);
    
  }// end loop over strips
  

  cChi->cd();
  
  g2->Draw("AB");
  
  cRes->cd();

  
  for (Int_t i=0; i<numPeaks; i++) {
    for (Int_t j=0; j<numDetectorElements; j++) {
      hResiduals->Fill(j+1,peaks[i],residualarray[j][i]);
      hResiduals->GetXaxis()->SetTitle("Strip");
      hResiduals->GetYaxis()->SetTitle("Ch");
      hResiduals->Draw("colz");
    }
  }
  
  
  fout.close();
  fin1.close();

  TFile * outfile = new TFile("linearfit.root", "RECREATE");
  outfile->cd();
  list->Write();
  outfile->Close();


}
