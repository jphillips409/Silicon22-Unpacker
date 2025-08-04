// g++ calib.cxx -std=c++0x  -o Calibrate
#include <iostream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <fstream>

using namespace std;
void calibrate(string sourcenames[], string rootfiles[], int NumSources ) {
  string sourcename_corr[NumSources];
  string histname = "CalHist.root";
  string temp;
  bool goodsources = true;
  bool goodroot = true;
  // Checks sources entered exist and changes the format if needed.
  for(int i = 0; i < NumSources; i++ ) {
    if( (sourcenames[i].compare("co60") == 0) || (sourcenames[i].compare("60co") == 0) || (sourcenames[i].compare("Co60") == 0) || (sourcenames[i].compare("60Co") == 0) ) sourcename_corr[i] = "co60";
    else if( (sourcenames[i].compare("co56") == 0) || (sourcenames[i].compare("56co") == 0) || (sourcenames[i].compare("Co56") == 0) || (sourcenames[i].compare("56Co") == 0) ) sourcename_corr[i] = "co56";
    else if( (sourcenames[i].compare("eu152") == 0) || (sourcenames[i].compare("152eu") == 0) || (sourcenames[i].compare("Eu152") == 0) || (sourcenames[i].compare("152Eu") == 0) ) sourcename_corr[i] = "eu152";
    else if( (sourcenames[i].compare("ba133") == 0) || (sourcenames[i].compare("133ba") == 0) || (sourcenames[i].compare("Ba133") == 0) || (sourcenames[i].compare("133Ba") == 0) ) sourcename_corr[i] = "ba133";
    else if( (sourcenames[i].compare("ambe") == 0) ) sourcename_corr[i] = "ambe";
    else if( (sourcenames[i].compare("y88") == 0) || (sourcenames[i].compare("88y") == 0) || (sourcenames[i].compare("Y88") == 0) || (sourcenames[i].compare("88Y") == 0) ) sourcename_corr[i] = "y88";
    else if( (sourcenames[i].compare("online") == 0) || (sourcenames[i].compare("Online") == 0)) sourcename_corr[i] = "Online";
    else {
      goodsources = false;
      cout << "Source: " << sourcenames[i] << " Not recognised nothing will happen" << endl;
    }
  }

  for(int i = 0; i < NumSources; i++ ) {
    ifstream f(rootfiles[i].c_str());
    if(!f.good()) {
      cout << "File: " << rootfiles[i] << " Does Not Exists" << endl;
      goodroot = false;
    }
  }
  if(goodsources && goodroot){
    ofstream outfile("sourcelist.dat");
    for(int i = 0; i < NumSources; i++ ) outfile << sourcename_corr[i] << "\t" << rootfiles[i] << endl;
    outfile.close();
    ifstream h(histname.c_str());
    if(h.good()){
      cout << "Calibrations histograms have already been sorted" << endl;
      cout << "Resort y/n?" << endl;
      cin >> temp;
      if(strcmp(temp.c_str(),"y") == 0 || strcmp(temp.c_str(),"Y") == 0) {
        cout << "Resorting calibration histograms" << endl;
        system("root -q histogram_make.c");
      } else {
        cout << "Calibrating existing histograms" << endl;
      }
    } else {
      system("root -q histogram_make.c");
    }
    system("root -q linear_energy.C");
    system("root -q quad_energy.C");
  }
}

int main(int argc, char **argv){

  string sourcenames[10];
  string rootfiles[10];
  int numsources = ((argc-1)/2);
  int check = ((argc-1)%2);
  int count = 0;
  if(argc > 2 && check==0) {
    for(int i = 1; i < argc; i++) {
      if(i%2 == 1){
        sourcenames[count] = argv[i];
      } else if(i%2 == 0){
	rootfiles[count] = argv[i];
	count++;
      }
    }
    calibrate(sourcenames, rootfiles, numsources);
  } else {
    cout << "Wrong Number of Arguments, need to be in form;" << "\n" << "Source1 AnalysisTree1 ... SourceX AnalysisTreeX" << endl;
  }
  return 0;
}
