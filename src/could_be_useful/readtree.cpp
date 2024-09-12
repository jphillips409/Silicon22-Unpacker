#include "event.h"
#include "histo_sort.h"
#include "histo_read.h"
#include "det.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include <memory>
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
  //  TFile f("forest.root");
  histo_read * Histo_read = new histo_read();
  //histo_sort * Histo_sort = new histo_sort();
  //forest * Forest = new forest();
 
  Event event;
  det Det(Histo_read);


  ifstream runFile;
  runFile.open("numbers.read");
  if (runFile.is_open() == 0)
    {
      cout << " couldn't open runfile " << endl;
      return 1;
    }

  ostringstream outstring;
  int number;
  int argcounter = 1;

  for(;;) //loop over run numbers
    {
      if (argc == 1)
	{
	  runFile >> number;
	  if (runFile.eof())break;
	  if (runFile.bad())break;
	}
      else if (argc > 1)
	{
	  if (argcounter == argc)
	    {
	      break;
	    }
	  else 
	    {
	      number = atoi(argv[argcounter]);
	      argcounter++;
	    }
	}
      outstring.str("");

      outstring << "forest_" << number<< ".root";
      
      string name = outstring.str();
      cout << name << endl;

      TFile file (name.c_str());
      TTree * tree = (TTree*)file.Get("T");



      tree->SetBranchAddress("nfront",&event.nfront);
      tree->SetBranchAddress("frontE",event.frontE);
      tree->SetBranchAddress("frontT",event.frontT);
      tree->SetBranchAddress("frontID",event.frontID);
      
      tree->SetBranchAddress("nback",&event.nback);
      tree->SetBranchAddress("backE",event.backE);
      tree->SetBranchAddress("backT",event.backT);
      tree->SetBranchAddress("backID",event.backID);
      
      tree->SetBranchAddress("ncsi",&event.ncsi);
      tree->SetBranchAddress("csiE",event.csiE);
      tree->SetBranchAddress("csiER",event.csiER);
      tree->SetBranchAddress("csiT",event.csiT);
      tree->SetBranchAddress("csiID",event.csiID);
      
      int N = tree->GetEntries();
      for(int i=0; i<N; i++)
	{
	  tree->GetEntry(i);
	  Det.loadTree(&event);
	  Det.analyze(i,0); //0 means nothing, just a placeholder since this is not being used at the moment and i want this to compile with the current version of the det object
	}
    }
  cout << "Finished loading" << endl;
  Histo_read->write();
}
