// C++ file to read MSU event files 
//file numbers.beam contains runs to sort
//uses class hira to unpack hira data
//write out spectra in file sort.root
//originally written by Kyle Brown + Robert Charity
//modifications from Nicolas Dronchi

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "det.h"
#include <ctime>
#include "histo_sort.h"
#include "histo_read.h"
//#include "forest.h"
#include "TString.h"


using namespace std;

int main(int argc, char* argv[])
{

  clock_t t;
  t = clock();

  string suffix = "beampurity";
  int setting = 0;
  if (argc > 1)
  {
    setting = 1;
    runnumber = stoi(argv[1]);
    suffix = argv[1];
  } 
  else
    abort();
    

  histo_sort * Histo_sort = new histo_sort(suffix);
  histo_read * Histo_read = new histo_read(suffix);

  //forest * Forest = new forest();
  unsigned short *point,*fragmentstart;
  int unsigned words;
  int unsigned type;
  
  int physicsEvent = 0;
  int physicsEventGood = 0;
  int physicsEventCounter = 0;
  int scalerBuffer = 0;
  int Npauses = 0;
  int Nresumes = 0;
  int runno = 0;
  int NS800 = 0;
  int NGoodS800 = 0;
  int NSecondary = 0;
  int NGoodSecondary = 0;
  int totalRuntime = 0;

  float transport_efficiency = 0.62;

  det Det(Histo_sort, Histo_read, setting);

  //Set the sourceID for the S800 (normally 2) and secondary DAQ
  //This will need to be revisited
  Det.SiID = 1;
  Det.S800ID = 2;
  
  ifstream scalarFile;
  
  //check if this file exists
  if (runFile.is_open() == 0 )
  {
    cout << " could open runfile " << endl;
    return 1;
  }
  
  ostringstream outstring;
  int number;
  int runnum, live_Rand, raw_Rand, RF_Scint; 
  ifstream evtfile;
  ifstream scalarfile;
  FILE *pFile;
  bool fileProblem = false;
  bool endOfRun = false;
  bool first = true; 
  int argcounter = 1;
  for (;;)  //loop over run numbers
  {
    number = runnumber;

/////////////////////////////////////////////////////////////....................
    //check to see if at end of file
    if (runFile.eof())break;
    if (runFile.bad())break;

    if (evtfile.is_open()) 
      cout << "problem previous file not closed" << endl;
    if (scalarfile.is_open())
      cout << "problem previous scalar file not closed" << endl;

    //read in the scaler file the corresponds to the run number
    outstring.str("");
    string scalardirec = Form("/data1/ca36_nov2020/scalers/");
    outstring << scalardirec << "run" << setfill('0') << setw(4) << number << ".csv";

    ostringstream ss;
    cout << endl;
    cout << "reading scalar file: " << outstring.str().c_str() << endl;
    scalarfile.open(outstring.str().c_str(), ifstream::in);
    if (!scalarfile.is_open())
      cout << "could not open scalar file" << endl;
    
    ss << scalarfile.rdbuf();
    string scalar_cont = ss.str();

    istringstream sstream(scalar_cont);
    string record;
    vector<string> items;
    char delim = ',';

    while (getline(sstream,record)) //seperates by /n to get each line
    {
      istringstream line(record);
      while (getline(line,record, delim)) //seperates by , to get each item
      {
        items.push_back(record);
      }
    }   

    int possec = 4;
    //int posscaler = 10;
    int timesec = 0;
    for(;;)
    {
      try{ timesec = stoi(items.at(possec));}
      catch(...) { possec++; continue;}
      break;
    }
    //for(;;)
    //{
    //  if (items.at(posscaler) == )
    //  {
    //    
    //  }
    //  else
    //    posscalar++;
    //}


    totalRuntime += timesec;
    cout << "for run " << number << " beam time is " << timesec << " for total of " <<totalRuntime << endl;
    scalarfile.close();
    //end reading the scaler file

    //Forest->newTree(number);
    for (int iExtra=0;iExtra<3;iExtra++) //loop over split evtfiles
    {
      
      //the following loop accounts for files that were split
      endOfRun=false;
      fileProblem = 0;
      outstring.str("");

      string directory = Form("/data1/ca36_nov2020/");

      outstring << directory << "run" << number <<"/run-" << setfill('0') << setw(4) << number;
//        outstring << directory << "/run-" << setfill('0') << setw(4) << number;

      if (iExtra == 0)
        outstring<<"-00.evt";
      else
        outstring<<"-"<<setfill('0') << setw(2) << iExtra<<".evt";

      string name = outstring.str();
      //cout << name << endl;

      //open evt file
      evtfile.clear();
      evtfile.open(name.c_str(),ios::binary);

      //check to see if there are extra files     
      if (iExtra>0 && !fileProblem && !evtfile) 
      {
        break;
      }

      cout << "reading event file: "<< name << endl;

      if (evtfile.bad()) cout << "bad" << endl;
      if (evtfile.fail()) cout << "fail" << endl;
      if (!evtfile)
      {
        cout << "could not open event file" << endl;
        return 1;
      }
      int test=0;

      /////////////////////////////////////////////////////////////////////////////////////////
      //Set the distance for the Si and Fibers for different settings
      float SiDistance = 40.05;
      //99.79 mg/cm2 -- 0.54 mm target
      //188 mg/cm2 -- 1mm target
      float TargetThick = 91.8456; //99.79; 
      // if(number <= 60)
      //   {
      //     SiDistance = 33.13; //cm
      //     TargetThick = 99.79;
      //   }
      // else
      //   {
      //     SiDistance = 33.13; //cm
      //     TargetThick = 99.79;
      //   }

      Det.Hira->RingCounter->SetDistance(SiDistance);
      Det.Hira->RingCounter->SetTargetThickness(TargetThick);
      Det.Hira->XY_mon->setDistance(SiDistance);
      
      ////////////////////////////////////////////////////////////////////////////////////////
      
      //for(int i=0;i<1000;i++)  // loop over items in a evtfile
      for(;;)      
      {

      /////////////////////////////////////////////////////////////////////////////////////
      //Ring items have the following structure
      //High-Level Description Lower-Level Description Size (bytes)
      //Header  Inclusive Size 4
      //        Type  4
      // Body Header  Size = 20 4
      //        Timestamp       8
      //        Source ID       4
      //        Barrier Type    4
      // Body Data... >=0
      //
      //I am modifing here to convert this into a sort code that will read built events
      //KB Oct 2018
      //The body of built evens has the following structure
      /////////////////////////////////////////////////
      //High-Level Description Lower-Level Description Size (bytes)
      //Body size              Number of bytes in body 4
      ///////////////////////////////////////////////////
      //Fragment #0              Fragment Header       20
      //                         Ring Item Header      8
      //                         Ring Item Body Header 4 or 20
      //                         Ring Item Body Determined by Readout program and user
      ///////////////////////////////////////////////////
      //Fragment #1              Fragment Heade        20
      //                         Ring Item Header      8
      //                         Ring Item Body Header 4 or 20
      //                         Ring Item Body Determinable by Readout program and user
      /////////////////////////////////////////////////////////////////////////////////////

        int const hBufferWords = 4;
        int const hBufferBytes = hBufferWords*2;
        unsigned short hBuffer[hBufferWords];
        evtfile.read((char*)hBuffer,hBufferBytes);

        if(evtfile.eof())
        {
          cout << "eof found" << endl;
          fileProblem = true;
          break;
        }
        if(evtfile.bad())
        {
          cout << " bad found" << endl;
          fileProblem = true;
          break;
        }

        point = hBuffer;
        int offset = 0; //Where we are in the buffer
        int nbytes = *point++; //inclusive ring item size
        int nbytes2 = *point++;
        int type = *point++; //ring item type
        int type2 = *point;

        offset = 4;
 
        int dBufferBytes = nbytes - 8; //skipping the inclusive size and data type
        int dBufferWords = dBufferBytes/2; //calculating 16 bit words from bytes

        unsigned short dBuffer[dBufferWords];
        evtfile.read((char*)dBuffer,dBufferBytes);
        point = dBuffer;

        int BHsize = *point++;
        int BHsize2 = *point++;

        offset +=2;
        int64_t eventTstamp;
        int eventsourceID = 0;
        int eventBarrierType = 0;

        int FragmentSize =0;
        //This assumes that we are dealing with NSCLDAQ-11.0 or later which
        //contains a body header of size 20.
        //If not we are dealing with something older and modifications are necessary
        // KB, Sept 2018

        if(BHsize==20)
        {
          point +=4; //skipping the timestamp for now -- need to come back to this later Kyle
          eventsourceID = *point++;
          point++;
          eventBarrierType = *point++;
          point++;
          offset +=8;
        }
        else
        {
          //Buffers with no body header, for use with NSCLDAQ older than 11.0
        }

        int64_t fragmentTimestamp[5]={0};
        int fragmentsourceID[5] ={-1};
        int fragmentsize[5] ={0};
        int fragmentBarrierType[5] = {0};

        int nFragment = 0;
        int fragmentcounter=0;
        int BuiltSize = 0;
        int BuiltSize2 = 0;

        bool foundS800 = false;
        bool foundSecondary = false;
        Det.Reset();

        while (offset != nbytes/2)
        {
          if (type == 1)
          {
            runno = *point;   
            cout << "run number = " << runno << endl;
            offset = nbytes/2;
          }
          else if (type == 30)
          {
            nFragment=fragmentcounter;
            fragmentcounter++;
            //Reads the body size, which is only present at the start of a built event
            if(nFragment ==0)
            {
              BuiltSize = *point++;
              BuiltSize2 = *point++;
              offset +=2;
            }      
            int nwordsring = BuiltSize/2;    
            //Now to read the fragment header
            fragmentstart=point;
            point+=4; //Skipping the time stamp for now KB
            fragmentsourceID[nFragment]=*point++;
            point++;
            fragmentsize[nFragment] = *point++;
            point++;
            fragmentBarrierType[nFragment]=*point++;
            point++;

            offset+=10;

            //We may want to modify this so that it checks the position of point
            //when the data returns from the unpackers
            //I will think about implementing this at a later date
            //KB Oct 2018
            offset+=fragmentsize[nFragment]/2; //Skipping the #words of the payload;
            point +=14;

            if (physicsEvent%1000 == 0) 
              cout << '\xd'<< physicsEvent << flush;
            physicsEvent++;
            
            //Det.S800ID == 2, Det.SiID == 1
            if(fragmentsourceID[nFragment] == Det.S800ID)
            {
              NS800++;
            }
            else if(fragmentsourceID[nFragment] == Det.SiID)
            {
              NSecondary++;
            }

            bool stat = Det.unpack(point,runno,fragmentsourceID[nFragment]);

            if(stat)
            {
              if(fragmentsourceID[nFragment] == Det.S800ID)
              {
                foundS800 = true;
                physicsEventGood++;
                NGoodS800++;
              }
              else if(fragmentsourceID[nFragment] == Det.SiID)
              {
                foundSecondary = true;
                NGoodSecondary++;
                physicsEventGood++;
              }
            }
            //Let move to the next fragment
            point = fragmentstart +(fragmentsize[nFragment]/2)+10;

            //Turning off trees for now, I will return to this later, not needed for online
            //KB Oct 2018
            //Det.treeGrow(); 
    
          }//end type == 30
          else if (type == 31)
          {
            physicsEventCounter++;
            offset = nbytes/2;
          }
          else if (type == 2)
          {
            endOfRun = true;
            break;
          }
          else if (type == 20)
          {
            scalerBuffer++;
            offset = nbytes/2;
          }
          else if (type == 3)
          {
            Npauses++;
            offset = nbytes/2;
          }
          else if (type == 4)
          {
            Nresumes++;
            offset = nbytes/2;
          }
          else
          {
            offset = nbytes/2;
          }

          nFragment++;

        } //end loop over ring item

        if (foundS800 || foundSecondary) 
        {
          //here we read trigger coin and add
          Det.analyze(physicsEvent, runno);
        }

      } //end loop over items in a evtfile
      evtfile.close();
      evtfile.clear(); // clear event status in case we had a bad file
  
    } //end loop over file subsections
    //Forest->writeTree();
  
  } //end loop of run file numbers
  
  cout << '\n'<<"physics Events = " << physicsEvent << endl;
  cout << "Good physics Events = " << physicsEventGood << endl;
  
  if (physicsEvent > 0) cout << "bad/total = " << 
          (1.-(double)physicsEventGood/(double)physicsEvent)*100.<< " %"<< endl;
  
  cout << "physics Event Counters = " << physicsEventCounter << endl;
  cout << "scaler buffers = " << scalerBuffer << endl;
  cout << "Numbers of pauses = " << Npauses << endl;
  cout << "Number of resumes = " << Nresumes << endl;

/*
  cout << "proton mult 1 = " << Det.Hira->RingCounter->protonYield_s800[1] << endl;
  cout << "proton mult 2 = " << Det.Hira->RingCounter->protonYield_s800[2] << endl;
  cout << " alpha mult 1 =" << Det.Hira->RingCounter->alphaYield_s800 << endl;

  cout << "Ar31 beam " << endl;


  cout << "proton mult 1 = " << Det.Hira->RingCounter->protonYield_Ar31[1] << endl;
  cout << "proton mult 2 = " << Det.Hira->RingCounter->protonYield_Ar31[2] << endl;

  cout << "alpha mult 1 = " << Det.Hira->RingCounter->alphaYield_Ar31 << endl;
  
  cout <<"Ar31 beam, S28 residue" << endl;
  cout << "proton mult 0 = " << Det.Hira->RingCounter->protonYield_Ar31_S28[0] << endl;  
  cout << "proton mult 1 = " << Det.Hira->RingCounter->protonYield_Ar31_S28[1] << endl;
  cout << "proton mult 2 = " << Det.Hira->RingCounter->protonYield_Ar31_S28[2] << endl;
  cout << "proton mult 3 = " << Det.Hira->RingCounter->protonYield_Ar31_S28[3] << endl;

  cout << "alpha mult 1 = " << Det.Hira->RingCounter->alphaYield_Ar31_S28 << endl;
  
  cout << " S29 beam " << endl;
  cout << "proton mult 1 = " << Det.Hira->RingCounter->protonYield_S29[1] << endl;
  cout << "proton mult 2 = " << Det.Hira->RingCounter->protonYield_S29[2] << endl;

  cout <<"S29 beam and P27 residue" << endl;
  cout << "proton mult 1 = " << Det.Hira->RingCounter->protonYield_S29_P27[1] << endl;
  cout << "proton mult 2 = " << Det.Hira->RingCounter->protonYield_S29_P27[2] << endl;

*/

  //Det.addfiberpos();
/*
  Histo_sort->tree->Fill();
  cout << "number of s800 singles = " << Det.N_s800_singles << endl;
  cout << "number of coincidences = " << Det.N_coin << endl;
  cout << "number of single event = " << Det.N_singles << endl;
  cout << "Num of 37Ca beam = " << Det.NCa37beam << endl;
  cout << "Num of 37Ca beam w/ 35K residue = " << Det.beam_residue << endl; 
  cout << "Num of 37Ca beam, 35K residue, fiber= " << Det.beam_residue_fiber << endl;
  cout << "Num of 37Ca beam, 35K residue, NOfiber= " << Det.beam_residue_nofiber << "  "
       << Det.beam_residue_nofiber/(Det.beam_residue_nofiber + Det.beam_residue_fiber)
       << "%" << endl;
  cout << "Num of 37Ca beam, 35K residue, proton= " << Det.beam_residue_RCsoln << endl;
  cout << "Num of 37Ca beam, 35K residue, NOproton= " << Det.beam_residue_noRCsoln
       << Det.beam_residue_noRCsoln/(Det.beam_residue_noRCsoln + Det.beam_residue_RCsoln)
       << "%" << endl;
  cout << "Num of 37Ca beam, 35K residue, proton, fiber= " << Det.beam_residue_RCsoln_fiber << endl;
  cout << "Num of 37Ca beam, 35K residue, proton, NOfiber= " << Det.beam_residue_RCsoln_nofiber << endl;


  cout << "Num of 35K residue (no requirement on beam) = " << Det.NK35residue << endl;
  cout << "Num of 35K with fiber (no requirement on beam) = " << Det.NK35_withfiber << endl;

  cout << "S800 events that we unpacked correctly = " << (float)NGoodS800/(float)NS800*100. << "%" <<endl;
  cout << "RingCounter Events good " << NGoodSecondary << " and unpacked"<< NSecondary << endl;
  cout << "RingCounter events that we unpacked correctly = " << (float)NGoodSecondary/(float)NSecondary*100. << "%" <<endl;

  cout << "Number of good HiRA events = " << Det.Hira->NHira << endl;
  cout << "Number of good fiber events = " << Det.Hira->Nfiber << endl;
  cout << "Number of good S800 events = " << Det.NS800 << endl;

  cout << "Number of events with neighbor pie issue: " << Det.Hira->RingCounter->K35multipiecounter << endl;

  //cout << "Number of S28 found = " << Det.NS28 <<endl;
  //cout << "Number of S28 found with a good fiber signal = " << Det.NS28_withfiber << endl;
*/

  Histo_sort->write(); // this forces the histrograms to be read out to file
  Histo_read->write();

  
  t = clock() - t;
  cout << "run time: " << (float)t/CLOCKS_PER_SEC/60 << " min" << endl;

}
