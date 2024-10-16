// C++ file to read MSU event files 
//file numbers.beam contains runs to sort
//uses class hira to unpack hira data
//write out spectra in file sort.root
//originally written by Kyle Brown + Robert Charity
//modifications from Nicolas Dronchi and Johnathan Phillips for the Si22 2024 experiment

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

//TODO clean up this code a bit
int main(int argc, char* argv[])
{

  clock_t t;
  t = clock();

  //string suffix = "temp2";
  string suffix = ""; //TODO no suffix for now
  int setting = 1;

  //TODO choose a S800 setting. Defines gobbi distance and target thickness
  /*if (argc > 1)
  {
    suffix = argv[1];
    if (suffix == "setting1")
      setting = 1;
    else if (suffix == "setting2")
      setting = 2;
    else if (suffix == "setting3")
      setting = 3;
    else
      abort();
  }*/

  setting = 0; //just default 0 for now

  //TODO also useful but needs to be commented out
  //open file with run numbers to sort
  ifstream runFile;
  ifstream scalarFile;
  ofstream crossSections;

  /*if (suffix == "setting1")
    runFile.open("numbers.s1");
  else if (suffix == "setting2")
    runFile.open("numbers.s2");
  else if (suffix == "setting3")
    runFile.open("numbers.s3");
  else
    runFile.open("numbers.beam");*/

  //runFile.open("numbers.beam");

  suffix += "";

  histo_sort * Histo_sort = new histo_sort(suffix);
  histo_read * Histo_read = new histo_read(suffix);

  targthick * TargThick;
  TargThick = targthick::instance();

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

  
  //check if this file exists
  //if (runFile.is_open() == 0 )
  //{
    //cout << " could open runfile " << endl;
    //return 1;
  //}
  
  ostringstream outstring;
  int number;
  int runnum, live_Rand, raw_Rand, RF_Scint; 
  ifstream evtfile;
  //ifstream scalarfile;
  FILE *pFile;
  bool fileProblem = false;
  bool endOfRun = false;
  bool first = true; 
  int argcounter = 1;


  int firstrunnum;

  // default case, get run numbers from numbers.beam
  if (argc == 1)
  {
    runFile.open("numbers.beam");
    //check if file exits
    if (!runFile.is_open())
    {
      cout << "could not open numbers.beam" << endl;
      return 1; 
    }
    runFile >> firstrunnum;
    runFile.clear(); //clear bad state
    runFile.seekg(0); //return to start
  }

  else if (argc == 2) // check to see if the argument if a file and if so get run numbers from it.
  {
    runFile.open(argv[argcounter]);
    if (runFile.is_open())
    {
      argc = 1;
      cout << "using runs from " << argv[argcounter] << endl;
    }     
    firstrunnum = atoi(argv[argcounter]);
  } 
  TargThick->SetThick(runnum);

  det Det(Histo_sort, Histo_read, setting);

  //Set the sourceID for the S800 (normally 2) and secondary DAQ
  //TODO need correct ID's for HINP, CAESAR, S800, and Janus
  Det.SiID = 1;
  Det.JanusID = 2;
  Det.S800ID = 3;

  int runcount = 0;
  for (;;)  //loop over run numbers
  {
    if (argc == 1) // get runnumbers from a file
    {
      runFile >> number;
      if (runFile.eof()) break;
      if (runFile.bad()) break;
    }
    else if (argc > 1)  // get run numbers from keyboard input
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


    //check to see if at end of file
    if (runFile.eof())break;
    if (runFile.bad())break;


    //TODO comment out scalarfile for now
    if (evtfile.is_open()) 
      cout << "problem previous file not closed" << endl;
    //if (scalarfile.is_open())
      //cout << "problem previous scalar file not closed" << endl;

    //read in the scaler file the corresponds to the run number
    outstring.str("");
    //string scalardirec = Form("/data1/ca36_nov2020/scalers/");
    //outstring << scalardirec << "run" << setfill('0') << setw(4) << number << ".csv";

    ostringstream ss;
    cout << endl;
    //cout << "reading scalar file: " << outstring.str().c_str() << endl;
    //scalarfile.open(outstring.str().c_str(), ifstream::in);
    //if (!scalarfile.is_open())
      //cout << "could not open scalar file" << endl;
    
    //ss << scalarfile.rdbuf();
    //string scalar_cont = ss.str();

    //istringstream sstream(scalar_cont);
    string record;
    vector<string> items;
    char delim = ',';

   // while (getline(sstream,record)) //seperates by /n to get each line
    //{
      //istringstream line(record);
      //while (getline(line,record, delim)) //seperates by , to get each item
      //{
        //items.push_back(record);
      //}
    //}   

    int possec = 4;
    //int posscaler = 10;
    int timesec = 0;
    //TODO don't know what this does, comment out for now
    /*for(;;)
    {
      try{ timesec = stoi(items.at(possec));}
      catch(...) { possec++; continue;}
      break;
    }*/
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
    //scalarfile.close();
    //end reading the scaler file

    //Forest->newTree(number);
    for (int iExtra=0;iExtra<3;iExtra++) //loop over split evtfiles
    {
      
      //the following loop accounts for files that were split
      endOfRun=false;
      fileProblem = 0;
      outstring.str("");

      string directory = Form("");

      outstring << directory << "data/run" << number <<"/run-" << setfill('0') << setw(4) << number;
//        outstring << directory << "/run-" << setfill('0') << setw(4) << number;

      if (iExtra == 0)
        outstring<<"-00.evt";
      else
        outstring<<"-"<<setfill('0') << setw(2) << iExtra<<".evt";

      string name = outstring.str();
        
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
      //TODO target distance and thickness will vary across S800 settings

      //Set the distance for the Si and Fibers for different settings
      float SiDistance = 40.05;//*0.99;
      //99.79 mg/cm2 -- 0.54 mm target
      //91.18456 from 0.497 mm thick Be target using rigidity measurement
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

      //TODO verify that targ dist and thick work. Also set fiber dist
      Det.Gobbi->SetTarget(SiDistance,TargetThick);

      //Det.Hira->XY_mon->setDistance(SiDistance);
      
      ////////////////////////////////////////////////////////////////////////////////////////
      
      //for(int i=0;i<5;i++)  // loop over items in a evtfile
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

        //TODO because janus had odd event lengths, we can read those events out to the buffer. We must pass the ifstream
        //     to the unpacker

        //TODO replace "point" with evtfile read statements and pass to unpack

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

        int nbytes = hBuffer[0];
        int nbytes2 = hBuffer[1];
        int type = hBuffer[2];
        int type2 = hBuffer[3];

        //cout << "here " << hex << nbytes << " " << nbytes2 << " " << type << " " << type2 << dec << endl;
        int offset = 0;
        offset = 8;
 
        int dBufferBytes = nbytes - 8; //skipping the inclusive size and data type
        int dBufferWords = dBufferBytes/2; //calculating 16 bit words from bytes

        unsigned short BHsize_arr[2];
        evtfile.read((char*)BHsize_arr,4);
        int BHsize = BHsize_arr[0];
        int BHsize2 = BHsize_arr[1];

        //cout << hex << BHsize << " " << BHsize2 << dec << endl;

        offset +=4;
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
          evtfile.ignore(8);
          //TODO skipping the timestamp for now -- need to come back to this later Johnathan
          unsigned short eventsource_arr[4];
          evtfile.read((char*)eventsource_arr,8);
          eventsourceID = eventsource_arr[0];

          eventBarrierType = eventsource_arr[2];

          offset +=16;

          //cout << hex << eventsourceID << " " << eventBarrierType << " " << dec << endl;
          //cout << "offset " << offset << endl;
        }
        //else
        //{
          //Buffers with no body header, for use with NSCLDAQ older than 11.0
        //}

        int64_t fragmentTimestamp[5]={0};
        int32_t fragmentsourceID[5] ={-1};
        int32_t fragmentsize[5] ={0};
        int32_t fragmentBarrierType[5] = {0};

        int nFragment = 0;
        int fragmentcounter=0;
        int BuiltSize = 0;
        int BuiltSize2 = 0;

        bool foundS800 = false;
        bool foundSecondary = false;
        Det.Reset();

        while (offset != nbytes)
        {
          if (type == 1)
          {
            evtfile.ignore(4);
            offset+=4;
            unsigned short runno_arr[1];
            evtfile.read((char*)runno_arr,2);
            runno = runno_arr[0];   
            cout << "run number = " << runno << endl;
            //Ignore the rest of this event
            offset+=2;
            evtfile.ignore(nbytes - offset);
            offset = nbytes;
          }
          else if (type == 12)
          {
            //Don't do anything, already skipped 8 bytes
            //cout << "ring format" << endl; 

            //Skip 4 bytes
            evtfile.ignore(4);
            offset +=4;
          }
          else if (type == 30)
          {
            nFragment=fragmentcounter;
            fragmentcounter++;
            //Reads the body size, which is only present at the start of a built event
            if(nFragment ==0)
            {
              unsigned short BuiltSize_arr[2];
              evtfile.read((char*)BuiltSize_arr,4);
              BuiltSize = BuiltSize_arr[0];
              BuiltSize2 = BuiltSize_arr[1];
              offset +=4;

              //cout << hex << BuiltSize << " " << BuiltSize2 << dec << endl;
            }      
            int nwordsring = BuiltSize;    
            //Now to read the fragment header
            evtfile.ignore(8); //Skipping the time stamp for now JP
            unsigned int fragment_arr[3];
            evtfile.read((char*)fragment_arr,12);
            fragmentsourceID[nFragment]=fragment_arr[0];
            
            fragmentsize[nFragment] = fragment_arr[1];
            int fragsize = fragmentsize[nFragment];

            fragmentBarrierType[nFragment] = fragment_arr[2];
            offset += 20;

            //cout << "here " << hex << fragmentsourceID[nFragment] << " " << fragmentsize[nFragment] << " " << fragmentBarrierType[nFragment] << dec << endl;

            //Ignore repeated header
            evtfile.ignore(28);
            offset+=28;
            //cout << "Source ID " << fragmentsourceID[nFragment] << endl;
            //We may want to modify this so that it checks the position of point
            //when the data returns from the unpackers
            //I will think about implementing this at a later date
            //KB Oct 2018
            //offset+=fragmentsize[nFragment]/2; //Skipping the #words of the payload;
            //point +=14;

            if (physicsEvent%5000 == 0) 
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
            //TODO turn off unpacking till header works
            bool stat = Det.unpack(&evtfile,runno,fragmentsourceID[nFragment], fragsize);

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
              else if(fragmentsourceID[nFragment] == Det.JanusID)
              {
                foundSecondary = true;
                NGoodSecondary++;
                physicsEventGood++;
              }
            }
            //Let move to the next fragment
            //TODO this will only work if unpacking is off
            //evtfile.ignore(fragsize - 28);
            offset += fragsize - 28;

            //Turning off trees for now, I will return to this later, not needed for online
            //KB Oct 2018
            //Det.treeGrow(); 
    
          }//end type == 30
          else if (type == 31)
          {
            physicsEventCounter++;
            evtfile.ignore(nbytes - offset);
            offset = nbytes;
          }
          else if (type == 2)
          {
            endOfRun = true;
            break;
          }
          else if (type == 20)
          {
            scalerBuffer++;
            evtfile.ignore(nbytes - offset);
            offset = nbytes;
          }
          else if (type == 3)
          {
            Npauses++;
            evtfile.ignore(nbytes - offset);
            offset = nbytes;
          }
          else if (type == 4)
          {
            Nresumes++;
            evtfile.ignore(nbytes - offset);
            offset = nbytes;
          }
          else if (type == 42)
          {
            evtfile.ignore(nbytes - offset);
            offset = nbytes;
          }
          else
          {
            evtfile.ignore(nbytes - offset);
            offset = nbytes;
          }

          nFragment++;
        } //end loop over ring item




        if (foundS800 || foundSecondary) 
        {
          //here we read trigger coin and add
          //TODO turn off analysis until unpacking works
          Det.analyze(physicsEvent, runno);
        }

      } //end loop over items in a evtfile
      evtfile.close();
      evtfile.clear(); // clear event status in case we had a bad file
  
    } //end loop over file subsections
    //Forest->writeTree();
    runcount++;
  } //end loop of run file numbers
  
  cout << '\n'<<"physics Events = " << physicsEvent << endl;
  cout << "Good physics Events = " << physicsEventGood << endl;
  crossSections.close();
  
  if (physicsEvent > 0) cout << "bad/total = " << 
          (1.-(double)physicsEventGood/(double)physicsEvent)*100.<< " %"<< endl;
  
  cout << "physics Event Counters = " << physicsEventCounter << endl;
  cout << "scaler buffers = " << scalerBuffer << endl;
  cout << "Numbers of pauses = " << Npauses << endl;
  cout << "Number of resumes = " << Nresumes << endl;

  cout << "good residues = " << Det.Nresidue << endl;
  cout << "Bad residues = " << Det.Nbadresidue << endl;

  //cout << "Ar33_37Cabeam = " << Det.Ar33_37Cabeam << endl;
  //cout << "Ar33_36Kbeam = " << Det.Ar33_36Kbeam << endl;

  cout << "number of s800 singles = " << Det.N_s800_singles << endl;
  cout << "number of coincidences = " << Det.N_coin << endl;
  cout << "number of single event = " << Det.N_singles << endl;

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

  int roll = rand() % (20 - 1 + 1) + 1;
  string line,line2;
  line = R"(          _-_.      )";
  cout << line << endl;
  line = R"(     _-',^. `-_.    )";
  cout << line << endl;
  line = R"( ._-' ,'   `.   `-_ )";
  cout << line << endl;
  line = R"(!`-_._________`-':::)";
  cout << line << endl;
  line = R"(!   /\        /\::::)";
  cout << line << endl;
  if (roll >= 10)
    line = R"(;  /  \  )";
  else
    line = R"(;  /  \   )";
  line2 = R"(  /..\:::)";
  cout << line << roll << line2 << endl;
  line = R"(! /    \    /....\::)";
  cout << line << endl;
  line = R"(!/      \  /......\:)";
  cout << line << endl;
  line = R"(;--.___. \/_.__.--;;)";
  cout << line << endl;
  line = R"( '-_    `:!;;;;;;;' )";
  cout << line << endl;
  line = R"(    `-_, :!;;;''    )";
  cout << line << endl;
  line = R"(        `-!'        )";
  cout << line << endl;
  if (roll == 20) { cout << "Nat 20!" << endl;}
  else if (roll == 1) { cout << "Nat 1 :(" << endl;}

  return 0;

}
