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
  
  int evt_counter = 0;
  int physicsEvent = 0;
  int physicsEventGood = 0;
  int physicsEventCounter = 0;
  int scalerBuffer = 0;
  int S800frag = 0;
  int VMEfrag = 0;
  int Janusfrag = 0;
  int onefrag = 0;
  int twofrag = 0;
  int threefrag = 0;
  int Npauses = 0;
  int Nresumes = 0;
  int runno = 0;
  int NS800 = 0;
  int NGoodS800 = 0;
  int NSecondary = 0;
  int NGoodSecondary = 0;
  int totalRuntime = 0;

  int S800multevents = 0;
  int Janusmultevents = 0;
  int VMEmultevents = 0;
  
  //Breaks down S800 coincidence evts into different parts
  int S800Coinc_evt = 0; //Total number of S800 coincidence events
  int S800Coinc_VME = 0; //S800 with VME, no fibers
  int S800Coinc_VME_Si = 0; //S800 with VME, no fibers, only silicons
  int S800Coinc_VME_CsI = 0; //S800 with VME, no fibers, only CsI
  int S800Coinc_VME_Both = 0; //S800 with VME and Si/CsI, not solution gates
  int S800Coinc_VME_Complete = 0; //S800 with VME, no fibers, Si and CsI
  int S800Coinc_1Fiber = 0; //S800 with 1 fiber
  int S800Coinc_2Fiber = 0; //S800 with 2 fibers, not matched
  int S800Coinc_2plusFiber = 0; //S800 with more than 2 boards, garbage
  int S800Coinc_Fibermatch = 0; //S800 with 2 matched fibers
  int S800Coinc_VME_1Fiber = 0; //S800, VME, 1 fiber
  int S800Coinc_VME_2Fiber = 0; //S800 VME, 2 fibers, not matched
  int S800Coinc_VME_2plusFiber = 0; //S800 with VME and more than 2 boards, garbage
  // Full event includes S800, VME, and matched fibers
  int S800Coinc_Complete = 0;
  int S800Coinc_Complete_Si23 = 0;
  int S800Coinc_Complete_Si23_Mg20_2p = 0;
  int S800Coinc_Complete_Si23_Si22_p = 0;
  int S800Coinc_Complete_Si23_Si23_p = 0;

  int S800Coinc_Complete_Si23_Ne17_3p = 0;
  int S800Coinc_Complete_Si23_Ne18_3p = 0;
  
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

  int setting = 1; //Default value of setting 1

  //For switching the distance of Gobbi to target for setting 1
  bool setting1_og = true;
  bool setting1_a = false;

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

      string numbersname = argv[argcounter];

      //Choose setting based on file
      if (numbersname.compare("numbers.set1") == 0) setting = 1;
      if (numbersname.compare("numbers.set2") == 0) setting = 2;
      if (numbersname.compare("numbers.set3") == 0) setting = 3;

    }     
    firstrunnum = atoi(argv[argcounter]);
  } 


  //Print out setting
  cout << "S800 Setting: " << setting << endl;


  TargThick->SetThick(runnum);

  det Det(Histo_sort, Histo_read, setting);

  //Set the sourceID for the S800 (normally 2) and secondary DAQ
  //TODO need correct ID's for HINP, CAESAR, S800, and Janus
  Det.SiID = 1;
  Det.JanusID = 3;
  Det.S800ID = 2;

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

      //outstring << directory << "/user/e21006/stagearea/merged/experiment/run" << number <<"/run-" << setfill('0') << setw(4) << number;
      //outstring << directory << "/user/e21006/stagearea/standalone/experiment/run" << number <<"/run-" << setfill('0') << setw(4) << number;
      //outstring << directory << "/mnt/rawdata/e21006/merged/experiment/run" << number <<"/run-" << setfill('0') << setw(4) << number;
//        outstring << directory << "/run-" << setfill('0') << setw(4) << number;

      //For data on Neutron
      outstring << directory << "/data2/Si22_nov2024/run" << number <<"/run-" << setfill('0') << setw(4) << number;

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
      float SiDistance = 55.2;//*0.99;
      float TargetThick = 0;
      if (setting == 0) //Gobbi further from target before 216, setting 1
      {
        TargetThick = 97.77768; 
        SiDistance = 55.2;
      }
      if (setting == 1 && number < 216) //Gobbi further from target before 216, setting 1
      {
        TargetThick = 97.77768;//99.79; 
        SiDistance = 55.383;//55.2;
      }
      if (setting == 1 && number >= 216) //Gobbi closer to target after 216, setting 1a
      {
        TargetThick = 97.77768; 
        SiDistance = 43.896;
      }
      if (setting == 2)
      {
        TargetThick = 99.79*2; 
        SiDistance = 34.8;
      }
      if (setting == 3)
      {
        TargetThick = 99.79*4; 
        SiDistance = 50.;
      }

			float FiberDistance = SiDistance + 10.557;//10.151+2.664; //Add some constant distance for Fibers, cm. Approximately assume that Si start at Gobbi cart.
      //99.79 mg/cm2 -- 0.54 mm target
      //91.18456 from 0.497 mm thick Be target using rigidity measurement
      //188 mg/cm2 -- 1mm target
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

      //print out distances for each run
      cout << "Run " << number << " Si dist (cm): " << SiDistance << " Fib dist (cm): " << FiberDistance << endl; 

      //TODO verify that targ dist and thick work. Also set fiber dist
      Det.Gobbi->SetTarget(SiDistance,TargetThick); //Uses cm
			Det.Janus->SetTarget(FiberDistance); //Janus uses cm for angle

      //Det.Hira->XY_mon->setDistance(SiDistance);
      
      ////////////////////////////////////////////////////////////////////////////////////////
  
      //Performs check of the target, aborts if not good
      if (TargetThick <= 0)
      {
        cout << "Target must be > 0 " << endl;
        abort();
      }

      //Performs check of Janus target, aborts if not good
      if (Det.Janus->dist_pub <= 0)
      {
        cout << "Janus target must be > 0 " << endl;
        abort();
      }
  
      //for(int i=0;i<500000;i++)  // loop over items in a evtfile
      for(;;)      
      {
      evt_counter++;
      //cout << "EVT START" << endl;
      //cout << evt_counter << endl;
      
      int S800_here = 0;
      int VME_here = 0;
      int Janus_here = 0;

      int S800_mult = 0;
      int Janus_mult = 0;
      int VME_mult = 0;

      unsigned long long FragHeadTStamp = 0;
      int64_t S800TStamp = 0;
      int64_t VMETStamp = 0;
      int64_t JanusTStamp = 0;

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

        unsigned short nbytes1 = hBuffer[0];
        unsigned short nbytes2 = hBuffer[1];
        unsigned short type1 = hBuffer[2];
        unsigned short type2 = hBuffer[3];

        //cout << "event header " << hex << nbytes1 << " " << nbytes2 << " " << type1 << " " << type2 << " " << dec << evt_counter << endl;
       
        unsigned int nbytes = nbytes1 + (nbytes2 << 16);
        
        unsigned int type = type1 + (type2 << 16);

         
        //Skip massive events, likely smashing multiple events together
        //if (nbytes > 4000)
        //{
        //  cout << "Skip massive event" << endl;
        //  evtfile.ignore(nbytes - hBufferBytes);
        //  continue;
        //}
       
        int offset = 0;
        offset = 8;
 
        if (type == 30) physicsEventCounter++;
 
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
          //int64_t evtTstamp_arr[1];
          //evtfile.read((char*)evtTstamp_arr,16);
          //eventTstamp = evtTstamp_arr[0];
          
          //TODO skipping the timestamp for now -- need to come back to this later Johnathan
          unsigned short eventsource_arr[4];
          evtfile.read((char*)eventsource_arr,8);
          eventsourceID = eventsource_arr[2];

          eventBarrierType = eventsource_arr[0];

          offset +=16;

          //cout << hex << eventsourceID << " " << eventBarrierType << " " << dec << endl;
          //cout << "offset " << offset << endl;
        }
        //else
        //{
          //Buffers with no body header, for use with NSCLDAQ older than 11.0
        //}

        unsigned short fragmentTimestamp[4]={0};
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
        //if (evt_counter > 3 && evt_counter < 3000)
        //{
          //evtfile.ignore(nbytes-offset);
          //continue;
        //} 
        while (offset != nbytes)
        {
          if (type == 1)
          {
            unsigned short runno_arr[1];
            evtfile.read((char*)runno_arr,2);
            runno = runno_arr[0]; 
            offset+=2;

            evtfile.ignore(4);
            offset+=4;
            cout << "run number = " << runno*1 << endl;
            //Ignore the rest of this event
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
            evtfile.read((char*)fragmentTimestamp,8);
            //for (int j=0;j<4;j++) cout << hex << fragmentTimestamp[j] << " ";
            //cout << endl;
            FragHeadTStamp = (long long)fragmentTimestamp[0] + ((long long)fragmentTimestamp[1] << 16) + ((long long)fragmentTimestamp[2] << 32) + ((long long)fragmentTimestamp[3] << 48);
            //cout << "Tstamp parts " << hex << fragmentTimestamp[0] << " " << fragmentTimestamp[1] << " " << fragmentTimestamp[2] << " " << fragmentTimestamp[3] << dec << endl;
            //cout << "FratHeadTStamp " << dec << FragHeadTStamp << endl;

            //evtfile.ignore(8); //Skipping the time stamp for now JP
            
            unsigned int fragment_arr[3];
            evtfile.read((char*)fragment_arr,12);
            fragmentsourceID[nFragment]=fragment_arr[0];
            
            fragmentsize[nFragment] = fragment_arr[1];
            int fragsize = fragmentsize[nFragment];

            fragmentBarrierType[nFragment] = fragment_arr[2];
            offset += 20;

           // cout << "here " << hex << fragmentsourceID[nFragment] << " " << fragmentsize[nFragment] << " " << fragmentBarrierType[nFragment] << dec << endl;

            //if (fragmentBarrierType[nFragment] != 0 || fragmentBarrierType[nFragment] != 0)
            //{
              //evtfile.ignore(fragsize+4);
              //offset+=fragsize+4;
              //continue;
            //}

            if (FragHeadTStamp > 104467440700)
            {
              cout << "Crazy large timestamp - skip " << endl;
              cout << "Timestamp: " << FragHeadTStamp << endl;
            }

            if (FragHeadTStamp == 0 || FragHeadTStamp < 100) //Also skip tiny timestamps that don't make sense
            {
              cout << "Janus delay causes strange first event with 0 t stamp - skip" << endl;
              evtfile.ignore(fragsize);
              offset += fragsize;
              continue;
            } 

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

            //if (physicsEvent%100000 == 0) 
              //cout << '\xd'<< physicsEvent << flush;
            physicsEvent++;
            
            //Det.S800ID == 2, Det.SiID == 1
            if(fragmentsourceID[nFragment] == Det.S800ID)
            {
              NS800++;
              S800frag++;
              S800_here = 1;
              S800TStamp = FragHeadTStamp;
              S800_mult++;
            }
            else if(fragmentsourceID[nFragment] == Det.SiID)
            {
              NSecondary++;
              VMEfrag++;
              VME_here = 1;
              VMETStamp = FragHeadTStamp;
              VME_mult++;
            }
            else if(fragmentsourceID[nFragment] == Det.JanusID)
            {
              if (Janus_here == 0) Janusfrag++;
              NSecondary++;
              Janus_here = 1;
              JanusTStamp = FragHeadTStamp;
              Janus_mult++;
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
            evtfile.ignore(nbytes - offset);
            offset = nbytes;
          }
          else if (type == 2)
          {
            endOfRun = true;
            cout << "End of run block" << endl;
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

        if (evt_counter%100000 == 0) 
          cout << '\xd'<< evt_counter << flush;

        //Different if statements
        //If sorting real runs, only analyze those with S800, VME, and 2 Janus boards
        //Must be turned off for calibrations to work
        //if (foundS800 || foundSecondary) 
        if (foundS800 && foundSecondary && Janus_mult == 2)
        {
          //here we read trigger coin and add
          //TODO turn off analysis until unpacking works
          Det.analyze(physicsEvent, runno);
        }
        else { //Make sure Janus gets cleared
          Det.Janus->janusevts.clear();
	        Det.Janus->Histo->clear();
        }
        
        
        if (S800_here == 0 && Janus_here == 0 && VME_here == 1)
        {
          onefrag += 1;
        }
        if (S800_here == 0 && Janus_here == 1 && VME_here == 0)
        {
          onefrag += 1;
        }
        if (S800_here == 1 && Janus_here == 0 && VME_here == 0)
        {
          onefrag += 1;
        }
        if (S800_here == 0 && Janus_here == 1 && VME_here == 1)
        {
          twofrag += 1;
        }
        if (S800_here == 1 && VME_here == 1 && Janus_here == 0)
        {
          twofrag += 1;
          if (Det.s800->Trig.registr & 2) S800Coinc_VME++;
          if (Det.Gobbi->S800coinc_CsI == true && Det.s800->Trig.registr & 2) S800Coinc_VME_CsI++;
          if (Det.Gobbi->S800coinc_Si == true && Det.s800->Trig.registr & 2) S800Coinc_VME_Si++;
          if (Det.Gobbi->S800coinc_both == true && Det.s800->Trig.registr & 2) S800Coinc_VME_Both++;
          if (Det.Gobbi->S800_complete == true && Det.s800->Trig.registr & 2) S800Coinc_VME_Complete++;
        }
        if (S800_here == 1 && Janus_here == 1 && VME_here == 0)
        {
          twofrag += 1;
          if (Janus_mult == 1 && Det.s800->Trig.registr & 2) S800Coinc_1Fiber++;
          if (Janus_mult == 2 && Det.Janus->fibmatch == false && Det.s800->Trig.registr & 2) S800Coinc_2Fiber++;
          if (Janus_mult > 2 && Det.s800->Trig.registr & 2) S800Coinc_2plusFiber++;
          if (Janus_mult == 2 && Det.Janus->fibmatch == true && Det.s800->Trig.registr & 2) S800Coinc_Fibermatch++;
        }
        if (S800_here == 1 && Janus_here == 1 && VME_here == 1)
        {
          threefrag += 1;
          if (Janus_mult == 1 && Det.s800->Trig.registr & 2) S800Coinc_VME_1Fiber++;
          if (Janus_mult == 2 && Det.Janus->fibmatch == false && Det.s800->Trig.registr & 2) S800Coinc_VME_2Fiber++;
          if (Janus_mult > 2 && Det.s800->Trig.registr & 2) S800Coinc_VME_2plusFiber++;

          if (Det.Gobbi->S800_complete == true  && Det.s800->Trig.registr & 2 && Det.Janus->fibmatch == true) S800Coinc_Complete++;
          if (Det.Gobbi->S800_complete == true  && Det.s800->Trig.registr & 2 && Det.Janus->fibmatch == true && Det.s800->beam_pid->Z == 14 && Det.s800->beam_pid->A == 23) S800Coinc_Complete_Si23++;

          if (Det.Gobbi->S800_complete == true  && Det.s800->Trig.registr & 2 && Det.Janus->fibmatch == true && Det.Si23_Mg20_2p_flag == true) S800Coinc_Complete_Si23_Mg20_2p++;

          if (Det.Gobbi->S800_complete == true  && Det.s800->Trig.registr & 2 && Det.Janus->fibmatch == true && Det.Si23_Si22_p_flag == true) S800Coinc_Complete_Si23_Si22_p++;


	        if (Det.Gobbi->S800_complete == true  && Det.s800->Trig.registr & 2 && Det.Janus->fibmatch == true && Det.Si23_Si23_p_flag == true) S800Coinc_Complete_Si23_Si23_p++;

	        if (Det.Gobbi->S800_complete == true  && Det.s800->Trig.registr & 2 && Det.Janus->fibmatch == true && Det.Ne17_3p_flag == true) S800Coinc_Complete_Si23_Ne17_3p++;

	        if (Det.Gobbi->S800_complete == true  && Det.s800->Trig.registr & 2 && Det.Janus->fibmatch == true && Det.Ne18_3p_flag == true) S800Coinc_Complete_Si23_Ne18_3p++;

        }


        
        //Plots the deltaT for each fragment relative to another
        if (VMETStamp != 0 && JanusTStamp != 0) Histo_sort->DeltaT_VME_Janus->Fill(VMETStamp - JanusTStamp);
        if (VMETStamp != 0 && S800TStamp != 0) Histo_sort->DeltaT_VME_S800->Fill(VMETStamp - S800TStamp);
        if (S800TStamp != 0 && JanusTStamp != 0) Histo_sort->DeltaT_S800_Janus->Fill(S800TStamp - JanusTStamp);

        if (VMETStamp != 0 && JanusTStamp != 0) Histo_sort->DeltaT_VME_Janus_vsEvtCnt->Fill(physicsEventCounter,VMETStamp - JanusTStamp);
        if (VMETStamp != 0 && S800TStamp != 0) Histo_sort->DeltaT_VME_S800_vsEvtCnt->Fill(physicsEventCounter,VMETStamp - S800TStamp);
        if (S800TStamp != 0 && JanusTStamp != 0) Histo_sort->DeltaT_S800_Janus_vsEvtCnt->Fill(physicsEventCounter,S800TStamp - JanusTStamp);

        if (S800_mult > 1) S800multevents++;
        if (VME_mult > 1) VMEmultevents++;
        if (Janus_mult > 1) Janusmultevents++;
        

        //Break down the coincidence events into source IDs and other efficiencies
        if (Det.s800->Trig.registr & 2)
        {
          S800Coinc_evt++;
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

  //cout << "number of s800 singles = " << Det.N_s800_singles << endl;
  //cout << "number of coincidences = " << Det.N_coin << endl;
  //cout << "number of single event = " << Det.N_singles << endl;


  cout << "Number of different sources" << endl;
  cout << "Num S800 " << S800frag << endl;
  cout << "Num VME " << VMEfrag << endl;
  cout << "Num Janus " << Janusfrag << endl;
  cout << "Number of sources in events " << "    Total events    " << evt_counter << endl;
  cout << "Num one source " << onefrag << endl;
  cout << "Num two source " << twofrag << endl;
  cout << "Num three source " << threefrag << endl;
  cout << "Add three sources " << onefrag + twofrag + threefrag << endl;

  cout << "Total SRS events: " << Det.s800->SRStot << endl;
  cout << "Good SRS events: " << Det.s800->SRSgood << endl;
  cout << "Bad SRS events: " << Det.s800->SRSbad << endl;

  cout << "Total num of S800 " << NS800 << endl;
  cout << "Num of huge S800 events: " << Det.s800->S800Huge << endl;
  cout << "Total num of protons " << Det.Gobbi->num_protons << endl;

  cout << "Matched Blue AND Red Janus: " << Det.Janus->evt_BRMatch << "  %" << 100*(float)Det.Janus->evt_BRMatch / (float)Janusfrag<< endl;
  cout << "Janus Blue NOT Red: " << Det.Janus->evt_BNOTR  << "  %" << 100*(float)Det.Janus->evt_BNOTR / (float)Janusfrag << endl;
  cout << "Matched Red NOT Blue Janus: " << Det.Janus->evt_RNOTB << "  %" << 100*(float)Det.Janus->evt_RNOTB / (float)Janusfrag << endl;

  cout << "Num events with mult S800 " << S800multevents << endl;
  cout << "Num events with mult VME " << VMEmultevents << endl;
  cout << "Num events with mult Janus " << Janusmultevents << endl;

  cout << "***************************************************************" << endl;
  cout << "               S800 coincidence event breakdown                " << endl;
  cout << "Total S800 coincidence evts: " << S800Coinc_evt << endl;
  cout << "S800 coincidence with just VME total: " << S800Coinc_VME << "   %   " << (float)S800Coinc_VME/S800Coinc_evt << endl;
  cout << "S800 coincidence with just VME Si: " << S800Coinc_VME_Si << "   %   " << (float)S800Coinc_VME_Si/S800Coinc_evt << endl;
  cout << "S800 coincidence with just VME CsI: " << S800Coinc_VME_CsI << "   %   " << (float)S800Coinc_VME_CsI/S800Coinc_evt << endl;
  cout << "S800 coincidence with just VME and Si/CsI: " << S800Coinc_VME_Both << "   %   " << (float)S800Coinc_VME_Both/S800Coinc_evt << endl;
  cout << "S800 coincidence with just complete VME: " << S800Coinc_VME_Complete << "   %   " << (float)S800Coinc_VME_Complete/S800Coinc_evt << endl;
  cout << "S800 coincidence with just 1 fiber board: " << S800Coinc_1Fiber << "   %   " << (float)S800Coinc_1Fiber/S800Coinc_evt << endl;
  cout << "S800 coincidence with just 2, unmatched fiber boards: " << S800Coinc_2Fiber << "   %   " << (float)S800Coinc_2Fiber/S800Coinc_evt << endl;
  cout << "S800 coincidence with only more than 2 fiber boards (junk): " << S800Coinc_2plusFiber << "   %   " << (float)S800Coinc_2plusFiber/S800Coinc_evt << endl;
  cout << "S800 coincidence with just matched fiber boards: " << S800Coinc_Fibermatch << "   %   " << (float)S800Coinc_Fibermatch/S800Coinc_evt << endl;
  cout << "S800 coincidence with VME and 1 fiber board: " << S800Coinc_VME_1Fiber << "   %   " << (float)S800Coinc_VME_1Fiber/S800Coinc_evt << endl;
  cout << "S800 coincidence with VME and 2 fiber boards, unmatched: " << S800Coinc_VME_2Fiber << "   %   " << (float)S800Coinc_VME_2Fiber/S800Coinc_evt << endl;
  cout << "S800 coincidence with VME and 3 fiber boards (junk): " << S800Coinc_VME_2plusFiber << "   %   " << (float)S800Coinc_VME_2plusFiber/S800Coinc_evt << endl;
  cout << "Complete S800 coincidence events : " << S800Coinc_Complete << "   %   " << (float)S800Coinc_Complete/S800Coinc_evt << endl;
  cout << "Complete S800 coincidence events with Si-23 : " << S800Coinc_Complete_Si23 << "   %   " << (float)S800Coinc_Complete_Si23/S800Coinc_evt << endl;
  cout << "Total number of Si23 beam : " << Det.N_coin_Si23 << "   %   " << endl;
  cout << "Total number of Si23 beam + 1p : " << Det.N_coin_Si23_1p << "   %   " << endl;
  cout << "Total number of Si23 beam + 2p : " << Det.N_coin_Si23_2p << "   %   " << endl;
  cout << "Total number of Si23 beam + 3p : " << Det.N_coin_Si23_3p << "   %   " << endl;
  cout << "Total number of Si23 beam + Mg20 : " << Det.N_coin_Si23_Mg20 << "   %   " << endl;
  cout << "Total number of Si23 beam + 1p + Mg20 : " << Det.N_coin_Si23_Mg20_1p << "   %   " << endl;
  cout << "Total number of Si23 beam + 2p + Mg20 : " << Det.N_coin_Si23_Mg20_2p << "   %   " << endl;
  cout << "Total number of Si23 beam + 2p + Mg20 + Janus : " << S800Coinc_Complete_Si23_Mg20_2p << "   %   " << (float)S800Coinc_Complete_Si23_Mg20_2p/S800Coinc_evt << endl; 
  cout << "Total number of Si23 beam + p + Si22 : " << Det.N_coin_Si22_1p << "   %   " << (float)Det.N_coin_Si22_1p/S800Coinc_evt << endl;
  cout << "Total number of Si23 beam + p + Si22 + Janus : " << S800Coinc_Complete_Si23_Si22_p << "   %   " << (float)S800Coinc_Complete_Si23_Si22_p/S800Coinc_evt << endl;
    cout << "Total number of Si23 beam + p + Si23 + Janus : " << S800Coinc_Complete_Si23_Si23_p << "   %   " << (float)S800Coinc_Complete_Si23_Si23_p/S800Coinc_evt << endl; 
  cout << "Total number of Si23 beam + 3p + Ne17 : " << Det.N_coin_Ne17_3p << "   %   " << (float)Det.N_coin_Ne17_3p/S800Coinc_evt << endl;
  cout << "Total number of Si23 beam + 3p + Ne17 + Janus : " << S800Coinc_Complete_Si23_Ne17_3p << "   %   " << (float)S800Coinc_Complete_Si23_Ne17_3p/S800Coinc_evt << endl;
  cout << "Total number of Si23 beam + 3p + Ne18 : " << Det.N_coin_Ne18_3p << "   %   " << (float)Det.N_coin_Ne18_3p/S800Coinc_evt << endl;
  cout << "Total number of Si23 beam + 3p + Ne18 + Janus : " << S800Coinc_Complete_Si23_Ne18_3p << "   %   " << (float)S800Coinc_Complete_Si23_Ne18_3p/S800Coinc_evt << endl;
  cout << "***************************************************************" << endl;

  
  //Histo_sort->tree->Fill();

  Histo_sort->write(); // this forces the histrograms to be read out to file
  Histo_read->write();

  
  t = clock() - t;
  cout << "run time: " << (float)t/CLOCKS_PER_SEC/60 << " min" << endl;

  /*int roll = rand() % (20 - 1 + 1) + 1;
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
  else if (roll == 1) { cout << "Nat 1 :(" << endl;}*/

	cout <<    "             ,,__					"	<< endl;
	cout <<    "    ..  ..   / o._) 	  ____       _     _     _      "	<< endl;                 
	cout <<    "   /--'/--\\  \\-'||  	 / ___| ___ | |__ | |__ (_)     "	<< endl;            
	cout <<    "  /        \\_/ / |  	| |  _ / _ \\|  _ \\|  _ \\| |  "	<< endl;            
	cout <<    ".'\\  \\__\\  __.'.' 	| |_| | (_) | |_) | |_) | |     "	<< endl;       
	cout <<    "  )\\ |  )\\ |        	 \\____|\\___/|____/|____/|_|   "	<< endl;
	cout <<    " // \\\\ // \\\\					"	<< endl;
	cout <<    "||_  \\\\|_  \\\\_      	"	<< endl;       
	cout <<    "'--' '--'' '--'						"	<< endl;

  return 0;

}
