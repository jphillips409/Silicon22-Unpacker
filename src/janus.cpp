// in the janus class, the data file is stored into vectors of events.
// the event class is detailed in eventCAEN.cpp
// created 12/1/2023 by Charlie Fallon
// modified 7/15/2024 by Henry Webb
// modified 9/10/24 for the Si22 code by Johnathan Phillips. Removed evt matching function, evts already merged.
//                                                           Added analyze() function

#include "janus.h"

#include "CAENd5202.h"
#include "fiber.h"
#include "histo_sort.h"

#include <iostream>
#include <ios>

//TODO do I need this config path? Seems needed for gain matching
#define CONFIGPATH "cal/fibers/"

//TODO need to link with correct histos
// Constructor
janus::janus(histo_sort * Histo1, double d) {
	Histo = Histo1;
	SIPMevent = new Event();
	Fiber = new fiber();
	distance = d;
  dist_pub = distance;

  //TODO need correct gain matching and directories
	ReadGains(string(CONFIGPATH) + "blue_gain_matching.txt", bluegains);
	ReadGains(string(CONFIGPATH) + "red_gain_matching.txt", redgains);
}

// Destructor
janus::~janus() {}

void janus::SetTarget(double Targetdist)
{
  //TODO Target position will need to change with S800 setting
	distance = Targetdist;
  dist_pub = distance;
}


// Read scaling values for gain matching from file
void janus::ReadGains(string ifname, double* arr) {
	ifstream ifile;
	ifile.open(ifname, ios::in);
	if (!ifile.is_open())
		throw invalid_argument("Supplied input file does not open properly");
	
	double data;
	for (int i = 0; i < 64; i++) {
		ifile >> data;
		if (ifile.eof())
			throw invalid_argument("Supplied input file shorter than expected length");
		arr[i] = data;
	}
}

// Unpack class handles the opened data file, unpacks each event
bool janus::unpack(ifstream *pevtfile) { 
	nevts = 0;
  long nbytes = 0; //different from nbytes 1 & 2
	int nbytes1 = 0;
  int nbytes2 = 0;
  int nbytes_temp = 0; //Need this to check for bad boards
  int tunit = 0; // time unit
  int acqmode = 0;
  int evtsize1 = 0;
  int evtsize2 = 0;
  int boardnum = -1;
  

	// Read event header
	// THIS MUST BE DONE ONCE BEFORE READING AN EVENT!
  unsigned char hbuffer[4];
  pevtfile->read((char*)hbuffer,4);
  nbytes1 = hbuffer[0]; // inclusive size of the event
  nbytes2 = hbuffer[1]; // inclusive size of the event
  //cout << "nbytes " << nbytes1 << " " << nbytes2 <<  " " << (nbytes2 << 8) << endl;
  int inclsize = 0;
  inclsize = nbytes1;
  inclsize += (nbytes2 << 8);
  //cout << hex << "inclsize " << inclsize << dec << endl;
  //cout << "inclsize " << inclsize << endl;
  //if (nbytes2 != 0) nbytes1 += nbytes2;
  //if (nbytes2 != 0) n
  tunit = hbuffer[2];
  acqmode = hbuffer[3];
 // cout << "t unit and acq mode " << tunit << " " <<acqmode << endl;

  int offset = 4; //tracking the event length in bytes

	// Declare common variables
	eventTiming hit;
	float tot;
	float toa;
	unsigned int pos;

	// Event loop (timing-only mode)
	bool val;
	
	
	nbytes = SIPMevent->ReadEventFromStream(pevtfile, redgains, bluegains); // reads next event
	//cout << "here nbytes from SIPM evt " << nbytes << endl;

  //cout << "Janus Timestamp " << SIPMevent->GetTimeStamp() << endl;
	//if (nbytes == -1) break; // stop at end of file

	//Event SIPMeventcur(SIPMevent);
	//cout << SIPMeventcur->Print(true) << endl;

	// Loop through hits in event
	unsigned char boardID = SIPMevent->GetBoardID();
  //cout << "boardID " << boardID*1 << endl;
	int nhits = (int)SIPMevent->GetNHits();
	for (int i = 0; i < nhits; i++) {
    if (boardID > 1)
    {
      cout << "BoardID > 1" << endl;
      break;
    }
		hit = SIPMevent->GetTimingEvent(i);
		tot = hit.ToTmatched;
		toa = hit.ToA;
		pos = hit.pos;
		// Fill histograms
		if (toa > -1)
			Histo->toa_hist->Fill(toa);
		if (tot > -1 && boardID == 0)
			Histo->tot_summary_red->Fill(pos, tot);
		else if (tot > -1 && boardID == 1)
			Histo->tot_summary_blue->Fill(pos, tot);
    
    //Fill 1d raw hists
    if (boardID == 0)
    {
      Histo->B0TOA_R[pos]->Fill(toa);
      Histo->B0TOT_R[pos]->Fill(tot);
    }
    if (boardID == 1)
    {
      Histo->B1TOA_R[pos]->Fill(toa);
      Histo->B1TOT_R[pos]->Fill(tot);
    }
  }
		
    //Fill the event 
  if (boardID < 2) janusevts.push_back(*SIPMevent);
	SIPMevent->clear(); //TODO does this need to be clear for each hit or for each event? I think event
		// Fill output tree
		//Histo->FillTree(*SIPMeventcur);

    //JANUS EVENTS NOW HAVE 7 FFs, changed right before Si22 experiment
    //If we hit event end, break, account for 4 FFs after event ******CHANGED*********
    //cout << dec << "nbytes " << nbytes << " offset " << offset << " inclsize " << inclsize << endl;
  //cout << "here2" << endl;
  if (nbytes !=  (inclsize - 7 - offset)) //Sometimes the trailer FFs won't be 7, skip for now and advance ifstream
  {
    pevtfile->ignore(inclsize - (nbytes + offset));
    //return false; //Works as long as you get the correct nbytes each time
    return true;
  }
  //cout << "here3" << endl;
  
  pevtfile->ignore(7);
  //cout << "END OF JANUS EVT" << endl;
	return true;
}

void janus::analyze(const s800_results& S800_results)
{

  fibmatch = false;

  if (janusevts.size() != 2 || janusevts.size() > 2)
  {
    if (janusevts.size() > 0)
    {
      if (janusevts[0].GetBoardID() == 0) evt_RNOTB++;    
      if (janusevts[0].GetBoardID() == 1) evt_BNOTR++;    
    }    
    return;
  }

  if (janusevts[0].GetBoardID() == janusevts[1].GetBoardID())
  { 
    cout << "SAME BOARD!!!" << endl;
    abort();
  }
  

  //cout << "Jan evt size " << janusevts.size() << endl;
  int idhorz = -1;
  int idvert = -1;
	eventTiming ev;

  //Forget why this is needed. Changes if board 1 or board 0 is the vert or horz?
  if (janusevts[0].GetBoardID() == 0)
  {
    idhorz = 1;
    idvert = 0;
  }
  else
  {
    idhorz = 0;
    idvert = 1;
  }

  Fiber->make_2d(&(janusevts[idhorz]), &(janusevts[idvert]), distance, S800_results);

  //Check for bad max fibers
  if (Fiber->badtx == true || Fiber->badty == true) return;

  fibmatch = true;

  evt_BRMatch++;

	//	Write histograms and tree here
  //Histo->FillMatchedTree(Fiber, redbuffevents[i], bluebuffevents[j]);
  
  //vert vs horz
  //The x value comes from the vertical fiber. So does vert mean y value (horz fiber) or the ver fiber (x value)
  
	Histo->Fiber_ixiy->Fill(Fiber->ix, Fiber->iy);
	Histo->Fiber_xy->Fill(Fiber->x, Fiber->y);

	ev = janusevts[idhorz].GetTimingEvent(Fiber->posmaxhorz);
  //if (janusevts[idhorz].GetBoardID() == 1) cout << "test here!!!!!!! " << endl;
	Histo->tot_summary_blue_matched->Fill(ev.pos, ev.ToTmatched);
	Histo->Fiber_toax->Fill(ev.ToA);

	ev = janusevts[idvert].GetTimingEvent(Fiber->posmaxvert);
	Histo->tot_summary_red_matched->Fill(ev.pos, ev.ToTmatched);
	Histo->Fiber_toay->Fill(ev.ToA);

  // Plot hit map for individual events
	double temppos;
	for (int k = 0; k<janusevts[idhorz].GetNHits(); k++)
  {
		ev = janusevts[idhorz].GetTimingEvent(k);
		temppos = -1*(ev.pos-0.5)*0.5 + 16; //mm
		Histo->Fiber_totx->AddBinContent(Histo->Fiber_totx->GetXaxis()->FindBin(temppos), ev.ToTmatched);
		Histo->Fiber_postotx->AddBinContent(Histo->Fiber_postotx->GetXaxis()->FindBin(ev.pos), ev.ToTmatched);
		Histo->Fiber_postoax->AddBinContent(Histo->Fiber_postoax->GetXaxis()->FindBin(ev.pos), ev.ToA);
	 }

	for (int k = 0; k<janusevts[idvert].GetNHits(); k++)
  {
	  ev = janusevts[idvert].GetTimingEvent(k);
		temppos = -1*(ev.pos-0.5)*0.5 + 16; //mm
		Histo->Fiber_toty->AddBinContent(Histo->Fiber_toty->GetXaxis()->FindBin(temppos), ev.ToTmatched);
		Histo->Fiber_postoty->AddBinContent(Histo->Fiber_postoty->GetXaxis()->FindBin(ev.pos), ev.ToTmatched);
		Histo->Fiber_postoay->AddBinContent(Histo->Fiber_postoay->GetXaxis()->FindBin(ev.pos), ev.ToA);
	}

  return;

}

