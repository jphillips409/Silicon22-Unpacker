// Modified on 15 July by Henry Webb
// Modified on 12 September 24 by Johnathan Phillips. Removed read header function (not needed for merged buffer).

#include "CAENd5202.h"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <time.h>

// File created to unpack the CAEN DT5202 into a class that stores all
// of the data for each event. This version of the unpacker is designed
// to work only with timing-only events (acqMode==0x02)

eventTiming::eventTiming() {}

void eventTiming::clear() {
	pos = 0;
	ToA = -1;
	ToT = -1;
	ToTmatched = -1;
}

/**********************************************************************************/

string Event::Print(bool b = false) const {
  ostringstream oss;
  oss << "board:" << (int)boardID << "\t timeStamp:" << fixed << timeStamp << "\t NHits:" << NHits << endl;

  // column headers
  if (b)
    oss << "Ev# | channel | pos | ToA | ToT" << endl;

  // hit loop
  eventTiming ev;
  for (int i = 0; i < NHits; i++) {
    ev = dataTiming[i];
    oss << i << " " << ev.getChan() << " " << ev.pos << " " << ev.ToA << " " << ev.ToT << endl;
  }

  return oss.str();
}

Event::Event() {
  clear();
}

Event::Event(Event* rhs) {
	// Copy event header info
	boardID = rhs->boardID;
	timeStamp = rhs->timeStamp;
	NHits = rhs->NHits;	

	// Copy event data
	dataTiming = rhs->dataTiming;
}

//set_vals in need Big Endian style
void Event::set_short(unsigned short &t, char*& p)
{
  t = (*p++ << 8);
  t = t | *p++;
}
void Event::set_24bit(unsigned int &t, char*& p)
{
  t = (*p++ << 8);
  t = (t | *p++) << 8;
  t = t | *p++;
}

// Reads one event from the stream and saves it to the private variables
long Event::ReadEventFromStream(ifstream *pfs, double* redgains, double* bluegains)
{
	pfs->peek();
  //if (!pfs->good()) //Don't think I need this
  //{
    //cout << "not good" << endl;
    //return -1;
  //}

	// Get initial position
	std::streampos initialPos = pfs->tellg();

  // Peak at the first part to deterime how large of a buffer to create
  size_t peaksize = 2;
  char peaker[peaksize];
  pfs->read((char*)peaker, peaksize);
  char* pbuf = peaker;
	unsigned short eventSize;
  set_val(eventSize, pbuf);

  // Create the buffer (size 2 less because we already read the first part)
  char buf[eventSize-2];
  pfs->read((char*)buf, eventSize-2);
  pbuf = buf;

  set_val(boardID, pbuf);
  set_val(timeStamp, pbuf);
  set_val(NHits, pbuf);
  //cout << "here " << eventSize << " " << boardID << " " << timeStamp << " " << NHits << endl;
  eventTiming Ev;
	unsigned char type; // 0x10, if only the ToA value is saved for that channel; 0x20, if only the ToT value is saved for that channel; 0x30, if both ToA and ToT values are saved
  for (int n=0; n<NHits; n++) {
    Ev.clear();
    
    set_val(Ev.chan, pbuf);
		Ev.pos = (((unsigned int)Ev.chan - ((unsigned int)Ev.chan % 2)) / 2) + (((unsigned int)Ev.chan % 2) * 32);
    set_val(type, pbuf);
    if (type == 0x10 || type == 0x30) set_val(Ev.ToA, pbuf);
    if (type == 0x20 || type == 0x30) set_val(Ev.ToT, pbuf);
		Ev.ToTmatched = ((double)Ev.ToT) * (((boardID == 0) * redgains[Ev.pos]) + ((boardID == 1) * bluegains[Ev.pos]));
    dataTiming.push_back(Ev);
  }

	// Get final position
  std::streampos finalPos = pfs->tellg();
  //cout << finalPos - initialPos << endl;
	// Return # of bytes read
  return long(finalPos - initialPos);
}

void Event::clear()
{
  // Event Header (Timing Mode)
  boardID = 0;
  timeStamp = 0;
  NHits = 0; // Number of recorded hits
  dataTiming.clear();
}

eventTiming Event::FindMax(function<bool(eventTiming, eventTiming)> comp) {
	return *std::max_element(dataTiming.begin(), dataTiming.end(), comp);
}



