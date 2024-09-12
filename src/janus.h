#ifndef janus_
#define janus_

#include <fstream>
#include <string>
#include <vector>
#include "histo_sort.h"
#include "histo_read.h"

class Event;
class fiber;

using namespace std;

class janus {

public:
	janus(histo_sort*, double);
	~janus();
	histo_sort* Histo;

	bool unpack(ifstream*);
	void analyze();
	
	Event* SIPMevent;
	fiber* Fiber;

	// Vector of most recent 20 events

  //Make a vector of events
  vector<Event> janusevts;

	vector<Event*> redbuffevents;
	vector<Event*> bluebuffevents;

	int Nsingles = 0;
	int Nmatched = 0;
	int Nskipped = 0;

	long nevts;

private:
	// Scale values for gain matching
	double bluegains[64];
	double redgains[64];

	float distance;

	void ReadGains(string, double*);
};

#endif
