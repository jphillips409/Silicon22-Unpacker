#ifndef janus_
#define janus_

#include <fstream>
#include <string>
#include <vector>
#include "histo_sort.h"
#include "histo_read.h"
#include "s800_results.h"

class Event;
class fiber;

using namespace std;

class janus {

public:
	janus(histo_sort*, double);
	~janus();
	histo_sort* Histo;

	bool unpack(ifstream*);
  void SetTarget(double Targetdist);
	void analyze(const s800_results& S800_results);
	
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

  //For event efficiency
  int evt_BRMatch = 0;
  int evt_BNOTR = 0;
  int evt_RNOTB = 0;

  //Flag for matched
  bool fibmatch = false;

	long nevts;

  float dist_pub;

private:
	// Scale values for gain matching
	double bluegains[64];
	double redgains[64];

	float distance;

	void ReadGains(string, double*);
};

#endif
