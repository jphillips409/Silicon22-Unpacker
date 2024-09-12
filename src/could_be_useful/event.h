#ifndef event_
#define event_

const int Nstrips = 14*32;
const int Ncsi = 14*4;

struct Event
{
  int nfront;
  float frontE[Nstrips];
  int frontT[Nstrips];
  int frontID[Nstrips];
  int nback;
  float backE[Nstrips];
  int backT[Nstrips];
  int backID[Nstrips];
  int ncsi;
  float csiE[Ncsi];
  int csiER[Ncsi];
  int csiID[Ncsi];
  int csiT[Ncsi];
};
#endif
