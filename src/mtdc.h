#ifndef MTDC
#define MTDC
#include <cstdlib>
/**
 * !\brief handles the readout from a MesytecTDC
 *
 * This class deals with the read out from a number of
 * MTDC VME modules
 */

class mtdc {
 public:
  unsigned short moduleid;
  unsigned short tdc_res;
  float		 res;
  unsigned short nwords;
  unsigned short channel[128];
  unsigned short data[128];
  unsigned short tstamplow;
  unsigned short tstamphigh;
  unsigned short tstampxtend;
  bool		 trigger;

  unsigned short number; // number of converted channels
  unsigned short crate; // crate number
  unsigned short geographic; // geagraphic address of module
  unsigned short underflow[128];
  unsigned short overflow[128];
  unsigned short * read(unsigned short *);
};
#endif
