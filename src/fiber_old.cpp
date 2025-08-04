// Created by Nicholas Dronchi
// Modified by Henry Webb
//                 _.-.
//          _.----'    `.
//         /             `.
//        /                `-._
// VK   .'     __.----.        `-._
//   .-'      /        `-._       =
//  //       /             `-.    =
//   //    .'                 `--.=
//    // .'
//     //

#include "fiber.h"

#include "CAENd5202.h"
#include "constants.h"

#include <algorithm>
#include <iostream>
#include <cmath>

fiber::fiber() {}

void fiber::clear() {
	posmaxhorz = -1;
	posmaxvert = -1;
	sumhorz = 0;
	sumvert = 0;
	multTrimmedX = 0;
	multTrimmedY = 0;
	tstampdiff = 0;
	tdiffx.clear();
	tdiffy.clear();
	badtx = false;
	badty = false;
	xindices.clear();
	yindices.clear();

  thetadeg = -999.;
  phideg = -999.;
  theta = -999.;
  phi = -999.;
}

/**
 * This function takes in a pair of matched fiber events and performs a center of gravity
 * calculation (first statistical moment), along with calculating other values for the
 * given event pair. The time gate values applied on a hit-wise basis were determined by
 * looking at histograms for the time difference of all the hits in both layers with the
 * "center hit" in the front layer (i.e. the front layer hit with the largest ToT). Note
 * also that the fiber layers are as such:
 * 
 * Vertical (blue) fiber gives x position
 * Horizontal (red) fiber gives y position
 */
bool fiber::make_2d(Event* horz, Event* vert, double distance,  s800_results S800_results) {
	// Set default values
	clear();

	// Apply front layer ToA gate, find center hit
	eventTiming ev;
	int mult = horz->GetNHits();
	double maxToT = 0;
	double maxToTalt1 = 0;
	int posmaxhorzalt = 0;
  //cout << "mult " << mult << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  //cout << "here1" << endl;
	for (int i = 0; i < mult; i++) {

		// Find alternate max ToT hit index
		ev = horz->GetTimingEvent(i);
		if (ev.ToTmatched > maxToTalt1) {
			maxToTalt1 = ev.ToTmatched;
			posmaxhorzalt = i;
		}
		//cout << "here" << endl;
    //TODO no time gates for now
		//if ((ev.ToA < 650) || (ev.ToA > 740) || (ev.ToTmatched < maxToT) continue;
		if (ev.ToTmatched < maxToT) continue;
		//cout << "here2" << endl;
		maxToT = ev.ToTmatched;
		posmaxhorz = i;
	}

	if (posmaxhorz == -1) {
		posmaxhorz = posmaxhorzalt;
		badtx = true;
	}
  //cout << "here2" << endl;
	// Advance declaration of variables
  double momhorz = 0;
  double momvert = 0;
	double threshhorz = 0;
	double threshvert = 0;
	double sumfibhorz = 0;
	double sumfibvert = 0;
	double PH;
	double tdiff;
	int tempi;

  if (posmaxhorz > 63 || posmaxhorz < -63)
  {
    cout << "Bad horizontal max position, ignore " << endl;
		badtx = true;
    return false;
  }

	// Horizontal fibers (front layer)
	eventTiming maxHitBlue = horz->GetTimingEvent(posmaxhorz);
	ix = maxHitBlue.pos;
  for (int i = 0; i < mult; i++) {

		// Apply hit-wise time gate
		ev = horz->GetTimingEvent(i);
		tdiff = (double)(ev.ToA - maxHitBlue.ToA) * 0.5; //ns
		tdiffx.push_back(tdiff);  
    //TODO no hit wise-time gate
		//if (tdiff < 0. || tdiff > 15.) continue;

		// Calculate values
		PH = max(ev.ToTmatched - threshhorz, 0.);
		sumhorz += PH;
    momhorz += PH * ev.pos;
		sumfibhorz += ev.pos;
		multTrimmedX++;
		tempi = i;
		xindices.push_back(tempi);
  }

	// Vertical fibers (back layer)
	tstampdiff = (vert->GetTimeStamp() - horz->GetTimeStamp()) * 1000.; //ns
	mult = vert->GetNHits();
	maxToT = 0;
	double maxToTalt2 = 0;
	int posmaxvertalt;
  for (int i = 0; i < mult; i++) {

		// Find alternate max ToT hit index
		ev = vert->GetTimingEvent(i);
		if (ev.ToTmatched > maxToTalt2) {
			maxToTalt2 = ev.ToTmatched;
			posmaxvertalt = i;
		}

		// Apply hit-wise time gate
		tdiff = (double)(ev.ToA - maxHitBlue.ToA) * 0.5 + tstampdiff; //ns
		tdiffy.push_back(tdiff);
    //TODO no hit wise-time gate
		//if (tdiff < 230. || tdiff > 270.) continue;

		// Get index of max ToT hit
		if (ev.ToTmatched > maxToT) {
			maxToT = ev.ToTmatched;
			posmaxvert = i;
		}

		// Calculate values
		PH = max(ev.ToTmatched - threshvert, 0.);
		sumvert += PH;
    momvert += PH * ev.pos;
		sumfibvert += ev.pos;
		multTrimmedY++;
		tempi = i;
		yindices.push_back(tempi);
  }

	// For testing, to analyze events that don't pass red fiber hit time gate
	if (posmaxvert == -1) {
    //cout << "maxtotalt " << maxToTalt2 << " alt " << posmaxvertalt << endl;
		posmaxvert = posmaxvertalt;
		badty = true;
	}

  if (posmaxvert > 63 || posmaxvert < -63)
  {
    cout << "Bad vertical max position, ignore " << endl;
		badty = true;
    return false;
  }
	eventTiming maxHitRed = vert->GetTimingEvent(posmaxvert);
	iy = maxHitRed.pos;

	// CHANNEL TO FIBER
	// fib = ((iCh - (iCh % 2)) / 2) + ((iCh % 2) * 32)

  // Calculate center of gravity
  double fibhorz = momhorz / sumhorz;
  double fibvert = momvert / sumvert;
	avgfibx = sumfibhorz / multTrimmedX;
	avgfiby = sumfibvert / multTrimmedY;
  //cout << "distance " << distance << endl;
  // Calculate physical parameters of matched event
  x = 1*(fibhorz+0.5)*0.5 - 16;              //mm
  y = 1*(fibvert+0.5)*0.5 - 16;              //mm
  //x = 1*(posmaxhorz+0.5)*0.5 - 16;              //mm //no center of grav calc
 // y = 1*(posmaxvert+0.5)*0.5 - 16;              //mm //no center of grav calc
  //x = -1*(fibhorz-0.5)*0.5 - 16;              //mm
  //y = -1*(fibvert-0.5)*0.5 - 16;              //mm
  //cout << x << " " << y << endl;
  //x = -1*(fibhorz-0.5)*0.5 + 16;              //mm
  //y = -1*(fibvert-0.5)*0.5 + 16;              //mm
  //r = sqrt((x * x) + (y * y)) * 0.1;          //cm
  //theta = atan(r / distance); //radians
  r = sqrt((x * x) + (y * y) + (distance*10*distance*10)) * 0.1;          //cm
  theta = acos(distance / r); //radians
  phi = atan2(y , x); //radians
  thetadeg = theta * rad_to_deg; //deg
  phideg = phi * rad_to_deg;         //deg
  //cout << " distance " << distance << " x " << x << " y " << y << " thetadeg " << thetadeg << " thetarad " << theta << endl;
	// Calculate additional values
	xdev = x - (-1*(ix-0.5)*0.5 + 16); //mm
	ydev = y - (-1*(iy-0.5)*0.5 + 16); //mm

  return true;
}
