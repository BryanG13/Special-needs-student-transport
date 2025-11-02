#ifndef OTHERS_H
#define OTHERS_H

#include "parameters/hardcoded.h"
#include <vector>

extern int estbus;         // Estimated passengers per bus
extern int Np;             // Cumulative sum counter (used to build the rp array)
extern int rp[nBuses + 1]; // Rank/Priority boundaries array (used for bus selection)

extern std::vector<double> BusPa; // arrival times of passengers in a bus

// location coordinates
extern double passengers[C][2];         // passenger locations
extern double mandatory[N][2];          // mandatory stop locations
extern double optional[(N - 1) * M][2]; // optional stop locations

// Arrival times of the passengers
extern double arrivals[C];

// Travel times
extern double traveltimep[C][Stations];        // travel times of people between passangers and stations
extern double traveltimes[Stations][Stations]; // travel times of buses between stations

// define sets----------------------------------------------------------------------------------------------
extern int I[nBuses];      // set for buses
extern int J[Stations];    // set for stops
extern int F[N];           // set for mandatory stops
extern int O[(N - 1) * M]; // set for optional stops
extern int P[C];           // set for passenger requests

// Cloesest stops
extern int closestS[Stations][Stations]; // closest stops to on another
extern int closestPS[C][Stations];       // closes stops to clients

// temp
extern int indexpt[C];         // passenger indices sorted by arrival time
extern double temparrivals[C]; // temporary copy of arrival times for sorting

// 2 opt
extern std::vector<int> route;                                 // route sequence for 2-opt optimization
extern int index1, index2, start, end_i, nextindex, mid, temp; // indices for 2-opt operations

// processing of the parameters: reading in data, calculating some useful metrics as input for optimization
void processingParameters();

#endif