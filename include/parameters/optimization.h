#ifndef OPT_H
#define OPT_H

#include "parameters/hardcoded.h"
#include <vector>

// Decision variables: passenger-to-station assignment (x)
extern bool ***xsol;  // passenger-to-station assignment: solution being evaluated
extern bool ***xk;    // passenger-to-station assignment: current best solution
extern bool ***xbest; // passenger-to-station assignment: overall best solution

// Decision variables: station-to-station routing (y)
extern bool yk[nBuses][Stations][Stations];    // station-to-station routing: current best routing
extern bool ysol[nBuses][Stations][Stations];  // station-to-station routing: solution being evaluated routing
extern bool ybest[nBuses][Stations][Stations]; // station-to-station routing: overall best routing

// Arrival times
extern double Ak[nBuses];    // current best arrival time for each bus
extern double Dk[nBuses];    // current best idle/delay time before departure for each bus
extern double Dsol[nBuses];  // solution idle/delay time
extern double Asol[nBuses];  // solution arrival time
extern double Abest[nBuses]; // overall best arrival time
extern double Dbest[nBuses]; // overall best idle/delay time

extern int bCapi[nBuses]; // used capacity of each bus in current solution

extern std::vector<double> rt;           // runtime at each improvement
extern std::vector<double> best_srt;     // best solution runtimes
extern std::vector<int> prog, best_prog; // cost progression over iterations

// Algorithm variables
extern double cost;                   // total objective function value
extern double currentcost, ncost, dE; // current cost, new cost, and cost difference for SA
extern double elapsed_time;           // total elapsed time for the run
extern int it;                        // iteration counter for LNS
extern double nextcost;               // cost of proposed solution

// Initialization phase variables
extern int pa;              // number of passengers assigned so far
extern int cap, countcap;   // current bus capacity usage and counter
extern int minDist;         // minimum distance value
extern int mini;            // index of minimum
extern int sumx;            // sum of x values (for checking if passenger is assigned)
extern double sum, sum_min; // utility sums

extern int pindex[C];                      // station assignment for each passenger during initialization
extern double earliestArr, latestArr, Arr; // earliest, latest, and chosen arrival times
extern int minps;                          // minimum passenger-to-station distance
extern int minpsi;                         // index of closest station for passenger

extern double timewindow; // time window between passengers
extern double firsttw;    // arrival time of first passenger on bus
extern int p_0, p_new;    // previous and new passenger indices for median calculation
extern int bestit;        // iteration number of best solution found

extern int destroyed[C]; // tracks which bus each destroyed passenger is reassigned to (-1 if not destroyed)
extern int di;           // passenger index to destroy/remove
extern int bi;           // new bus for the removed passenger
extern int pbi;          // previous bus the removed passenger was on
extern int pst;          // previous station the removed passenger went to

extern double firstmin, secmin, thirdmin; // first, second, third minimum distance values
extern int firstminj, secminj, thirdminj; // indices of first, second, third minimum
extern int newj;                          // new station for reassigned passenger
extern bool infeas;                       // flag for infeasible solution
extern bool inGraph, inSuccessor;         // flags for checking if station is in route

// Bus selection variables
extern double value_buses[nBuses];  // ranking values for buses (based on avg arrival time difference)
extern int rank_buses[nBuses];      // sorted bus indices by ranking value
extern int it_bi;                   // selected bus index from weighted distribution
extern double earlyTW[nBuses];      // earliest passenger arrival time on each bus
extern double lateTW[nBuses];       // latest passenger arrival time on each bus
extern int trybuses, decide, which; // counters and decisions for bus/station selection

// Simulated Annealing parameters and 2-opt
extern int maxdist;  // maximum distance between nearest neighbors
extern double T;     // initial temperature
extern double alpha; // cooling rate
extern double T_end; // end temperature
extern int L;        // iterations per temperature
extern int m;        // number of nearest neighbors to consider

// Best solution tracking
extern double best_cost;      // best cost found across all runs
extern double best_rt;        // runtime of best solution
extern double last, bestlast; // time of last improvement

// delete raw pointers
void deletePointers();

// initialize pointers
void initilalizeX();

#endif