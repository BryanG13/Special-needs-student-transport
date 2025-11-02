#ifndef HARDCODED_H
#define HARDCODED_H

// Instance/run configuration ----------------------------------------------------------------
constexpr int instance = 14; // Instance number for indexing purposes (compile-time constant)

constexpr int nRUN = 10; // Number of runs (compile-time constant)

// WEIGHT FACTORS (compile-time constants)--------------------------------------------------
constexpr float c1 = 0.25f;              // (For travel-time of the buses)
constexpr float c2 = 0.35f;              // (For walking time of the passengers)
constexpr float c3 = 1.0f - c1 - c2;     // (For the absolute difference in desired arrival time and actual arrival time)

// Define parameters-----------------------------------------------------------------------
constexpr int nBuses = 5;                 // number of buses available
constexpr int N = 10;                     // number of mandatory stations
constexpr int M = 3;                      // number of stations in cluster
constexpr int Stations = (N - 1) * M + N; // amount of Stations
constexpr int C = 40;                     // number of clients in opt horizon

constexpr int bCapacity = 15;          // Bus capacity
constexpr float pspeed = 1.0f;         // passengers speed in meter per second
constexpr float bspeed = 40.0f / 3.6f; // bus speed in m/s
constexpr int delta = 30;              // acceleration and deceleration time  in seconds
constexpr int tau = 5;                 // dwell time coefficient in seconds
constexpr int d = 20 * 60;             // threshold of individual walking time in seconds
constexpr int d_time1 = 15 * 60;       // Maximum amount of time a passenger can arrive too early
constexpr int d_time2 = 5 * 60;        // Maximum amount of time a passenger can arrive too late

// LNS and destroy/repair parameters (compile-time constants)
constexpr int f_IT = 5000;        // max iterations for the LNS cycles

// `ndestroy` may be changed at runtime by the heuristic (randomized destroy size).
// Use an inline variable in the header so it has external linkage but remains
// defined once across translation units while still being mutable.
inline int ndestroy = 3;       // default number of passengers to remove

constexpr int secondS = 100 - 25; // Second station selection threshold (percentage)
constexpr int thirdS = 5;         // Third station selection threshold (percentage)

#endif // HARDCODED_H
