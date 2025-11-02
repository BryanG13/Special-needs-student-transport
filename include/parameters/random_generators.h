#ifndef RANDOM_GENERATORS_H
#define RANDOM_GENERATORS_H

#include <random>

// Random number generator
extern std::default_random_engine rng_engine;

// Random distributions (initialized after parameters are loaded)
extern std::uniform_int_distribution<int> rand_passenger;     // random passenger selection (0 to C-1)
extern std::uniform_int_distribution<int> rand_station;       // random station selection (0 to Stations-1)
extern std::uniform_int_distribution<int> rand_percentage;    // random percentage (1 to 100)
extern std::uniform_int_distribution<int> rand_bus_weighted;  // weighted bus selection (1 to Np)
extern std::uniform_int_distribution<int> rand_destroy_count; // random number of passengers to destroy (3 to 5)

// Function to initialize distributions after parameters are loaded
void initializeRandomDistributions();

#endif // RANDOM_GENERATORS_H
