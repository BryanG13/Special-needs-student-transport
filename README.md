# Special Needs Student Transport Optimization

An optimization heuristic for routing vehicles to transport students with special needs to their respective schools.

## Overview

This project implements a metaheuristic to solve the vehicle routing problem for transporting students with special educational needs. The algorithm considers various constraints including vehicle capacity, wheelchair accommodations, time windows, and transfer points.

## Features

- **Vehicle Capacity Management**: Handles different capacity requirements for wheelchair and non-wheelchair students
- **Time Window Constraints**: Ensures students arrive at schools within specified time windows
- **Transfer Support**: Allows students to transfer between vehicles at designated transfer points
- **Multiple Schools**: Routes students to their assigned schools efficiently
- **Travel Time Limits**: Enforces maximum travel time constraints for student comfort

## Building the Project

### Prerequisites
- CMake 3.22 or higher
- C++17 compatible compiler
- Threading support

### Build Instructions

```bash
# Create and navigate to build directory
mkdir -p build && cd build

# Configure with CMake
cmake ..

# Build the project
cmake --build . --config Release
```

Or use the provided VS Code task:
```bash
# In VS Code, run the "build" task
```

## Project Structure

```
├── include/
│   ├── optimization/      # LNS and initial solution algorithms
│   ├── parameters/        # Problem parameters and configuration
│   └── sideFunctions/     # Helper functions for reporting and sorting
├── src/
│   └── main.cpp          # Main implementation
├── Data/
│   ├── input/            # Input data files
│   └── results/          # Output results
└── CMakeLists.txt        # Build configuration
```

## Input Data

The program expects:
- Travel time matrices between locations
- Student information (wheelchair requirements, assigned schools)
- School time windows (earliest/latest arrival times)
- Vehicle capacity and parameters

## Algorithm

The solution approach uses:
1. **Initial Solution**: Nearest neighbor heuristic with simulated annealing
2. **LNS Operators**: Destroy and repair methods to improve solutions
3. **Feasibility Checks**: Ensures all constraints are satisfied

## License

See [LICENSE](LICENSE) file for details.

## Author

Bryan Galarza
