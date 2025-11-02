

// Initialize the LNS algorithm
void initializeLNS();

//  Destroy part of the solution
void destroyOperator();

// Rebuild part of the destroyed solution
void rebuildOperator();

// Acceptance criteria
void acceptanceCriteria(std::clock_t start_time);

// 2opt improvement after a destroy-repair cycle
void improvementLNS();