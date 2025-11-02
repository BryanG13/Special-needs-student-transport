#include <iostream>
//#include <istream>
#include <fstream>
#include <sstream>
#include <stdlib.h> 
#include <random>
//#include <time.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

/*
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC_NEW
#include <crtdbg.h>
#include <assert.h>
#endif
*/

#include <array>

//--------------------------------------------------- PARAMETERS -----------------------------------------------------------------------
const int tau_1 = 30; //dwell time at boarding for non-wheelchair students, in seconds
const int tau_2 = 120; //dwell time at boarding for wheelchair students, in seconds
const int tau_3 = 5; //dwell time at aligting and transfer for non-wheelchair students, in seconds
const int tau_4 = 70; //dwell time at aligting and transfer for wheelchair students, in seconds
const int D_t = 90 * 60; //Max time students can be travelling
const int C_b = 36; // Bus capcity
bool* W_ch; // Indicator to see if a student p in P needs a wheel chair, needs to be read
int* St_Sc; // Indicator to see to which school a student p in P is assigned to, needs to be read
int* t_ea; // Earliest arrival times for student p in P
int* t_la; // Latest arrival times for student p in P

int P = 0; // Number of students
int S = 0; // Number of schools
int L = 0; // Total number of locations
int T = -2; // number ID of the location for a tra nsfer

int** TT;// travel times between locations
int** closestS;

int** yk;
int** ysol;

int* best_route;


std::default_random_engine generator;
std::uniform_real_distribution<double> one_to_zero(0, 1);

using namespace std;
std::piecewise_linear_distribution<double> triangular_distribution(double min, double peak, double max){
	std::array<double, 3> i{ min, peak, max };
	std::array<double, 3> w{ 1, 0.65, 0.25 };
	return std::piecewise_linear_distribution<double>{i.begin(), i.end(), w.begin()};
}

inline void read_travel_matrix(string path, int **& TT){
	ifstream  data(path);
	string line;
	vector<vector<string>> parsedCsv;
	while (getline(data, line)){
		stringstream lineStream(line);
		string cell;
		vector<string> parsedRow;
		while (getline(lineStream, cell, ',')){
			parsedRow.push_back(cell);
		}

		parsedCsv.push_back(parsedRow);
	}

	//cout << " DATA\n";
	for (int i = 1; i <= L; i++) {
		for (int j = 1; j <= L; j++) {
			TT[i - 1][j - 1] = stoi(parsedCsv[i][j]);
			//cout << TT[i - 1][j - 1] << "\t";
		}
		//cout << endl;
	}
	data.close();
	//exit(0);
}

inline int giveIndex(vector<int>& cities, const int edge, const int n) {
	int i;
	for (i = 0; i < n; i++) {
		if (cities[i] == edge) break;
	}
	return i;
}

inline int partition(int*& index, int*& dist,  const int low, const int high)
{
	double pivot = dist[high]; // pivot  
	int i = (low - 1); // Index of smaller element  
	int t, t0;
	for (int j = low; j <= high - 1; j++) {
		// If current element is smaller than the pivot  
		if (dist[j] < pivot) {
			i++; // increment index of smaller element  
			t = index[j];
			index[j] = index[i];
			index[i] = t;

			t0 = dist[j];
			dist[j] = dist[i];
			dist[i] = t0;

		}
	}
	t = index[high];
	index[high] = index[i + 1];
	index[i + 1] = t;

	t0 = dist[high];
	dist[high] = dist[i + 1];
	dist[i + 1] = t0;

	return (i + 1);
}

inline void quickSort(int*& index, int*& dist, const int low, const int high) {
	if (low < high) {
		/* pi is partitioning index, arr[p] is now
		at right place */
		int pi = partition(index, dist, low, high);

		// Separately sort elements before  
		// partition and after partition  
		quickSort(index, dist, low, pi - 1);
		quickSort(index, dist, pi + 1, high);
	}
}

inline void make_feasible2(vector<int>& cities, const int b, const  int n, vector<int>& schools, const vector<int>& l_stu) {
	int N_s = n, N_c = 0;
	//vector<int> schools;
	int i, j;

	for (i = 0; i < N_s; i++) { // delete schools in te route first;
		if (cities[i] >= P) {
			cities.erase(cities.begin() + i);
			N_s--;
			i--;
		}
	}

	N_c = schools.size();
	int** schools_last = new int* [N_c];
	for (j = 0; j < N_c; j++) {
		schools_last[j] = new int[2];
		schools_last[j][0] = schools[j];
		schools_last[j][1] = cities.back();
		for (i = N_s - 1; i >=0 ; i--) {
			if (schools_last[j][0] == St_Sc[cities[i]]|| cities[i] == l_stu[j]) {
				schools_last[j][1] = cities[i];
				break;
			}
		}
	}
	
	//check if there are any transfers on this bus b
	bool trans_b = false;
	vector<int> tf_schools; // list of potential schools that need to be visited after T
	for (i = 0; i < P; i++) {
		if (ysol[i][1] == b) { //if this bus is a second bus in a transfer
			//make sure that T is visited before their assigned school
			trans_b = true;
			if (!std::count(tf_schools.begin(), tf_schools.end(), St_Sc[i])) {
				tf_schools.push_back(St_Sc[i]);
			}
			//cout << " /////  second bus in transfer detected for student " << i << " with assigned school S_" << St_Sc[i] - P << ", first bus is " << b_ysol[i][0] << endl;
		}
	}
	int i2;
	for (i = 0; i < N_s; i++) {
		if (ysol[cities[i]][1] != -1) { //if this is the first bus in the transfer
			// make sure they are visisdted before T
			//cout << " /////  first bus in transfer detected for student " << cities[i] << " second bus is " << b_ysol[cities[i]][1] << endl;
			for (j = 0; j < N_c; j++) {
				if (schools_last[j][0] == T) {
					auto it = find(cities.begin(), cities.end(), schools_last[j][1]);
					i2 = it - cities.begin();
					if (i > i2) {
						schools_last[j][1] = cities[i];
					}
					break;
				}
			}
		}
	}

	//insert the schools in the right position now
	for (i = 0; i < N_c; i++) {
		// cout << " School " << schools_last[i][0] << " last student " << schools_last[i][1] << endl;
		auto it = find(cities.begin(), cities.end(), schools_last[i][1]);
		if (it != cities.end()) {
			cities.insert(it + 1, schools_last[i][0]);
			N_s++;
		}
		else {
			cout << "SOMETHING WENT WRONG, CANNOT FIND STUDENT in route anymore" << endl;
		}
	}
	int i_s, i_t, N_tf;
	if (trans_b) { // if b was the second bus in the transfer
		// cout << " /////  second bus in transfer detected!!\n";
		auto it = find(cities.begin(), cities.end(), T);
		i_t = it - cities.begin();
		i_s;
		N_tf = tf_schools.size();
		for (i = 0; i < N_tf; i++) {
			//cout << " check \n";
			//cout << " School S_" << tf_schools[i] -P << endl;
			auto it2 = find(cities.begin(), cities.end(), tf_schools[i]);
			i_s = it2 - cities.begin();
			// cout << " route size: " << cities.size() << " i_s " << i_s << " i_t " << i_t << " N_c: " << N_c << endl;
			cout << endl;
			if (i_t > i_s) { //if transfer visisted after school swap schools
				cities.insert(it2 + 1, T); //insert infront of school
				cities.erase(i_s + cities.begin()); // erase school
				cities.insert(i_t + cities.begin() + 1, tf_schools[i]);
				cities.erase(i_t + cities.begin());
				i_t = i_s; //update index of transfer
			}
		}
	}

	for (i = 0; i < N_c; i++) {
		delete schools_last[i];
	}
	delete schools_last;
}

inline int make_feasible(vector<int>& cities, const int b,  const int n, const bool initial, vector<int>& schools, const int N_c, vector<int>& tf_schools, const int N_tf, const bool trans_b) {
	int N_s = n;
	int i, j, p, cost;
	if (initial) {
		if(cities.back() < P) cities.push_back(St_Sc[cities[0]]);
		for (i = 0; i < N_s - 1; i++) {
			if (cities[i] >= P) {
				cities.erase(i + cities.begin());
				break;
			}
		}
		
		cost = 0;
		N_s = cities.size();
		for (i = 1; i < N_s; i++) {
			cost += TT[cities[i - 1]][cities[i]];
			if (cities[i] < P)cost += tau_1 * (1 - W_ch[cities[i]]) + tau_2 * W_ch[cities[i]];
			else {
				for (p = 0; p < P; p++) {
					if ((ysol[p][0] == b && St_Sc[p] == cities[i]) || (ysol[p][1] == b && cities[i] == T)) cost += tau_3 * (1 - W_ch[p]) + tau_4 * W_ch[p];
				}
			}
		}
		//if (N_s != n) cout << " whattt\n";
		return cost;
		
	}
	
	for (i = 0; i < N_s; i++) { // delete schools in te route first;
		if (cities[i] >= P) {
			cities.erase(cities.begin() + i);
			N_s--;
			i--;
		}
	}

	//vector<vector<int>> schools_last;
	int** schools_last = new int* [N_c];
	for (j = 0; j < N_c; j++) {
		schools_last[j] = new int[2];
		schools_last[j][0] = schools[j];
		schools_last[j][1] = cities.back();
		for (i = N_s - 1; i >= 0; i--) {
			if (schools_last[j][0] == St_Sc[cities[i]]) {
				schools_last[j][1] = cities[i];
				break;
			}
		}
		//cout << "school: " << schools_last[j][0] << " last student " << schools_last[j][1] << endl;
	}

	//check if there are any transfers on this bus b
	
	int i2 = 0;
	for (i = 0; i < N_s; i++) {
		if (ysol[cities[i]][1] != -1) { //if this is the first bus in the transfer
			// make sure they are visisdted before T
			//cout << " /////  first bus in transfer detected for student " << cities[i] << " second bus is " << b_ysol[cities[i]][1] << endl;
			for (j = 0; j < N_c; j++) {
				if (schools_last[j][0] == T) {
					auto it = find(cities.begin(), cities.end(), schools_last[j][1]);
					i2 = it - cities.begin();
					if (i > i2) {
						schools_last[j][1] = cities[i];
					}
					break;
				}
			}
		}
	}
	
	//insert the schools in the right position now
	for (i = 0; i < N_c; i++) {
		//cout << " School " << schools_last[i][0] << " last student " << schools_last[i][1] << endl;
		auto it = find(cities.begin(), cities.end(), schools_last[i][1]);
		if (it != cities.end()){
			cities.insert(it + 1, schools_last[i][0]);
			N_s++;
		}
		else {
			cout << "SOMETHING WENT WRONG, CANNOT FIND STUDENT in route anymore"  << endl;
		}
	}
	int i_t, i_s;
	if (trans_b) { // if b was the second bus in the transfer
		//cout << " /////  second bus in transfer detected!!\n";
		auto it = find(cities.begin(), cities.end(), T);
		i_t = it - cities.begin();
		i_s;
		for (i = 0; i < N_tf; i++) {	
			//cout << " check \n";
			//cout << " School S_" << tf_schools[i] -P << endl;
			auto it2 = find(cities.begin(), cities.end(), tf_schools[i]);
			i_s = it2 - cities.begin();
			//cout << " route size: " << cities.size() << " i_s " << i_s << " i_t " << i_t << " N_c: " << N_c << endl;
			if (i_t > i_s) { //if transfer visisted after school swap schools
				cities.insert(it2 + 1, T); //insert infront of school
				cities.erase(i_s + cities.begin()); // erase school
				cities.insert(i_t + cities.begin() + 1, tf_schools[i]);
				cities.erase(i_t + cities.begin());

				i_t = i_s; //update index of transfer
			}
		}
	}
	//cout << " calculate cost now \n";
	cost = 0;
	N_s = n;
	for (i = 1; i < N_s; i++) {
		cost += TT[cities[i - 1]][cities[i]];
		if (cities[i] < P)cost += tau_1 * (1 - W_ch[cities[i]]) + tau_2 * W_ch[cities[i]];
		else {
			for (p = 0; p < P; p++) {
				if((ysol[p][0] == b && St_Sc[p] == cities[i]) || (ysol[p][1] == b && cities[i] == T)) cost += tau_3 * (1 - W_ch[p]) + tau_4 * W_ch[p];
			}
		}
	}
	//if (N_s != n) cout << " whattt\n";
	//delete mem
	for (i = 0; i < N_c; i++) {
		delete schools_last[i];
	}
	delete schools_last;
	return cost;
}

inline bool give_max_shifts(const vector<int>& route, const vector<int>& times, const int N_r,  int& past, int& future, int& mpast, int& mfuture) {
	past = 100000;
	future = 100000;
	mpast = 0;
	mfuture = 0;
	int temp, temp2;
	bool FEAS = true;
	for (int i = 0; i < N_r; i++) {
		if (route[i] >= P) {  //if school
			temp = times[i] - t_ea[route[i] - P]; // diff with LB
			temp2 = t_la[route[i] - P] - times[i]; // diff with UB (all assuming its a feasible Dsol !!!!! --> not anymore?)
			if (temp < 0) {
				FEAS = false;
				if (-temp > mfuture) mfuture = -temp;
			}
			else if (temp < past) past = temp;
			if (temp2 < 0) {
				FEAS = false;
				if (-temp2 > mpast) mpast = -temp2;
			}
			if (temp2 < future) future = temp2;
		}
	}
	return FEAS;
}

inline bool time_table2(const vector<int>& xsol, vector<int>& b_Dsol, const int n_v, int* alights, const int N_sc) {

	int N_r = xsol.size();
	int fsch, l, j;
	//cout << "N_r " << N_r << endl;

	for (l = 0; l < N_r; l++) {
		if (xsol[l] >= P) { // find first school visited
			fsch = l;
			break;
		}
	}

	if (N_sc == 1) {
		b_Dsol.clear();
		//cout << "first school " << xsol[fsch] - P << endl;
		b_Dsol.insert(b_Dsol.begin(), t_ea[xsol[fsch] - P]); // schedule first school to arrive asap 
		b_Dsol.insert(b_Dsol.begin(), b_Dsol[0] - TT[xsol[fsch]][xsol[fsch - 1]] - alights[0]);
		//cout << "check\n";
		for (l = fsch - 2; l >= 0; l--) { // determine timetable of previous stops
			b_Dsol.insert(b_Dsol.begin(), b_Dsol[0] - TT[xsol[l]][xsol[l + 1]] - tau_1 * (1 - W_ch[xsol[l + 1]]) - tau_2 * W_ch[xsol[l + 1]]);
		}
		//cout << "check\n";
		return true;
	}
	else {
		//cout << "first school " << xsol[fsch] - P << " at position " << fsch << endl;
		vector<int> Dsol;

		Dsol.insert(Dsol.begin(), t_ea[xsol[fsch] - P]); // schedule first school to arrive asap 
		Dsol.insert(Dsol.begin(), Dsol[0] - TT[xsol[fsch]][xsol[fsch - 1]] - alights[0]);
		for (l = fsch - 2; l >= 0; l--) { // determine timetable of previous stops
			Dsol.insert(Dsol.begin(), Dsol[0] - TT[xsol[l]][xsol[l + 1]] - tau_1 * (1 - W_ch[xsol[l + 1]]) - tau_2 * W_ch[xsol[l + 1]]);
		}

		int s = 1;
		int travel_time = 0;
	

		for (l = fsch + 1; l < N_r; l++) { // determine the next stops if fsch is not the last stop
			travel_time = Dsol[l - 1] + TT[xsol[l - 1]][xsol[l]];
			if (xsol[l] < P) {
				travel_time += tau_1 * (1 - W_ch[xsol[l]]) + tau_2 * W_ch[xsol[l]];
			}
			else {
				travel_time += alights[s];
				s++;
			}
			Dsol.push_back(travel_time);
		}

		int maxFS = 0, maxBS = 0, minBS = 0, minFS = 0, FW = 0, BW = 0;
		bool FS = false, FS2 = false;

		FS = give_max_shifts(xsol, Dsol, N_r, maxBS, maxFS, minBS, minFS);
		if (FS == true) {
			b_Dsol.clear();
			for (int i = 0; i < N_r; i++) {
				b_Dsol.push_back(Dsol[i]);
			}
			return true;
		}
		else {
			if (minBS == 0 && minFS != 0) { //if no need to go back but you need to go forwards
				if (minFS > maxFS) return false;
				else {
					for (int i = 0; i < N_r; i++) {
						Dsol[i] += minFS;
					}
					FS = true;
				}
			}
			else if (minBS != 0 && minFS == 0) { // if no need to go forwards but we need to go back
				if (minBS > maxBS) return false;
				else {
					for (int i = 0; i < N_r; i++) {
						Dsol[i] -= minBS;
					}
					FS = true;
				}
			}
			else { // both back and forwards needed ???
				if (minFS > maxFS || minBS > maxBS) return false;
				else {
					s = 1;
					for (int i = fsch + 1; i < N_r; i++) { // start after the first school
						if (xsol[i] >= P) {
							if (t_ea[xsol[i] - P] > Dsol[i]) { // if too early
								//move FW
								FW = t_ea[xsol[i] - P] - Dsol[i];
								if (FW <= maxFS) {
									for (int y = 0; y < N_r; y++) {
										Dsol[y] += FW;
									}
									FS2 = give_max_shifts(xsol, Dsol, N_r,  maxBS, maxFS, minBS, minFS);
									if (FS2) {
										FS = true;
										break;
									}
								}
								else {
									// try to wait
									int s2 = s;
									Dsol[i] = t_ea[xsol[i] - P];
									for (int y = i + 1; y < N_r; y++) {
										travel_time = Dsol[i - 1] + TT[xsol[i - 1]][xsol[i]];
										if (xsol[i] < P) {
											travel_time += tau_1 * (1 - W_ch[xsol[i]]) + tau_2 * W_ch[xsol[i]];
											//cout << xsol[i] << endl;
										}
										else {
											travel_time += alights[s2];
											s2++;
										}
									}
									FS2 = give_max_shifts(xsol, Dsol, N_r, maxBS, maxFS, minBS, minFS);
									if (FS2) {
										FS = true;
										break;
									}
								}
							}
							else if (t_la[xsol[i] - P] < Dsol[i]) { // if too late
								// move BW
								BW = Dsol[i] - t_la[xsol[i] - P];
								if (BW < maxBS) {
									for (int y = 0; y < N_r; y++) {
										Dsol[y] -= BW;
									}
									FS2 = give_max_shifts(xsol, Dsol, N_r, maxBS, maxFS, minBS, minFS);
									if (FS2) {
										FS = true;
										break;
									}
								}
								else {
									return false;
								}
							}
							s++;
						}
					}
				}

			}
		}

		if (FS) {
			b_Dsol.clear();
			for (int i = 0; i < N_r; i++) {
				b_Dsol.push_back(Dsol[i]);
			}
			return FS;
		}
		else return false;
	}
}

inline int max_TT(vector<vector<int>>& b_routes, const int location, int indexr, const bool school_in, const int alt_school, const int bus, vector<vector<int>>& b_times, const bool p_in) {
	//cout << " Student " << location << endl;
	if (!p_in) { // if there is no added student  --> in case we deal with the second bus in transfer
		int tt = 0, max_tt = 0, N_r = 0;
		for (int p = 0; p < P; p++) {
			tt = 0;
			//if bus is second bus in transfer
			int start = -1, end = -1, alight = 0;
			int bus1 = ysol[p][0];
			if (bus1 != -1 && ysol[p][1] == bus) { // if bus is the second bus in a transfer
				N_r = b_routes[bus1].size();
				for (int i = 0; i < N_r; i++) { // time on previous bus 
					if (b_routes[bus1][i] == p) start = b_times[bus1][i];
					if (b_routes[bus1][i] == T) {
						end = b_times[bus1][i];
						break;
					}
				}
				tt = (end - start) + tau_3 * (1 - W_ch[location]) + tau_4 * W_ch[location]; // plus extra alighting of new student 
				start = -1, end = -1, alight = 0;
				N_r = b_routes[bus].size();
				for (int i = 0; i < N_r; i++) { // time on second bus
					if (b_routes[bus][i] == T) start = b_times[bus][i];
					if (b_routes[bus][i] == St_Sc[p]) {
						end = b_times[bus][i];
						break;
					}
					if (b_routes[bus][i] < P)if (start != -1 && St_Sc[p] == St_Sc[b_routes[bus][i]])alight += tau_3 * (1 - W_ch[b_routes[bus][i]]) + tau_4 * W_ch[b_routes[bus][i]];
				}
				tt += (end - start - alight);
			}
			else if (bus1 == bus && ysol[p][1] != -1) { // if bus is the first  bus in transfer
				start = -1, end = -1, alight = 0;
				N_r = b_routes[bus].size();
				for (int i = 0; i < N_r; i++) {
					if (b_routes[bus][i] == p) start = b_times[bus][i];
					if (b_routes[bus][i] == T) {
						end = b_times[bus][i];
						break;
					}
				}
				tt = (end - start);

				int bus2 = ysol[p][1];
				start = -1, end = -1;
				N_r = b_routes[bus2].size();
				for (int i = 0; i < N_r; i++) {
					if (b_routes[bus2][i] == T) start = b_times[bus2][i];
					if (b_routes[bus2][i] == St_Sc[p]) {
						end = b_times[bus2][i];
						break;
					}
					if (b_routes[bus2][i] < P)if (start != -1 && St_Sc[p] == St_Sc[b_routes[bus2][i]])alight += tau_3 * (1 - W_ch[b_routes[bus2][i]]) + tau_4 * W_ch[b_routes[bus2][i]];
				}
				tt += (end - start - alight) + tau_3 * (1 - W_ch[location]) + tau_4 * W_ch[location];

			}
			else if (bus1 == bus && ysol[p][1] == -1) { // only one bus 
				start = -1, end = -1, alight = 0;
				N_r = b_routes[bus].size();
				for (int i = 0; i < N_r; i++) {
					if (b_routes[bus][i] == p) start = b_times[bus][i];
					if (b_routes[bus][i] == St_Sc[p]) {
						end = b_times[bus][i];
						break;
					}
					if (b_routes[bus][i] < P)if (start != -1 && St_Sc[p] == St_Sc[b_routes[bus][i]])alight += tau_3 * (1 - W_ch[b_routes[bus][i]]) + tau_4 * W_ch[b_routes[bus][i]];
				}
				tt = (end - start - alight) + tau_3 * (1 - W_ch[location]) + tau_4 * W_ch[location];
			}
			

			if (tt > max_tt) max_tt = tt;
		}
		if (max_tt <= D_t) {
			auto itr0 = find(b_routes[bus].begin(), b_routes[bus].end(), T);
			int j = itr0 - b_routes[bus].begin();
			b_times[bus][j] += tau_3 * (1 - W_ch[location]) + tau_4 * W_ch[location];

			return max_tt;
		}
		else return 100000;

	}
	int N_r = b_routes[bus].size();
	//cout << "N_r: " << N_r << endl << " route:\n";
	
	vector<int> route;
	vector<int> l_stu;
	vector<int> schools;
	for (int i = 0; i < N_r; i++) {
		if (b_routes[bus][i] < P) route.push_back(b_routes[bus][i]);
		else {
			schools.push_back(b_routes[bus][i]);
			if (i != 0) {
				l_stu.push_back(b_routes[bus][i - 1]);
			}
			else l_stu.push_back(-1);
		}
		//cout << b_route[i] << "\t";
	}
	if (school_in) {
		if (!std::count(schools.begin(), schools.end(), St_Sc[location])) {
			schools.push_back(St_Sc[location]);
			l_stu.push_back(-1);
		}
	}
	N_r = route.size();
	

	if (indexr >= N_r) indexr = N_r - 1;
	route.insert(route.begin() + indexr, location);
	N_r = route.size();
	
	//cout << "check "<<  endl;
	int N_sc = schools.size();
	for (int s = 0; s < N_sc; s++) {
		//insert_school(route, schools[s], P, St_Sc, b_ysol, T, bus, sec_tr);
		route.push_back(schools[s]);
		//cout << "S_" << schools[s] - P << endl;
	}
	N_r = route.size();
	int tt = 0;

	//cout << "Before make feasible\n";
	make_feasible2(route, bus, N_r, schools, l_stu);
	//cout << " after make feasible\n";
	schools.clear();
	for (int i = 0; i < N_r; i++) {
		if (route[i] >= P) {
			schools.push_back(route[i]);
		}
		//if (bus == 14) cout << route[i] << "\t";
	}
	//if (bus == 14) cout << endl;

	if (!school_in) {
		//cout << " check \n";
		int is = -1, ip = -1;
		for (int i = 0; i < N_r; i++) {
			if (route[i] == alt_school) is = i;
			else if (route[i] == location) ip = i;
		}
		//cout << " ip " << ip << " is " << is << endl;
		if (is < 0) {
			return 1000000;
		}
		else if (ip > is) {
			route.erase(route.begin() + ip);
			route.insert(route.begin() + is - 1, location);
		}
	}
	//cout << " check 2\n";

	int i = 0, p = 0, fs = 0, ls = 0;
	//cout << " Schools size  " << N_sc << endl;
	int ssmaxtime = 0, ssmintime = 0, sspottime = 0, ssrtime = 0, alight = 0, alight2 = 0, n_schools = 0;
	bool ss_feas = true;

	int* alights = new int[N_sc];
	for (int s = 0; s < N_sc; s++) {
		alight = 0;
		for (int i = 0; i < N_r; i++) {
			if (route[i] < P) {
				if (St_Sc[route[i]] == schools[s] || (schools[s] == T && ysol[route[i]][1] != -1) || schools[s] == alt_school) {
					alight += tau_3 * (1 - W_ch[route[i]]) + tau_4 * W_ch[route[i]];
				}
			}
		}
		if (schools[s] == T) { // boarding time
			for (int p = 0; p < P; p++) {
				if (ysol[p][1] == bus) {
					alight += tau_3 * (1 - W_ch[p]) + tau_4 * W_ch[p];
				}
			}
		}
		alights[s] = alight;
	}

	//calculate the timetable
	vector<int> time;
	bool tfeas = time_table2( route, time, bus, alights, N_sc);
	if (tfeas) {
		// calculate the max tt
		int max_tt = 0;
		N_r = route.size();
		for (p = 0; p < P; p++) {
			tt = 0;
			if ((ysol[p][0] == bus || p == location) && ysol[p][1] == -1) { // if bus is a vehicle and no transfer
				int start = -1, end = -1, alight = 0;
				for (int i = 0; i < N_r; i++) {
					if (route[i] == p) start = time[i];
					if (route[i] == St_Sc[p]) {
						end = time[i];
						break;
					}
					if (route[i] < P)if (start != -1 && St_Sc[p] == St_Sc[route[i]]) alight += tau_3 * (1 - W_ch[route[i]]) + tau_4 * W_ch[route[i]];
				}
				tt = end - start - alight;
				//cout << " travel time " << tt / 60 << " of student " << p << " on vehicle " << bus << endl;
				
			}
			else if ((ysol[p][0] == bus && ysol[p][1] != -1 ) || (p == location && !school_in)) { // if bus is the first bus in transfer
				int start = -1, end = -1, alight = 0;
				for (int i = 0; i < N_r; i++) {
					if (route[i] == p) start = time[i];
					if (route[i] == T) {
						end = time[i];
						break;
					}
				}
				tt = end - start;
				int bus2 = ysol[p][1];
				if (bus2 != -1) {
					start = -1, end = -1, alight = 0;
					int N_r2 = b_routes[bus2].size();
					for (int i = 0; i < N_r2; i++) {
						if (b_routes[bus2][i] == T) start = b_times[bus2][i];
						if (b_routes[bus2][i] == St_Sc[p]) {
							end = b_times[bus2][i];
							break;
						}
						if (b_routes[bus2][i] < P)if (start != -1 && St_Sc[p] == St_Sc[b_routes[bus2][i]])alight += tau_3 * (1 - W_ch[b_routes[bus2][i]]) + tau_4 * W_ch[b_routes[bus2][i]];
					}
					tt += (end - start - alight);
				}
			}

			if (tt > max_tt) max_tt = tt;
		}

		if (max_tt <= D_t) {
			//update route
			N_r = route.size();
			b_routes[bus].clear();
			b_times[bus].clear();
			for (int i = 0; i < N_r; i++) {
				b_routes[bus].push_back(route[i]);
				b_times[bus].push_back(time[i]);
				//cout << b_times[bus][i] << "\t";
			}
			//cout << endl;
			
		}
		//cout << endl;
		delete alights;
		return max_tt;
	}
	else {
		//cout << endl;
		delete alights;
		return 100000;
	}
}

inline int improve_route(const int n, vector<int>& cities, const int  N_it, const bool I_S, const int b, const bool initial) {
	// still need to work in the precendence rule
	//clock_t start_t = clock();
	int i, j, l;
	int newcost, currentcost;
	int dE;
	//double Time = 0, Timeit = 0;

	//int* tempdist = new int[n];
	//int* index = new int[n];
	int** closestn = new int* [L];
	if (n == L) {
		for (i = 0; i < L; i++) {
			closestn[i] = new int[n];
			for (j = 0; j < n; j++) {
				closestn[i][j] = closestS[i][j];
			}
			
		}
	}
	else {
		for (i = 0; i < L; i++) {
			closestn[i] = new int[n];
		}
		for (i = 0; i < n; i++) {
			//initialize and copy distance of cities to temp array
			l = 0;
			for (j = 0; j < L; j++) {
				if (std::count(cities.begin(), cities.end(), closestS[cities[i]][j])) { // if location j is in the list 
					closestn[cities[i]][l] = closestS[cities[i]][j];
					l++;
				}
			}
			/*
			for (j = 0; j < n; j++) {
				if (cities[j] != i) tempdist[j] = dist[i][cities[j]];
				else tempdist[j] = 100000;
				index[j] = cities[j];
			}

			//sort according to dist
			quickSort(index, tempdist, 0, n - 1);

			//keep track of best neighbors
			for (l = 0; l < n; l++) {
				closestn[i][l] = index[l]; //index
				//std::cout << index[l] << " ";
			}
			//std::cout << "\n";
			*/
		}
	}
	
	//cout << "check after closestN\n";
	/*
	if (!I_S) {
		cout << " route: " << endl;
		for (i = 0; i < n; i++) {
			cout << cities[i] << "\t";
		}
		cout << "\nNumber of stops " << n << endl;
		for (i = 0; i < L; i++) {
			if (closestn[i][0] != 100000) {
				cout << "i " << i << endl;
				for (j = 0; j < n; j++) {
					cout << closestn[i][j] << "\t";
				}
				cout << endl;
			}
		}
	}
	//*/
	if (I_S) { // If we need an intial feasible solution
		//Make feasible solution -------------------------------
		cities[0] = 0;
		bool in = false;
		i = cities[0]; //city number 
		j = 0; //index of cities List

		int tryi;
		//std::cout << "city: " << i << " ";
		//std::cout << "j: " << j << "\n";
		while (j + 1 < n) {
			//std::cout << "location: " << i << " ";
			//std::cout << "j: " << j + 1 << "\n";
			in = false;
			//try the nearest neighbor
			//check if this neighbor is already in list
			for (l = j - 1; l >= 0; l--) if (cities[l] == closestn[i][0]) {
				in = true;
				break;
			}
			if (!in) {
				i = closestn[i][0];
				cities[j + 1] = i;
			}
			//if already in, try next ones
			else {
				tryi = 1;
				while (in) {
					//std::cout << " not nearest \n";
					//try again nearest node
					//check if its already in 
					in = false;
					for (l = j - 1; l >= 0; l--) if (cities[l] == closestn[i][tryi]) {
						in = true;
						break;
					}
					//if (tryi == n) break;
					tryi++;
				}
				tryi--;
				i = closestn[i][tryi];
				cities[j + 1] = i;
			}
			
			j++;
		}
	}
	
	currentcost = 0;
	for (int i = 1; i < n; i++) {
		currentcost += TT[cities[i - 1]][cities[i]];
		if (cities[i] < P)currentcost += tau_1 * (1 - W_ch[cities[i]]) + tau_2 * W_ch[cities[i]];
	}
	//cout << endl;
	//SA PARAMETERS -----------------------------------------------------------------------------
	//number of nearest neighbors 
	int m = 20; if (m > n) m = n;
	double T_max = 100, Temp = T_max;
	double lam = 10;
	int sig = 0;


	//INITIATE PARAMETERS
	//std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, m - 1);
	std::uniform_int_distribution<int> distribution1(0, n / 3);
	//std::uniform_real_distribution<double> distribution2(0, 1);
	int edges[2], newedges[2]; // edges are the positions in the list cities 
	edges[0] = cities[n / 3];
	edges[1] = closestn[edges[0]][n / 2];
	int it = 0, itbad = 0;
	int start, end, temp, mid, nextindex, index1, index2;

	//to keep track of the best found yet
	int bestsol = 0;
	int* bestcities = new int[n];
	for (i = 0; i < n; i++) {
		bestcities[i] = cities[i];
	}
	bestsol = currentcost;
	vector<int> temp_cities;//temp cities
	for (i = 0; i < n; i++) {
		temp_cities.push_back(cities[i]);
	}
	
	//-------------------------------------------------------------------- Start Simulated Annealing  -------------------------------
	int iter = 0;
	vector<int> schools;
	for (i = 0; i < n; i++) { // delete schools in te route first;
		if (cities[i] >= P) {
			schools.push_back(cities[i]);
		}
	}

	int N_c = schools.size();

	bool trans_b = false;
	vector<int> tf_schools; // list of potential schools that need to be visited after T
	for (i = 0; i < P; i++) {
		if (ysol[i][1] == b) { //if this bus is a second bus in a transfer
			//make sure that T is visited before their assigned school
			trans_b = true;
			if (!std::count(tf_schools.begin(), tf_schools.end(), St_Sc[i])) {
				tf_schools.push_back(St_Sc[i]);
			}
		}
	}
	int N_tf = tf_schools.size();
	//cout << " check before iterations\n";
	while (iter < N_it) {
		iter++;
		//cout << "iteration " << iter << endl;
		//std::cout <<" Best cost "  << bestcost << "\n";

		//randomly select NEW edges for new iteration
		//if(!I_S)std::cout << " current edges " << edges[0] << " " << edges[1] << "\n";
		newedges[0] = closestn[edges[0]][distribution1(generator)];
		//if (!I_S)std::cout << " new edges " << newedges[0] ;
		newedges[1] = closestn[newedges[0]][distribution(generator)];
		//if (!I_S)std::cout  << " " << newedges[1] << "\n";
		//std::cout << " rand " << newi << " \n";

		if (newedges[0] == newedges[1]) continue;

		//new solution ---------------------------------------
		index1 = giveIndex(cities, newedges[0], n);
		index2 = giveIndex(cities, newedges[1], n);

		if (index1 > index2) {
			start = index2;
			end = index1;
		}
		else {
			start = index1;
			end = index2;
		}

		if (end == n) {
			end = n - 1;
			nextindex = n;
		}
		else nextindex = end + 1;

		temp_cities.clear();
		for (i = 0; i < n; i++) {
			temp_cities.push_back(cities[i]);
		}
		//2opt
		
		mid = (end - start) / 2;
		/*
		if (!I_S) {
			if (mid >= n) cout << " mid too large=" << mid << endl;
			if (start + mid >= n) cout << " start too large=" << start << endl;
			if (end >= n) cout << " end too large=" << end << endl;
			if (end + 1 - mid < 0) cout << " end too small=" << end << endl;
		}
		*/

		for (i = 1; i <= mid; i++) {
			temp = temp_cities[start + i];
			temp_cities[start + i] = temp_cities[end + 1 - i];
			temp_cities[end + 1 - i] = temp;
		}
		/*
		if (!initial) {
			cout << " before make feasible\n";
			bool stop_n = false;
			for (int y = 0; y < temp_cities.size(); y++) {
				for (int z = 0; z < temp_cities.size(); z++) {
					if (z != y && temp_cities[z] == temp_cities[y]) {
						stop_n = true;
						break;
					}
				}
				if (temp_cities[y] < P)cout << temp_cities[y] << "\t";
				else cout << "S_" << temp_cities[y] - P << "\t";
			}
			cout << endl;

			if (stop_n) {
				cout << "//////////////////////////////////////////////////////////////////////////// DOUBLESSSSSSSSSSSSSS\n";
				exit(0);
			}
		}
		//*/
		newcost = make_feasible(temp_cities, b, n, initial, schools, N_c, tf_schools, N_tf, trans_b);
		/*
		if (!initial) {
			cout << " after make feasible\n";
			bool stop_n = false;
			for (int y = 0; y < temp_cities.size(); y++) {
				for (int z = 0; z < temp_cities.size(); z++) {
					if (z != y && temp_cities[z] == temp_cities[y]) {
						stop_n = true;
						break;
					}
				}
				if (temp_cities[y] < P)cout << temp_cities[y] << "\t";
				else cout << "S_" << temp_cities[y] - P << "\t";
			}
			cout << endl;

			if (stop_n) {
				cout << "//////////////////////////////////////////////////////////////////////////// DOUBLESSSSSSSSSSSSSS\n";
				exit(0);
			}
		}
		//*/
		//if (!I_S)cout << "check 1\n";
		dE = newcost - currentcost;

		if (dE < 0) {

			//std::cout << " Found better dE =" << dE<< "   +++++++++++++++++++++++\n";
			currentcost = newcost;

			for (i = 0; i < n; i++) {
				cities[i] = temp_cities[i];
			}

			edges[0] = newedges[0];
			edges[1] = newedges[1];

			/*
			std::cout << "city " << newedges[0] << " and city " << newedges[1] << "\n";
			for (i = 0; i <= n; i++) {
				std::cout << cities[i] << " ";
			}
			std::cout << " \n";
			*/
			if (currentcost < bestsol) {
				bestsol = currentcost;
				for (i = 0; i < n; i++) {
					bestcities[i] = cities[i];
				}
			}
		}
		else if (exp(-dE / Temp) > one_to_zero(generator)) {
			//std::cout << " Found better dE =" << dE << "  -------------------------\n";
			itbad++;
			currentcost = newcost;

			for (i = 0; i < n; i++) {
				cities[i] = temp_cities[i];
			}

			edges[0] = newedges[0];
			edges[1] = newedges[1];

			/*
			std::cout << "city " << newedges[0] << " and city " << newedges[1] << "\n";
			for (i = 0; i <= n; i++) {
				std::cout << cities[i] << " ";
			}
			std::cout << " \n";
			*/
		}

		//temperature
		if (dE > 0)sig++;
		else if (dE < 0)sig = 0;
		Temp = T_max + lam * log(1 + sig);
		//std::cout << " nAccepted: " << (double)itbad / it * 100 << "% \n";

		/*
		//after temp is too low, just do a random search for remaining  computation time
		if (T < T_end) {
			T = 0.000000000000000000001;
			std::uniform_int_distribution<int> distribution(0, n-2);
			std::uniform_int_distribution<int> distribution1(0, n-2);
		}
		*/

		//next it
		it++;
		//Time = (float)(clock() - start_t) / CLK_TCK;

	}
	//cout << " check after itereations\n";

	for (i = 0; i < n; i++) {
		cities[i] = bestcities[i];
	}


	//std::cout << " nAccepted: " << (double) itbad/it*100 << "% \n";
	//std::cout << " T: " << T << "\n";
	//std::cout << " iterations: " << it << "\n";
	//clear memory
	//delete tempdist;
	//delete index;
	for (i = 0; i < L; i++) {
		delete closestn[i];
	}
	delete closestn;
	delete bestcities;
	return bestsol;
}

inline bool insert_location(vector<vector<int>>& routes, const int location, const int bus, vector<vector<int>>& times) {
	int curr_school;
	int indexr = routes[bus].size();
	//cout << "++++++ Trying to insert ";
	//cout << "student " << location << endl;
	int routeS = routes[bus].size();
	if (routeS == 0) {
		//cout << "      first location in route\n";
		routes[bus].insert(routes[bus].begin(), location);
		routes[bus].push_back(St_Sc[location]);

		times[bus].clear();
		//cout << "first school " << xsol[fsch] - P << endl;
		times[bus].insert(times[bus].begin(), t_ea[St_Sc[location] - P]); // schedule first school to arrive asap 
		times[bus].insert(times[bus].begin(), times[bus][0] - TT[St_Sc[location]][location] - tau_1 * (1 - W_ch[location]) - tau_2 * W_ch[location]);
		return true;
	}
	curr_school = St_Sc[location];
	auto itr = find(best_route, best_route + L, location);
	int current = distance(best_route, itr);
	// find best place to insert location
	bool in = false;
	int  ttemp;
	//cout << "current " << current << endl;
	for (int j = 0; j < routeS; j++) {
		if (routes[bus][j] != location) {
			auto itr = find(best_route, best_route + L, routes[bus][j]);
			ttemp = distance(best_route, itr);
			//cout << "ttemp " << ttemp << endl;
			if (ttemp > current) {
				indexr = j;
				break;
			}
		}
		else {
			return false;
		}
	}
	
	
	//cout << "indexr " << indexr << " route length " << routeS << endl;
	bool school_in = true, p_in = true;
	int alt_school = -1;
	int tt = max_TT(routes, location, indexr, school_in, alt_school, bus, times, p_in);
	//cout << " max tt : " << tt / 60 << endl;
	
	if (tt > D_t) {
		//cout << " NOT okay for travel time\n";
		return false;
	}
	else {
		//cout << " OKay for travel time\n";
		return true;
	}
}

inline int transfer_possible(vector<vector<int>>& times, vector<vector<int>>& routes, const int n_v, const int location, const int b_old, const int b_old2) { // returns transfer bus and adjust the routes
	
	int dummy1, dummy2;
	//if (!count(routes[n_v].begin(), routes[n_v].end(), T)) return -1; // transfer stop not visited by the first bus, no transfer!!!!
	int currentschool = St_Sc[location];
	// determine where to insert passenger p 
	int indexr = routes[n_v].size();
	//cout << "++++++ Trying to insert ";
	//cout << "student " << location << endl;
	int routeS = routes[n_v].size();

	auto itr = find(best_route, best_route + L, location);
	int current = distance(best_route, itr);
	// find best place to insert location
	int  ttemp = 0;
	//cout << "current " << current << endl;
	for (int j = 0; j < routeS; j++) {
		if (routes[n_v][j] != location) {
			auto itr = find(best_route, best_route + L, routes[n_v][j]);
			ttemp = distance(best_route, itr);
			//cout << "ttemp " << ttemp << endl;
			if (ttemp > current) {
				indexr = j;
				break;
			}
		}
	}
	vector<int> time;
	int N_t = times[n_v].size();
	for (int i = 0; i < N_t; i++) {
		time.push_back(times[n_v][i]);
	}

	int B = routes.size(), N_r = 0, i_transfer = -1, i_target = -1, tt1 = 0, tt2 = 0, curr_cap = 0, l = 0;
	//cout << "check\n";
	// Check if you can add the passenger without adding the school
	bool school_in1 = false, school_in2 = true, p_in1 = true, p_in2 = false;
	int alt_school2 = -1;
	tt1 = max_TT(routes, location, indexr, school_in1, T, n_v, times, p_in1);
	cout << "MAX travel time on first vehicle " << tt1 / 60.0 << endl;

	int tts1 = 0, tts2 = 0, tts = 0;
	int it_trans1 = -1, tt_trans1 = -1, tt_trans2 = -1;
	//int alight1 = 0, alight2 = 0, boarding1 = 0, boarding2 = 0;
	int past1 = 0, past2 = 0, future1 = 0, future2 = 0, shiftP = 0, shiftF = 0;

	

	if (tt1 <= D_t) {
		routeS = routes[n_v].size();
		// travel time of student p on fist vehicle
		tts1 = 0;
		routeS = routes[n_v].size();
		for (int j = 0; j < routeS; j++) {
			if (routes[n_v][j] == location) tts = times[n_v][j];
			if (routes[n_v][j] == T) {
				tts1 = times[n_v][j];
				break;
			}
		}
		tts1 -= tts;

		
		auto itr = find(routes[n_v].begin(), routes[n_v].end(), T);
		it_trans1 = itr - routes[n_v].begin();
		tt_trans1 = times[n_v][it_trans1]; // current arrival time at transfer 

		// max time to shift back or forwards
		give_max_shifts(routes[n_v], times[n_v], routeS,  past1, future1, dummy1, dummy2);
		// see if other buses with enough capacity arrive at schools that are visited by n_v
		for (int b = 0; b < B; b++) {
			if (b != n_v && b != b_old && b != b_old2) {
				cout << "~~~~~ possible second vehicle " << b << endl;
				i_target = -1;
				i_transfer = -1;
				N_r = routes[b].size();
				// determine if the bus visists target AND transfer school and also determin current cap
				for (int i = 0; i < N_r; i++) {
					if (currentschool == routes[b][i]) {
						i_target = i;
					}
					else if (T == routes[b][i]) {
						i_transfer = i;
					}
				}
				if (i_target != -1 && i_transfer != -1 && i_transfer < i_target) { // if it does and if it also visits other schools of n_v and before the target school
					cout << " second route OK\n";
					curr_cap = 0;

					// students picked up before transfer and that not alight at transfer
					for (int i = 0; i <= i_transfer; i++) {
						if (routes[b][i] < P)curr_cap += (3 * W_ch[routes[b][i]] + 1);
						else {
							for (int y = 0; y < P; y++) { // passengers alighting before the transfer and at the transfer
								if ((ysol[y][0] == b) && St_Sc[y] == routes[b][i]) curr_cap -= (3 * W_ch[y] + 1);
							}
						}
					}
					// students picked up after transfer
					for (int i = i_transfer + 1; i < N_r; i++) {
						if (routes[b][i] < P)curr_cap += (3 * W_ch[routes[b][i]] + 1); 
					}
					// students transferring over to this vehicle
					for (int y = 0; y < P; y++) {
						if (ysol[y][1] == b) curr_cap += (3 * W_ch[y] + 1);
					}
					if (curr_cap + (3 * W_ch[location] + 1) <= C_b) {
						cout << " capacity OK\n";
						//check if u can add student to transfer bus b
						// find best place to insert location
						indexr = routes[b].size();
						for (int j = 0; j < N_r; j++) {
							if (routes[b][j] != location) {
								auto itr = find(best_route, best_route +L, routes[b][j]);
								ttemp = distance(best_route, itr);
								//cout << "ttemp " << ttemp << endl;
								if (ttemp > current) {
									indexr = j;
									break;
								}
							}
						}

						tt2 = max_TT( routes, location, indexr, school_in2, alt_school2, b, times, p_in2);
						cout << "MAX travel time on second vehicle " << tt2 / 60.0 << endl;
						if (tt2 <= D_t) { // if yes
							// check if you student p isnt travelling for too long
							N_r = routes[b].size();
							tts2 = - 1, tts = -1;
							int alight2 = 0;
							for (int j = 0; j < N_r; j++) {
								if (routes[b][j] == T) tts = times[b][j];
								if (routes[b][j] == St_Sc[location]) {
									tts2 = times[b][j];
									break;
								}
								if (routes[b][j] < P)if (tts != -1 && St_Sc[location] == St_Sc[routes[b][j]])alight2 += tau_3 * (1 - W_ch[routes[b][j]]) + tau_4 * W_ch[routes[b][j]];
							}
							tts2 -= tts;
							tts2 -= alight2;
							
							tts = tts1 + tts2;
							//cout << " Check 4\n";
							if (tts <= D_t) {
								cout << "Travel time of passenger " << tts / 60.0 << " OK\n";
								
								// in that case, see if they can transport student with feasible timetable: n_v needs to (be able to) arrive before arrival of second bus
								tt_trans2 = times[b][i_transfer]; // current arrival time at transfer
								bool feasible = give_max_shifts(routes[b], times[b], N_r,  past2, future2, dummy1, dummy2);
								if (!feasible) {
									cout << "INFEASIBLE timetable for second vehicle " << b << " in transfer when trying to insert student " << location << " with assigned school S_" << St_Sc[location] - P << endl;
									for (int y = 0; y < routes[b].size(); y++) {
										int seconds = times[b][y];
										int minutes = seconds / 60;
										minutes = minutes % 60;
										string smin = "";
										if (minutes >= 10) {
											smin = to_string(minutes);
										}
										else smin = "0" + to_string(minutes);
										int hours = seconds / 3600;
										if (routes[b][y] < P)cout << hours << ":" << smin << " (" << routes[b][y] << ")\t";
										else cout << hours << ":" << smin << " (S_" << routes[b][y] - P << ")\t";
									}
									exit(0);
								}
								
								if (tt_trans2 >= tt_trans1) {// if the bues arrives before transfer, good, just see we can minimise travel times
									cout << "Vehicle " << b << " arrives at " << tt_trans2/3600.0 << ", after vehicle " << n_v << " at " << tt_trans1/3600.0 << endl;
									//cout << "v_" << n_v << "--> max FW : " << future1/60.0 << endl;
									// see we can minimise travel times
									shiftF = min(tt_trans2 - tt_trans1, future1);
									shiftP = min(tt_trans2 - tt_trans1, past2);

									if (shiftF == tt_trans2 - tt_trans1) {	// by only moving n_v to the future
										cout << " move tt of vehicle " << n_v << " to the future, by " << shiftF/60.0 << " min\n";
										//make shift to tt of n_v
										routeS = times[n_v].size();
										for (int l = 0; l < routeS; l++) {
											times[n_v][l] += shiftF;
										}
										
										return b; // if feasible transfer found then stop search
									}
									else if (shiftP == tt_trans2 - tt_trans1) {	// by only moving b to the past
										cout << " move tt of vehicle " << b << " to the past, by " << shiftP / 60.0 << " min\n";
										//make shift to tt of b
										N_r = times[b].size();
										for (int l = 0; l < N_r; l++) {
											times[b][l] -= shiftP;
										}
										return b; // if feasible transfer found then stop search
									}
									else  { // by moving both
										shiftP = min(shiftP, D_t - (tts + tt_trans2 - tt_trans1 - shiftF)); // make sure u dont overcorrect
										cout << " TRY TO move tt of vehicle " << b << " to the past by " << shiftP/60.0 << " min and vehicle " << n_v << " to the future, by " << shiftF / 60.0 << " min\n";
										if (max(max(tts, tt2), tt1) + tt_trans2 - tt_trans1 - shiftP - shiftF <= D_t) {
											//make shift to tt of b and n_v
											routeS = times[n_v].size();
											for (int l = 0; l < routeS; l++) {
												times[n_v][l] += shiftF;
											}
											N_r = times[b].size();
											for (int l = 0; l < N_r; l++) {
												times[b][l] -= shiftP;
											}
											return b; // if feasible transfer found then stop search
										}
										else {
											cout << "total timetable not feasible\n";
											auto itr0 = find(routes[b].begin(), routes[b].end(), T);
											int j = itr0 - routes[b].begin();
											times[b][j] -= (tau_3 * (1 - W_ch[location]) + tau_4 * W_ch[location]);
										}
									}
									
												
								}
								else { // else if its still possible to shift n_v to the past or b to the future or both
									// see we can minimise travel times
									cout << "Vehicle " << b << " arrives at " << tt_trans2/3600.0 << " so it DOES NOT arrive after vehicle " << n_v << " at " << tt_trans2/3600.0<<  endl;
									shiftF = min(tt_trans1 - tt_trans2, future2);
									shiftP = min(tt_trans1 - tt_trans2, past1);

									/*
									cout << " vehicle " << n_v << endl;
									for (int y = 0; y < times[n_v].size(); y++) {
										int seconds = times[n_v][y];
										int minutes = seconds / 60;
										minutes = minutes % 60;
										string smin = "";
										if (minutes >= 10) {
											smin = to_string(minutes);
										}
										else smin = "0" + to_string(minutes);
										int hours = seconds / 3600;
										if (routes[n_v][y] < P)cout << hours << ":" << smin << " (" << routes[n_v][y] << ")\t";
										else cout << hours << ":" << smin << " (S_" << routes[n_v][y] - P << ")\t";
									}
									cout << "\n veicle " << b << endl;
									for (int y = 0; y < times[b].size(); y++) {
										int seconds = times[b][y];
										int minutes = seconds / 60;
										minutes = minutes % 60;
										string smin = "";
										if (minutes >= 10) {
											smin = to_string(minutes);
										}
										else smin = "0" + to_string(minutes);
										int hours = seconds / 3600;
										if (routes[b][y] < P)cout << hours << ":" << smin << " (" << routes[b][y] << ")\t";
										else cout << hours << ":" << smin << " (S_" << routes[b][y] - P << ")\t";
									}
									//*/
									if (shiftP == tt_trans1 - tt_trans2) {	// by only moving n_v to the past
										//make shift to tt of n_v 
										cout << " move tt of vehicle " << n_v << " to the past, by " << shiftP / 60.0 << " min\n";
										
										routeS = times[n_v].size();
										for (int l = 0; l < routeS; l++) {
											times[n_v][l] -= shiftP;
										}
										return b; // if feasible transfer found then stop search
									}
									else if (shiftF == tt_trans1 - tt_trans2) {	// by only moving b to the future
										cout << " move tt of vehicle " << b << " to the future, by " << shiftF / 60.0 << " min\n";
										//make shift to tt of b
								
										N_r = times[b].size();
										for (int l = 0; l < N_r; l++) {
											times[b][l] += shiftF;
										}
										return b; // if feasible transfer found then stop search
									}
									else { // by moving both
										shiftF = min(shiftP, D_t - (tts + tt_trans1 - tt_trans2 - shiftP)); // make sure u dont overcorrect
										cout << " TRY TO move tt of vehicle " << n_v << " to the past, by " << shiftP / 60.0 << " min and vehicle " << b << " to the future, by " << shiftF / 60.0 << " min\n";
										cout << "waiting time: " << (tt_trans1 - tt_trans2) / 60.0 - (shiftP + shiftF) / 60.0 << " min\n";
										if ( tt_trans1 - tt_trans2 - shiftP - shiftF == 0) { // only if there is no waiting time, otherwaise the second transfer departs before the first transfer
											//make shift to tt of b and n_v
											
											routeS = times[n_v].size();
											for (int l = 0; l < routeS; l++) {
												times[n_v][l] -= shiftP;
											}
											
											N_r = times[b].size();
											for (int l = 0; l < N_r; l++) {
												times[b][l] += shiftF;
											}
											return b; // if feasible transfer found then stop search
										}
										else {
											cout << "total timetable not feasible\n";
											auto itr0 = find(routes[b].begin(), routes[b].end(), T);
											int j = itr0 - routes[b].begin();
											times[b][j] -= (tau_3 * (1 - W_ch[location]) + tau_4 * W_ch[location]);
										}
									}
								}
								

							}
							else {
								cout << "total travel time too high: " << tts/60.0 << "\n";
								//cout << " check\n";
								auto itr0 = find(routes[b].begin(), routes[b].end(), T);
								int j = itr0 - routes[b].begin();
								times[b][j] -= (tau_3 * (1 - W_ch[location]) + tau_4 * W_ch[location]);
							}

						}
					}
				}
			}
		}
	}
	else return -1;
	// if not feasible --> go back to vehicle 
	auto itr0 = find(routes[n_v].begin(), routes[n_v].end(), location);
	routes[n_v].erase(itr0);

	N_t = time.size();
	times[n_v].clear();
	for (int i = 0; i < N_t; i++) {
		times[n_v].push_back(time[i]);
	}
	return -1;

	
}

inline int max_TT0(const int bus, const int bus_t, const int wait, vector<vector<int>>& b_routes, vector<vector<int>>& b_times) {
	// calculate the max tt
	int max_tt = 0, tt = 0;
	int N_r = b_routes[bus].size();
	int p, i, start, end, alight, bus2, N_r2;
	for (p = 0; p < P; p++) {
		tt = 0;
		if ((ysol[p][0] == bus) && ysol[p][1] == -1) { // if bus is a vehicle and no transfer
			start = -1, end = -1, alight = 0;
			for (i = 0; i < N_r; i++) {
				if (b_routes[bus][i] == p) start = b_times[bus][i];
				if (b_routes[bus][i] == St_Sc[p]) {
					end = b_times[bus][i];
					break;
				}
				if (b_routes[bus][i] < P)if (start != -1 && St_Sc[p] == St_Sc[b_routes[bus][i]]) alight += tau_3 * (1 - W_ch[b_routes[bus][i]]) + tau_4 * W_ch[b_routes[bus][i]];
			}
			tt = end - start - alight;
			//cout << " travel time " << tt / 60 << " of student " << p << " on vehicle " << bus << endl;

		}
		else if ((ysol[p][0] == bus && ysol[p][1] != -1)) { // if bus is the first bus in transfer
			start = -1, end = -1, alight = 0;
			for (i = 0; i < N_r; i++) {
				if (b_routes[bus][i] == p) start = b_times[bus][i];
				if (b_routes[bus][i] == T) {
					end = b_times[bus][i];
					break;
				}
			}
			tt = end - start;
			bus2 = ysol[p][1];
			if (bus2 != -1) {
				start = -1, end = -1, alight = 0;
				N_r2 = b_routes[bus2].size();
				for (i = 0; i < N_r2; i++) {
					if (b_routes[bus2][i] == T) start = b_times[bus2][i];
					if (b_routes[bus2][i] == St_Sc[p]) {
						end = b_times[bus2][i];
						break;
					}
					if (b_routes[bus2][i] < P)if (start != -1 && St_Sc[p] == St_Sc[b_routes[bus2][i]])alight += tau_3 * (1 - W_ch[b_routes[bus2][i]]) + tau_4 * W_ch[b_routes[bus2][i]];
				}
				tt += (end - start - alight);
			}
		}
		if (ysol[p][0] == bus && ysol[p][1] == bus_t) tt += wait;

		if (tt > max_tt) max_tt = tt;
	}
	return max_tt;
}

inline bool correct_tt_transfers(const int b, const int b_tf2, vector<vector<int>>& xsol, vector<vector<int>>& Dsol) {
	auto it1 = find(xsol[b].begin(), xsol[b].end(), T);
	int i_t1 = it1 - Dsol[b].begin();
	int tt_trans1 = Dsol[b][i_t1]; // dep of the first transfer

	auto it2 = find(xsol[b_tf2].begin(), xsol[b_tf2].end(), T);
	int i_t2 = it1 - Dsol[b_tf2].begin();
	int tt_trans2 = Dsol[b_tf2][i_t2]; // dep of the first transfer

	int past1, past2, future1, future2, dummy1, dummy2;
	int routeS = xsol[b].size();
	give_max_shifts(xsol[b], Dsol[b], routeS,  past1, future1, dummy1, dummy2);
	int N_r = xsol[b_tf2].size();
	give_max_shifts(xsol[b_tf2], Dsol[b_tf2], N_r,  past2, future2, dummy1, dummy2);

	int shiftP, shiftF, mtt1, wait, l;

	if (tt_trans2 >= tt_trans1) {// if the bus arrives before transfer, good, just see we can minimise travel times
		//cout << "    v_" << b_tf2 << " departs after v_" << b << endl;
		shiftF = min(tt_trans2 - tt_trans1, future1);
		shiftP = min(tt_trans2 - tt_trans1, past2);

		if (shiftF == tt_trans2 - tt_trans1) {	// by only moving n_v to the future
			//make shift to tt of n_v
			//cout << " move tt of vehicle " << b << " to the future, by " << shiftF / 60.0 << " min\n";
			routeS = Dsol[b].size();
			for (l = 0; l < routeS; l++) {
				Dsol[b][l] += shiftF;
			}

			return true; // if feasible transfer found then stop search
		}
		else if (shiftP == tt_trans2 - tt_trans1) {	// by only moving b to the past
			//cout << " move tt of vehicle " << b_tf2 << " to the past, by " << shiftP / 60.0 << " min\n";
			//make shift to tt of b
			N_r = Dsol[b_tf2].size();
			for (l = 0; l < N_r; l++) {
				Dsol[b_tf2][l] -= shiftP;
			}
			return true; // if feasible transfer found then stop search
		}
		else { // by moving both
			if (shiftP + shiftF >= (tt_trans2 - tt_trans1)) { // then feasible
				if (shiftP > shiftF) {
					shiftP = (tt_trans2 - tt_trans1) - shiftF;
				}
				else {
					shiftF = (tt_trans2 - tt_trans1) - shiftP;
				}
				routeS = Dsol[b].size();
				for (l = 0; l < routeS; l++) {
					Dsol[b][l] += shiftF;
				}
				N_r = Dsol[b_tf2].size();
				for (l = 0; l < N_r; l++) {
					Dsol[b_tf2][l] -= shiftP;
				}
				return true; // if feasible transfer found then stop search
			}
			else { // there will be waiting time
				//cout << " move tt of vehicle " << b << " to the future, by " << shiftF / 60.0 << " min and vehicle " << b_tf2 << " to the past, by " << shiftP / 60.0 << " min\n";
				wait = (tt_trans2 - tt_trans1) - (shiftP + shiftF);
				mtt1 = max_TT0(b, b_tf2, wait, xsol, Dsol);
				if (mtt1 <= D_t) {
					routeS = Dsol[b].size();
					for (l = 0; l < routeS; l++) {
						Dsol[b][l] += shiftF;
					}
					N_r = Dsol[b_tf2].size();
					for (l = 0; l < N_r; l++) {
						Dsol[b_tf2][l] -= shiftP;
					}
					return false;
				}
				else return false;
			}
		}
	}
	else { // else if its still possible to shift n_v to the past or b to the future or both
		shiftF = min(tt_trans1 - tt_trans2, future2);
		shiftP = min(tt_trans1 - tt_trans2, past1);
		//cout << "    v_" << b_tf2 << " departs before v_" << b << endl;

		if (shiftP == tt_trans1 - tt_trans2) {	// by only moving n_v to the past
			//make shift to tt of n_v 
			//cout << " move tt of vehicle " << b << " to the past, by " << shiftP / 60.0 << " min\n";
			routeS = Dsol[b].size();
			for (l = 0; l < routeS; l++) {
				Dsol[b][l] -= shiftP;
			}
			return true; // if feasible transfer found then stop search
		}
		else if (shiftF == tt_trans1 - tt_trans2) {	// by only moving b to the future
			//cout << " move tt of vehicle " << b_tf2 << " to the future, by " << shiftF / 60.0 << " min\n";
			//make shift to tt of b

			N_r = Dsol[b_tf2].size();
			for (l = 0; l < N_r; l++) {
				Dsol[b_tf2][l] += shiftF;
			}
			return true; // if feasible transfer found then stop search
		}
		else { // by moving both
			//cout << " move tt of vehicle " << b_tf2 << " to the future, by " << shiftF / 60.0 << " min and vehicle " << b << " to the past, by " << shiftP / 60.0 << " min\n";
			if (shiftP + shiftF >= (tt_trans2 - tt_trans1)) { // then feasible
				if (shiftP > shiftF) {
					shiftP = (tt_trans2 - tt_trans1) - shiftF;
				}
				else {
					shiftF = (tt_trans2 - tt_trans1) - shiftP;
				}
				routeS = Dsol[b_tf2].size();
				for (l = 0; l < routeS; l++) {
					Dsol[b_tf2][l] += shiftF;
				}
				N_r = Dsol[b].size();
				for (l = 0; l < N_r; l++) {
					Dsol[b][l] -= shiftP;
				}
				return true; // if feasible transfer found then stop search
			}
			else { // there will be waiting time
				wait = (tt_trans2 - tt_trans1) - (shiftP + shiftF);
				mtt1 = max_TT0(b, b_tf2, wait, xsol, Dsol);
				if (mtt1 <= D_t) {
					routeS = Dsol[b].size();
					for (l = 0; l < routeS; l++) {
						Dsol[b][l] -= shiftP;
					}
					N_r = Dsol[b_tf2].size();
					for (l = 0; l < N_r; l++) {
						Dsol[b_tf2][l] += shiftF;
					}
					return false;
				}
				else return false;
			}
		}
	}
}

inline int isFEAS(vector<vector<int>>& b_xsol, vector<vector<int>>& b_Dsol) {
	int second_cost = 0;
	int b, t, s, p;
	int B = b_xsol.size(), N_r = 0;
	int tt = 0;
	int c_p = 0, c_a = 0, c_p2 = 0;

	for (b = 0; b < B; b++) {
		c_p = 0;
		N_r = b_xsol[b].size();
		for (s = 0; s < N_r; s++) {
			for (t = 0; t < N_r; t++) {
				if (s != t && b_xsol[b][s] == b_xsol[b][t]) {
					return -1;
				}
			}
			if (b_xsol[b][s] < P) { // Location is student
				c_p += (3 * W_ch[b_xsol[b][s]] + 1);
			}
			else { // location is school now
				if (b_Dsol[b][s] < t_ea[b_xsol[b][s] - P]) {
					return -1;
				}
				else if (b_Dsol[b][s] > t_la[b_xsol[b][s] - P]) {
					return -1;
				}
				if (c_p > C_b) {
					return -1;
				}
				c_a = 0;
				for (t = s - 1; t >= 0; t--) { // correct capacity of studnets assigned to this school 
					if (b_xsol[b][t] < P) if (St_Sc[b_xsol[b][t]] == b_xsol[b][s]) c_a += (3 * W_ch[b_xsol[b][t]] + 1); //students get off
				}
				c_p -= c_a;
				if (b_xsol[b][s] == T) {
					c_p2 = 0;
					for (p = 0; p < P; p++) {
						if (ysol[p][1] == b) {
							c_p2 += (3 * W_ch[p] + 1);  // students that will get on
						}
					}
					if (c_p + c_p2 > C_b) {
						return -1;
					}
				}
			}
		}
	}
	
	//Objective function value
	int* scosts = new int[B];
	for (b = 0; b < B; b++) {
		scosts[b] = -1;
	}
	int tt2 = 0;
	int start, end, alight, start2, end2, alight2, N_r2, b2, i;
	bool dess;
	for (p = 0; p < P; p++) {
		if (ysol[p][0] > -1) {
			b = ysol[p][0];
			b2 = ysol[p][1];
			N_r = b_xsol[b].size();
			//if not transfer
			if (b2 == -1) {
				start = -1, end = -1, alight = 0;
				for (i = 0; i < N_r; i++) {
					if (b_xsol[b][i] == p) start = b_Dsol[b][i];
					if (b_xsol[b][i] == St_Sc[p]) {
						end = b_Dsol[b][i];
						break;
					}
					if (b_xsol[b][i] < P)if (start != -1 && St_Sc[p] == St_Sc[b_xsol[b][i]])alight += tau_3 * (1 - W_ch[b_xsol[b][i]]) + tau_4 * W_ch[b_xsol[b][i]];
				}
				tt = end - start - alight;
				if (tt > scosts[b]) scosts[b] = tt;
				if (tt > D_t) {
					delete scosts;
					return -1;
				}
				else if (end < start) {
					delete scosts;
					return -1;
				}
			}
			else { // if transfer
				start = -1, end = -1, alight = 0;
				for (i = 0; i < N_r; i++) {
					if (b_xsol[b][i] == p) start = b_Dsol[b][i];
					if (b_xsol[b][i] == T) {
						end = b_Dsol[b][i];
						break;
					}
				}
				tt = end - start;
				if (tt > scosts[b]) scosts[b] = tt;
				if (tt < 0) {
					delete scosts;
					return -1;
				}
				N_r2 = b_xsol[b2].size();
				start2 = -1, end2 = -1, alight2 = 0;
				for (i = 0; i < N_r2; i++) {
					if (b_xsol[b2][i] == T) start2 = b_Dsol[b2][i];
					if (b_xsol[b2][i] == St_Sc[p]) {
						end2 = b_Dsol[b2][i];
						break;
					}
					if (b_xsol[b2][i] < P)if (start2 != -1 && St_Sc[p] == St_Sc[b_xsol[b2][i]])alight2 += tau_3 * (1 - W_ch[b_xsol[b2][i]]) + tau_4 * W_ch[b_xsol[b2][i]];
				}
				tt2 = (end2 - start2 - alight2);

				if (tt2 > scosts[b2]) scosts[b2] = tt2;
				if (tt + tt2 + abs(end - start2) > D_t) {
					dess = false;
					//cout << " too long travel time but lets try to fix it\n";
					dess = correct_tt_transfers(b, b2, b_xsol, b_Dsol);
					
					if (!dess) {
						delete scosts;
						return -1;
					}
				}
				if (start2 > end2) {
					delete scosts;
					return -1;
				}
				tt = tt + tt2 + abs(end - start2);
			}
		}
		else {
			delete scosts;
			return -1;
		}
	}
	
	second_cost = 0;
	for (b = 0; b < B; b++) {
		second_cost += scosts[b];
	}

	delete scosts;
	return second_cost;
}

inline void destroy_solution(const int b, const int b2, const int p, vector<vector<int>>& xsol, vector<vector<int>>& Dsol) {
	bool dess = true;
	//destroy student p
	ysol[p][0] = -1; // assign old vehicle to student
	ysol[p][1] = -1;

	auto it = find(xsol[b].begin(), xsol[b].end(), p);
	int l = it - xsol[b].begin();
	Dsol[b].erase(l + Dsol[b].begin()); //remove from timetable
	xsol[b].erase(it); //remove from route
	int j = l - 1;
	//cout << "l: " << l <<" --> place of the student " << endl;
	int N_b = 0, N_r = 0, N_sc = 0;

	//cout << " start destroying!!!\n";
	// check to see if its the last student needing to go to the assigned school
	N_r = xsol[b].size();
	N_b = 0;
	int trans = 0;
	vector<int> schools;
	for (int i = 0; i < N_r; i++) {
		if (xsol[b][i] < P) {
			if (St_Sc[xsol[b][i]] == St_Sc[p])N_b++; // for students wanting to arrive at their school 
			if (St_Sc[p] == T && ysol[xsol[b][i]][1] != -1) trans++; // for students wanting to use the bus to transfer at T
		}
		else {
			N_sc++;
			schools.push_back(xsol[b][i]);
		}
	}
	//int b_tf1 = -2, b_tf2 = -2;
	for (int s = 0; s < P; s++) {
		if (s != p) {
			if (ysol[s][1] == b) { // if this is the second vehicle in transfer for s
				//b_tf1 = ysol[s][0];
				trans++;
				if (St_Sc[p] == St_Sc[s]) N_b++;
			}
		}
	}
	
	if (N_b + trans == 0) {
		//cout << " no more students visiting S_" << St_Sc[p] - P << " or transferring at S_" << T - P << endl;
		auto it2 = find(xsol[b].begin(), xsol[b].end(), St_Sc[p]);
		if (it2 != xsol[b].end()) { // if school in the route
			int l2 = it2 - xsol[b].begin();
			Dsol[b].erase(l2 + Dsol[b].begin()); //remove school from timetable
			xsol[b].erase(it2); //remove school from route
			auto it3 = find(schools.begin(), schools.end(), St_Sc[p]);
			schools.erase(it3);
			N_sc--;
			if (l2 - 1 > j) {
				j = l2 - 1;
			}
		}
	}
	N_r = xsol[b].size();
	if (j >= 0 && N_r != 0) { // if its the fisrst stop, no need to ajust timemtable
		//cout << " lets start the descruction, j = " << j << " lenght of route " << xsol[b].size() << "\n";
		int* alights = new int[N_sc];
		int alight = 0;
		N_r = xsol[b].size();
		for (int s = 0; s < N_sc; s++) {
			alight = 0;
			for (int i = 0; i < N_r; i++) {
				if (xsol[b][i] < P) {
					if (St_Sc[xsol[b][i]] == schools[s] || (schools[s] == T && ysol[xsol[b][i]][1] != -1)) {
						alight += tau_3 * (1 - W_ch[xsol[b][i]]) + tau_4 * W_ch[xsol[b][i]];
					}
				}
			}
			if (schools[s] == T) { // boarding time
				for (int p = 0; p < P; p++) {
					if (ysol[p][1] == b) {
						alight += tau_3 * (1 - W_ch[p]) + tau_4 * W_ch[p];
					}
				}
			}
			alights[s] = alight;
		}
		time_table2(xsol[b], Dsol[b], b, alights, N_sc);
		delete alights;
	}
	//cout << " route destroyed\n";
	

	//undo transfer as well
	if (b2 != -1) {
		//cout << " also destroy transfer\n";
		// check to see if its the last student needing to go to the assigned school
		N_r = xsol[b2].size();
		N_b = 0;
		trans = 0;
		schools.clear();
		N_sc = 0;
		for (int i = 0; i < N_r; i++) {
			if (xsol[b2][i] < P) {
				if (St_Sc[xsol[b2][i]] == St_Sc[p])N_b++; // for students wanting to arrive at their school 
				if (St_Sc[xsol[b2][i]] == T) trans++; // for students wanting to use the bus to transfer at T
			}
			else {
				N_sc++;
				schools.push_back(xsol[b][i]);
			}
		}
		//b_tf1 = -2, b_tf2 = -2;
		for (int s = 0; s < P; s++) {
			if (s != p) {
				if (ysol[s][1] == b2) { // if this is the second vehicle in transfer for s
					//b_tf1 = ysol[s][0];
					trans++;
					if (St_Sc[p] == St_Sc[s]) N_b++;
				}
			}
		}

		if (N_b == 0) {
			//cout << " no more students visiting " << St_Sc[p] << endl;
			auto it2 = find(xsol[b2].begin(), xsol[b2].end(), St_Sc[p]);
			if (it2 != xsol[b2].end()) { // if school in the route
				int l2 = it2 - xsol[b2].begin();
				Dsol[b2].erase(l2 + Dsol[b2].begin()); //remove school from timetable
				xsol[b2].erase(it2); //remove school from route
				//remove school from schools 
				auto it3 = find(schools.begin(), schools.end(), St_Sc[p]);
				schools.erase(it3);
				N_sc--;
				if (l2 - 1 > j) j = l2 - 1;
			}
		}

		if (trans == 0) {
			//cout << " no more students visiting " << St_Sc[p] << endl;
			auto it2 = find(xsol[b2].begin(), xsol[b2].end(), T);
			if (it2 != xsol[b2].end()) { // if school in the route
				int l2 = it2 - xsol[b2].begin();
				Dsol[b2].erase(l2 + Dsol[b2].begin()); //remove school from timetable
				xsol[b2].erase(it2); //remove school from route
				//remove school from schools 
				auto it3 = find(schools.begin(), schools.end(), T);
				schools.erase(it3);
				N_sc--;
				if (l2 - 1 > j) j = l2 - 1;
			}
		}
		//cout << " route destroyed\n";
		N_r = xsol[b2].size();
		if (N_r != 0) {
			int alight = 0;
			int* alights = new int[N_sc];
			for (int s = 0; s < N_sc; s++) {
				alight = 0;
				for (int i = 0; i < N_r; i++) {
					if (xsol[b2][i] < P) {
						if (St_Sc[xsol[b2][i]] == schools[s] || (schools[s] == T && ysol[xsol[b2][i]][1] != -1)) {
							alight += tau_3 * (1 - W_ch[xsol[b2][i]]) + tau_4 * W_ch[xsol[b2][i]];
						}
					}
				}
				if (schools[s] == T) { // boarding time
					for (int p = 0; p < P; p++) {
						if (ysol[p][1] == b2) {
							alight += tau_3 * (1 - W_ch[p]) + tau_4 * W_ch[p];
						}
					}
				}
				alights[s] = alight;
			}
			time_table2(xsol[b2], Dsol[b2], b2, alights, N_sc);
			delete alights;
		}
	}
	//cout << "End destruction\n";
}

inline bool repair_solution(const int b, const int b2, const int p, const int n_v, const int s, vector<vector<int>>& xsol, vector<vector<int>>& Dsol) {
	int curr_capT = 0; // current cap before transfer
	int l = -1, x, z, y, N_imp = 5000;
	bool I_S = false, initial = false;

	int i_t = -1, b_trans = -1;
	int N_r = xsol[n_v].size();
	for (l = N_r - 1; l >= 0; l--) {
		if (xsol[n_v][l] == s) break;
	}
	for (i_t = N_r - 1; i_t >= 0; i_t--) {
		if (xsol[n_v][i_t] == T) {
			break;
		}
	}

	int curr_cap = 0; // current cap max
	int cap = 0;
	for (y = 0; y < N_r - 1; y++) {
		if (xsol[n_v][y] < P)cap += (3 * W_ch[xsol[n_v][y]] + 1); // total current capacity
		else {
			for (z = y - 1; z >= 0; z--) {
				if (xsol[n_v][z] < P && St_Sc[xsol[n_v][z]] == xsol[n_v][y]) cap -= (3 * W_ch[xsol[n_v][z]] + 1);
			}
			if (xsol[n_v][y] == T) {
				for (x = 0; x < P; x++) {
					if (ysol[x][1] == n_v) { //vehicle n_v is the second transfer 
						cap += (3 * W_ch[x] + 1);
					}
				}
			}
		}
		//cout << " stop: " << xsol[n_v][y] << " cap " << cap << endl;
		if (cap > curr_cap) curr_cap = cap; // if the capacity is bigger now
	}
	//cout << "max current capacity: " << curr_cap << "\n";
	
	// add transfer students
	if (i_t != -1) { //there are possible transfer students only if T is present
		for (l = 0; l < i_t; l++) {
			if (xsol[n_v][l] < P) curr_capT += (3 * W_ch[xsol[n_v][l]] + 1); // capacity before the transfer
		}
		//cout << " capacity: " << curr_capT << " until transfer\n";
	}

	
	if (curr_cap + (3 * W_ch[p] + 1) <= C_b) {
		// if vehicle n_v already visits school of p --> insert and try to see if its feasible
		if (l != -1) {
			//cout << " S_" << St_Sc[p] - P << " already exisits in v_" << n_v << endl;
			// if feasible --> accept
			if (insert_location(xsol, p,  n_v, Dsol)) {
				//cout << " ------------------------------------------> new vehicle accepted (exisiting school)\n";
				ysol[p][0] = n_v; // assign vehicle to student
				ysol[p][1] = -1;
				return true;
			}
			else {
				// try to improve route first
				//cout << " try improving route\n";
				improve_route( N_r, xsol[n_v], N_imp, I_S, n_v, initial); // try to optimize route to see if add student is possible
				// calculate max capacity again
				curr_cap = 0; // current cap max
				cap = 0;
				for (y = 0; y < N_r - 1; y++) {
					if (xsol[n_v][y] < P)cap += (3 * W_ch[xsol[n_v][y]] + 1); // total current capacity
					else {
						for (z = y - 1; z >= 0; z--) {
							if (xsol[n_v][z] < P && St_Sc[xsol[n_v][z]] == xsol[n_v][y]) cap -= (3 * W_ch[xsol[n_v][z]] + 1);
						}
						if (xsol[n_v][y] == T) {
							for (x = 0; x < P; x++) {
								if (ysol[x][1] == n_v) { //vehicle n_v is the second transfer 
									cap += (3 * W_ch[x] + 1);
								}
							}
						}
					}
					//cout << " stop: " << xsol[n_v][y] << " cap " << cap << endl;
					if (cap > curr_cap) curr_cap = cap; // if the capacity is beigger now
				}
				if (curr_cap + (3 * W_ch[p] + 1) <= C_b) {
					if (insert_location(xsol, p, n_v, Dsol)) {
						//cout << " ------------------------------------> new vehicle accepted (exisiting school after improvement) \n";
						ysol[p][0] = n_v; // assign vehicle to student
						ysol[p][1] = -1;
						return true;
					}
					else {
						// if not feasible --> go back to vehicle b
						//cout << "  new vehicle rejected\n";
						ysol[p][0] = b; // assign old vehicle to student
						ysol[p][1] = b2;
						return false;
					}
				}
				else {
					ysol[p][0] = b; // assign old vehicle to student
					ysol[p][1] = b2;
					return false;
				}
			}

		}
		// if vehice n_v not part of the route --> see if transfer is possible 
		else {
			//cout << " S_" << St_Sc[p] -P << " does NOT exisit in v_" << n_v << endl;
			// if feasible --> accept 
			if (insert_location(xsol, p, n_v, Dsol)) {
				//cout << " -------------------------------------------->  new vehicle accepted (school added too)\n";
				ysol[p][0] = n_v; // assign vehicle to student
				ysol[p][1] = -1;
				return true;
			}
			// if not feasible  --> try to improve route first
			else {
				/*
				for (int y = 0; y < xsol[n_v].size(); y++) {
					if(xsol[n_v][y]<P) cout << xsol[n_v][y] << "(S_"<< St_Sc[xsol[n_v][y]] - P << ")\t";
					else cout << "S_" << xsol[n_v][y] - P << "\t";
				}
				cout << endl << " size " << xsol[n_v].size() << " N_r " << N_r << endl;
				//*/

				//cout << " try to improve route now\n";
				improve_route( N_r, xsol[n_v], N_imp, I_S,  n_v, initial); // try to optimize route to see if add student is possible
				// calculate max capacity again
				curr_cap = 0; // current cap max
				cap = 0;
				for (y = 0; y < N_r - 1; y++) {
					if (xsol[n_v][y] < P)cap += (3 * W_ch[xsol[n_v][y]] + 1); // total current capacity
					else {
						for (z = y - 1; z >= 0; z--) {
							if (xsol[n_v][z] < P && St_Sc[xsol[n_v][z]] == xsol[n_v][y]) cap -= (3 * W_ch[xsol[n_v][z]] + 1);
						}
						if (xsol[n_v][y] == T) {
							for (x = 0; x < P; x++) {
								if (ysol[x][1] == n_v) { //vehicle n_v is the second transfer 
									cap += (3 * W_ch[x] + 1);
								}
							}
						}
					}
					//cout << " stop: " << xsol[n_v][y] << " cap " << cap << endl;
					if (cap > curr_cap) curr_cap = cap; // if the capacity is beigger now
				}
				if (curr_cap + (3 * W_ch[p] + 1) <= C_b) {
					if (insert_location(xsol, p, n_v, Dsol)) {
						//cout << " ---------------------------------------------> new vehicle accepted (school added too, after improvement)\n";
						ysol[p][0] = n_v; // assign vehicle to student
						ysol[p][1] = -1;
						return true;
					}
					else if (T != -1 && i_t != -1) {
						// if not feasible  --> try TRANSFER: Check if you can add the passenger without adding the school
						//cout << " try TRANSFER\n";
						b_trans = transfer_possible(Dsol, xsol, n_v, p,  b, b2);
						if (b_trans != -1) {
							//cout << " --------------------------------------------> new vehicle accepted with TRANSFER v_" << b_trans << endl;
							ysol[p][0] = n_v; // assign vehicle to student
							ysol[p][1] = b_trans;
							return true;
						}
						else {// if not feasible --> go back to vehicle b
							//cout << "  new vehicle rejected\n";
							ysol[p][0] = b; // assign old vehicle to student
							ysol[p][1] = b2;
							return false;
						}
					}
					else {
						//cout << "  new vehicle rejected\n";
						ysol[p][0] = b; // assign old vehicle to student
						ysol[p][1] = b2;
						return false;
					}
				}
				else {
					ysol[p][0] = b; // assign old vehicle to student
					ysol[p][1] = b2;
					return false;
				}
			}
		}
	}
	else if (curr_capT + (3 * W_ch[p] + 1) <= C_b && i_t != -1 && T != -1) {
		//cout << " try TRANSFER (due to lower capacity)\n";
		b_trans = transfer_possible(Dsol, xsol, n_v, p,  b, b2);
		if (b_trans != -1) {
			//cout << " --------------------------------------------> new vehicle accepted with TRANSFER v_" << b_trans << endl;
			ysol[p][0] = n_v; // assign vehicle to student
			ysol[p][1] = b_trans;
			return true;
		}
		else {// if not feasible --> go back to vehicle b
			//cout << "  new vehicle rejected\n";
			ysol[p][0] = b; // assign old vehicle to student
			ysol[p][1] = b2;
			return false;
		}
	}
	else {
		//cout << " TOO MUCH CAPACITY\n";
		ysol[p][0] = b; // assign old vehicle to student
		ysol[p][1] = b2;
		return false;
	}
	return true;
}

inline int printpluscost(const vector<vector<int>>& b_xsol, const vector<vector<int>>& b_Dsol, int ** ysol, bool print, bool write, string region, string tt_type) {
	int second_cost = 0;
	int b, t, s, p;
	int c_rejec = 0;
	int B = b_xsol.size(), N_r = 0;
	int tt = 0;
	int c_p = 0;
	int n_taxi = 0, m_c = 0;

	ofstream sol_p("C:/Users/bgalarza/Desktop/BG/School bus routing for special needs/Solutions/" + region + "/" + tt_type + "/formal_solution.txt");

	if (print) {
		int seconds, minutes, hours;
		string smin;
		cout << "++++++++++++++++++++++++++++++++++++++++ TIMETABLE\n";
		for (b = 0; b < B; b++) {
			cout << " VEHICLE " << b << endl;
			N_r = b_Dsol[b].size();
			for (s = 0; s < N_r; s++) {
				seconds = b_Dsol[b][s];
				minutes = seconds / 60;
				minutes = minutes % 60;
				smin = "";
				if (minutes >= 10) {
					smin = to_string(minutes);
				}
				else smin = "0" + to_string(minutes);
				hours = seconds / 3600;
				if (b_xsol[b][s] < P)cout << hours << ":" << smin << " (" << b_xsol[b][s] << ")\t";
				else cout << hours << ":" << smin << " (S_" << b_xsol[b][s] - P << ")\t";
			}
			cout << "\n";
		}
	}
	if (write) {
		int seconds, minutes, hours;
		string smin, ssec;
		sol_p << "---------------------------- Route and timetable:\n";
		for (b = 0; b < B; b++) {
			sol_p << "VEHICLE " << b + 1 << endl;
			N_r = b_Dsol[b].size();
			for (s = 0; s < N_r; s++) {
				minutes = b_Dsol[b][s] / 60;
				minutes = minutes % 60;
				seconds = b_Dsol[b][s] % 60;
				smin = "";
				ssec = "";
				if (minutes >= 10) {
					smin = to_string(minutes);
				}
				else smin = "0" + to_string(minutes);

				if (seconds >= 10) {
					ssec = to_string(seconds);
				}
				else ssec = "0" + to_string(seconds);

				hours = b_Dsol[b][s] / 3600;
				if (b_xsol[b][s] < P)sol_p << hours << ":" << smin << ":" << ssec << " (" << b_xsol[b][s] << ")\t";
				else sol_p << hours << ":" << smin << ":" << ssec << " (S_" << b_xsol[b][s] - P << ")\t";
			}
			sol_p << "\n";
		}
	}

	if (print) cout << "++++++++++++++++++++ ROUTE\n";
	if (write) sol_p << "------------------- Capacity utilisation:\n";
	for (b = 0; b < B; b++) {
		c_p = 0;
		m_c = 0;
		if (print) cout << " VEHICLE " << b << endl;
		if (write) sol_p << " VEHICLE " << b + 1 << endl;
		N_r = b_xsol[b].size();
		for (s = 0; s < N_r; s++) {
			if (b_xsol[b][s] < P) { // Location is student
				c_p += (3 * W_ch[b_xsol[b][s]] + 1);
			}
			else { // location is school now
				if (m_c < c_p) m_c = c_p;
				int c_a = 0;
				for (t = s - 1; t >= 0; t--) { // correct capacity of studnets assigned to this school 
					if (b_xsol[b][t] < P) if (St_Sc[b_xsol[b][t]] == b_xsol[b][s]) c_a += (3 * W_ch[b_xsol[b][t]] + 1); //students get off
				}
				c_p -= c_a;
				if (print) cout << "There are " << c_a << " students alighting at S_" << b_xsol[b][s] - P << ",   " << c_p << " students left\n";
				if (write) sol_p << "At school " << b_xsol[b][s] - P << ", there are " << c_a + c_p << " students onboard the vehicle, " << c_a << " students alight here and " << c_p << " students are still onboard\n";
				if (b_xsol[b][s] == T) {
					bool transfer = false;
					int c_p2 = 0;
					for (p = 0; p < P; p++) {
						if (ysol[p][1] == b) {
							c_p2 += (3 * W_ch[p] + 1);  // students that will get on
							transfer = true;
						}
					}
					if (c_p + c_p2 > m_c) m_c = c_p + c_p2;
					if (print && transfer) cout << "\n==> " << c_p2 << " students board at S_" << b_xsol[b][s] - P << " for transfer " << " ==> new total students: " << c_p + c_p2 << endl;
					if (write && transfer) sol_p << "\n==> " << c_p2 << " students board at school " << b_xsol[b][s] - P << " during a transfer " << " changing the total occupancy to " << c_p + c_p2 << endl;
				}
			}
		}
		if (m_c < 20)n_taxi++;
	}
	if (print) cout << endl;

	if (print)cout << "+++++++++++++++++++++++++++++ ASSIGNMENTS \n";
	if (write)sol_p << "----------------------------- Student assignments:\n";
	//Objective function value
	int* scosts = new int[B];
	for (b = 0; b < B; b++) {
		scosts[b] = -1;
	}
	int tt2 = 0;
	int start, end, alight, start2, end2, alight2, N_r2, b2, i;
	bool dess;
	for (p = 0; p < P; p++) {
		if (ysol[p][0] > -1) {
			b = ysol[p][0];
			b2 = ysol[p][1];
			N_r = b_xsol[b].size();
			//if not transfer
			if (b2 == -1) {
				start = -1, end = -1, alight = 0;
				for (i = 0; i < N_r; i++) {
					if (b_xsol[b][i] == p) start = b_Dsol[b][i];
					if (b_xsol[b][i] == St_Sc[p]) {
						end = b_Dsol[b][i];
						break;
					}
					if (b_xsol[b][i] < P)if (start != -1 && St_Sc[p] == St_Sc[b_xsol[b][i]])alight += tau_3 * (1 - W_ch[b_xsol[b][i]]) + tau_4 * W_ch[b_xsol[b][i]];
				}
				tt = end - start - alight;
				if (tt > scosts[b]) scosts[b] = tt;
			}
			else { // if transfer
				start = -1, end = -1, alight = 0;
				for (i = 0; i < N_r; i++) {
					if (b_xsol[b][i] == p) start = b_Dsol[b][i];
					if (b_xsol[b][i] == T) {
						end = b_Dsol[b][i];
						break;
					}
				}
				tt = end - start;
				if (tt > scosts[b]) scosts[b] = tt;
				N_r2 = b_xsol[b2].size();
				start2 = -1, end2 = -1, alight2 = 0;
				for (i = 0; i < N_r2; i++) {
					if (b_xsol[b2][i] == T) start2 = b_Dsol[b2][i];
					if (b_xsol[b2][i] == St_Sc[p]) {
						end2 = b_Dsol[b2][i];
						break;
					}
					if (b_xsol[b2][i] < P)if (start2 != -1 && St_Sc[p] == St_Sc[b_xsol[b2][i]])alight2 += tau_3 * (1 - W_ch[b_xsol[b2][i]]) + tau_4 * W_ch[b_xsol[b2][i]];
				}
				tt2 = (end2 - start2 - alight2);

				if (tt2 > scosts[b2]) scosts[b2] = tt2;
				tt = tt + tt2 + abs(end - start2);
			}
			if (print) {
				cout << " Student p_" << p << "\tTravelled " << round(tt / 0.6) / 100 << " minutes on\tv_" << b;
				if (ysol[p][1] != -1) cout << " and v_" << ysol[p][1];
				cout << "\t -->  to S_" << St_Sc[p] - P << endl;
			}
			if (write) {
				int minutes = start/ 60;
				minutes = minutes % 60;
				int seconds = start % 60;
				string smin = "";
				string ssec = "";
				if (minutes >= 10) {
					smin = to_string(minutes);
				}
				else smin = "0" + to_string(minutes);

				if (seconds >= 10) {
					ssec = to_string(seconds);
				}
				else ssec = "0" + to_string(seconds);

				int hours = start / 3600;
				sol_p << " Student " << p << " got picked up at " << hours << ":" << smin << ":" << ssec << ",\ttravelled " << round(tt / 0.6) / 100 << " minutes on vehicle " << b + 1;
				if (ysol[p][1] != -1) {
					minutes = end2 / 60;
					minutes = minutes % 60;
					seconds = start % 60;
					smin = "";
					ssec = "";
					if (minutes >= 10) {
						smin = to_string(minutes);
					}
					else smin = "0" + to_string(minutes);

					if (seconds >= 10) {
						ssec = to_string(seconds);
					}
					else ssec = "0" + to_string(seconds);

					hours = end2 / 3600;
					sol_p << " and vehicle " << ysol[p][1] + 1;
					sol_p << "\tand got dropped off at " << hours << ":" << smin << ":" << ssec << " at school " << St_Sc[p] - P << endl;
				}
				else {
					minutes = end / 60;
					minutes = minutes % 60;
					seconds = start % 60;
					smin = "";
					ssec = "";
					if (minutes >= 10) {
						smin = to_string(minutes);
					}
					else smin = "0" + to_string(minutes);

					if (seconds >= 10) {
						ssec = to_string(seconds);
					}
					else ssec = "0" + to_string(seconds);

					hours = end / 3600;
					sol_p << "\tand got dropped off at " << hours << ":" << smin << ":" << ssec << " at school " << St_Sc[p] - P << endl;
				}
			}
		}
	}

	second_cost = 0;
	for (b = 0; b < B; b++) {
		second_cost += scosts[b];
	}

	if (print) {
		cout << "-------------------- COST (nr. vehicles): " << B << " (" << n_taxi << " can be taxis)\n\t\tSum of max travel times : " << round(second_cost / 0.6) / 100 << " min " << endl;
	}
	if (write) {
		int sec_pb = second_cost / B;
		sol_p << "\n\n++++++++++++++++++++ COST (nr. vehicles): " << B << " vehicles, out of which " << n_taxi << " can be taxis and " << B - n_taxi << " have to be buses\n\t\t Sum of vehicle travel times : " << second_cost / 60 << " minutes and " << second_cost % 60 << " seconds\n";
		sol_p << "\t\t Average vehicle travel time: " << sec_pb / 60 << " minutes and " << sec_pb % 60 << " seconds\n";
	}
	delete scosts;
	sol_p.close();
	return second_cost;
}

int main() {
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

	int p, k, l, c, i, j, b, s;
	string region = "Diest"; //"Diest" or "Gavere"
	string tt_type = "avg"; // "min", "max" or "avg"
	
	bool BigRoute_OPT = false;
	if (region == "Gavere") {
		P = 100;
		S = 1;
		L = P + S;
	}
	else if (region == "Diest") {
		P = 360;
		S = 3;
		L = P + S;
		T = L - 1;
	}
	else {
		cout << " WRONG REGION NAME --> abort\n";
		exit(0);
	}

	
	//--------------------------------------------------- PARAMETERS -----------------------------------------------------------------------
	W_ch = new bool[P]; // Indicator to see if a student p in P needs a wheel chair, needs to be read
	St_Sc = new int[P]; // Indicator to see to which school a student p in P is assigned to, needs to be read
	t_ea = new int[S]; // Earliest arrival times for student p in P
	t_la= new int[S]; // Latest arrival times for student p in P


	// Early arrival times
	ifstream filea("C:/Users/bgalarza/Desktop/BG/School bus routing for special needs/Data/" + region + "/e_arrivals.txt");
	i = 0;
	while (i < S) {
		filea >> t_ea[i];
		i++;
	}
	filea.close();

	// Early arrival times
	ifstream fileal("C:/Users/bgalarza/Desktop/BG/School bus routing for special needs/Data/" + region + "/l_arrivals.txt");
	i = 0;
	while (i < S) {
		fileal >> t_la[i];
		i++;
	}
	fileal.close();

	// Wheelchair use
	ifstream filew("C:/Users/bgalarza/Desktop/BG/School bus routing for special needs/Data/" + region + "/wheelchair.txt");
	i = 0;
	while (i < P) {
		filew >> W_ch[i];
		i++;
	}
	filew.close();

	// Assigned schools
	ifstream filesa("C:/Users/bgalarza/Desktop/BG/School bus routing for special needs/Data/" + region + "/School_assignment.txt");
	i = 0;
	while (i < P) {
		filesa >> St_Sc[i];
		St_Sc[i] += P - 1; // to enumerate in entire list
		i++;
	}
	filesa.close();



	// Read in travel times
	TT = new int* [L];// travel times between locations
	for (i = 0; i < L; i++) {
		TT[i] = new int[L];
	}

	string path = "C:/Users/bgalarza/Desktop/BG/School bus routing for special needs/Data/" + region + "/travel_time_matrix_" + tt_type + ".csv";

	read_travel_matrix(path, TT);
	
	// -------------------------------------------------- VARIABLES ------------------------------------------------------------------------
	vector<vector<int>> xk; //routes of each vehicle b
	vector<vector<int>> xsol;

	vector<vector<int>> Dk; // timetable of each vehicle b
	vector<vector<int>> Dsol;

	yk = new int* [P]; // assignment of vehicles to stundents, each elements are the vehicles assigned to the student p in P, if only one vehicle assignmedn then the second is -1
	ysol =  new int* [P];
	for (p = 0; p < P; p++) {
		ysol[p] = new int[2];
		yk[p] = new int[2];
		yk[p][0] = -1;
		yk[p][1] = -1;
		ysol[p][0] = -1;
		ysol[p][1] = -1;
	}

	// ------------------------------------------------ PRE PROCESSING --------------------------------------------------------------------
	int* tempdist = new int[L];
	int* index = new int[L];
	closestS = new int* [L];
	for (i = 0; i < L; i++) {
		closestS[i] = new int[L];
	}
	
	for (i = 0; i < L; i++) {
		//initialize and copy distance of cities to temp array
		for (j = 0; j < L; j++) {
			if (j != i) tempdist[j] = TT[i][j];
			else tempdist[j] = 100000;
			index[j] = j;
		}

		//sort according to dist
		quickSort(index, tempdist, 0, L - 1);

		//keep track of best neighbors 
		for (l = 0; l < L; l++) {
			closestS[i][l] = index[l]; //index
		}
	}

	// get best route
	int N_it = 150000;
	best_route = new int[L];
	for (i = 0; i < L; i++) {
		best_route[i] = i;
	}
	
	if (BigRoute_OPT) {
		
		N_it = 30000000;
		cout << "BIG route to be improved (will stop the run)\n";
		vector<string> types = { "min", "avg", "max" };
		cout << region << endl;
		for (int z = 0; z < 3; z++) {
			cout << types[z] << endl;
			path = "C:/Users/bgalarza/Desktop/BG/School bus routing for special needs/Data/" + region + "/travel_time_matrix_" + types[z] + ".csv";
				
			read_travel_matrix(path, TT);

			for (i = 0; i < L; i++) {
				best_route[i] = i;
			}
			int best_sol = INT32_MAX;
			int current_sol;
			for (int x = 0; x < 10; x++) {
				cout << " itration " << x + 1 << endl;
				vector<int> best_route_curr;
				for (i = 0; i < L; i++) {
					best_route_curr.push_back(i);
				}
				current_sol = improve_route(L, best_route_curr, N_it, 1,  -1, 0);

				if (current_sol < best_sol) {
					cout << "--> better found\n";
					best_sol = current_sol;
					for (i = 0; i < L; i++) {
						best_route[i] = best_route_curr[i];
					}
				}
			}
			//write down the big route 
			ofstream BR_p("C:/Users/bgalarza/Desktop/BG/School bus routing for special needs/Data/" + region + "/BR_" + types[z] + ".txt");
			for (i = 0; i < L; i++) {
				BR_p << best_route[i] << endl;
			}
			BR_p.close();
		}
			
		exit(0);
	}
	else {
		cout << "BIG route improved --> read solution\n";
		// read optimized
		ifstream BR_r("C:/Users/bgalarza/Desktop/BG/School bus routing for special needs/Data/" + region + "/BR_" + tt_type + ".txt");
		i = 0;
		while (i < L) {
			BR_r >> best_route[i];
			i++;
		}
		BR_r.close();
	}
	cout << " end improvement\n";
	
	//remove memory
	delete index;
	delete tempdist;
	

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++ INITIAL SOLUTION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Use only buses and only students of the same school in the same bus, no transfers
	int N_b = 0, tt = 0;
	int curr_cap = 0;
	int curr_school = 0;
	int students_left = P;
	int vehicle = -1;
	int alight = 0;
	bool* ind_students = new bool[P]; //indicator array to see which stundents have been assigned (1) and which haven't (0)
	for (p = 0; p < P; p++) {
		ind_students[p] = 0;
	}

	int N_improv = 60000;
	bool I_S = false, initial = true;
	while (students_left > 0) {
		curr_cap = 0; // restart current capacity
		alight = 0; // alighting time
		vector<int> x1;
		vector<int> D1;
		xsol.push_back(x1); // initialize by adding first vehicle
		Dsol.push_back(D1);
		vehicle++; // index of current vehicle
		//decide which school is next
		for (int l = 0; l < L; l++) {
			p = best_route[l];
			if (p < P) {
				if (ind_students[p] == 0) {
					curr_school = St_Sc[p];
					break;
				}
			}
		}
		// cout << "++++++++++++++++++++++++++++++++++++++++++++ Vehicle " << vehicle <<" current school S_" <<curr_school-P << endl;
		for (l = 0; l < L; l++) {
			p = best_route[l];
			if (p < P) {
				if (ind_students[p] == 0 && curr_school == St_Sc[p]) { // not assigned stundents of the same school
					if (curr_cap + (3 * W_ch[p] + 1) <= C_b) { // add as many students of the same school to a vehicle (consider capacity + max travel time)
						//cout << " capacity left\n";
						if (insert_location(xsol, p, vehicle, Dsol)) { // add passenger location to the route in the best possible place
							ysol[p][0] = vehicle; // assign first vehicle to student
							curr_cap = curr_cap + (3 * W_ch[p] + 1); // add capacity, considering the wheelchairs
							ind_students[p] = 1; // update assigned students
							students_left--;
							alight += tau_3 * (1 - W_ch[p]) + tau_4 * W_ch[p];
							//cout << "\t-->Accepted\n";
						}
						else {
							//cout << "-------------------------> try to improve the route to add student " << p << endl;					
							N_b = xsol[vehicle].size();
							improve_route(N_b, xsol[vehicle], N_improv, I_S, vehicle, initial); // try to optimize route to see if add student is possible
							if (insert_location(xsol, p, vehicle, Dsol)) { // add passenger location to the route in the best possible place
								ysol[p][0] = vehicle; // assign first vehicle to student
								curr_cap = curr_cap + (3 * W_ch[p] + 1); // add capacity, considering the wheelchairs
								ind_students[p] = 1; // update assigned students
								students_left--;
								alight += tau_3 * (1 - W_ch[p]) + tau_4 * W_ch[p];
								//cout << "\t-->Accepted (after impr route)\n";
							}
							else {// when bus is full, optimize (route and) timetable, and school to the route
								N_b = xsol[vehicle].size();
								improve_route(N_b, xsol[vehicle], N_improv, I_S, vehicle, initial); // try to optimize route to see if add student is possible
								break;
							}
						}
					}
					else { // when bus is full, optimize (route and) timetable, and school to the route
						N_b = xsol[vehicle].size();
						improve_route(N_b, xsol[vehicle], N_improv, I_S, vehicle, initial); // try to optimize route to see if add student is possible
						break;
					}
				}
			}
		}
		if (l == L && curr_cap != 0) {
			N_b = xsol[vehicle].size();
			improve_route(N_b, xsol[vehicle], N_improv, I_S, vehicle, initial); // try to optimize route to see if add student is possible
		}
		if (curr_cap == 0) {
			vehicle--;
			Dsol.pop_back();
			xsol.pop_back();
		}
	}

	// Display solution and cheack feasibility 
	bool print = true, write = false;
	//printpluscost(xsol, Dsol, ysol, print, write, region, tt_type);
	int scost_current = isFEAS(xsol, Dsol);
	if (scost_current == -1) {
		cout << "  +++++++++++++++++++++++++++++++++++ INITAL SOLUTION INFEASIBLE!!! :(((\n";
		printpluscost(xsol, Dsol, ysol, print, write, region, tt_type);
		exit(0);
	}

	cout << "++++++ Initial solution with " << xsol.size() << " vehicles and  avg max travel time per bus: " << scost_current / xsol.size() / 60.0 << " min" << endl;
	
	// +++++++++++++++++++++++++++++++++++++++++++++++++++ HEURISTIC ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// NOTE:	adding a transfer should never add another vehicle, modify existing routes to make an assignment possible
	//			transfers only make sense when a vehicle is too full to transport everyone directly?? --> why not prev bus bring the passengers all the way???
	// NOTE:	there are students from the same location --> add them to the same bus if possible? especially if the same school too?
	// NOTE:	for the routing, add students to the route according to how far they are from their school (or BR?), then if the max travel time is exceeeded, 
	//			improve the route (2opt?) to see if still possible

	// ++++++++++++++++++ LNS: in each iteration, undo a number of student assignments; 
	//						* destroying is simply removing a student, 
	//						* repair is reassigning and making sure everything is feasible 
	//						  --> if current bus does visit St_Sc[p] then either add bus to trip OR do a transfer (if transfer place already visited) AND if not possible to reassign add bus
	int b2 = 0;
	piecewise_linear_distribution<double> des_st = triangular_distribution(0, P / 10, P - 1);
	normal_distribution<double> new_assign(0.0, 15.0); // used to 0-10
	vector<int> old_xsol0, old_xsol1, old_xsol2, old_Dsol0, old_Dsol1, old_Dsol2;
	N_it = 15000;

	int n_des = int(P * 0.06); // used to be 30 for Diest
	int* destroyed = new int[n_des];
	int n_v = 0;
	int itr = 0;
	bool create = false, destr = false;

	int N_r = 0, N_r2 = 0, b_trans = -1, fsch = -1, travel_time = 0;
	int B_current = xsol.size(), B_new = 0, scost_new = 0;

	//current solution
	xk.clear();
	Dk.clear();
	vector<int> x1, D1;

	for (b = 0; b < B_current; b++) {
		x1.clear();
		D1.clear();
		N_r = xsol[b].size();
		for (i = 0; i < N_r; i++) {
			x1.push_back(xsol[b][i]);
			D1.push_back(Dsol[b][i]);
		}
		xk.push_back(x1);
		Dk.push_back(D1);
	}
	for (p = 0; p < P; p++) {
		yk[p][0] = ysol[p][0];
		yk[p][1] = ysol[p][1];
	}
	
	int* student_rank = new int[P]; // rank to determined which student is destroyed
	int* student_score = new int[P]; // score based on how many students are onboard their bus 
	for (p = 0; p < P; p++) {
		student_score[p] = xsol[ysol[p][0]].size();
		student_rank[p] = p;
	}
	int ttemp = 0;
	print = false;

	double T_max = 1000, Temp = T_max;
	double lam = 500;
	int sig = 0;
	int dE = 0;

	while (itr < N_it) { // enter loop
		//*
		xsol.clear();
		Dsol.clear();
		for (b = 0; b < B_current; b++) {
			x1.clear();
			D1.clear();
			N_r = xk[b].size();
			for (i = 0; i < N_r; i++) {
				x1.push_back(xk[b][i]);
				D1.push_back(Dk[b][i]);
			}
			xsol.push_back(x1);
			Dsol.push_back(D1);
		}
		for (p = 0; p < P; p++) {
			ysol[p][0] = yk[p][0];
			ysol[p][1] = yk[p][1];
		}
		//*/
		if ((itr + 1) % 1000 == 0)cout << "  ******************************************* Interation " << itr + 1 << "\tvehicles: " << B_current << "\tmax travel time: " << scost_current / B_current / 60.0 << " min/vehicle" << endl;

		for (i = 0; i < n_des; i++) {
			destroyed[i] = -1;
		}

		// ------------------ DESTROY:  determine the destroyed students
		//cout << " destruction: \n";
		for (p = 0; p < P; p++) {
			student_score[p] = xsol[ysol[p][0]].size();
			student_rank[p] = p;
		}
		quickSort(student_rank, student_score, 0, P - 1); // sort from smallest to biggest
		for (i = 0; i < n_des; i++) {
			ttemp = int(des_st(generator));
			j = student_rank[ttemp];
			while (find(destroyed, destroyed + n_des, j) != destroyed + n_des) {
				ttemp = int(des_st(generator));
				j = student_rank[ttemp];
			}
			destroyed[i] = j;
			//cout << " j: " << j << "\n";
		}
		//cout << endl;

		// +++++++++++++++++ REPAIR: determine assignments, route and timetable 
		for (i = 0; i < n_des; i++) {
			p = destroyed[i];
			b = ysol[p][0];
			b2 = ysol[p][1];
			s = St_Sc[p];
			n_v = b;
			auto it = find(best_route, best_route + L, p);
			j = it - best_route;
			while (n_v == b || n_v == b2) {
				//determine new assignments
				l = int(new_assign(generator));
				if (j + l >= L || j + l < 0) {
					l *= -1;
				}
				k = best_route[j + l];
				while (k == p || k >= P) {
					l = int(new_assign(generator));
					if (j + l >= L || j + l < 0) {
						l *= -1;
					}
					k = best_route[j + l];
				}
				//new vehicle assiginment

				n_v = ysol[k][0];
			}
			/*
			if (b2 != -1) {
				if(n_v!=b2) cout << "---------- destroyed p_" << p << " with assigned S_" << St_Sc[p] - P << "  now assiged to v_" << n_v << " used to be v_" << ysol[p][0] << " and v_" << b2 << endl;
				else cout << "---------- destroyed p_" << p << " with assigned S_" << St_Sc[p] - P << "  now assiged to v_" << n_v << " used to be v_" << ysol[p][0] << " and v_" << b2 << " --> SAME\n";
			}
			else {
				if (n_v == b) cout << "---------- destroyed p_" << p << " with assigned S_" << St_Sc[p] - P << "  now assiged to v_" << n_v << " used to be v_" << ysol[p][0] << " --> SAME\n";
				else cout << "----------- destroyed p_" << p << " with assigned S_" << St_Sc[p] - P << " now assiged to v_" << n_v << " used to be v_" << ysol[p][0] << endl;
			}
			//*/

			//continue;
			//determine new routes and timetables of both destroyed and to be rebuild 
			
			//cout << N_r << endl;
			// remember the OGs
			N_r = xsol[n_v].size();
			old_xsol0.clear();
			old_Dsol0.clear();
			for (j = 0; j < N_r; j++) {
				old_xsol0.push_back(xsol[n_v][j]);
				old_Dsol0.push_back(Dsol[n_v][j]);
			}

			N_r = xsol[b].size();
			old_xsol1.clear();
			old_Dsol1.clear();
			for (j = 0; j < N_r; j++) {
				old_xsol1.push_back(xsol[b][j]);
				old_Dsol1.push_back(Dsol[b][j]);
			}
			if (b2 != -1) {
				N_r2 = xsol[b2].size();
				old_xsol2.clear();
				old_Dsol2.clear();
				for (j = 0; j < N_r2; j++) {
					old_xsol2.push_back(xsol[b2][j]);
					old_Dsol2.push_back(Dsol[b2][j]);
				}
			}

			// ---------------------------- First comes destruction
			destroy_solution(b, b2, p, xsol, Dsol);
			// ++++++++++++++++++++++++++++++ Now comes creation 
			create = repair_solution(b, b2, p, n_v, s, xsol, Dsol);
			if (!create) {
				N_b = old_xsol0.size();
				xsol[n_v].clear();
				Dsol[n_v].clear();
				for (j = 0; j < N_b; j++) {
					xsol[n_v].push_back(old_xsol0[j]);
					Dsol[n_v].push_back(old_Dsol0[j]);
				}
				N_b = old_xsol1.size();
				xsol[b].clear();
				Dsol[b].clear();
				for (j = 0; j < N_b; j++) {
					xsol[b].push_back(old_xsol1[j]);
					Dsol[b].push_back(old_Dsol1[j]);
				}
				if (b2 != -1) {
					N_b = old_xsol2.size();
					xsol[b2].clear();
					Dsol[b2].clear();
					for (j = 0; j < N_b; j++) {
						xsol[b2].push_back(old_xsol2[j]);
						Dsol[b2].push_back(old_Dsol2[j]);
					}
				}
			}
			//else student_score[p] = xsol[n_v].size(); // update student score
			/*
			else {
				//if (itr + 1 == 14) printpluscost(xsol, Dsol, ysol, TT, tau_1, tau_2, tau_3, tau_4, P, S, D_t, C_b, t_ea, t_la, W_ch, St_Sc, true, FEAS, T);
				//cout << " ---------------------------------------------> ACCEPTED\n";
			}
			
			cout << "AFTER:  vehicle " << n_v << endl;
			for (int y = 0; y < xsol[n_v].size(); y++) {
				cout << xsol[n_v][y] << "\t";
			}
				cout << endl;
			*/
			//cout << "--- end repair-destroy\n";
		}
		// eliminate "empty buses" 
		B_new  = xsol.size();
		for (b = 0; b < B_new; b++) {
			if (xsol[b].size() == 0) {
				// adjust ysol first  --> everything larger than b needs to go down by one
				for (p = 0; p < P; p++) {
					if (ysol[p][0] > b) ysol[p][0]--;
					if (ysol[p][1] > b) ysol[p][1]--;
				}
				//cout << " yay one vehicle (v_" << b << ") less   --> on iter: " << itr + 1 << endl;
				xsol.erase(xsol.begin() + b);
				Dsol.erase(Dsol.begin() + b);
				b--;
				B_new--;
				
			}
		}
		//cout << " check FEASIBILITY\n";
		//see if feasible and determine new second cost (accept if better cost of equal cost and better second cost)
		//if( itr + 1 == 1562) printpluscost(xsol, Dsol, ysol, TT, tau_1, tau_2, tau_3, tau_4, P, S, D_t, C_b, t_ea, t_la, W_ch, St_Sc, true, FEAS, T);
		scost_new = isFEAS(xsol, Dsol);
		if (scost_new != -1) {
			if (B_new < B_current) {
				cout << " solution with less vehicles: " << B_new << ", with max tt : " << scost_new/60.0 << " min --> on iter : " << itr + 1 << endl;
				// accept new solution
				B_current = B_new;

				xk.clear();
				Dk.clear();
				for (b = 0; b < B_current; b++) {
					x1.clear();
					D1.clear();
					N_r = xsol[b].size();
					for (i = 0; i < N_r; i++) {
						x1.push_back(xsol[b][i]);
						D1.push_back(Dsol[b][i]);
					}
					xk.push_back(x1);
					Dk.push_back(D1);
				}
				for (p = 0; p < P; p++) {
					yk[p][0] = ysol[p][0];
					yk[p][1] = ysol[p][1];
				}
			}
			else if (B_new == B_current) {
				dE = scost_new - scost_current;
				if (dE < 0) {
					//accept new solution
					//cout << " solution with smaller max tt: " << scost_new/60.0 << " min   -- > on iter : " << itr + 1 << endl;
					scost_current = scost_new;
					xk.clear();
					Dk.clear();
					for (b = 0; b < B_current; b++) {
						x1.clear();
						D1.clear();
						N_r = xsol[b].size();
						for (i = 0; i < N_r; i++) {
							x1.push_back(xsol[b][i]);
							D1.push_back(Dsol[b][i]);
						}
						xk.push_back(x1);
						Dk.push_back(D1);
					}
					for (p = 0; p < P; p++) {
						yk[p][0] = ysol[p][0];
						yk[p][1] = ysol[p][1];
					}
				}
				else if (exp(-dE / Temp) > one_to_zero(generator)) {
					scost_current = scost_new;
					xk.clear();
					Dk.clear();
					for (b = 0; b < B_current; b++) {
						x1.clear();
						D1.clear();
						N_r = xsol[b].size();
						for (i = 0; i < N_r; i++) {
							x1.push_back(xsol[b][i]);
							D1.push_back(Dsol[b][i]);
						}
						xk.push_back(x1);
						Dk.push_back(D1);
					}
					for (p = 0; p < P; p++) {
						yk[p][0] = ysol[p][0];
						yk[p][1] = ysol[p][1];
					}
				}
				//temperature
				if (dE > 0)sig++;
				else if (dE < 0)sig = 0;
				Temp = T_max + lam * log(1 + sig);
				//else cout << " solution with max tt: " << scost_new / 60.0 << " min   -- > on iter : " << itr + 1 << endl;
			}
		}
		/*
		else {
			//cout << "WOT.DA.HELL. infeasible ...\n";
			printpluscost(xk, Dk, yk, print, write, region, tt_type);
			cout << "INFEAS on itration " << itr + 1 << endl;
			//exit(0);
		}
		//exit(0);
		//cout << "next iteration\n";
		//*/
		itr++;
	}
	print = true, write = true;
	printpluscost(xk, Dk, yk, print, write, region, tt_type);


	// write to file
	ofstream yk_p("C:/Users/bgalarza/Desktop/BG/School bus routing for special needs/Solutions/" + region + "/" + tt_type + "/ysol.txt");
	ofstream xk_p("C:/Users/bgalarza/Desktop/BG/School bus routing for special needs/Solutions/" + region + "/" + tt_type + "/xsol.txt");
	ofstream Dk_p("C:/Users/bgalarza/Desktop/BG/School bus routing for special needs/Solutions/" + region + "/" + tt_type + "/Dsol.txt");
	for (i = 0; i < P; i++) {
		yk_p << yk[i][0] << "\t" << yk[i][0] << endl;
	}
	yk_p.close();
	int B = xk.size();
	for (b = 0; b < B; b++) {
		xk_p << " Vehicle " << b + 1 << endl;
		Dk_p << " Vehicle " << b + 1 << endl;
		N_b = xk[b].size();
		for (i = 0; i < N_b; i++) {
			xk_p << xk[b][i] << "\t";
			Dk_p << Dk[b][i] << "\t";
		}
		xk_p << endl;
		Dk_p << endl;
	}
	Dk_p.close();
	yk_p.close();

	// DELETE data
	for (i = 0; i < L; i++) {
		delete TT[i];
		delete closestS[i];
	}
	delete TT;
	delete closestS;

	for (i = 0; i < P; i++) {
		delete ysol[i];
		delete yk[i];
	}
	delete ysol;
	delete yk;


	delete ind_students;
	delete t_ea;
	delete t_la;
	delete W_ch;
	delete St_Sc;
	delete destroyed;
	delete student_rank;
	delete student_score;
	delete best_route;
	

	
	//_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
	//_CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
	//_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
	//_CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDOUT);
	//_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
	//_CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);
	//_CrtDumpMemoryLeaks();
}