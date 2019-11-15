#include <iostream>
#include "tcs_ssn_h/vertex.h"
#include<fstream>
#include"tcs_ssn_h/readFiles.hpp"
#include <list>
#include"tcs_ssn_h/heap.h"
#include "tcs_ssn_h/dij.hpp"
#include"tcs_ssn_h/parameterSettings.h"
#include <algorithm>
#include "tcs_ssn_h/gendef.h"
#include <iterator>
#include <unordered_map>
#include "buildIndex.h"

// define arrays






// the main function
void main() {

	std::unordered_map<pair, bool, pair_hash> map;
	map[std::make_pair(1, 2)] = true;

	std::cerr << map[std::make_pair(2, 2)];

	read_vertices(sn_vrtx, 10, sn_V_file);
	sn_read_edges(snGraph, sn_E_file);
	rn_read_edges(rnGraph, rnEdges, No_rn_E, 10, rn_E_file);

	

	std::cerr<< "\n"<<  rn_Dij(0, 200 )<< "\n";
	std::cerr << "\n"<< sn_Dij(0, 2764)<< "\n";
};



