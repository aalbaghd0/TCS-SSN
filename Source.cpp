#include <iostream>
#include "tcs_ssn_h/vertex.h"
#include<fstream>
#include"tcs_ssn_h/readFiles.hpp"
#include <list>
#include"tcs_ssn_h/heap.h"
#include "tcs_ssn_h/dij.hpp"
#include"tcs_ssn_h/parameterSettings.h"



// define arrays
vertex *sn_vrtx = new vertex[No_sn_V];
vertex* rn_vrtx = new vertex[No_rn_V];

std::list<int> * snGraph = new std::list<int>[No_sn_V];
std::list<int>* rnGraph = new std::list<int>[No_rn_V];

edges *snEdges = new edges[No_sn_E];
edges* rnEdges = new edges[No_rn_E];


// the main function
void main(int args, char arg[]) {
	read_vertices(sn_vrtx, 10, sn_V_file);
	sn_read_edges(snGraph, sn_E_file);
	rn_read_edges(rnGraph, rnEdges, No_rn_E,10, rn_E_file);
	std::cerr<< "\n"<< rn_Dij(rnGraph, rnEdges, No_rn_V, 0, 200)<< "\n";
	std::cerr << "\n"<< sn_Dij(snGraph, No_sn_V, 0, 2764)<< "\n";
};



