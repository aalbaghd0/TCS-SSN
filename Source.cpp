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

// define arrays



// the main function
void main(int args, char arg[]) {
	read_vertices(sn_vrtx, 10, sn_V_file);
	sn_read_edges(snGraph, sn_E_file);
	rn_read_edges(rnGraph, rnEdges, No_rn_E,10, rn_E_file);

	Heap* hp = new Heap();
	HeapEntry* he = new HeapEntry();
	hp->init(2);
	he->son1 = 5;
	he->key = 5;
	hp->insert(he);
	he->son1 = 20;
	he->key = 20;
	hp->insert(he);
	
	he->son1 = 19;
	he->key = 19;
	hp->insert(he);

	he->son1 = 30;
	he->key = 30;
	hp->insert(he);

	he->son1 = 1;
	he->key = 1;
	hp->insert(he);

	he->son1 = 0;
	he->key = 0;
	hp->insert(he);
	
	he->son1 = 11;
	he->key = 11;
	hp->insert(he);

	he->son1 = 2;
	he->key = 2;
	hp->insert(he);
	hp->deleteEntry(2);
	while (hp->used > 0) {
		hp->remove(he);
		std::cerr << he->son1 << "\n";
	}
	//hp->deleteEntry(19);
	std::cerr<< hp->check();
	std::cerr<< "\n"<<  rn_Dij(0, 200 )<< "\n";
	std::cerr << "\n"<< sn_Dij(0, 2764)<< "\n";
};



