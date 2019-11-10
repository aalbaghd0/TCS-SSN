#include <iostream>
#include "vertex.h"
#include<fstream>
#include"readFiles.hpp"
#include <list>
#include"heap.h"
#include "dij.hpp"


#define min					(a, b) (((a) < (b))? (a) : (b)  )
#define max					(a, b) (((a) > (b))? (a) : (b)  )

// file sizes
#define	No_sn_V				40000
#define No_sn_E				604304 * 2
#define No_rn_V				6105
#define No_rn_E				7035 * 2
// file names
char sn_V_file[20] =		"../SN_vertices.txt";
char sn_E_file[20] =		"../SN_edges.txt";
char rn_V_file[20] =		"../RN_vertices.txt";
char rn_E_file[20] =		"../RN_edges.txt";

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
	double cost = rn_Dij(rnGraph, rnEdges, No_rn_V, 0, 9);
};



