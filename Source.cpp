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
	
	std::cerr<< number_nodes(6, 2);
	read_vertices(sn_vrtx, 10, sn_V_file);
	sn_read_edges(snGraph, sn_E_file);
	rn_read_edges(rnGraph, rnEdges, No_rn_E, 10, rn_E_file);

	std::cerr << inf_score(5, 1);
	
	truss_decomposition();
	sn_piv_select();
	std::unordered_map<int, std::set<int>> GGG;
	int* f_piv = new int[100];
	int no_new_piv = 2;
	int no_prev_piv = No_index_piv;
	int* layer_piv = new int[No_index_piv];
	std::memcpy(layer_piv, index_piv, sizeof(index_piv[0]) * No_index_piv);


	while (no_new_piv > 0) {
		GGG = Index_piv_select(no_new_piv, layer_piv, no_prev_piv, f_piv);
		
		no_prev_piv = no_new_piv;
		no_new_piv = no_new_piv / 2;
		std::memcpy(layer_piv, f_piv, sizeof(f_piv[0]) * no_prev_piv);


		for (int i = 0; i < no_prev_piv; ++i) {
			std::cerr << "Index_Piv " << f_piv[i] << " --->> ";
			for (std::set<int>::iterator it = GGG[f_piv[i]].begin(); it != GGG[f_piv[i]].end(); ++it) {
				std::cerr << *it << " ";
			}
			std::cerr << "\n" << "\n";
		}
	}
};



