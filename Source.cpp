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
#include "tcs_ssn_h/IndexSummurizing.h"
#include "tcs_ssn_h/IndexTraversal.hpp"

// define arrays


// this is a good code
// the main function
void main() {

	read_vertices(sn_vrtx, 10, sn_V_file);
	sn_read_edges(snGraph, sn_E_file);
	rn_read_edges(rnGraph, rnEdges, No_rn_E, 10, rn_E_file);

	
	
/***********  Build Index Functions  *************/
	//find_road_network_pivots();
	//find_social_network_pivots();

	//get_socialNetwork_connected();
	//get_roadNetwork_connected();
	
	int* cand_piv = new int[No_index_piv];

	std::cerr << "\n pivots------ ";
	for (int i = 0; i < No_index_piv; ++i) {
	labelA:
		int git = uniform(0, No_rn_V - 1);
		if (!isInTheArray(cand_piv, No_index_piv, git)) {
			cand_piv[i] = git;
			//std::cerr << git << " ";
		}

		else
			goto labelA;
	}

	gen_subgraphs_update(cand_piv);
	
	
	indexing();
	setParentOfNodes();


/**********  summurizin the index functions  *************/
	summraize_truss();
	summraize_RN_lb_dist();
	summraize_SN_lb_dist();
	summurize_keywords();
	summ_ub_in_influence_score();
	summ_ub_out_influence_score();
	computeIn_inf_Out_inf_for_vertices();


/***********  INDEX TRAVERSAL  *************************/
	indexTrav(4);




}



