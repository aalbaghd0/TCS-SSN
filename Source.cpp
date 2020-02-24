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
	/*
	Heap* hp = new Heap();
	hp->init(1);
	{
		HeapEntry* he = new HeapEntry();
		he->key = 0;
		he->son1 = 0;
		hp->insert(he);
		delete he;
	}
	{
		HeapEntry* he = new HeapEntry();
		he->key = 5;
		he->son1 = 5;
		hp->insert(he);
		delete he;
	}
	{
		HeapEntry* he = new HeapEntry();
		he->key = 7;
		he->son1 = 7;
		hp->insert(he);
		delete he;
	}
	HeapEntry* he = new HeapEntry();
	he->key = 3;
	he->son1 = 3;
	hp->insert(he);
	delete he;

	{
		HeapEntry* he = new HeapEntry();
		he->key = 1;
		he->son1 = 1;
		hp->insert(he);
		delete he;
	}
	std::cerr << hp->check();
	hp->deleteEntry(0);
	while (hp->used > 0){
		hp->remove(he);
		std::cerr << he->son1<<" ";
	}
	*/
	



	//char user_file[100] = "data/200K_Gow_Cali/RT_groups_social_network.txt";
	//char road_file[100] = "data/200K_Gow_Cali/RT_groups_road_network.txt";
	//char output_file[100] = "data/200K_Gow_Cali/New_mapping_By_groups.txt";
	//char rn_Edges_file[100] = "data/200K_Gow_Cali/RN_cali_edges.txt";
	// this is what you will use BiDataset
	//mapping_social_to_road(user_file, road_file, output_file, rn_Edges_file);

	//char the_sn_edges_file[100] = "data/real_Gowalla+San Joaquin/gowalla_edges.txt";
	//char read_the_sn_edges_file[100]  = "data/real_Gowalla+San Joaquin/Gow_edges_ready.txt";
	//gen_edge_probability(the_sn_edges_file, read_the_sn_edges_file);

	read_vertices(sn_vrtx, 10, sn_V_file);
	sn_read_edges(snGraph, sn_E_file);
	rn_read_edges(rnGraph, rnEdges, No_rn_E, 10, rn_E_file);

/***********  Build Index Functions  *************/
	//get_socialNetwork_connected();
	//get_roadNetwork_connected();
	
	//find_road_network_pivots();
	//find_social_network_pivots();

	
	

	//gen_subgraphs_update(cand_piv);
	indexing();
	//setParentOfNodes(); 100


/**********  summurizin the index functions  *************/
	summraize_truss();
	summraize_RN_lb_dist();
	summraize_SN_lb_dist();
	summurize_keywords();
	summ_ub_in_influence_score();
	summ_ub_out_influence_score();
	computeIn_inf_Out_inf_for_vertices();


/***********  INDEX TRAVERSAL  *************************/
	
	indexTrav();

}



