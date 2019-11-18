#pragma once
#ifndef INDEX_HPP
#define INDEX_HPP

#include <iostream>
#include <list>
#include "tcs_ssn_h/parameterSettings.h"
#include "tcs_ssn_h/vertex.h"
#include <set>
#include "tcs_ssn_h/heap.h"
#include <set>
#include <unordered_map>
#include "tcs_ssn_h/dij.hpp"
#include"tcs_ssn_h/vertex.h"
#include "tcs_ssn_h/heap.h"
#include"tcs_ssn_h/gendef.h"
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>

/*
	functions prototypes
*/
double sn_dist(int src, int dst);
double rn_dist(int src, int dst);
double rn_dist_for_users(int src, int dst);
double X_sc(std::set<int> G[], int number_subgraphs);
double X_st(std::set<int> G[], int number_subgraphs);
double X_inf(std::set<int> G[], int number_subgraphs);
double evaluate_subgraphs(std::set<int> G[], int number_subgraphs);
double evaluate_Indexsubgraphs(std::set<int> G[], int no_of_subgraphs);


double quality(int v, int piv) {
	return (rn_dist_for_users(v, piv) + sn_dist(v, piv));
}

void gen_subgraphs(int cand_index_piv[], std::set<int> G[]) {
	double qual_rslt;
	int best_quality;
	int assign;
	
	for (int v = 0; v < No_sn_V; v++) {
		best_quality = INT_MAX;
		qual_rslt = 0.0;
		for (int piv = 0; piv < No_index_piv; piv++) {
			
			qual_rslt = quality(v, cand_index_piv[piv]);

			if (qual_rslt < best_quality) {
				assign = piv;
				best_quality = qual_rslt;
			}
		}
		G[assign].insert(v);
	}
}

void sn_piv_select(std::list<int> rnGraph, std::list<int> snGraph) {
	double global_cost = INT_MAX;
	int* S_p = new int[No_index_piv];
	int* new_S_p = new int[No_index_piv];
	std::set<int>* G = new std::set<int>[No_index_piv];
	std::set<int>* new_G = new std::set<int>[No_index_piv];
	
	double cost_G = 0.0;
	double local_cost = 0.0;
	int new_cost = 0.0;

	int global_iter = 5;
	int swap_iter = 10;

	int get_piv = 0;
	int new_piv = 0;

	bool cost_was_updated = false;
	for (int a = 1; a < global_iter; a++) {
		
		// select random pivots at first
		for (int i = 0; i < No_index_piv; ++i) {
			S_p[i] = uniform(0, No_sn_V);
		}
		// get subgraphs based on pivots
		gen_subgraphs(S_p, G);

		// evaluate the cost function
		local_cost = evaluate_subgraphs(G, No_subgraphs);
		
		for (int  b = 1; b < swap_iter; b++) {
			get_piv = uniform(0, No_index_piv);
			new_piv = uniform(0, No_sn_V);

			memcpy(new_S_p, S_p, No_index_piv);
			new_S_p[get_piv] = new_piv;
			
			gen_subgraphs(new_S_p, new_G);
			new_cost = evaluate_subgraphs(new_G, No_subgraphs);

			if (new_cost > local_cost) {
				local_cost = new_cost;
				memcpy(S_p, new_S_p, No_index_piv);
			}
		}
		if (local_cost > global_cost) {
			memcpy(index_piv, S_p, No_index_piv);
			global_cost = local_cost;
		}
	}

	/*
	for (int i = 0; i < No_index_piv; ++i) {
		for (std::set<int>::iterator it = G[i].begin(); it != G[i].end(); ++it) {
			index[GlobalIndexIter].insert(*it);
			GlobalIndexIter++;
		}
	}
	*/
}


double evaluate_subgraphs(std::set<int> G[], int no_of_subgraphs) {
	double rslt1 = W1 *  X_sc(G, no_of_subgraphs);
	double rslt2 = W2 * ( 1 - X_st(G, no_of_subgraphs)  );
	double rslt3 = W1 * ( 1 - X_inf(G, no_of_subgraphs) );

	return rslt1 + rslt2 + rslt3;
}



double X_sc(std::set<int> G[], int no_of_subgraphs) {
	std::unordered_map<pair, bool, pair_hash> map;

	double sub_rslt = 0.0;
	double rslt = 0.0;

	for (int g = 0; g < no_of_subgraphs; g++) {
		for (std::set<int>::iterator it = G[g].begin(); it != G[g].end(); ++it) {
			for (std::set<int>::iterator it2 = G[g].begin(); it2 != G[g].end(); ++it2) {
				if (*it != *it2) {

					if (!map[std::make_pair(*it, *it2)] && !map[std::make_pair(*it2, *it)]) {
						sub_rslt = sub_rslt + rn_dist_for_users(*it, *it2);
						map[std::make_pair(*it, *it2)] = true;
						map[std::make_pair(*it2, *it)] = true;
					}
				}
			}
		}
		rslt = rslt + sub_rslt;
		sub_rslt = 0.0;
	}
	return rslt;
}

double X_st(std::set<int> G[], int no_of_subgraphs) {
	std::unordered_map<pair, bool, pair_hash> map;

	double sub_rslt = 0.0;
	double rslt = 0.0;

	for (int g = 0; g < no_of_subgraphs; g++) {
		for (std::set<int>::iterator it = G[g].begin(); it != G[g].end(); ++it) {
			for (std::set<int>::iterator it2 = G[g].begin(); it2 != G[g].end(); ++it2) {
				if (*it != *it2) {

					if (!map[std::make_pair(*it, *it2)] && !map[std::make_pair(*it2, *it)]) {

						sub_rslt = sub_rslt + ((sn_vrtx[*it].truss + sn_vrtx[*it2].truss) / sn_dist(*it, *it2));

						map[std::make_pair(*it, *it2)] = true;
						map[std::make_pair(*it2, *it)] = true;
					}
				}
			}
		}
		rslt = rslt + sub_rslt;
		sub_rslt = 0.0;
	}
	return rslt;
}

double X_inf(std::set<int> G[], int no_of_subgraphs) {
	std::unordered_map<pair, bool, pair_hash> map;

	double sub_rslt = 0.0;
	double rslt = 0.0;

	for (int g = 0; g < no_of_subgraphs; g++) {
		for (std::set<int>::iterator it = G[g].begin(); it != G[g].end(); ++it) {
			for (std::set<int>::iterator it2 = G[g].begin(); it2 != G[g].end(); ++it2) {
				if (*it != *it2) {

					if (!map[std::make_pair(*it, *it2)] && !map[std::make_pair(*it2, *it)]) {
						
						sub_rslt = sub_rslt + (inf_score(*it, *it2));

						map[std::make_pair(*it, *it2)] = true;
						map[std::make_pair(*it2, *it)] = true;
					}
				}
			}
		}
		rslt = rslt + sub_rslt;
		sub_rslt = 0.0;
	}
	return rslt;
}

/*
Here we compute distances functions
-- we start with the road network shortest path distance
*/


/*
GIVEN::     RN, src, and dst
ENSURE::    shortest path distance between src and dst 
*/
double rn_dist_for_users(int src, int dst) {
	double temp_rslt = 0.0;

	for (int i = 0; i < No_CKINs; ++i) {// for all user locations
		for ( int j = 0; j < No_CKINs; j++) {
			// map the user location to the from the social to the road network, match the 
			// user location to a vertex on the road network
			int src_map = sn_vrtx[src].ckins[i];
			int dst_map = sn_vrtx[dst].ckins[j];

			temp_rslt = temp_rslt + rn_Dij(src_map, dst_map);
		}
	}
	return temp_rslt / (No_CKINs * 2);
}
////////////////////////////////////////////////////////
/*
GIVEN::     SN, src, and dst
ENSURE::    shortest path distance between src and dst
*/
double sn_dist(int src, int dst) {
	return sn_Dij(src, dst);
}

////////////////////////////////////////////////////////
/*
GIVEN::     RN, src, and dst
ENSURE::    shortest path distance between src and dst
*/
double rn_dist(int src, int dst) {
	return rn_Dij(src, dst);
}

//////////////
/*
REQUIRE::		two vertices a and b
ENSURE ::		rslt--> is the itersection set of edges with the edge a--b
				we return the number of intersected vertices
*/


int intersect(std::set<int>* common_edges, int a, int b) {

	double rslt = 0.0;
	std::set<int> intersect;

	set_intersection(sn_vrtx[snEdges[a].from].nbrs.begin(), sn_vrtx[snEdges[a].from].nbrs.end(),
		sn_vrtx[snEdges[b].to].nbrs.begin(), sn_vrtx[snEdges[b].to].nbrs.end(),
		std::inserter(intersect, intersect.begin()));
	rslt = intersect.size();
	// v -- a, v -- b, a -- v, b -- v

	for (std::set<int>::iterator it = intersect.begin();
		it != intersect.end(); ++it) {

		common_edges->insert(hash_edge[std::make_pair(a, *it)]);
		common_edges->insert(hash_edge[std::make_pair(*it, a)]);
	}
	return rslt;
}
//////////////////////////////////////////////////////
/*
GIVEN	::	social network graph
ENSURE	::	the maximum support of each edge
*/

void truss_decomposition() {
	int k = 3;
	Heap* hp = new Heap();
	hp->init(2);
	std::set<int>* pool = new std::set<int>;
	int e;
	int min_sup = INT_MAX;
	int min_temp = INT_MAX;

	for (int e = 0; e < No_sn_E; e++) {
		snEdges[e].sup = intersect(pool, snEdges[e].from, snEdges[e].to);
		HeapEntry* he = new HeapEntry();
		he->son1 = e;

		if (snEdges[e].sup < min_sup)
			min_sup = snEdges[e].sup;

		he->key = snEdges[e].sup;

		hp->insert(he);
		delete he;
	}

label2:
	while (true) {
		HeapEntry* he = new HeapEntry();
		hp->remove(he);
		e = he->son1;
		delete he;
		if (snEdges[e].sup >= (k - 2))
			break;

		intersect(pool, snEdges[e].from, snEdges[e].to);
		for (std::set<int>::iterator it = pool->begin(); it != pool->end(); ++it) {
			snEdges[*it].sup = snEdges[*it].sup - 1;

			hp->deleteEntry(*it);
			snEdges[*it].sup = snEdges[*it].sup - 1;

			HeapEntry* he = new HeapEntry();
			he->son1 = *it;
			he->key = snEdges[*it].sup;
			hp->insert(he);
			delete he;
		}
	}

	if (hp->used > 0) {
		k = k + 1;
		goto label2;
	}
	// assign sup of truss values
	hp->~Heap();
	
	for (int i = 0; i < No_sn_V; ++i) {
		Heap* hp = new Heap();
		hp->init(2);

		for (std::set<int>::iterator it = sn_vrtx[i].myedges.begin();
			it != sn_vrtx[i].myedges.end(); ++it) {
			HeapEntry* he = new HeapEntry();
			he->son1 = *it;
			he->key = -snEdges[*it].sup;
			hp->insert(he);
			delete he;
		}
		HeapEntry* he = new HeapEntry();
		hp->remove(he);
		sn_vrtx[i].truss = he->son1;
		delete he;
		hp->~Heap();
	}
	
}
 
///////////////////////////////

/*
REQUIRES :: spgraphs
ENSURES  :: build the hybrid index
*/

//امرر كراف كامل, وامرر 
void get_index_subgraphs(std::set<int> G[], int piv_set[], int no_piv_set,int new_pivots[], int no_new_piv) {
	double qual_rslt;
	int best_quality;
	int assign;
	
	for (int v = 0; v < no_piv_set; v++) {
		best_quality = INT_MAX;
		qual_rslt = 0.0;
		for (int piv = 0; piv < no_new_piv; piv++) {

			qual_rslt = quality(piv_set[v], new_pivots[piv]);

			if (qual_rslt < best_quality) {
				assign = piv;
				best_quality = qual_rslt;
			}
		}
		G[assign].insert(v);
	}
}
//////
/*

	ENSURES :: find level_index pivots 
*/
int* Index_piv_select(int no_new_piv, int prev_piv[], int no_prev_piv) {
	double global_cost = INT_MAX;

	int* S_p = new int[no_new_piv];
	int* new_S_p = new int[no_new_piv];
	int* final_S_p = new int[no_new_piv];



	int* temp_arry = new int[no_new_piv];
	
	std::set<int>* G = new std::set<int>[no_new_piv];
	std::set<int>* new_G = new std::set<int>[no_new_piv];

	double cost_G = 0.0;
	double local_cost = 0.0;
	int new_cost = 0.0;

	int global_iter = 5;
	int swap_iter = 10;

	int get_piv = 0;
	int new_piv = 0;

	int git;
	for (int a = 1; a < global_iter; ++a) {
		// select random pivots at first
		for (int i = 0; i < no_new_piv; ++i) {
			labelA:
			git = uniform(0, no_prev_piv);
			int git_val = prev_piv[git];
			if (!isInTheArray(S_p, no_prev_piv, git_val))
				S_p[i] = git_val;
			else
				goto labelA;
		}
		// get subgraphs based on pivots
		get_index_subgraphs(G, prev_piv, no_prev_piv, S_p, no_new_piv);
		
		// evaluate the cost function
		local_cost = evaluate_Indexsubgraphs(G, no_new_piv);
		for (int b = 1; b < swap_iter; b++) {
			
			get_piv = uniform(0, no_new_piv);
			new_piv = uniform(0, no_prev_piv);

			memcpy(new_S_p, S_p, no_new_piv);
			new_S_p[get_piv] = new_piv;
			//get the neew subgraphs
			get_index_subgraphs(new_G, prev_piv, no_prev_piv, new_S_p, no_new_piv);
			//evaluate the new subgraphs
			new_cost = evaluate_subgraphs(new_G, no_new_piv);

			if (new_cost > local_cost) {
				local_cost = new_cost;
				memcpy(S_p, new_S_p, no_new_piv);
			}
		}
		if (local_cost > global_cost) {
			memcpy(final_S_p, S_p, no_new_piv);
			global_cost = local_cost;
		}
	}



	return final_S_p;
}

double evaluate_Indexsubgraphs(std::set<int> G[], int no_of_subgraphs) {
	
	double rslt1 = W1 * X_sc(G, no_of_subgraphs);

	return rslt1 ;
}


void generateTheIndex(int prev_piv[], int no_prev_pivots, int divBy) {
	
	int no_new_piv = no_prev_pivots / divBy;
	
	while (no_new_piv > 0) {
		int* temp_piv_arr = new int[no_new_piv];
		temp_piv_arr = Index_piv_select(no_new_piv, prev_piv, no_prev_pivots);
		
		no_prev_pivots = no_new_piv;
		no_new_piv = no_new_piv / divBy;
		
		memset(prev_piv, -1, no_prev_pivots);

		memcpy(prev_piv, temp_piv_arr, no_prev_pivots);
	}
}
#endif // !INDEX_HPP
