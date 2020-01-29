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

#include "tcs_ssn_h/dij.hpp"
#include"tcs_ssn_h/vertex.h"
#include "tcs_ssn_h/heap.h"
#include"tcs_ssn_h/gendef.h"
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>
#include <limits>
#include <fstream>


/*
	functions prototypes
*/
double sn_dist(int src, int dst);
double rn_dist(int src, int dst);
double rn_dist_for_users(int src, int dst);
double X_sc(std::unordered_map<int, std::unordered_set<int>> G, int pivots[], int number_subgraphs);
double X_st(std::unordered_map<int, std::unordered_set<int>> G, int pivots[], int number_subgraphs);
double X_inf(std::unordered_map<int, std::unordered_set<int>> G, int pivots[], int number_subgraphs);
double evaluate_subgraphs(std::unordered_map<int, std::unordered_set<int>> G, int pivots[], int number_subgraphs);
double evaluate_Indexsubgraphs(std::unordered_set<int> G[], int no_of_subgraphs);
double inf_score_users(int src, int dst);
double evaluate_SN_pivots(int S_p[], int no_pivs);
double evaluate_RN_pivots(int S_p[], int no_pivs);
std::unordered_map<int, std::unordered_set<int>> gen_subgraphs_update(int cand_piv[]);
void expand(int v, std::unordered_map<int, Heap>& hp, int assign[], int& condtition,
	Heap& pivHeap, std::unordered_map<pair, bool, pair_hash>& isSeen, int cand_piv[]);
double rn_dist_for_users_NoT_BFS(int src, int dst);
double rn_dist_rnPiv_to_user(int r_piv, int user);
float gaussian(float mean, float sigma);


double quality(int v, int piv) {
	double sn_dist_rslt = 0.0, rn_dist_rslt = 0.0;

/*
	if (!check_hash_sn_dist[std::make_pair(piv, v)]) { // if we don't have the distance, the we will compute it
		sn_Dij_to_all_vertices(piv); // find distance to all other vertices
									 // at the same time store distances to all other vertices
		sn_dist_rslt = hash_sn_dist[std::make_pair(piv, v)];
	}
	else {
		sn_dist_rslt = hash_sn_dist[std::make_pair(piv, v)];
	}

	*/
	return (rn_dist_for_users(piv, v));// + sn_dist_rslt);
}

std::unordered_map<int, std::unordered_set<int>> gen_subgraphs(int cand_index_piv[]) {
	double qual_rslt;
	int best_quality;
	int assign;

	std::unordered_map<int, std::unordered_set<int>> GG;

	for (int v = 0; v < No_sn_V; v++) {
		best_quality = INT_MAX;
		qual_rslt = 0.0;
		for (int piv = 0; piv < No_index_piv; piv++) {
			qual_rslt = quality(v, cand_index_piv[piv]);

			if ((qual_rslt < best_quality) && (GG[cand_index_piv[piv]].size() < 100)) {
				assign = cand_index_piv[piv];
				best_quality = qual_rslt;
			}
		}
		GG[assign].insert(v);
		//std::cerr << sizeof(GG[assign]);
	}

	/*
	for (int i = 0; i < No_index_piv; ++i) {
		std::cerr << "Index_Piv " << cand_index_piv[i] << " --->> ";
		for (std::unordered_set<int>::iterator it = GG[cand_index_piv[i]].begin(); it != GG[cand_index_piv[i]].end(); ++it) {
			std::cerr << *it << " ";
		}
		std::cerr << "\n" << "\n";
	}
	*/
	std::cerr << "\n" << "------- size of each node -------"<<"\n";
	for (int i = 0; i < No_index_piv; ++i) {
		std::cerr<< GG[cand_index_piv[i]].size() << " ";
	}
	return GG;
}

std::unordered_map<int, std::unordered_set<int>> sn_piv_select() {
	double global_cost = -INT_MAX;
	int* S_p = new int[No_index_piv];
	int* new_S_p = new int[No_index_piv];
	double final_cost;
	std::unordered_map<int, std::unordered_set<int>> GG;
	std::unordered_map<int, std::unordered_set<int>> new_GG;
	std::unordered_map<int, std::unordered_set<int>> G;


	double cost_G = 0.0;
	double local_cost = 0.0;
	double new_cost = 0.0;

	int global_iter = 2;
	int swap_iter = 2;

	int get_piv = 0;
	int new_piv = 0;
	bool readInitialPivots = true;

	FILE* fi = fopen("data/sn_init_pivots.txt", "r");

	bool cost_was_updated = false;
	for (int a = 1; a < global_iter; a++) {

		// for the first round, we read pivots from the file
		if (readInitialPivots) {
			if (fi == NULL) {
				std::cout << "The edge file cannot be open";
			}
			else {
				int i = 0;
				int getRide;
				fscanf(fi, "%d", &getRide);

				while (!feof(fi)) {
					fscanf(fi, "%d", &S_p[i]);
					++i;
				}
			}
		}
		else { // for the second round, we generate pivots rendomly
			// select random pivots at first
			for (int i = 0; i < No_index_piv; ++i) {
			labelA:
				int git = uniform(0, No_rn_V - 1);
				if (!isInTheArray(S_p, No_index_piv, git))
					S_p[i] = git;
				else
					goto labelA;
			}
		}
		readInitialPivots = false;
		// get subgraphs based on pivots
		//for (int i = 0; i < No_index_piv; ++i)
		//	std::cerr << S_p[i] << " ";

		GG = gen_subgraphs_update(S_p);

		// evaluate the cost function
		local_cost = evaluate_subgraphs(GG, S_p, No_subgraphs);
		//std::cerr << "THe Evaluation :: " << local_cost << "\n \n";


		for (int b = 1; b < swap_iter; b++) {
			get_piv = uniform(0, No_index_piv - 1);

		labelB:
			int git = uniform(0, No_rn_V - 1);
			if (!isInTheArray(S_p, No_index_piv, git))
				new_piv = git;
			else
				goto labelB;


			std::memcpy(new_S_p, S_p, sizeof(S_p[0]) * No_index_piv);

			new_S_p[get_piv] = new_piv;

			new_GG = gen_subgraphs_update(new_S_p);



			new_cost = evaluate_subgraphs(new_GG, new_S_p, No_subgraphs);
			//std::cerr << "THe Evaluation :: " << new_cost << "\n \n";
			if (new_cost > local_cost) {
				local_cost = new_cost;
				std::memcpy(S_p, new_S_p, sizeof(new_S_p[0]) * No_index_piv);

				// get the final subgraph
				//memcpy(final_G, new_G, sizeof(new_G) * No_index_piv);

				GG = new_GG;
			}
		}
		if (local_cost > global_cost) {
			memcpy(index_piv, S_p, sizeof(S_p[0]) * No_index_piv);
			global_cost = local_cost;
			G = GG;
		}
	}


	/*
	for (int i = 0; i < No_index_piv; ++i) {
		std::cerr << "Index_Piv " << index_piv[i] << " --->> ";
		for (std::unordered_set<int>::iterator it = G[index_piv[i]].begin(); it != G[index_piv[i]].end(); ++it) {
			std::cerr << *it << " ";
		}
		std::cerr << "\n" << "\n";
	}
	std::cerr << "THe Evaluation :: " << global_cost << "\n \n";
	*/
	return G;
}


double evaluate_subgraphs(std::unordered_map<int, std::unordered_set<int>> G, int pivots[], int no_of_subgraphs) {
	double rslt1 = 0;//W1 * X_sc(G, pivots, no_of_subgraphs);
	double rslt2 = 0;//W2 * (1 - X_st(G, pivots, no_of_subgraphs));
	double rslt3 = 0;//W1 * (1 - X_inf(G, pivots, no_of_subgraphs));

	return rslt1 + rslt2 + rslt3;
}



double X_sc(std::unordered_map<int, std::unordered_set<int>> G, int pivots[], int no_of_subgraphs) {
	std::unordered_map<pair, bool, pair_hash> map;

	double sub_rslt = 0.0;
	double rslt = 0.0;

	for (int g = 0; g < no_of_subgraphs; g++) {
		std::cerr << g << "\n";
		for (std::unordered_set<int>::iterator it = G[pivots[g]].begin(); it != G[pivots[g]].end(); ++it) {
			for (std::unordered_set<int>::iterator it2 = G[pivots[g]].begin(); it2 != G[pivots[g]].end(); ++it2) {
				if (*it != *it2) {
					if (!map[std::make_pair(*it, *it2)] && !map[std::make_pair(*it2, *it)]) { // to make sure not to compute distance more than once

							sub_rslt = rn_dist_for_users_NoT_BFS(*it, *it2);
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

double X_st(std::unordered_map<int, std::unordered_set<int>> G, int pivots[], int no_of_subgraphs) {
	std::unordered_map<pair, bool, pair_hash> map;

	double sub_rslt = 0.0;
	double rslt = 0.0;

	for (int g = 0; g < no_of_subgraphs; g++) {
		for (std::unordered_set<int>::iterator it = G[pivots[g]].begin(); it != G[pivots[g]].end(); ++it) {
			for (std::unordered_set<int>::iterator it2 = G[pivots[g]].begin(); it2 != G[pivots[g]].end(); ++it2) {
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

double X_inf(std::unordered_map<int, std::unordered_set<int>> G, int pivots[], int no_of_subgraphs) {
	std::unordered_map<pair, bool, pair_hash> map;

	double sub_rslt = 0.0;
	double rslt = 0.0;

	for (int g = 0; g < no_of_subgraphs; g++) {
		for (std::unordered_set<int>::iterator it = G[g].begin(); it != G[g].end(); ++it) {
			for (std::unordered_set<int>::iterator it2 = G[pivots[g]].begin(); it2 != G[pivots[g]].end(); ++it2) {
				if (*it != *it2) {

					if (!map[std::make_pair(*it, *it2)] && !map[std::make_pair(*it2, *it)]) {

						sub_rslt = sub_rslt + inf_score_users(*it, *it2);

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
	compute the influence score between two users
*/

double inf_score_users(int src, int dst) {
	double rslt = 0.0;

	if (check_hash_infScore[std::make_pair(src, dst)]) { // if it is already coomputed, we just return it
		rslt = hash_infScore[std::make_pair(src, dst)];
	}
	else { // if nott computed, we run bfs to compute to all vertices, and then we retuen the value
		inf_score_to_all_vertices(src);
		rslt = hash_infScore[std::make_pair(src, dst)];
	}

	return rslt;
}





/*
GIVEN::     RN, src, and dst
ENSURE::    shortest path distance between src and dst
*/
double rn_dist_for_users(int src, int dst) {
	double temp_rslt = 0.0;
	double temp = 0.0;

	if (check_hash_rnToUser_dist[std::make_pair(src, dst)]) { // if we already computed those two users
		return hash_rnToUser_dist[std::make_pair(src, dst)];
	}
	else { // else, compute them and store the result

		for (int i = 0; i < No_CKINs ; ++i) {// for all user locations
			for (int j = 0; j < No_CKINs ; j++) {
				// map the user location to the from the social to the road network, match the 
				// user location to a vertex on the road network
				int src_map = sn_vrtx[src].ckins[i];
				int dst_map = sn_vrtx[dst].ckins[j];

				if (!check_hash_rn_dist[std::make_pair(src_map, dst_map)]) { // if we don't have the distance, the we will compute it
					rn_Dij_to_all_vertices(src_map); // find distance to all other vertices
													 // at the same time store distances to all other vertices
					temp = hash_rn_dist[std::make_pair(src_map, dst_map)];
				}
				else {
					temp = hash_rn_dist[std::make_pair(src_map, dst_map)];
				}

				temp_rslt = temp_rslt + temp;
			}
		}

		check_hash_rnToUser_dist[std::make_pair(src, dst)] = true;
		hash_rnToUser_dist[std::make_pair(src, dst)] = temp_rslt / (No_CKINs * No_CKINs);

		check_hash_rnToUser_dist[std::make_pair(dst, src)] = true;
		hash_rnToUser_dist[std::make_pair(dst, src)] = temp_rslt / (No_CKINs * No_CKINs);

		return temp_rslt / (No_CKINs * No_CKINs);
	}
}

/*
GIVEN::     RN, src, and dst
ENSURE::    shortest path distance between src and dst
*/
double rn_dist_for_users_NoT_BFS(int src, int dst) {
	double temp_rslt = 0.0;
	double temp = 0.0;

	if (check_hash_rnToUser_dist[std::make_pair(src, dst)]) { // if we already computed those two users
		return hash_rnToUser_dist[std::make_pair(src, dst)];
	}
	else { // else, compute them and store the result

		for (int i = 0; i < No_CKINs; ++i) {// for all user locations
			for (int j = 0; j < No_CKINs; j++) {
				// map the user location to the from the social to the road network, match the 
				// user location to a vertex on the road network
				int src_map = sn_vrtx[src].ckins[i];
				int dst_map = sn_vrtx[dst].ckins[j];

				if (!check_hash_rn_dist[std::make_pair(src_map, dst_map)]) { // if we don't have the distance, the we will compute it
					rn_Dij(src_map, dst_map); // find distance to all other vertices
													 // at the same time store distances to all other vertices
					temp = hash_rn_dist[std::make_pair(src_map, dst_map)];
				}
				else {
					temp = hash_rn_dist[std::make_pair(src_map, dst_map)];
				}

				temp_rslt = temp_rslt + temp;
			}
		}

		check_hash_rnToUser_dist[std::make_pair(src, dst)] = true;
		hash_rnToUser_dist[std::make_pair(src, dst)] = temp_rslt / (No_CKINs * No_CKINs);

		check_hash_rnToUser_dist[std::make_pair(dst, src)] = true;
		hash_rnToUser_dist[std::make_pair(dst, src)] = temp_rslt / (No_CKINs * No_CKINs);

		return temp_rslt / (No_CKINs * No_CKINs);
	}
}
////////////////////////////////////////////////////////
/*
GIVEN::     SN, src, and dst
ENSURE::    shortest path distance between src and dst
*/
double sn_dist(int src, int dst) {
	if (!check_hash_sn_dist[std::make_pair(src, dst)]) { // if we don't have the distance, the we will compute it
		sn_Dij_to_all_vertices(src); // find distance to all other vertices
									 // at the same time store distances to all other vertices
		return hash_sn_dist[std::make_pair(src, dst)];
	}
	else {
		return  hash_sn_dist[std::make_pair(src, dst)];
	}
}

////////////////////////////////////////////////////////
/*
GIVEN::     RN, src, and dst
ENSURE::    shortest path distance between src and dst
*/
double rn_dist(int src, int dst) {
	return rn_Dij(src, dst);
}

////////////////////////////////////////////////////////////////////////
/*
REQUIRE::		two vertices a and b
ENSURE ::		rslt--> is the itersection set of edges with the edge a--b
				we return the number of intersected vertices
*/
///////////////////////////////////////////////////////////////////////
int intersect(std::unordered_set<int>* common_edges, int a, int b) {

	double rslt = 0.0;
	std::set<int> intersect;
	common_edges->clear();
	set_intersection(sn_vrtx[a].nbrs_set.begin(), sn_vrtx[a].nbrs_set.end(),
		sn_vrtx[b].nbrs_set.begin(), sn_vrtx[b].nbrs_set.end(),
		std::inserter(intersect, intersect.begin()));
	
	rslt = intersect.size();
	// v -- a, v -- b, a -- v, b -- v

	for (std::set<int>::iterator it = intersect.begin();
		it != intersect.end(); ++it) {

		common_edges->insert(hash_edge[std::make_pair(a, *it)]);
		common_edges->insert(hash_edge[std::make_pair(*it, a)]);

	}

	//std::cerr << rslt << "\n";
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
	std::unordered_set<int>* pool = new std::unordered_set<int>;
	int e;

	for (int e = 0; e < No_sn_E; e++) {
		snEdges[e].sup = intersect(pool, snEdges[e].from, snEdges[e].to);

		//std::cerr << snEdges[e].sup * 2 << pool->size() << "\n";

		HeapEntry* he = new HeapEntry();

		he->son1 = e;
		he->key = snEdges[e].sup;

		hp->insert(he);

		delete he;
	}

label2:
	while (hp->used > 0) {
		HeapEntry* he = new HeapEntry();
		hp->remove(he);
		e = he->son1;

		intersect(pool, snEdges[e].from, snEdges[e].to);

		std::unordered_set<int>::iterator it2;

		it2 = sn_vrtx[snEdges[e].from].nbrs.find(snEdges[e].to);
		if (it2 != sn_vrtx[snEdges[e].from].nbrs.end())
			sn_vrtx[snEdges[e].from].nbrs.erase(it2);



		it2 = sn_vrtx[snEdges[e].to].nbrs.find(snEdges[e].from);

		if (it2 != sn_vrtx[snEdges[e].to].nbrs.end())
			sn_vrtx[snEdges[e].to].nbrs.erase(it2);


		std::set<int>::iterator it22;

		it22 = sn_vrtx[snEdges[e].from].nbrs_set.find(snEdges[e].to);
		if (it22 != sn_vrtx[snEdges[e].from].nbrs_set.end())
			sn_vrtx[snEdges[e].from].nbrs_set.erase(it22);



		it22 = sn_vrtx[snEdges[e].to].nbrs_set.find(snEdges[e].from);

		if (it22 != sn_vrtx[snEdges[e].to].nbrs_set.end())
			sn_vrtx[snEdges[e].to].nbrs_set.erase(it22);


		delete he;
		if (snEdges[e].sup < (k - 2)) {

			for (std::unordered_set<int>::iterator it = pool->begin(); it != pool->end(); ++it) {
				//std::cerr << *it << " ";
				hp->deleteEntry(*it);

				snEdges[*it].sup = snEdges[*it].sup - 1;

				HeapEntry* he = new HeapEntry();
				he->son1 = *it;
				he->key = snEdges[*it].sup;
				hp->insert(he);
				delete he;
			}
		}
		else {
			HeapEntry* he = new HeapEntry();
			he->son1 = e;
			he->key = snEdges[e].sup;
			hp->insert(he);
			delete he;
			break;
		}
	}

	if (hp->used > 0) {
		k = k + 1;
		goto label2;
	}
	hp->~Heap();

	// assign sup of truss values
	
	for (int i = 0; i < No_sn_V; ++i) {
		int max = -INT_MAX;
		
		for (std::unordered_set<int>::iterator it = sn_vrtx[i].myedges.begin();
						it != sn_vrtx[i].myedges.end(); ++it) {

			if (max < snEdges[*it].sup) {
				max = snEdges[*it].sup;
			}

		}

		sn_vrtx[i].truss = max;

	}
}


///////////////////////////////

/*
REQUIRES :: spgraphs
ENSURES  :: build the hybrid index
*/

//امرر كراف كامل, وامرر 
std::unordered_map<int, std::unordered_set<int>> get_index_subgraphs(int piv_set[], int no_piv_set, int
	new_pivots[], int no_new_piv) {
	double qual_rslt;
	int best_quality;
	int assign;

	std::unordered_map<int, std::unordered_set<int>> Gr;
	for (int v = 0; v < no_piv_set; v++) {
		best_quality = INT_MAX;
		qual_rslt = 0.0;
		for (int piv = 0; piv < no_new_piv; piv++) {

			qual_rslt = quality(piv_set[v], new_pivots[piv]);

			if (qual_rslt < best_quality) {
				assign = new_pivots[piv];
				best_quality = qual_rslt;
			}
		}
		Gr[assign].insert(piv_set[v]);
	}
	/*
	for (int i = 0; i < no_new_piv; ++i) {
		std::cerr << "Index_Piv " << new_pivots[i] << " --->> ";
		for (std::unordered_set<int>::iterator it = Gr[new_pivots[i]].begin(); it != Gr[new_pivots[i]].end(); ++it) {
			std::cerr << *it << " ";
		}
		std::cerr << "\n" << "\n";
	}

	*/
	return Gr;
}


//////
/*
	GIVEN	:: a lower index pivot array
	ENSURES :: find level_index pivots
*/

std::unordered_map<int, std::unordered_set<int>> Index_piv_select(int no_new_piv, int prev_piv[], int no_prev_piv, int f_piv[]) {
	// define variables
	double global_cost = -INT_MAX;
	int* S_p = new int[No_index_piv];
	int* new_S_p = new int[No_index_piv];
	double final_cost;
	std::unordered_map<int, std::unordered_set<int>> GG;
	std::unordered_map<int, std::unordered_set<int>> new_GG;
	std::unordered_map<int, std::unordered_set<int>> GGG;


	memset(f_piv, -1, sizeof(f_piv[0]) * No_index_piv);

	double cost_G = 0.0;
	double local_cost = 0.0;
	double new_cost = 0.0;

	int global_iter = 2;
	int swap_iter = 2;

	int get_piv = 0;
	int new_piv = 0;

	bool cost_was_updated = false;
	for (int a = 1; a < global_iter; a++) { // set global iterator

		int git;
		// select random pivots at first
		for (int i = 0; i < no_new_piv; ++i) {
		labelA:
			git = uniform(0, no_prev_piv - 1);
			int git_val = prev_piv[git];  // get the values from the array itselt
			if (!isInTheArray(S_p, no_prev_piv, git_val))
				S_p[i] = git_val;
			else
				goto labelA;
		}

		// get subgraphs based on pivots
		GG = get_index_subgraphs(prev_piv, no_prev_piv, S_p, no_new_piv);


		// evaluate the cost function
		local_cost = evaluate_subgraphs(GG, S_p, no_new_piv);
		//std::cerr << "THe Evaluation :: " << local_cost << "\n \n";


		for (int b = 1; b < swap_iter; b++) {
			get_piv = uniform(0, no_new_piv - 1);

		labelB:
			int git = prev_piv[uniform(0, no_prev_piv - 1)];

			if (!isInTheArray(S_p, no_new_piv, git))
				new_piv = git;
			else
				goto labelB;

			memcpy(new_S_p, S_p, sizeof(S_p[0]) * no_new_piv);
			new_S_p[get_piv] = new_piv;


			new_GG = get_index_subgraphs(prev_piv, no_prev_piv, new_S_p, no_new_piv);


			new_cost = evaluate_subgraphs(new_GG, new_S_p, no_new_piv);

			//std::cerr << "THe Evaluation :: " << new_cost << "\n \n";
			if (new_cost > local_cost) {
				local_cost = new_cost;
				std::memcpy(S_p, new_S_p, sizeof(new_S_p[0]) * no_new_piv);

				// get the final subgraph
				//memcpy(final_G, new_G, sizeof(new_G) * No_index_piv);

				GG = new_GG;
			}
		}
		if (local_cost > global_cost) {
			memcpy(f_piv, S_p, sizeof(S_p[0]) * no_new_piv);
			global_cost = local_cost;
			GGG = GG;
		}
	}


	/*
	// get the new assigned fathers in a hashmap
	for (int i = 0; i < no_new_piv; ++i) {
		std::cerr << "Index_Piv " << f_piv[i] << " --->> ";
		for (std::unordered_set<int>::iterator it = GGG[f_piv[i]].begin(); it != GGG[f_piv[i]].end(); ++it) {
			std::cerr << *it << " ";
		}
		std::cerr << "\n" << "\n";
	}
	std::cerr << "THe Evaluation :: " << global_cost << "\n \n";
	*/
	return GGG;
}



/*
	GIVEN		:: a social network and a number ??, of pivots
	ENSURES		:: the social network pivot set SN_piv_set
*/
void find_social_network_pivots() {

	const int h = No_SN_piv;
	int globelIter = 2, max = No_sn_V, min = -INT_MAX;
	double cost = 0;
	double d1, d2, diff, maxEva = 0;
	double localCost, newCost = 0;
	int swapIter = 2;
	int rand_piv, npiv;
	double globaCost = -INT_MAX;

	int S_p[h], new_S_p[h];
	int new_piv = -1;
	for (int a = 1; a <= globelIter; ++a) {

		// select random pivots at first
		for (int i = 0; i < No_SN_piv; ++i) {
		labelA:
			int git = uniform(0, No_sn_V - 1);
			if (!isInTheArray(S_p, No_SN_piv, git))
				S_p[i] = git;
			else
				goto labelA;
		}

		//evaluate the pivot set
		localCost = evaluate_SN_pivots(S_p, No_SN_piv);

		for (int b = 1; b < swapIter; b++) {

			int get_piv = uniform(0, No_SN_piv - 1);

		labelB:
			int git = uniform(0, No_sn_V - 1);
			if (!isInTheArray(S_p, No_SN_piv, git))
				new_piv = git;
			else
				goto labelB;

			std::memcpy(new_S_p, S_p, sizeof(S_p[0]) * No_SN_piv);

			newCost = evaluate_SN_pivots(new_S_p, No_SN_piv);

			if (newCost > localCost) {
				localCost = newCost;
				std::memcpy(S_p, new_S_p, sizeof(S_p[0]) * No_SN_piv);

			}
		}

		if (localCost > globaCost) {
			memcpy(SN_piv_set, S_p, sizeof(S_p[0]) * No_SN_piv);
			globaCost = localCost;
		}
	}

	// get the distance from the social network pivot to every social network vrtx
	for (int v = 0; v < No_sn_V; ++v) {
		for (int piv = 0; piv < No_SN_piv; piv++) {
			sn_vrtx[v].sn_distToPiv[std::make_pair(v, SN_piv_set[piv])] = sn_dist(v, SN_piv_set[piv]);
		}
	}

	std::cerr << "\n ------ sn pivots ------ \n";
	for (int i = 0; i < No_SN_piv; ++i) {
		std::cerr << SN_piv_set[i] << " --- ";
	}
}


double evaluate_SN_pivots(int S_p[], int no_pivs) {

	double rslt = 0.0;
	double diff = -INT_MAX;
	double temp = -INT_MAX;

	for (int i = 0; i < 100; i++) {

		for (int j = 100; j < 200; j++) {

			for (int k = 0; k < no_pivs; k++) {

				temp = abs(sn_dist(S_p[k], i) - sn_dist(S_p[k], j));

				if (diff < temp)
					diff = temp;
			}

			rslt = rslt + diff;

		}

	}

	return rslt;

}


void find_road_network_pivots() {

	int globelIter = 2;
	int swapIter = 2;

	double cost = 0;
	double localCost, newCost = 0;
	int rand_piv, npiv;
	double globaCost = -INT_MAX;

	int S_p[No_RN_piv], new_S_p[No_RN_piv];
	int new_piv = -1;
	for (int a = 1; a <= globelIter; ++a) {

		// select random pivots at first
		for (int i = 0; i < No_RN_piv; ++i) {
		labelA:
			int git = uniform(0, No_rn_V - 1);
			if (!isInTheArray(S_p, No_RN_piv, git))
				S_p[i] = git;
			else
				goto labelA;
		}

		//evaluate the pivot set
		localCost = evaluate_RN_pivots(S_p, No_RN_piv);

		for (int b = 1; b < swapIter; b++) {

			int get_piv = uniform(0, No_RN_piv - 1);

		labelB:
			int git = uniform(0, No_rn_V - 1);
			if (!isInTheArray(S_p, No_RN_piv, git))
				new_piv = git;
			else
				goto labelB;
			
			std::memcpy(new_S_p, S_p, sizeof(S_p[0]) * No_RN_piv);

			newCost = evaluate_RN_pivots(new_S_p, No_RN_piv);

			if (newCost > localCost) {
				localCost = newCost;
				std::memcpy(S_p, new_S_p, sizeof(S_p[0]) * No_RN_piv);

			}
		}

		if (localCost > globaCost) {
			memcpy(RN_piv_set, S_p, sizeof(S_p[0]) * No_RN_piv);
			globaCost = localCost;
		}
	}

	// save the distance from pivots to all other vertices
	for (int v = 0; v < No_sn_V; ++v) {
		for (int piv = 0; piv < No_RN_piv; piv++) {
			sn_vrtx[v].rn_distToPiv[std::make_pair(v, RN_piv_set[piv])] = rn_dist_rnPiv_to_user(RN_piv_set[piv], v);
		}
	}

	/* REQUIRES :: reoad network and social network users, each social network user 
				   has multiple locations over the road network.
	//ENSURES :: find the average road network distance (with ckins) from each user check-in locations to 
				 the road network pivot.
	

	for (int piv = 0; piv < No_RN_piv; piv++) {

		// run bfs to get distance to all other road network vertices
		rn_Dij_to_all_vertices(RN_piv_set[piv]);

		// compute the average distance from the road network pivot to all
		// social netwok users
		for (int s = 0; s < No_sn_V; ++s) {
			for (int c = 0; c < No_CKINs; c++) {

			}
		}
	}
	*/
	std::cerr << "\n ------ rn pivots ------ \n";
	for (int i = 0; i < No_RN_piv; ++i) {
		std::cerr << RN_piv_set[i] << " --- ";
	}
}


double rn_dist_rnPiv_to_user(int r_piv, int user) {
	double rslt = 0.0;
	double temp = 0.0;

	for (int i = 0; i < No_CKINs - 1; ++i) {

		int src_map = sn_vrtx[user].ckins[i];

		if (!check_hash_rn_dist[std::make_pair(r_piv, src_map)]) { // if we don't have the distance, the we will compute it
			rn_Dij_to_all_vertices(r_piv); // find distance to all other vertices
											 // at the same time store distances to all other vertices
			temp = hash_rn_dist[std::make_pair(r_piv, src_map)];
		}
		else {
			temp = hash_rn_dist[std::make_pair(r_piv, src_map)];
		}

		rslt = rslt + temp;

	}

	return rslt/No_CKINs;
}

double evaluate_RN_pivots(int S_p[], int no_pivs) {

	double rslt = 0.0;
	double diff = -INT_MAX;
	double temp = -INT_MAX;

	for (int i = 0; i < 10; i++) { // for each social network vrtx

		for (int j = 10; j < 20; j++) { // for each social network vrtx

			for (int k = 0; k < no_pivs; k++) { // for each road network pivot

				temp = abs(rn_dist_rnPiv_to_user(S_p[k], i) - rn_dist_rnPiv_to_user(S_p[k], j));

				if (diff < temp)
					diff = temp;
			}
			
			rslt = rslt + diff;

		}

	}
	

	return rslt;
}


//////////////////////////////////////////////////////////////////////
/*
 GIVEN		:: a first layer subgraphs with their corresponding pivots
 ENSURES	:: a tree index as an array of objects
*/
//////////////////////////////////////////////////////////////////////
void indexing() {
	// compute the truss for each user in the social network
	
	truss_decomposition();

	//define a graph
	std::unordered_map<int, std::unordered_set<int>> GGG;

	// get the index pivots set (index_piv[])
	// and devide the subgraph based on those pivots and return the subgraphs in GGG
	GGG = sn_piv_select();


	int assign_counter = INDEXSIZE - 1;

	// for subgraph, # of subgraphs == # of index pivots
	for (int j = 0; j < No_index_piv; ++j) {
		// for each subgraph
		for (std::unordered_set<int>::iterator it = GGG[index_piv[j]].begin(); it != GGG[index_piv[j]].end(); ++it) {
			//store all children of the node in the tree, position assign_counter
			tree[assign_counter].child.insert(*it);

			// store the tree position of each node
			hash_my_position_in_tree[*it] = assign_counter;

			// assign the levelto each node
			tree[assign_counter].level = 0;

		}
		assign_counter--;
	}



	int* f_piv = new int[No_index_piv];
	
	int PivsAtLevelCou = 1;
	int no_new_piv = PivsInLevelTree[PivsAtLevelCou];
	PivsAtLevelCou++;



	int no_prev_piv = No_index_piv;
	int* layer_piv = new int[No_index_piv];
	std::memcpy(layer_piv, index_piv, sizeof(index_piv[0]) * No_index_piv);


	int level = 1;
	while (no_new_piv > 0) {
		// layer the graph and get 
		GGG = Index_piv_select(no_new_piv, layer_piv, no_prev_piv, f_piv);

		//printing
		//*
		for (int i = 0; i < no_new_piv; ++i) {
			std::cerr << "Index_Piv " << f_piv[i] << " --->> ";
			for (std::unordered_set<int>::iterator it = GGG[f_piv[i]].begin(); it != GGG[f_piv[i]].end(); ++it) {
				std::cerr << *it << " ";
			}
			std::cerr << "\n" << "\n";
		}
		//*/


		for (int j = 0; j < no_new_piv; ++j) {
			//trying to assign the children to the node
			for (std::unordered_set<int>::iterator it = GGG[f_piv[j]].begin(); it != GGG[f_piv[j]].end(); ++it) {

				tree[assign_counter].child.insert(hash_my_position_in_tree[*it]);
				
				// set the level
				tree[assign_counter].level = level;

				hash_my_position_in_tree[*it] = assign_counter;
				

			}

			assign_counter--;

		}

		// get a new level
		level = level + 1;

		if (PivsAtLevelCou < LengthOfPivLevels) {
			no_prev_piv = no_new_piv;
			no_new_piv = PivsInLevelTree[PivsAtLevelCou];
			PivsAtLevelCou++;
		}
		else {
			no_new_piv = 0;
		}
		

		std::memcpy(layer_piv, f_piv, sizeof(f_piv[0]) * no_prev_piv);

	}

	
	std::cerr << " --------------------------- " << "\n";
	for (int i = 0; i < INDEXSIZE - No_index_piv; i++) {
		std::cerr << "the node is " << i << " --> ";
		for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) {
			std::cerr << *it << " ";
		}
		std::cerr << std::endl;
	}
}

/*
GIVEN	 :: tree node (index)
ENSURES  :: assigns the parent to each node, parent of root = -1;
*/
void setParentOfNodes() {
	
	
	// set the parent to the root
	tree[0].parent = -1;

	int i = 1;
	while (!(tree[i].level == 0)) {
		for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) {
			
			tree[*it].parent = i;

		}
		i++;
	}
}


/*
	GIVEN	:: a social network
	ENSURES :: this social netowk is connected, if not, then add more edges to make it connected

*/
void get_socialNetwork_connected(){

	std::unordered_set<int>* sett = new std::unordered_set<int>[100000];
	int* havingSet = new int[No_sn_V];
	bool* check_havingSet = new bool[No_sn_V];
	int setCou = 0;
	
	for (int i = 0; i < No_sn_V; ++i)
		check_havingSet[i] = false;

	for (int i = 0; i < No_sn_V; ++i) {
		if ( !check_havingSet[i]) { // if the vrtx have no group

			Heap* hp = new Heap();
			hp->init(2);

			HeapEntry* he = new HeapEntry();
			he->son1 = i;
			he->key = 0;
			hp->insert(he);
			delete he;

			
			while (hp->used > 0) {
				HeapEntry* he = new HeapEntry();
				hp->remove(he);

				int vr = he->son1;

				delete he;
				sett[setCou].insert(vr);
				havingSet[vr] = setCou;
				check_havingSet[vr] = true;
				
				for (std::list<int>::iterator it = snGraph[vr].begin(); it != snGraph[vr].end(); ++it) {
					
					if (!check_havingSet[*it]) { // we put it into the heap if it was not explored before
						HeapEntry* he = new HeapEntry();
						he->son1 = *it;
						he->key = 1;
						hp->insert(he);
						delete he;
					
					}
				}
		
			}

			setCou++;
			hp->~Heap();

		}


	}

	
	std::cerr << setCou;
	
	std::ofstream fout;
	fout.open("e_____.txt");
	
	std::unordered_set<int>::iterator it;
	std::unordered_set<int>::iterator it2;
	
	for (int i = 0; i < setCou - 1; i++) {
			it = sett[i].begin();
			it2 = sett[i + 1].begin();
			fout << *it << " " << *it2 <<" " <<((double)rand() / (RAND_MAX)) << " " << ((double)rand() / (RAND_MAX)) << " " << ((double)rand() / (RAND_MAX)) << "\n";
	}
	
	fout.close();

	delete[] sett;
	delete[] havingSet;
	delete[] check_havingSet;
}

/*
	GIVEN	:: a road network
	ENSURES :: this road netowk is connected, if not, then add more edges to make it connected

*/
void get_roadNetwork_connected() {

	std::unordered_set<int>* sett = new std::unordered_set<int>[100000];
	int* havingSet = new int[No_rn_V];
	bool* check_havingSet = new bool[No_rn_V];
	int setCou = 0;

	for (int i = 0; i < No_rn_V; ++i)
		check_havingSet[i] = false;

	for (int i = 0; i < No_rn_V; ++i) {
		if (!check_havingSet[i]) { // if the vrtx have no group

			Heap* hp = new Heap();
			hp->init(2);

			HeapEntry* he = new HeapEntry();
			he->son1 = i;
			he->key = 0;
			hp->insert(he);
			delete he;


			while (hp->used > 0) {
				HeapEntry* he = new HeapEntry();
				hp->remove(he);

				int vr = he->son1;

				delete he;
				sett[setCou].insert(vr);
				havingSet[vr] = setCou;
				check_havingSet[vr] = true;

				for (std::list<int>::iterator it = rnGraph[vr].begin(); it != rnGraph[vr].end(); ++it) {

					if (!check_havingSet[*it]) { // we put it into the heap if it was not explored before
						HeapEntry* he = new HeapEntry();
						he->son1 = *it;
						he->key = 1;
						hp->insert(he);
						delete he;

					}
				}

			}

			setCou++;
			hp->~Heap();

		}


	}


	std::cerr << setCou;

	std::ofstream fout;
	fout.open("rne_____.txt");

	std::unordered_set<int>::iterator it;
	std::unordered_set<int>::iterator it2;

	for (int i = 0; i < setCou - 1; i++) {
		it = sett[i].begin();
		it2 = sett[i + 1].begin();
		fout << *it << " " << *it2 << "\n";
	}

	fout.close();

	delete[] sett;
	delete[] havingSet;
	delete[] check_havingSet;
}


void expand(int v, std::unordered_map<int, Heap>& hp, int assign[], int& condtition,
		Heap& pivHeap, std::unordered_map<pair, bool, pair_hash>& isSeen, int cand_piv[]) {

	int tt = 0;
	while (hp[v].used > 0) {
		HeapEntry* he = new HeapEntry();
		hp[v].remove(he);

		int cand = he->son1;
		int w = he->key;
		delete he;

		int min_key = INT_MAX;

		if (assign[cand] == -1) {
			assign[cand] = v;
			tt = 1;
			condtition++;
			for (std::list<int>::iterator it = snGraph[cand].begin(); it != snGraph[cand].end(); it++) {


				if (assign[*it] == -1) {
					if (!(isSeen[std::make_pair(v, *it)])) {
						if (!isInTheArray(cand_piv, No_index_piv, *it)) {

							HeapEntry* he = new HeapEntry();
							he->son1 = *it;
							//he->key = w + rn_edge_info[std::make_pair(cand, *it)].weight;
							he->key = w + 1;
							if (min_key > he->key)
								min_key = he->key;

							hp[v].insert(he);
							delete he;

							isSeen[std::make_pair(v, *it)] = true;
						}
					}
				}

			}
		}
		if (min_key< INT_MAX) {
			HeapEntry* he2 = new HeapEntry();
			he2->son1 = v;
			he2->key = min_key;
			pivHeap.insert(he2);
			delete he2;
			break;
		}
	}

}

/*

Gen. subgraphs function
Given a pivot set and a network,
The gen_subgraph_update, will 

*/

std::unordered_map<int, std::unordered_set<int>> gen_subgraphs_update(int cand_piv[]) {

	std::unordered_map<int, Heap> hpp;
	Heap pivHeap;
	pivHeap.init(2);
	std::unordered_map<pair, bool, pair_hash> isSeen;

	for (int piv = 0; piv < No_index_piv; piv++) {
		hpp[cand_piv[piv]].init(2);
		HeapEntry* he = new HeapEntry();
		he->son1 = cand_piv[piv];
		he->key = 0;
		hpp[cand_piv[piv]].insert(he);

		HeapEntry* he2 = new HeapEntry();
		he2->son1 = cand_piv[piv];
		he2->key = 0;
		pivHeap.insert(he2);
		delete he2;
	}


	// create an array to assign vrtices to pivots
	// initiate it with -1
	int* assign = new int[No_sn_V];
	std::memset(assign, -1, sizeof(assign[0]) * No_sn_V);


	int piv = 0;
	int condtition = 0;
	while (condtition < No_sn_V) {
		
		HeapEntry* he = new HeapEntry();
		pivHeap.remove(he);
		piv = he->son1;
		delete he;

		expand(piv, hpp, assign, condtition, pivHeap, isSeen, cand_piv);

	}


	std::unordered_map<int, std::unordered_set<int>> GG;
	for (int i = 0; i < No_sn_V; i++) {
		GG[assign[i]].insert(i);
	}

	///*

	std::cerr << "\n" << "---- size of parts ----" << "\n";
	for (int p = 0; p < No_index_piv; p++) {
		//std::cerr << cand_piv[p] << " -->> ";
		//for (std::unordered_set<int>::iterator it = GG[cand_piv[p]].begin(); it != GG[cand_piv[p]].end(); ++it) {
		//	std::cerr << *it << " ";
		//}
		//std::cerr << "\n";

		std::cerr << GG[cand_piv[p]].size() << " ";
	}
	//*/

	/*
	std::cerr <<"\n" <<"----all assignments----" << "\n";
	for (int i = 0; i < No_rn_V; i++) {
		std::cerr << i << " " << assign[i]<< "\n" ;
	}
	*/

	return GG;
}

void print_sn_vrtx_coordinates() {
	for (int i = 0; i < No_sn_V; ++i) {

	}
}


////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////

double Ecl(double x1, double y1, double x2, double y2) {
	double x = x1 - x2; //calculating number to square in next step
	double y = y1 - y2;
	double dist;

	dist = pow(x, 2) + pow(y, 2);       //calculating Euclidean distance
	dist = sqrt(dist);

	return dist;
}

void mapping_socialVertices_to_roads(char* sn_users, char* rn_users, char* output) {
	std::ofstream out_;
	out_.open(output);

	rn_vertecies* rn_v = new rn_vertecies[No_rn_V];
	rn_vertecies* sn_v = new rn_vertecies[No_sn_V];

	FILE* snFile = fopen(sn_users, "r");
	FILE* rnFile = fopen(rn_users, "r");
	
	if (snFile == NULL || rnFile == NULL) {
		std::cout << "could not open the file to read social network vertices";
		return;
	}

	int getRide;
	fscanf(snFile, "%d", &getRide);
	int cont = 0;
	while (!feof(snFile)) {
		fscanf(snFile, "%d %lf %lf", &sn_v[cont].id, &sn_v[cont].x, &sn_v[cont].y);
		cont++;
	}

	fscanf(rnFile, "%d", &getRide);
	cont = 0;
	while (!feof(rnFile)) {
		fscanf(rnFile, "%d %lf %lf", &rn_v[cont].id, &rn_v[cont].x, &rn_v[cont].y);
		cont++;
	}



	int* ass = new int[No_sn_V];
	
	

	for (int s = 0; s < No_sn_V; s++) {
		double min = +INT_MAX;
		Heap* hp = new Heap();
		hp->init(2);
		for (int r = 0; r < No_rn_V; r++) {

			double E_dist = Ecl(sn_v[s].x, sn_v[s].y, rn_v[r].x, rn_v[r].y);

			if (hp->used < 6) {		
				HeapEntry* he = new HeapEntry();
				he->son1 = r;
				he->key = - E_dist;
				hp->insert(he);
				delete he;
			}
			else {
				HeapEntry* he = new HeapEntry();
				hp->remove(he);

				if ((-he->key) > E_dist) {
					delete he;
					HeapEntry* he = new HeapEntry();
					he->son1 = r;
					he->key = -E_dist;
					hp->insert(he);
				}
				else {
					hp->insert(he);
					delete he;
				}
			}
		}
		
		out_ << s;
		int a[6];
		int i = 5;
		HeapEntry* he = new HeapEntry();
		while (hp->used> 0) {
			hp->remove(he);
			a[i] = he->son1;
			--i;
		}

		for (int i = 0; i < 6; i++)
			out_ << " " << a[i];
		for (int i = 0; i < 10; i++)
			out_ << " " << uniform (0, 9);

		out_ << "\n";
		delete he;
		delete hp;
	}
	
	delete[] sn_v;
	delete[] rn_v;
}


float gaussian(float mean, float sigma) {
	float v1, v2;
	float s;
	float x;

	do {
		v1 = 2 * uniform_dou(0, 1) - 1;
		v2 = 2 * uniform_dou(0, 1) - 1;
		s = v1 * v1 + v2 * v2;
	} while (s >= 1.);

	x = v1 * sqrt(-2. * log(s) / s);

	/*  x is normally distributed with mean 0 and sigma 1.  */
	x = x * sigma + mean;

	return (x);
}
#endif // !INDEX_HPP
