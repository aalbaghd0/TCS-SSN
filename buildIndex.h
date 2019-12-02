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
double X_sc(std::unordered_map<int, std::set<int>> G, int pivots[], int number_subgraphs);
double X_st(std::unordered_map<int, std::set<int>> G, int pivots[], int number_subgraphs);
double X_inf(std::unordered_map<int, std::set<int>> G, int pivots[],int number_subgraphs);
double evaluate_subgraphs(std::unordered_map<int, std::set<int>> G, int pivots[], int number_subgraphs);
double evaluate_Indexsubgraphs(std::set<int> G[], int no_of_subgraphs);
double inf_score_users(int src, int dst);


double quality(int v, int piv) {
	double sn_dist_rslt = 0.0, rn_dist_rslt = 0.0;
	

	if ( ! check_hash_sn_dist[std::make_pair(piv, v)]) { // if we don't have the distance, the we will compute it
		sn_Dij_to_all_vertices(piv); // find distance to all other vertices
							         // at the same time store distances to all other vertices
		sn_dist_rslt = hash_sn_dist[std::make_pair(piv, v)];
	}
	else {
		sn_dist_rslt = hash_sn_dist[std::make_pair(piv, v)];
	}

	return (rn_dist_for_users(piv, v) + sn_dist_rslt);
}

std::unordered_map<int, std::set<int>> gen_subgraphs(int cand_index_piv[]) {
	double qual_rslt;
	int best_quality;
	int assign;

	std::unordered_map<int, std::set<int>> GG;

	for (int v = 0; v < No_sn_V; v++) {
		best_quality = INT_MAX;
		qual_rslt = 0.0;
		for (int piv = 0; piv < No_index_piv; piv++) {
			qual_rslt = quality(v, cand_index_piv[piv]);
			
			if (qual_rslt < best_quality) {
				assign = cand_index_piv[piv];
				best_quality = qual_rslt;
			}
		}
		GG[assign].insert(v);
	}

	/*
	for (int i = 0; i < No_index_piv; ++i) {
		std::cerr << "Index_Piv " << cand_index_piv[i] << " --->> ";
		for (std::set<int>::iterator it = GG[cand_index_piv[i]].begin(); it != GG[cand_index_piv[i]].end(); ++it) {
			std::cerr << *it << " ";
		}
		std::cerr << "\n" << "\n";
	}
	*/
	return GG;
}

std::unordered_map<int, std::set<int>> sn_piv_select() {
	double global_cost = -INT_MAX;
	int* S_p = new int[No_index_piv];
	int* new_S_p = new int[No_index_piv];
	double final_cost;
	std::unordered_map<int, std::set<int>> GG;
	std::unordered_map<int, std::set<int>> new_GG;
	std::unordered_map<int, std::set<int>> G;
	

	double cost_G = 0.0;
	double local_cost = 0.0;
	double new_cost = 0.0;

	int global_iter = 3;
	int swap_iter = 4;

	int get_piv = 0;
	int new_piv = 0;

	bool cost_was_updated = false;
	for (int a = 1; a < global_iter; a++) {
		
		// select random pivots at first
		for (int i = 0; i < No_index_piv; ++i) {
			labelA:
			int git = uniform(0, No_sn_V);
			if (!isInTheArray(S_p, No_index_piv, git))
				S_p[i] = git;
			else
				goto labelA;
		}
		// get subgraphs based on pivots
		GG = gen_subgraphs(S_p);

		// evaluate the cost function
		local_cost = evaluate_subgraphs(GG, S_p, No_subgraphs);
		//std::cerr << "THe Evaluation :: " << local_cost << "\n \n";
		
		
		for (int  b = 1; b < swap_iter; b++) {
			get_piv = uniform(0, No_index_piv);
			
		labelB:
			int git = uniform(0, No_sn_V);
			if (!isInTheArray(S_p, No_index_piv, git))
				new_piv = git;
			else
				goto labelB;
			

			std::memcpy(new_S_p, S_p, sizeof(S_p[0]) * No_index_piv);
			
			new_S_p[get_piv] = new_piv;
			
			new_GG = gen_subgraphs(new_S_p);
			
			

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


		for (int i = 0; i < No_index_piv; ++i) {
			std::cerr << "Index_Piv " << index_piv[i] << " --->> ";
			for (std::set<int>::iterator it = G[index_piv[i]].begin(); it != G[index_piv[i]].end(); ++it) {
				std::cerr << *it << " ";
			}
			std::cerr << "\n" << "\n";
		}
		std::cerr << "THe Evaluation :: " << global_cost << "\n \n";

		return G;
}


double evaluate_subgraphs(std::unordered_map<int, std::set<int>> G, int pivots[], int no_of_subgraphs) {
	double rslt1 = W1 *  X_sc(G, pivots, no_of_subgraphs);
	double rslt2 = W2 * ( 1 - X_st(G, pivots,no_of_subgraphs)  );
	double rslt3 = W1 * ( 1 - X_inf(G, pivots,no_of_subgraphs) );

	return rslt1 + rslt2 + rslt3;
}



double X_sc(std::unordered_map<int, std::set<int>> G, int pivots[], int no_of_subgraphs) {
	std::unordered_map<pair, bool, pair_hash> map;

	double sub_rslt = 0.0;
	double rslt = 0.0;

	for (int g = 0; g < no_of_subgraphs; g++) {
		for (std::set<int>::iterator it = G[pivots[g]].begin(); it != G[pivots[g]].end(); ++it) {
			for (std::set<int>::iterator it2 = G[pivots[g]].begin(); it2 != G[pivots[g]].end(); ++it2) {
				if (*it != *it2) {
					if (!map[std::make_pair(*it, *it2)] && !map[std::make_pair(*it2, *it)]) { // to make sure not to compute distance more than once
						sub_rslt = rn_dist_for_users(*it, *it2);						
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

double X_st(std::unordered_map<int, std::set<int>> G, int pivots[], int no_of_subgraphs) {
	std::unordered_map<pair, bool, pair_hash> map;

	double sub_rslt = 0.0;
	double rslt = 0.0;

	for (int g = 0; g < no_of_subgraphs; g++) {
		for (std::set<int>::iterator it = G[pivots[g]].begin(); it != G[pivots[g]].end(); ++it) {
			for (std::set<int>::iterator it2 = G[pivots[g]].begin(); it2 != G[pivots[g]].end(); ++it2) {
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

double X_inf(std::unordered_map<int, std::set<int>> G, int pivots[], int no_of_subgraphs) {
	std::unordered_map<pair, bool, pair_hash> map;

	double sub_rslt = 0.0;
	double rslt = 0.0;

	for (int g = 0; g < no_of_subgraphs; g++) {
		for (std::set<int>::iterator it = G[g].begin(); it != G[g].end(); ++it) {
			for (std::set<int>::iterator it2 = G[pivots[g]].begin(); it2 != G[pivots[g]].end(); ++it2) {
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
		
		for (int i = 0; i < No_CKINs; ++i) {// for all user locations
			for (int j = 0; j < No_CKINs; j++) {
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
int intersect(std::set<int>* common_edges, int a, int b) {

	double rslt = 0.0;
	std::set<int> intersect;
	common_edges->clear();
	set_intersection(sn_vrtx[a].nbrs.begin(), sn_vrtx[a].nbrs.end(),
		sn_vrtx[b].nbrs.begin(), sn_vrtx[b].nbrs.end(),
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
	std::set<int>* pool = new std::set<int>;
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

		std::set<int>::iterator it2;

		it2 = sn_vrtx[snEdges[e].from].nbrs.find(snEdges[e].to);
		if(it2 != sn_vrtx[snEdges[e].from].nbrs.end())
			sn_vrtx[snEdges[e].from].nbrs.erase(it2);


		it2 = sn_vrtx[snEdges[e].to].nbrs.find(snEdges[e].from);
		
		if (it2 != sn_vrtx[snEdges[e].to].nbrs.end())
			sn_vrtx[snEdges[e].to].nbrs.erase(it2);


		delete he;
		if (snEdges[e].sup < (k - 2)) {

			for (std::set<int>::iterator it = pool->begin(); it != pool->end(); ++it) {
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
std::unordered_map<int, std::set<int>> get_index_subgraphs(int piv_set[], int no_piv_set,int 
									new_pivots[], int no_new_piv) {
	double qual_rslt;
	int best_quality;
	int assign;

	std::unordered_map<int, std::set<int>> Gr;
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
		for (std::set<int>::iterator it = Gr[new_pivots[i]].begin(); it != Gr[new_pivots[i]].end(); ++it) {
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

std::unordered_map<int, std::set<int>> Index_piv_select(int no_new_piv, int prev_piv[], int no_prev_piv, int f_piv[]) {
	// define variables
	double global_cost = -INT_MAX;
	int* S_p = new int[No_index_piv];
	int* new_S_p = new int[No_index_piv];
	double final_cost;
	std::unordered_map<int, std::set<int>> GG;
	std::unordered_map<int, std::set<int>> new_GG;
	std::unordered_map<int, std::set<int>> GGG;

	
	memset(f_piv, -1, sizeof(f_piv[0]) * No_index_piv);

	double cost_G = 0.0;
	double local_cost = 0.0;
	double new_cost = 0.0;

	int global_iter = 5;
	int swap_iter = 10;

	int get_piv = 0;
	int new_piv = 0;

	bool cost_was_updated = false;
	for (int a = 1; a < global_iter; a++) { // set global iterator

		int git;
		// select random pivots at first
		for (int i = 0; i < no_new_piv; ++i) {
		labelA:
			git = uniform(0, no_prev_piv);
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
			get_piv = uniform(0, no_new_piv);

		labelB:
			int git = prev_piv[uniform(0, no_prev_piv)];
			
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
	
	
	///*
	// get the new assigned fathers in a hashmap
	for (int i = 0; i < no_new_piv; ++i) {
		std::cerr << "Index_Piv " << f_piv[i] << " --->> ";
		for (std::set<int>::iterator it = GGG[f_piv[i]].begin(); it != GGG[f_piv[i]].end(); ++it) {
			std::cerr << *it << " ";
		}
		std::cerr << "\n" << "\n";
	}
	std::cerr << "THe Evaluation :: " << global_cost << "\n \n";
	//*/
	return GGG;
}


/*
 GIVEN		:: a first layer subgraphs with their corresponding pivots
 ENSURES	:: a tree index as an array of objects
*/
void indexing() {

	int c = 0;
	std::cerr << number_nodes(No_index_piv, 2, c) << " " << c;


	truss_decomposition();
	std::unordered_map<int, std::set<int>> GGG;

	GGG = sn_piv_select();
	
	int assign_counter = c - 1;
	for (int j = 0; j < No_index_piv; ++j) {
		for (std::set<int>::iterator it = GGG[index_piv[j]].begin(); it != GGG[index_piv[j]].end(); ++it) {
			tree[assign_counter].child.insert(*it);
			hash_father_list[*it] = assign_counter;
			tree[assign_counter].level = 0;
			std::cerr << "the vertex " << *it << " its position in the tree " << assign_counter << "\n";
		}
		assign_counter--;
	}

	std::cerr << "---- the index pivots final -------";
	for (int i = 0; i < No_index_piv; i++) {
		std::cerr << index_piv[i] << " ";
	}
	std::cerr << "----------";

	int* f_piv = new int[100];
	int no_new_piv = 2;

	int no_prev_piv = No_index_piv;
	int* layer_piv = new int[No_index_piv];
	std::memcpy(layer_piv, index_piv, sizeof(index_piv[0]) * No_index_piv);

	std::cerr << "\n" << "---- the index pivots final -------";
	for (int i = 0; i < No_index_piv; i++) {
		std::cerr << layer_piv[i] << " ";
	}
	std::cerr << "----------";

	while (no_new_piv > 0) {
		GGG = Index_piv_select(no_new_piv, layer_piv, no_prev_piv, f_piv);

		for (int j = 0; j < no_new_piv; ++j) {
			//trying to assign the children to the node
			for (std::set<int>::iterator it = GGG[f_piv[j]].begin(); it != GGG[f_piv[j]].end(); ++it) {
				//tree[assign_counter].child.insert(*it);
				tree[assign_counter].child.insert(hash_father_list[*it]);
				
				//tree[assign_counter].ptr.insert(hash_father_list[*it]);
				hash_father_list[*it] = assign_counter;
			}
			assign_counter--;
		}

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
	for (int i = 0; i < c; i++) {
		std::cerr << "the node is " << i << " --> ";
		for (std::set<int>::iterator it = tree[i].ptr.begin(); it != tree[i].ptr.end(); ++it) {
			std::cerr << *it;
			
		}
		std::cerr << " the parent is:: " << tree[i].parent;
		std::cerr << std::endl;
	}
	std::cerr << " --------------------------- " << "\nl";
	for (int i = 0; i < c; i++) {
		std::cerr << "the node is " << i << " --> ";
		for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) {
			std::cerr << *it;
		}
		std::cerr << std::endl;
	}
}

#endif // !INDEX_HPP
