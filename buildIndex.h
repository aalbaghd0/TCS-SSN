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

double quality(int v, int piv) {
	return (rn_dist(v, piv) + sn_dist(v, piv));
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
	for (int a = 1; a < global_iter; a++) {
		
		// select random pivots at first
		for (int i = 0; i < No_index_piv; ++i) {
			S_p[i] = uniform(0, No_sn_V);
		}
		// get subgraphs based on pivots
		gen_subgraphs(S_p, G);

		// evaluate the cost function
		local_cost = evaluate_subgraphs(G);
		
		for (int  b = 1; b < swap_iter; b++) {
			get_piv = uniform(0, No_index_piv);
			new_piv = uniform(0, No_sn_V);

			memcpy(new_S_p, S_p, No_index_piv);
			new_S_p[get_piv] = new_piv;
			
			gen_subgraphs(new_S_p, new_G);
			new_cost = evaluate_subgraphs(new_G);

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
}


double evaluate_subgraphs(std::set<int> G[]) {
	double rslt1 = W1 *  X_sc(G);
	double rslt2 = W2 * ( 1 - X_st(G)  );
	double rslt3 = W1 * ( 1 - X_inf(G) );

	return rslt1 + rslt2 + rslt3;
}

double X_sc(std::set<int> G[]) {
	std::pair <int, int> p1, p2, p3;
	std::unordered_map<std::pair<int, int>, bool> map;

	double sub_rslt = 0.0;
	double rslt = 0.0;

	for (int g = 0; g < No_subgraphs; g++) {
		for (std::set<int>::iterator it = G[g].begin(); it != G[g].end(); ++it) {
			for (std::set<int>::iterator it2 = G[g].begin(); it2 != G[g].end(); ++it2) {
				if (*it != *it2) {
					p1 = std::make_pair(*it,  *it2);
					p2 = std::make_pair(*it2, *it);

					if (!map[p1] && !map[p2]) {
						sub_rslt = sub_rslt + rn_dist(*it, *it2);
						map[p1] = true;
						map[p2] = true;
					}
				}
			}
		}
		rslt = rslt + sub_rslt;
		sub_rslt = 0.0;
	}
	return rslt;
}

double X_st(std::set<int> G[]) {
	std::pair <int, int> p1, p2, p3;
	std::unordered_map<std::pair<int, int>, bool> map;

	double sub_rslt = 0.0;
	double rslt = 0.0;

	for (int g = 0; g < No_subgraphs; g++) {
		for (std::set<int>::iterator it = G[g].begin(); it != G[g].end(); ++it) {
			for (std::set<int>::iterator it2 = G[g].begin(); it2 != G[g].end(); ++it2) {
				if (*it != *it2) {
					p1 = std::make_pair(*it, *it2);
					p2 = std::make_pair(*it2, *it);

					if (!map[p1] && !map[p2]) {
						sub_rslt = sub_rslt + ((truss(*it) + truss(*it2)) / sn_dist(*it, *it2));
						map[p1] = true;
						map[p2] = true;
					}
				}
			}
		}
		rslt = rslt + sub_rslt;
		sub_rslt = 0.0;
	}
	return rslt;
}

double X_inf(std::set<int> G[]) {
	std::pair <int, int> p1, p2, p3;
	std::unordered_map<std::pair<int, int>, bool> map;

	double sub_rslt = 0.0;
	double rslt = 0.0;

	for (int g = 0; g < No_subgraphs; g++) {
		for (std::set<int>::iterator it = G[g].begin(); it != G[g].end(); ++it) {
			for (std::set<int>::iterator it2 = G[g].begin(); it2 != G[g].end(); ++it2) {
				if (*it != *it2) {
					p1 = std::make_pair(*it, *it2);
					p2 = std::make_pair(*it2, *it);

					if (!map[p1] && !map[p2]) {
						sub_rslt = sub_rslt + (inf_score(*it, *it2));
						map[p1] = true;
						map[p2] = true;
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
double rn_dist(int src, int dst) {
	double temp_rslt = 0.0;

	for (int i = 0; i < No_CKINs; ++i) {
		for ( int j = 0; j < No_CKINs; j++) {
			// map the user location to the from the social to the road network, match the 
			// user location to a vertex on the road network
			int src_map = mapping(src, sn_vrtx[src].ckins[i].first, sn_vrtx[src].ckins[i].second);
			int dst_map = mapping(dst, sn_vrtx[dst].ckins[j].first, sn_vrtx[dst].ckins[j].second);

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

///

/*
GIVEN	::	social network graph
ENSURE	::	the maximum support of each edge
*/
void truss_decomposition() {
	int k = 3;
	int get_edge;
	std::pair<int, int> pr;
	Heap* hp = new Heap();
	std::unordered_map<std::pair<int, int>, int> is_there;
	
	for (int e = 0; e < No_sn_E; e++) {
		snEdges[e].sup = intersect(sn_vrtx[snEdges[e].from].nbrs, sn_vrtx[snEdges[e].to].nbrs);
		HeapEntry* he = new HeapEntry();
		he->son1 = e;
		he->key = snEdges[e].sup;
		delete he;
	}

	while (true) {
		HeapEntry* he = new HeapEntry();
		hp->remove(he);
		if (he->key <= (k - 2)) {
			delete he;
			hp->~Heap();
			break;
		}

		get_edge = he->son1;
		delete he;
		int mx = max(snEdges[get_edge].from, snEdges[get_edge].to);
		for (std::set<int>::iterator it = sn_vrtx[mx].nbrs.begin();
			it != sn_vrtx[mx].nbrs.end(); ++it) {
			pr = std::make_pair(mx, *it);
			if (is_there[pr] != 100) {
				is_there[std::make_pair(mx, *it)] = 100;
			}
		}
	}
	

}


double intersect(std::set<int> a, std::set<int> b) {
	double rslt = 0.0;

	return rslt;
}


void bucketSort(float arr[], int n) {
	// 1) Create n empty buckets 
	vector<float> b[n];

	// 2) Put array elements in different buckets 
	for (int i = 0; i < n; i++) {
		int bi = n * arr[i]; // Index in bucket 
		b[bi].push_back(arr[i]);
	}

	// 3) Sort individual buckets 
	for (int i = 0; i < n; i++)
		sort(b[i].begin(), b[i].end());

	// 4) Concatenate all buckets into arr[] 
	int index = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < b[i].size(); j++)
			arr[index++] = b[i][j];
}

//////
#endif // !INDEX_HPP
