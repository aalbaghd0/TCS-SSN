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


double quality(int v, int piv) {
	return (??);
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
	int* P = new int[No_index_piv];
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
			memcpy(P, S_p, No_index_piv);
			global_cost = local_cost;
		}
	}
	return P;
}


double evaluate_subgraphs(std::set<int> G[]) {

}


#endif // !INDEX_HPP
