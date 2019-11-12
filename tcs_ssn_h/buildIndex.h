#pragma once
#ifndef INDEX_HPP
#define INDEX_HPP

#include <iostream>
#include <list>
#include "tcs_ssn_h/parameterSettings.h"
#include "tcs_ssn_h/vertex.h"
#include <set>
#include "tcs_ssn_h/heap.h"

void gen_subgraphs(std::list<int> rnGraph, std::list<int> snGraph, edges rnEdges, edges snEdges, 
	int sn_piv[], int sn_piv_no ) {
	int G[index_piv];
	
	for (int i = 0; i < No_sn_V; i++) {
		int best_quality = INT_MAX;
		for (int j = 0; j < index_piv; j++) {
			
		}
	}
}

void sn_piv_select(int indx_piv_num, std::list<int> snGraph, edges snEdges) {

}



#endif // !INDEX_HPP
