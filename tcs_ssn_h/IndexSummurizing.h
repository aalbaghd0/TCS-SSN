#ifndef INDEX_SUMM
#define INDEX_SUMM
#include "parameterSettings.h"
#include "vertex.h"
#include "../buildIndex.h"
/*
	GIVEN	:: an array of dnodes
	ENSURES :: summurize the data required on the index for the traversal
*/


//////////////////////////////////////////////////////////////////////////////////
/*
	GIVEN	:: an array of dnodes
	ENSURES :: summurize the maximum truss value in each node
*/
void summraize_truss() {
	for (int i = INDEXSIZE - 1; i >= 0; i--) { // for all nodes from bottom to up

		int max_t = -INT_MAX;

		for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) { // for all children in this node
			if (tree[i].level == 0) { // we deal with objects
				if (max_t < sn_vrtx[*it].truss)
					max_t = sn_vrtx[*it].truss;
			}
			else { // not a leaf node
				if (max_t < tree[*it].truss)
					max_t = tree[*it].truss;
			}
		}
		tree[i].truss = max_t;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/*
	GIVEN	:: an array of dnodes
	ENSURES :: summurize the maximum and the minimum road network distance to the road network pivots in each node
*/
void summraize_RN_lb_dist() {

	for (int i = INDEXSIZE - 1; i >= 0; i--) { // for all nodes from the bottom up

		if (tree[i].level == 0) { // if the node is a leaf

			int min_dist = INT_MAX;
			int min_piv = -1;

			int max_dist = -INT_MAX;
			int max_piv = -1;

			for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) { // for all children in the node

				for (int r_piv = 0; r_piv > No_RN_piv; r_piv++) { // for all RN pivots

					int val = rn_dist_for_users(*it, RN_piv_set[r_piv]); // get the distance from vertices to pivots

					if (min_dist > val) { // grt the min
						min_dist = val;
						min_piv = RN_piv_set[r_piv];
					}

					if (max_dist < val) { // get the maximum distance
						max_dist = val;
						max_piv = RN_piv_set[r_piv];
					}

				}

			}

			// here we save the distances
			tree[i].rn_min_dist_to_piv = min_dist;
			tree[i].rn_minimum_piv = min_piv;

			tree[i].rn_max_dist_to_piv = max_dist;
			tree[i].rn_maximum_piv = max_piv;


		}
		else { // if it is an intermediate node

			int min_dist = INT_MAX;
			int min_piv = -1;

			int max_dist = -INT_MAX;
			int max_piv = -1;

			for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) { // for all node children in the node

				if (min_dist > tree[i].rn_min_dist_to_piv) { // grt the min
					min_dist = tree[i].rn_min_dist_to_piv;
					min_piv = tree[i].rn_minimum_piv;
				}

				if (max_dist < tree[i].rn_max_dist_to_piv) { // grt the min
					min_dist = tree[i].rn_max_dist_to_piv;
					min_piv = tree[i].rn_maximum_piv;
				}

			}

			// here we save the distances
			tree[i].rn_min_dist_to_piv = min_dist;
			tree[i].rn_minimum_piv = min_piv;

			tree[i].rn_max_dist_to_piv = max_dist;
			tree[i].rn_maximum_piv = max_piv;

		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////
/*
	GIVEN	:: an array of dnodes
	ENSURES :: summurize the maximum and the minimum social distance to the social network pivots in each node
*/
void summraize_SN_lb_dist() {

	for (int i = INDEXSIZE - 1; i >= 0; i--) { // for all nodes from the bottom up

		if (tree[i].level == 0) { // if the node is a leaf

			int min_dist = INT_MAX;
			int min_piv = -1;

			int max_dist = -INT_MAX;
			int max_piv = -1;

			for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) { // for all children in the node

				for (int s_piv = 0; s_piv > No_RN_piv; s_piv++) { // for all RN pivots

					int val = sn_dist(*it, RN_piv_set[s_piv]); // get the distance from vertices to pivots

					if (min_dist > val) { // grt the min
						min_dist = val;
						min_piv = RN_piv_set[s_piv];
					}

					if (max_dist < val) { // get the maximum distance
						max_dist = val;
						max_piv = RN_piv_set[s_piv];
					}

				}

			}

			// here we save the distances
			tree[i].sn_min_dist_to_piv = min_dist;
			tree[i].sn_minimum_piv = min_piv;

			tree[i].sn_max_dist_to_piv = max_dist;
			tree[i].sn_maximum_piv = max_piv;


		}
		else { // if it is an intermediate node

			int min_dist = INT_MAX;
			int min_piv = -1;

			int max_dist = -INT_MAX;
			int max_piv = -1;

			for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) { // for all node children in the node

				if (min_dist > tree[i].sn_min_dist_to_piv) { // grt the min
					min_dist = tree[i].sn_min_dist_to_piv;
					min_piv = tree[i].sn_minimum_piv;
				}

				if (max_dist < tree[i].sn_max_dist_to_piv) { // grt the min
					min_dist = tree[i].sn_max_dist_to_piv;
					min_piv = tree[i].sn_maximum_piv;
				}

			}

			// here we save the distances
			tree[i].sn_min_dist_to_piv = min_dist;
			tree[i].sn_minimum_piv = min_piv;

			tree[i].sn_max_dist_to_piv = max_dist;
			tree[i].sn_maximum_piv = max_piv;

		}
	}
}

//////////////////////////////////////////////////////////////////////////
/*
	GIVEN	:: an array of dnodes
	ENSURES :: summurize the keyword vector in each node
*/
void summurize_keywords() {
	for (int i = INDEXSIZE - 1; i >= 0; i--) {

		for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) {

			if (tree[i].level == 0) {

				tree[i].keys |= sn_vrtx[*it].key;

			}
			else {

				tree[i].keys |= tree[*it].keys;

			}

		}
	}
}

////////////////////////////////////////////////////////////////////////////////////
/*
	GIVEN	:: an index tree
	ENSURES	:: summurize the influence score -- in
*/
void summ_ub_in_influence_score() {

	for (int i = INDEXSIZE - 1; i >= 0; i--) { // for all nodes from the bottom up

		if (tree[i].level == 0) { // if the node is a leaf

			// for all u /in e_i
			for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) { // for all children in the node

				// for all v neighbours of u
				// in_nbrs are in the in neighbours of the vertex 
				for (std::unordered_set<int>::iterator it2 = sn_vrtx[*it].in_nbrs.begin(); it2 != sn_vrtx[*it].in_nbrs.end(); ++it2) {

					for (int k = 0; k < No_of_TOPICS; k++) {

						if (tree[i].in_topics_prob[k] < sn_edge_info[std::make_pair(*it, *it2)].topics[k]) {

							tree[i].in_topics_prob[k] = sn_edge_info[std::make_pair(*it, *it2)].topics[k];
						}

					}

				}

			}

		}
		else { // if it is an intermediate node

			// for all u /in e_i
			for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) { // for all children in the node

				// for all v neighbours of u
				// in_nbrs are in the in neighbours of the vertex 

				for (int k = 0; k < No_of_TOPICS; k++) {

					if (tree[i].in_topics_prob[k] < tree[*it].in_topics_prob[k]) {

						tree[i].in_topics_prob[k] = tree[*it].in_topics_prob[k];
					}

				}

			}

		}
	}
}
////////////////////////////////////////////////////////////////
/*
	GIVEN	:: an index tree
	ENSURES	:: summurize the influence score -- out
*/
void summ_ub_out_influence_score() {

	for (int i = INDEXSIZE - 1; i >= 0; i--) { // for all nodes from the bottom up

		if (tree[i].level == 0) { // if the node is a leaf

			// for all u /in e_i
			for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) { // for all children in the node

				// for all v neighbours of u
				// in_nbrs are in the in neighbours of the vertex 
				for (std::unordered_set<int>::iterator it2 = sn_vrtx[*it].out_nbrs.begin(); it2 != sn_vrtx[*it].out_nbrs.end(); ++it2) {

					for (int k = 0; k < No_of_TOPICS; k++) {

						if (tree[i].out_topics_prob[k] < sn_edge_info[std::make_pair(*it, *it2)].topics[k]) {

							tree[i].out_topics_prob[k] = sn_edge_info[std::make_pair(*it, *it2)].topics[k];

						}

					}

				}

			}

		}
		else { // if it is an intermediate node

			// for all u /in e_i
			for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) { // for all children in the node

				// for all v neighbours of u
				// in_nbrs are in the in neighbours of the vertex 

				for (int k = 0; k < No_of_TOPICS; k++) {

					if (tree[i].out_topics_prob[k] < tree[*it].out_topics_prob[k]) {

						tree[i].out_topics_prob[k] = tree[*it].out_topics_prob[k];

					}

				}

			}

		}
	}
}

#endif // !INDEXTRAVERSE

