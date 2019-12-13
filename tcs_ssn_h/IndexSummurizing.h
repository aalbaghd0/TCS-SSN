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

			double min_dist = INT_MAX;
			int min_piv = -1;

			int max_dist = -INT_MAX;
			int max_piv = -1;

			for (int r_piv = 0; r_piv > No_RN_piv; r_piv++) { // for all RN pivots

				for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) { // for all children in the node


					double val = rn_dist_for_users(*it, RN_piv_set[r_piv]); // get the distance from vertices to pivots

					if (min_dist > val) { // grt the min
						min_dist = val;
						//min_piv = RN_piv_set[r_piv];
					}

					if (max_dist < val) { // get the maximum distance
						max_dist = val;
						//max_piv = RN_piv_set[r_piv];
					}
				}

				tree[i].rn_min_dist_to_piv[std::make_pair(i, min_piv)] = min_dist;
				tree[i].rn_max_dist_to_piv[std::make_pair(i, max_piv)] = max_dist;

			}

			// here we save the distances
			//tree[i].rn_min_dist_to_piv = min_dist;
			//tree[i].rn_minimum_piv = min_piv;

			//tree[i].rn_max_dist_to_piv = max_dist;
			//tree[i].rn_maximum_piv = max_piv;


		}
		else { // if it is an intermediate node

			double min_dist = INT_MAX;
			int min_piv = -1;

			double max_dist = -INT_MAX;
			int max_piv = -1;

			for (int r_piv = 0; r_piv > No_RN_piv; r_piv++) { // for all RN pivots

				for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); ++it) { // for all node children in the node

					if (min_dist > tree[*it].rn_min_dist_to_piv[std::make_pair(*it, RN_piv_set[r_piv])]) { // grt the min
						min_dist = tree[*it].rn_min_dist_to_piv[std::make_pair(*it, RN_piv_set[r_piv])];
						//min_piv = tree[i].rn_minimum_piv;
					}

					if (max_dist < tree[*it].rn_max_dist_to_piv[std::make_pair(*it, RN_piv_set[r_piv])]) { // grt the min
						max_dist = tree[*it].rn_max_dist_to_piv[std::make_pair(*it, RN_piv_set[r_piv])];
						//min_piv = tree[i].rn_minimum_piv;
					}
				}
				tree[i].rn_min_dist_to_piv[std::make_pair(i, RN_piv_set[r_piv])] = min_dist;
				tree[i].rn_max_dist_to_piv[std::make_pair(i, RN_piv_set[r_piv])] = max_dist;
			}

			// here we save the distances
			//tree[i].rn_min_dist_to_piv[std::make_pair(RN_piv_set[r_piv], *it)]


			//tree[i].rn_min_dist_to_piv = min_dist;
			//tree[i].rn_minimum_piv = min_piv;

			//tree[i].rn_max_dist_to_piv = max_dist;
			//tree[i].rn_maximum_piv = max_piv;

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


///////////////////////////////////////////////////////
/*
REQUIRES :: a query vertex, v
ENSURES  ::
*/
int getTreeHight() {
	int height = 0;
	int ptr = INT_MAX;
	int i = INDEXSIZE - 1;

	while (i >= 0) {
		height++;
		i = tree[i].parent;
	}

	height = height + 1; // this is to add the root level

	return height;
}


///////////////////////////////////////////////////////
/*
REQUIRES :: a query vertex, v
ENSURES  ::
*/
std::unordered_map<pair, int, pair_hash> get_queryNode(int q) {
	int git;
	int treeHeight = getTreeHight();
	int* queryNodeIndex = new int[treeHeight];

	// first, we find the query node at leaf level (level == 0)
	bool found = false;
	for (int i = INDEXSIZE - 1; i >= 0; ++i) {

		for (std::unordered_set<int>::iterator it = tree[i].child.begin(); it != tree[i].child.end(); it++) {
			if (q == *it) {
				git = i;
				found = true;
				break;
			}
		}
		if (found)
			break;
	}
	// we build an array of the query vertex
	std::unordered_map<pair, int, pair_hash> queryNodeLevel;
	queryNodeIndex[treeHeight - 1] = git;
	queryNodeLevel[std::make_pair(q, 0)] = git;


	for (int i = treeHeight - 2; i >= 0; --i) {
		queryNodeIndex[i] = tree[queryNodeIndex[i - 1]].parent;
		queryNodeLevel[std::make_pair(q, i)] = tree[queryNodeLevel[std::make_pair(q, i - 1)]].parent;
	}

	return queryNodeLevel;
}




/*
	GIVEN		:: social network
	ENSURES		:: in_inf, and out_info for each social network vrtx
*/
void computeIn_inf_Out_inf_for_vertices() {

	for (int v = 0; v < No_sn_V; v++) {

		for (std::unordered_set<int>::iterator it = sn_vrtx[v].out_nbrs.begin(); it != sn_vrtx[v].out_nbrs.end(); it++) {

			for (int i = 0; i < No_of_TOPICS; ++i) {

				if (sn_vrtx[v].out_inf[i] < sn_edge_info[std::make_pair(v, *it)].topics[i])
					sn_vrtx[v].out_inf[i] = sn_edge_info[std::make_pair(v, *it)].topics[i];

				if (sn_vrtx[v].in_inf[i] < sn_edge_info[std::make_pair(*it, v)].topics[i])
					sn_vrtx[v].in_inf[i] = sn_edge_info[std::make_pair(*it, v)].topics[i];

			}

		}

	}
}



#endif // !INDEXTRAVERSE
