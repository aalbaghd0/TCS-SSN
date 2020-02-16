#ifndef INDEX_TRAV
#define INDEX_TRAV
#include "parameterSettings.h"
#include "vertex.h"
#include"heap.h"
#include"../buildIndex.h"
#include "IndexSummurizing.h"
#include <ctime>

/*
	This pice of the code will traverse the index to get the solution
*/
bool  Index_InfluenceScorePruning(int q, int theNode, int Qtopic[], int SizeQtopic);
double lb_SpatialDistance(int q, int candV);
bool SocialDistancePruning(int q, int candV);
double lb_SocialDistance(int q, int candV);
bool Index_SpatialDistancePruning(int q, int theNode);
bool  Index_InfluenceScorePruning(int q, int theNode, int Qtopic[], int SizeQtopic);
double Index_lb_dist_SN(int qVrtx, int theNode);
bool SocialDistancePruning(int q, int candV);
double Index_lb_infScore_in(int q, int theNode, int Qtopic[], int SizeQtopic);
double Index_lb_infScore_out(int q, int theNode, int Qtopic[], int SizeQtopic);
double Index_lb_dist_RN(int qVrtx, int theNode);
double lb_SocialDistance(int q, int candV);
bool SpatialDistancePruning(int q, int candV);
bool Index_SocialDistancePruning(int q, int theNode);
bool IndexStructuralCohesivenessPruning(int q, int theNode);
bool IndexKeywordbasedPruning (int q, int it, std::bitset<No_K> tt);
bool KeywordbasedPruning(int v, std::bitset<No_K> k_set);
bool StructuralCohesivenessPruning(int v);
bool InfluenceScorePruning(int q, int v, int Qtopic[], int SizeQtopic);



void indexTrav(int q, std::bitset<No_K> k_set) {
	// here we set the query topics set, they are just an index 0, 1, 2, ..
	int SizeQtopic = 2;
	int* Qtopic = new int[SizeQtopic];
	Qtopic[0] = 0;
	Qtopic[1] = 1;


	// set of counters
	int no_nodes_not_pruned = 0;
	int no_nodes_all = 0;

	int no_obj = 0;



	Heap* hp = new Heap();
	hp->init(2);

	std::unordered_set<int> S;

	//int treeHeight = getTreeHight();
	//int* queryNode = new int[treeHeight];

	std::unordered_map<pair, int, pair_hash> queryNodeLevel;
	queryNodeLevel = get_queryNode(q);

	int root = 0;

	HeapEntry* he = new HeapEntry();

	he->son1 = root;
	he->key = 0;
	hp->insert(he);
	delete he;
	int no_nodes_notPruned = 0;
	Heap* hp_p = new Heap();
	hp_p->init(2);


	//start recording the time
	std::clock_t c_start = std::clock();

	while (hp->used > 0) {
		HeapEntry* he = new HeapEntry();
		hp->remove(he);

		int theNode = he->son1;
		int key = he->key;
		delete he;
		// if key greater than threshold, treminate the loop
		if (key > SIGMA)
			break;

		if (tree[theNode].level == 0) { // if a leaf node

			for (std::unordered_set<int>::iterator it = tree[theNode].child.begin(); it != tree[theNode].child.end(); ++it) {
				no_obj++;
				if (!InfluenceScorePruning(q, *it, Qtopic, SizeQtopic) && !StructuralCohesivenessPruning(*it) && !KeywordbasedPruning(*it, k_set) && !SocialDistancePruning(q, *it)
					&& !SpatialDistancePruning(q, *it)) { // if the object cannot be pruned
					no_nodes_notPruned++;
					S.insert(*it);

				}

			}

		}
		else { // non-leaf node

			// optain the tree node that contains the query vertex
			//int queryNode = queryNodeLevel[std::make_pair(q, tree[theNode].level)];

			for (std::unordered_set<int>::iterator it = tree[theNode].child.begin(); it != tree[theNode].child.end(); ++it) {
				no_nodes_all++;
				if (!Index_InfluenceScorePruning(q, *it, Qtopic, SizeQtopic) && !IndexKeywordbasedPruning(q, *it, k_set) 
					&& !IndexStructuralCohesivenessPruning(q, *it) && !Index_SocialDistancePruning(q, *it) 
					&& !Index_SpatialDistancePruning(q, *it)) {

					HeapEntry* he = new HeapEntry();

					he->son1 = *it;
					he->key = Index_lb_dist_RN(q, *it);
					no_nodes_not_pruned ++;
					hp->insert(he);
					delete he;

				}
			}
		}



		//S = Refine(S);

		//return S;
	} 
	std::clock_t c_end = std::clock();
	long double time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
	std::cout << "CPU time used: " << time_elapsed_ms / 1000.0 << " s\n";

	std::cerr << "\n the size of ths candidate set"<< S.size()<<"\n";
	std::cerr << " \n # all nodes = " << no_nodes_all << " \n # all nodes = " << no_nodes_not_pruned << " \n";

	std::cerr << "nomber of objects survived the object pruning: " << no_obj << "\n";
}


/////////////////////////////////////////

//			OBJECT LEVEL PRUNING

////////////////////////////////////////

/*
	object level spatial distance pruning
	GIVEN	:: q, and a candidate vertex candV
	ENSURES :: true is it can be pruned, false, otherwise
*/
bool SpatialDistancePruning(int q, int candV) {

	double dist = lb_SpatialDistance(q, candV);

	if (dist > SIGMA)
		return true;
	else
		return false;

}

/*
	GIVEN	:: two social network vertices
	ENSURES	:: the lower bound of the spatial road network distance between them
*/
double lb_SpatialDistance(int q, int candV) {
	double dist = 0.0;
	double max = 0;
	for (int i = 0; i < No_RN_piv; i++) {

		dist = abs(sn_vrtx[q].rn_distToPiv[std::make_pair(q, RN_piv_set[i])] - sn_vrtx[candV].rn_distToPiv[std::make_pair(candV, RN_piv_set[i])]);

		if (max < dist)
			max = dist;

	}

	return max;
}

/*
	object level spatial social pruning
	GIVEN	:: q, and a candidate vertex candV
	ENSURES :: true is it can be pruned, false, otherwise
*/
bool SocialDistancePruning(int q, int candV) {

	double dist = lb_SocialDistance(q, candV);

	if (dist > No_Hops)
		return true;
	else
		return false;

}

/*
	GIVEN	:: two social network vertices
	ENSURES	:: the lower bound of the social network distance between them
*/
double lb_SocialDistance(int q, int candV) {
	double lb_dist = 0.0;
	double max = - INT_MAX;
	
	for (int i = 0; i < No_SN_piv; i++) {

		lb_dist = abs(sn_vrtx[q].sn_distToPiv[std::make_pair(q, SN_piv_set[i])] - 
					sn_vrtx[candV].sn_distToPiv[std::make_pair(candV, SN_piv_set[i])]);

		if (max < lb_dist)
			max = lb_dist;

	}

	return max;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/*
	::::::  KEYWORD PRUNING  ::::::
	GIVEN	:: social network vrtx and a query set of keywords
	ENSURES :: there exists at least one keyword (element) incommon between the candidate vertex and the query keyword set

*/
bool KeywordbasedPruning(int v, std::bitset<No_K> k_set) {


	k_set &= sn_vrtx[v].key;

	if (k_set.count() > 0)
		return false;
	else
		return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////

/*
	::::::  Stractural Pruning  ::::::
	GIVEN	:: social network vrtx
	ENSURES :: there exists at least one keyword (element) incommon between the candidate vertex and the query keyword set

*/
bool StructuralCohesivenessPruning(int v) {

	if (sn_vrtx[v].truss < Ktruss - 2)
		return true;
	else
		return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////////

/*
	::::::  Influence Score Pruning  ::::::
	GIVEN	:: social network vrtx, a query vrtx, and a set of topics
	ENSURES :: there exists at least one keyword (element) incommon between the candidate vertex and the query keyword set

*/

double lb_infScore_in(int q, int v, int Qtopic[], int SizeQtopic) {
	
	double lb_inf_in = 0.0;

	for (int i = 0; i < No_of_TOPICS; ++i) {

		if (isInTheArray(Qtopic, SizeQtopic, i)) {

			lb_inf_in = lb_inf_in + sn_vrtx[q].out_inf[i] * sn_vrtx[v].in_inf[i];

		}
	}
	return lb_inf_in;
}

double lb_infScore_out(int q, int v, int Qtopic[], int SizeQtopic) {

	double lb_inf_in = 0.0;

	for (int i = 0; i < No_of_TOPICS; ++i) {

		if (isInTheArray(Qtopic, SizeQtopic, i)) {

			lb_inf_in = lb_inf_in + sn_vrtx[q].in_inf[i] * sn_vrtx[v].out_inf[i];

		}
	}
	return lb_inf_in;
}

bool InfluenceScorePruning(int q, int v, int Qtopic[], int SizeQtopic) {

	//double lb_inf_in = lb_infScore_in(q, v, Qtopic, SizeQtopic);
	double lb_inf_out = lb_infScore_out(q, v, Qtopic, SizeQtopic);


	//if ((lb_inf_in < THETA) && (lb_inf_out < THETA))
	if (lb_inf_out < THETA)
		return true;
	else
		return false;

}
////////////////////////////////////////////////////////


//				INDEX LEVEL PRUNING

////////////////////////////////////////////////////////
/*
	lemma 7, Index spatial distance pruning
*/
bool Index_SpatialDistancePruning(int q, int theNode) {

	double lb_dist = Index_lb_dist_RN(q, theNode);

	if (lb_dist > SIGMA)
		return true;
	else
		return false;

}

/*
	lemma 8,  Influence  Score  Pruning 
*/
bool  Index_InfluenceScorePruning(int q, int theNode, int Qtopic[], int SizeQtopic) {

	double lb_inf_in = Index_lb_infScore_in(q, theNode, Qtopic, SizeQtopic);
	double lb_inf_out = Index_lb_infScore_out(q, theNode, Qtopic, SizeQtopic);


	if ( (lb_inf_in < THETA) && (lb_inf_out < THETA) )
		return true;
	else
		return false;
}
double Index_lb_infScore_in(int q, int theNode, int Qtopic[], int SizeQtopic) {

	double lb_inf_in = 0.0;

	for (int i = 0; i < No_of_TOPICS; ++i) {
		
		if (isInTheArray(Qtopic, SizeQtopic, i)) {

			lb_inf_in = lb_inf_in + sn_vrtx[q].in_inf[i] * tree[theNode].out_topics_prob[i];

		}
	}
	return lb_inf_in;
}

double Index_lb_infScore_out(int q, int theNode, int Qtopic[], int SizeQtopic) {

	double lb_inf_out = 0.0;

	for (int i = 0; i < No_of_TOPICS; ++i) {

		if (isInTheArray(Qtopic, SizeQtopic, i)) {

			lb_inf_out = lb_inf_out + sn_vrtx[q].out_inf[i] * tree[theNode].in_topics_prob[i];

		}
	}

	return lb_inf_out;
}


double Index_lb_dist_RN(int qVrtx, int theNode) {
	double git_max = - INT_MAX;

	for (int r_piv = 0; r_piv < No_RN_piv; ++ r_piv) {
		double dst = 0;
		/*
		std::cerr << "sn_vrtx[qVrtx].rn_distToPiv: "<< sn_vrtx[qVrtx].rn_distToPiv[std::make_pair(qVrtx, RN_piv_set[r_piv])]
			<< "\n tree[theNode].rn_min_dist_to_piv"<< tree[theNode].rn_min_dist_to_piv[std::make_pair(theNode, RN_piv_set[r_piv])]<<
			"\n sn_vrtx[qVrtx].rn_distToPiv - tree[theNode].rn_min_dist_to_piv: " << abs(sn_vrtx[qVrtx].rn_distToPiv[std::make_pair(qVrtx, RN_piv_set[r_piv])]
			- tree[theNode].rn_min_dist_to_piv[std::make_pair(theNode, RN_piv_set[r_piv])]) << "\nl -------------------- \n";
*/
		///*
		if (sn_vrtx[qVrtx].rn_distToPiv[std::make_pair(qVrtx, RN_piv_set[r_piv])]
							< tree[theNode].rn_min_dist_to_piv[std::make_pair(theNode, RN_piv_set[r_piv])]) {
			
			dst = abs(sn_vrtx[qVrtx].rn_distToPiv[std::make_pair(qVrtx, RN_piv_set[r_piv])]
							- tree[theNode].rn_min_dist_to_piv[std::make_pair(theNode, RN_piv_set[r_piv])]);
		
		}
		else if(sn_vrtx[qVrtx].rn_distToPiv[std::make_pair(qVrtx, RN_piv_set[r_piv])]
									> tree[theNode].rn_max_dist_to_piv[std::make_pair(theNode, RN_piv_set[r_piv])]){

			 dst = abs(sn_vrtx[qVrtx].rn_distToPiv[std::make_pair(qVrtx, RN_piv_set[r_piv])]
							- tree[theNode].rn_max_dist_to_piv[std::make_pair(theNode, RN_piv_set[r_piv])]);

		}else{

			dst = 0;
			//std::cerr << "\n sn_vrtx[qVrtx].rn_distToPiv - tree[theNode].rn_min_dist_to_piv" << abs(sn_vrtx[qVrtx].rn_distToPiv[std::make_pair(qVrtx, RN_piv_set[r_piv])]
				//- tree[theNode].rn_min_dist_to_piv[std::make_pair(theNode, RN_piv_set[r_piv])]) << "\nl";
			//dst = abs(sn_vrtx[qVrtx].rn_distToPiv[std::make_pair(qVrtx, RN_piv_set[r_piv])]
				//- tree[theNode].rn_min_dist_to_piv[std::make_pair(theNode, RN_piv_set[r_piv])]);
		}
		
		
		//*/

		//dst = abs(sn_vrtx[qVrtx].rn_distToPiv[std::make_pair(qVrtx, RN_piv_set[r_piv])]
			///- tree[theNode].rn_min_dist_to_piv[std::make_pair(theNode, RN_piv_set[r_piv])]);
		if (git_max < dst) {
			git_max = dst;
		}
		

	}
	//std::cerr << "\n"<< git_max;
	
	return git_max;

}


/*

	Lemma 9, Social Distance Based Prining -- INDEX LEVEL


*/

bool Index_SocialDistancePruning(int q, int theNode) {

	double lb_dist = Index_lb_dist_SN(q, theNode);

	if (lb_dist > No_Hops)
		return true;
	else
		return false;

}

double Index_lb_dist_SN(int qVrtx, int theNode) {

	double git_max = -INT_MAX;

	for (int s_piv = 0; s_piv < No_SN_piv; ++s_piv) {
		double dst;
		if (sn_vrtx[qVrtx].sn_distToPiv[std::make_pair(qVrtx, SN_piv_set[s_piv])]
			< tree[theNode].sn_min_dist_to_piv[std::make_pair(theNode, SN_piv_set[s_piv])]) {

			dst = abs(sn_vrtx[qVrtx].sn_distToPiv[std::make_pair(qVrtx, SN_piv_set[s_piv])]
				- tree[theNode].sn_min_dist_to_piv[std::make_pair(theNode, SN_piv_set[s_piv])]);

		}
		else if (sn_vrtx[qVrtx].sn_distToPiv[std::make_pair(qVrtx, SN_piv_set[s_piv])]
			> tree[theNode].sn_max_dist_to_piv[std::make_pair(theNode, SN_piv_set[s_piv])]) {
			
			//std::cerr <<"max dist " <<tree[theNode].sn_max_dist_to_piv[std::make_pair(theNode, SN_piv_set[s_piv])];
			dst = abs(sn_vrtx[qVrtx].sn_distToPiv[std::make_pair(qVrtx, SN_piv_set[s_piv])]
				- tree[theNode].sn_max_dist_to_piv[std::make_pair(theNode, SN_piv_set[s_piv])]);

		}
		else {
			//std::cerr << "max dist " << tree[theNode].sn_max_dist_to_piv[std::make_pair(theNode, SN_piv_set[s_piv])];

			dst = 0;

		}

		if (git_max < dst) {
			git_max = dst;
		}

	}

	return git_max;
}


/*

  Lemma 10 Index Structural Cohesiveness Pruning

*/

bool IndexStructuralCohesivenessPruning(int q, int theNode) {
	
	double lb_w = tree[theNode].truss;

	if (lb_w < Ktruss)
		return true;
	else
		return false;
}


/*
	
	Lemma 1-- Keyword based pruning   -- INDEX LEVEL

*/

bool IndexKeywordbasedPruning(int q, int theNode, std::bitset<No_K> k_set) {
	
	k_set &= tree[theNode].keys;

	if (k_set.count() > 0)
		return false;
	else
		return true;
}
#endif // !INDEX_TRAV

