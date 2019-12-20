#ifndef INDEX_TRAV
#define INDEX_TRAV
#include "parameterSettings.h"
#include "vertex.h"
#include"heap.h"
#include"../buildIndex.h"
#include "IndexSummurizing.h"


/*
	This pice of the code will traverse the index to get the solution
*/
bool  Index_InfluenceScorePruning(int q, int theNode, int Qtopic[], int SizeQtopic);
double ub_SpatialDistance(int q, int candV);
bool SocialDistancePruning(int q, int candV);
double ub_SocialDistance(int q, int candV);
bool Index_SpatialDistancePruning(int q, int theNode);
bool  Index_InfluenceScorePruning(int q, int theNode, int Qtopic[], int SizeQtopic);
double Index_lb_dist_RN(int qVrtx, int theNode);
bool SpatialDistancePruning(int q, int candV);
bool SpatialDistancePruning(int q, int candV);
double Index_lb_infScore_in(int q, int theNode, int Qtopic[], int SizeQtopic);
double Index_lb_infScore_out(int q, int theNode, int Qtopic[], int SizeQtopic);



void indexTrav(int q) {

	// here we set the query topics set, they are just an index 0, 1, 2, ..
	int SizeQtopic = 2;
	int* Qtopic = new int[SizeQtopic];
	Qtopic[0] = 0;
	Qtopic[1] = 1;


	Heap* hp = new Heap();
	hp->init(2);
	
	std::unordered_set<int> S;

	int treeHeight = getTreeHight();
	int* queryNode = new int[treeHeight];
	
	std::unordered_map<pair, int, pair_hash> queryNodeLevel;
	queryNodeLevel = get_queryNode(q);

	int root = 0;

	HeapEntry* he = new HeapEntry();

	he->son1 = root;
	he->key = 0;
	hp->insert(he);
	delete he;

	Heap* hp_p = new Heap();
	hp_p->init(2);
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

				if (!(SpatialDistancePruning(q, *it) && SocialDistancePruning(q, *it))) { // if the object cannot be pruned

					S.insert(*it);

				}

			}

		}
		else { // non-leaf node

			// optain the tree node that contains the query vertex
			int queryNode = queryNodeLevel[std::make_pair(q, tree[theNode].level)];
			
			for (std::unordered_set<int>::iterator it = tree[theNode].child.begin(); it != tree[theNode].child.end(); ++it) {

				if (!Index_SpatialDistancePruning(q, *it) && !Index_InfluenceScorePruning(q, *it, Qtopic, SizeQtopic)) {

					HeapEntry* he = new HeapEntry();

					he->son1 = *it;
					he->key = Index_lb_dist_RN(q, *it);

					hp->insert(he);
					delete he;
				}
			}

		}

	}

	//S = Refine(S);

	//return S;
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

	double dist = ub_SpatialDistance(q, candV);

	if (dist > SIGMA)
		return true;
	else
		return false;

}

/*
	GIVEN	:: two social network vertices
	ENSURES	:: the lower bound of the spatial road network distance between them
*/
double ub_SpatialDistance(int q, int candV) {
	double dist = 0.0;
	double min = INT_MAX;
	for (int i = 0; i < No_RN_piv; i++) {

		dist = sn_vrtx[q].rn_distToPiv[std::make_pair(q, RN_piv_set[i])] + sn_vrtx[candV].rn_distToPiv[std::make_pair(candV, RN_piv_set[i])];

		if (min > dist)
			min = dist;

	}

	return dist;
}

/*
	object level spatial social pruning
	GIVEN	:: q, and a candidate vertex candV
	ENSURES :: true is it can be pruned, false, otherwise
*/
bool SocialDistancePruning(int q, int candV) {

	double dist = ub_SocialDistance(q, candV);

	if (dist > No_Hops)
		return true;
	else
		return false;

}

/*
	GIVEN	:: two social network vertices
	ENSURES	:: the lower bound of the social network distance between them
*/
double ub_SocialDistance(int q, int candV) {
	double dist = 0.0;
	double min = INT_MAX;
	
	for (int i = 0; i < No_SN_piv; i++) {

		dist = sn_vrtx[q].sn_distToPiv[std::make_pair(q, SN_piv_set[i])] + sn_vrtx[candV].sn_distToPiv[std::make_pair(candV, SN_piv_set[i])];

		if (min > dist)
			min = dist;

	}

	return dist;
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


	if ( (lb_inf_in > THETA) && (lb_inf_out > THETA) )
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
		double dst;
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

		}
		 
		if (git_max < dst) {
			git_max = dst;
		}

	}

	return git_max;

}





#endif // !INDEX_TRAV

