#ifndef INDEX_TRAV
#define INDEX_TRAV
#include "parameterSettings.h"
#include "vertex.h"
#include"heap.h"
#include"../buildIndex.h"
#include "IndexSummurizing.h"
#include <ctime>
#include"../buildIndex.h"

/*
	This pice of the code will traverse the index to get the solution
*/

double lb_SpatialDistance(int q, int candV);
bool SocialDistancePruning(int q, int candV, int noOfHops);
double lb_SocialDistance(int q, int candV);
bool Index_SpatialDistancePruning(int q, int theNode, int sigma);
bool  Index_InfluenceScorePruning(int q, int theNode, int Qtopic[], int SizeQtopic, int theta);
double Index_lb_dist_SN(int qVrtx, int theNode);
bool SocialDistancePruning(int q, int candV);
double Index_lb_infScore_in(int q, int theNode, int Qtopic[], int SizeQtopic);
double Index_lb_infScore_out(int q, int theNode, int Qtopic[], int SizeQtopic);
double Index_lb_dist_RN(int qVrtx, int theNode);
double lb_SocialDistance(int q, int candV);
bool SpatialDistancePruning(int q, int candV, int sigma);
bool Index_SocialDistancePruning(int q, int theNode, int noOfHops);
bool IndexStructuralCohesivenessPruning(int q, int theNode, int truss_val);
bool IndexKeywordbasedPruning (int q, int it, std::bitset<No_K> tt);
bool KeywordbasedPruning(int v, std::bitset<No_K> k_set);
bool StructuralCohesivenessPruning(int v, int truss_val);
bool InfluenceScorePruning(int q, int v, int Qtopic[], int SizeQtopic, int theta);
void Refine(int q, int NoOfHops, double sigma, int truss_val,std::unordered_set<int> S, long double time_elapsed_s);
std::unordered_set<int> GreedyBaseLine(int q, int noOfHops, double sigma, int truss_val, double theta, int Qtopic[], std::bitset<No_K> k_set, int SizeQtopic, long double &time);
int intersect(int a, int b);



void indexTrav() {
	// here we set the query topics set, they are just an index 0, 1, 2, ..
	int max = 0;
	for (int i = 0; i < No_sn_V; i++) {
		if (max < snGraph[i].size())
			max = i;
	}
	std::cerr << "the vertex with maximun degree: " << max << std::endl;
	
	while (TRUE) {
		int noOfHops;
		double sigma;
		double theta;
		int SizeQtopic = 2;
		int* Qtopic = new int[SizeQtopic];
		int length_querySet;
		int key_val;
		std::bitset<No_K> k_set;
		int q;
		int truss_val = 0;
		std::cerr << "input query vertex (q): ";
		std::cin >> q;
		std::cerr << "input number of social network distace (no of hops): ";
		std::cin >> noOfHops;
		std::cerr << "input the truss value (k-truss): ";
		std::cin >> truss_val;
		std::cerr << "input road network distance threshold (sigma): ";
		std::cin >> sigma;
		std::cerr << "input the influence threshold (theta): ";
		std::cin >> theta;
		std::cerr << "input query topics 2 vales between 1 and 0: ";
		std::cin >> Qtopic[0];
		std::cerr << "input again vale between 1 and 0: ";
		std::cin >> Qtopic[1];
		std::cerr << "the length of the query set: ";
		std::cin >> length_querySet;
		for (int i = 0; i < length_querySet; i++) {
			std::cerr << "pick a random value between 0 and 9: ";
			std::cin >> key_val;
			k_set.set(key_val);
		}


		// set of counters
		int no_nodes_not_pruned = 0;
		int no_nodes_all = 0;

		int no_obj = 0;



		Heap* hp = new Heap();
		hp->init(2);

		std::unordered_set<int> S;

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
					if (!InfluenceScorePruning(q, *it, Qtopic, SizeQtopic, theta) && !StructuralCohesivenessPruning(*it, truss_val) && !KeywordbasedPruning(*it, k_set)
						&& !SocialDistancePruning(q, *it, noOfHops) && !SpatialDistancePruning(q, *it, sigma)) { // if the object cannot be pruned
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
					if (!Index_InfluenceScorePruning(q, *it, Qtopic, SizeQtopic, theta) && !IndexKeywordbasedPruning(q, *it, k_set)
						&& !IndexStructuralCohesivenessPruning(q, *it, truss_val) && !Index_SocialDistancePruning(q, *it, noOfHops)
						&& !Index_SpatialDistancePruning(q, *it, sigma)) {

						HeapEntry* he = new HeapEntry();

						he->son1 = *it;
						he->key = Index_lb_dist_RN(q, *it);
						no_nodes_not_pruned++;
						hp->insert(he);
						delete he;

					}
				}
			}



			//S = Refine(S);

			//return S;
		}
		std::clock_t c_end = std::clock();
		long double time_elapsed_s = (c_end - c_start) / 1000.0;
		std::cout << "CPU time used in Tree Traversal: " << time_elapsed_s << " s\n";

		std::cerr << "\n the size of ths candidate set" << S.size() << "\n";
		std::cerr << " \n # all nodes = " << no_nodes_all << " \n # all nodes = " << no_nodes_not_pruned << " \n";

		std::cerr << "nomber of objects survived the object pruning: " << no_obj << "\n";

		std::cerr << "\n---------------------------------------------------------------\n";

		//S.clear();
		long double time;
		
		Refine(q, noOfHops, sigma, truss_val, S, time_elapsed_s);
		GreedyBaseLine(q, noOfHops, sigma, truss_val, theta, Qtopic, k_set, SizeQtopic, time);

		S.clear();
		hp->~Heap();
		delete[] Qtopic;
		check_hash_rn_dist.clear();

	}
}


/////////////////////////////////////////

//			OBJECT LEVEL PRUNING

////////////////////////////////////////

/*
	object level spatial distance pruning
	GIVEN	:: q, and a candidate vertex candV
	ENSURES :: true is it can be pruned, false, otherwise
*/
bool SpatialDistancePruning(int q, int candV, int sigma) {

	double dist = lb_SpatialDistance(q, candV);

	if (dist > sigma)
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
bool SocialDistancePruning(int q, int candV, int noOfHops) {

	double dist = lb_SocialDistance(q, candV);

	if (dist > noOfHops)
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
bool StructuralCohesivenessPruning(int v, int truss_val) {

	if (sn_vrtx[v].truss < truss_val - 2)
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

bool InfluenceScorePruning(int q, int v, int Qtopic[], int SizeQtopic, int theta) {

	//double lb_inf_in = lb_infScore_in(q, v, Qtopic, SizeQtopic);
	double lb_inf_out = lb_infScore_out(q, v, Qtopic, SizeQtopic);


	//if ((lb_inf_in < THETA) && (lb_inf_out < THETA))
	if (lb_inf_out < theta)
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
bool Index_SpatialDistancePruning(int q, int theNode, int sigma) {

	double lb_dist = Index_lb_dist_RN(q, theNode);

	if (lb_dist > sigma)
		return true;
	else
		return false;

}

/*
	lemma 8,  Influence  Score  Pruning 
*/
bool  Index_InfluenceScorePruning(int q, int theNode, int Qtopic[], int SizeQtopic, int theta) {

	double lb_inf_in = Index_lb_infScore_in(q, theNode, Qtopic, SizeQtopic);
	double lb_inf_out = Index_lb_infScore_out(q, theNode, Qtopic, SizeQtopic);


	if ( (lb_inf_in < theta) && (lb_inf_out < theta) )
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

bool Index_SocialDistancePruning(int q, int theNode, int noOfHops) {

	double lb_dist = Index_lb_dist_SN(q, theNode);

	if (lb_dist > noOfHops)
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

bool IndexStructuralCohesivenessPruning(int q, int theNode, int truss_val) {
	
	double lb_w = tree[theNode].truss;

	if (lb_w < truss_val)
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


void Refine(int q, int noOfHops, double sigma, int truss_val, std::unordered_set<int> S, long double treeTime) {
	// get the maximum subgraph using BFS

	bool *pruned = new bool [No_sn_V];
	bool *visited = new bool [No_sn_V];

	for (int i = 0; i < No_sn_V; ++i) {
		pruned[i] = true;
		visited[i] = false;
	}

	for (std::unordered_set<int>::iterator it = S.begin(); it != S.end(); ++it) {
		pruned[*it] = false;
	}

	
	
	// we clean the social network graph by removing pruned vertices
	for (int i = 0; i < No_sn_V; ++i) {
		if (pruned[i] == false) {
			for (std::list<int>::iterator it = snGraph[i].begin(); it != snGraph[i].end(); it++) {
				if (pruned[*it] == false) {
					if (!f_edges[std::make_pair(i, *it)] && !f_edges[std::make_pair(*it, i)]) {

					snGraph_fnl[*it].push_back(i);
					snGraph_fnl[i].push_back(*it);
					f_edges[std::make_pair(i, *it)] = true;
					f_edges[std::make_pair(*it, i)] = true;

					}
				}
			}
		}
	}

	std::unordered_set<int> finalSet;

	int cand, cand2;
	double weight;
	
	int* dist = new int[No_sn_V];
	
	int pos = sn_vrtx[q].ckins[0];
	//start recording the time
	
	int test = 0;

	int firstTime = 0;
	// do BFS to get the road network distance starting from the query vertex
	long double rnTime = Refine_rn_Dij_to_all_vertices(pos, sigma);

	std::clock_t c_del_1_s, c_del_1_e;
	long double del1 = 0.0, del2 = 0.0, del3 = 0.0, del4 = 0.0, del5 = 0.0, del6 = 0.0;

	int counter = 0;
	std::clock_t c_start = std::clock();

	while (test != -1) {
		counter++;
		c_del_1_s = std::clock();
		test = -1;
		std::memset(dist, 0, sizeof(dist[0]) * No_sn_V);
		Heap* hp = new Heap();
		hp->init(2);
		HeapEntry* he = new HeapEntry();
		he->son1 = q;
		he->key = 0;
		hp->insert(he);
		delete he;
		finalSet.clear();
		c_del_1_e = std::clock();
		del4 = del4 + ((c_del_1_e - c_del_1_s));

		finalSet.insert(q);
		
		while (hp->used > 0) {

			c_del_1_s = std::clock();
			HeapEntry* he = new HeapEntry();
			hp->remove(he);
			cand = he->son1;
			weight = he->key;
			delete he;
			c_del_1_e = std::clock();
			del3 = del3 + ((c_del_1_e - c_del_1_s));

			for (std::list<int>::iterator it = snGraph_fnl[cand].begin(); it != snGraph_fnl[cand].end(); ++it) {
				
				if (!visited[*it]) {
					visited[*it] = true;

					c_del_1_s = std::clock();
					int posDist = sn_vrtx[*it].ckins[0];
					HeapEntry* he = new HeapEntry();
					//std::cerr << "rn_dist_for_users(q, *it) < sigma: " << rn_dist_for_users(q, *it) << std::endl;
					c_del_1_e = std::clock();
					del5 = del5 + ((c_del_1_e - c_del_1_s));
					
					he->son1 = *it;
					he->key = weight + 1;

					
		
					// social distance pruning
					if (noOfHops >= he->key
						&& //rn_dist_for_users(q, *it) < sigma
						check_hash_rn_dist[std::make_pair(pos, posDist)] && hash_rn_dist[std::make_pair(pos, posDist)] < sigma
						) {
						
						hp->insert(he);
						finalSet.insert(*it);
					}

					c_del_1_s = std::clock();
					if (firstTime == 0) {
						std::cerr << "\n**********************************************************\n";
						std::cerr << "The size of candidate set after indexing and one BFS:: " << finalSet.size() << std::endl;
					}
					c_del_1_e = std::clock();
					del1 = del1 + ((c_del_1_e - c_del_1_s));
					if(intersect(cand, *it) < truss_val -2 ){
						c_del_1_s = std::clock();
						snGraph_fnl[*it].remove(cand);
						snGraph_fnl[cand].erase(it);
						test = 1;
						c_del_1_e = std::clock();
						del1 = del1 + ((c_del_1_e - c_del_1_s));
					}
					delete he;
				}
			}
		}
		c_del_1_s = std::clock();
		hp->~Heap();
		c_del_1_e = std::clock();
		del6 = del6 + ((c_del_1_e - c_del_1_s));
	}

	std::clock_t c_end = std::clock();
	long double time_elapsed_ms = (c_end - c_start) / 1000.0;

	
	std::cout << "CPU time used for Refinment " << time_elapsed_ms << " s\n";
	std::cerr << "the Total CPU Time: "  <<time_elapsed_ms + rnTime +treeTime - (del1/ 1000.0) - (del2/1000.0) - (del3 / 1000.0) - (del4 / 1000.0) - (del5 / 1000.0) - (del6 / 1000.0) << "s\n";
	std::cerr << "the final set size:: " << finalSet.size() << std::endl;
	std::cerr << "No Of times we refine " << counter << std::endl;

	std::cerr << "\n**********************************************************\n";
	 
	
	delete[] visited;
	delete[] pruned;
	f_edges.clear();
	
}




int intersect(int a, int b) {
	//std::cerr << "the edge truss: " << snEdges[ hash_edge[std::make_pair(a, b)] ].sup << std::endl;
	double rslt = 0.0;
	std::unordered_set<int> intersect;
	set_intersection(sn_vrtx[a].nbrs_set.begin(), sn_vrtx[a].nbrs_set.end(),
		sn_vrtx[b].nbrs_set.begin(), sn_vrtx[b].nbrs_set.end(),
		std::inserter(intersect, intersect.begin()));

	rslt = intersect.size();
	return rslt;
}

/*
///////////////////////////////////////////////////////////////////////

					THE GREEDY BASE ALGORITHM

///////////////////////////////////////////////////////////////////////
*/
bool KeywordbasedPruning_GREEDY(int v, std::bitset<No_K> kk_set) {
	bool val = true;;
	for (int i = 0; i < No_K; i++) {
		if (kk_set[i] == sn_vrtx[v].key[i])
			val = false;
	}
	return val;
}
std::unordered_set<int> GreedyBaseLine(int q, int noOfHops, double sigma, int truss_val, double theta, int Qtopic[], std::bitset<No_K> k_set, int SizeQtopic, long double &time) {

	// filter by social network
	// filter by road network
	// filter by

	// get the maximum subgraph using BFS
	std::list<int>* snGraph_copy = new std::list<int>[No_sn_V];
	for (int i = 0; i < No_sn_V; i++) {
		for (std::list<int>::iterator it = snGraph[i].begin(); it != snGraph[i].end(); it++) {
			snGraph_copy[i].push_back(*it);
		}
	}


	bool* pruned = new bool[No_sn_V];
	bool* visited = new bool[No_sn_V];

	for (int i = 0; i < No_sn_V; ++i) {
		pruned[i] = false;
		visited[i] = false;
	}

	std::unordered_set<int> finalSet;

	std::unordered_set<int> remainingEdges;
	

	int cand, cand2;
	double weight;

	int* dist = new int[No_sn_V];

	int pos = sn_vrtx[q].ckins[0];
	//start recording the time

	long double del1= 0.0;
	std::clock_t c_s, c_e;

	long double rnTime = Refine_rn_Dij_to_all_vertices(pos, sigma);
	std::clock_t c_start = std::clock();
	int test = 0;
	while (test != -1) {
		test = -1;
		Heap* hp = new Heap();
		hp->init(2);
		HeapEntry* he = new HeapEntry();
		he->son1 = q;
		he->key = 0;
		hp->insert(he);
		delete he;
		finalSet.clear();
		finalSet.insert(q);
		remainingEdges.clear();
	
		while (hp->used > 0) {

			c_s = std::clock();
			HeapEntry* he = new HeapEntry();
			hp->remove(he);

			cand = he->son1;
			weight = he->key;
			delete he;
			c_e = std::clock();
			del1 = c_e - c_s;

			for (std::list<int>::iterator it = snGraph_copy[cand].begin(); it != snGraph_copy[cand].end(); it++) {
					if (!visited[*it]) {
						visited[*it] = true;
						HeapEntry* he = new HeapEntry();
						he->son1 = *it;
						he->key = weight + 1;

						int posDist = sn_vrtx[*it].ckins[0];

						// social distance pruning
						if (noOfHops >= he->key) {

							// structural pruning
							//if (!intersect(*it, truss_val)) {
								// topic based pruning
								//if (!InfluenceScorePruning(q, *it, Qtopic, SizeQtopic, theta)) {
									// keyword based pruning
							if (!KeywordbasedPruning_GREEDY(*it, k_set)) {
								//spatial distance pruning
								if (rn_dist_for_users(q, *it) < sigma
									//check_hash_rn_dist[std::make_pair(pos, posDist)] && hash_rn_dist[std::make_pair(pos, posDist)] < sigma
									) {
									remainingEdges.insert(hash_edge[std::make_pair(cand, *it)]);
									hp->insert(he);
									finalSet.insert(*it);
								}
								//}
							//}
							}
						}
						delete he;
					}
			}
		}
		for (std::unordered_set<int>::iterator it2 = remainingEdges.begin(); it2 != remainingEdges.end(); it2++) {
			int from, to;
			from = snEdges[*it2].from;
			to = snEdges[*it2].to;
			if (intersect(from, to) < truss_val - 2) {
				snGraph_copy[to].remove(from);
				snGraph_copy[from].remove(to);
				test = 1;
			}
		}
		
		std::cerr << "finalSet.size(): " << finalSet.size() << std::endl;
		
		for (std::unordered_set<int>::iterator it2 = finalSet.begin(); it2 != finalSet.end(); it2++) {
			for (std::unordered_set<int>::iterator it3 = finalSet.begin(); it3 != finalSet.end(); it3++) {
				if (*it2 != *it3 && !pruned[*it3] && !pruned[*it2]) {
					double dist = rn_dist_for_users(*it2, *it3);
					if (dist > sigma) {
						finalSet.erase(*it3);
						pruned[*it3] = true;
					}
				}

			}
		}
			hp->~Heap();
	}

	std::clock_t c_end = std::clock();
	long double time_elapsed_ms = (c_end - c_start) / 1000.0;
	time = time_elapsed_ms + rnTime;
	std::cerr << "\n**********************************************************\n";
	std::cout << "CPU time used for GREEDY Baseline " << time - (del1 / 1000.0) << " s\n";
	std::cerr << "the final set size:: " << finalSet.size() << std::endl;
	std::cerr << "\n**********************************************************\n";

	
	delete[] visited;
	delete[] pruned;

	return finalSet;
}


#endif // !INDEX_TRAV












/*
///////////////////////////////////////////////////////////////////////

					THE GREEDY BASE ALGORITHM

///////////////////////////////////////////////////////////////////////


void GreedyBaseLine(int q, int noOfHops, double sigma, int truss_val, double theta, int Qtopic[], std::bitset<No_K> k_set, int SizeQtopic) {

	std::clock_t c_start = std::clock();
	// filter by social network
	// filter by road network
	// filter by  

	// get the maximum subgraph using BFS

	bool* pruned = new bool[No_sn_V];
	bool* visited = new bool[No_sn_V];

	for (int i = 0; i < No_sn_V; ++i) {
		pruned[i] = false;
		visited[i] = false;
	}

	std::unordered_set<int> finalSet;

	Heap* hp = new Heap();
	hp->init(2);
	HeapEntry* he = new HeapEntry();
	he->son1 = q;
	he->key = 0;
	hp->insert(he);
	delete he;

	finalSet.insert(q);

	int cand, cand2;
	double weight;

	int* dist = new int[No_sn_V];

	int pos = sn_vrtx[q].ckins[0];
	//start recording the time


	long double rnTime = Refine_rn_Dij_to_all_vertices(pos, sigma);
	while (hp->used > 0) {

		HeapEntry* he = new HeapEntry();
		hp->remove(he);

		cand = he->son1;
		weight = he->key;
		delete he;

		for (std::list<int>::iterator it = snGraph[cand].begin(); it != snGraph[cand].end(); ++it) {
			if (!pruned[*it]) {
				if (!visited[*it]) {
					visited[*it] = true;
					HeapEntry* he = new HeapEntry();
					he->son1 = *it;
					he->key = weight + 1;

					int posDist = sn_vrtx[*it].ckins[0];

					// social distance pruning
					if (noOfHops >= he->key) {
						//spatial distance pruning
						if (check_hash_rn_dist[std::make_pair(pos, posDist)] && hash_rn_dist[std::make_pair(pos, posDist)] < sigma) {
							// structural pruning
							if (!StructuralCohesivenessPruning(*it, truss_val)) {
								// topic based pruning
								if (!InfluenceScorePruning(q, *it, Qtopic, SizeQtopic, theta)) {
									// keyword based pruning
									if (!KeywordbasedPruning(*it, k_set))

										hp->insert(he);
									finalSet.insert(*it);
								}
							}
						}
					}
					delete he;
				}
			}

		}
	}

	std::clock_t c_end = std::clock();
	long double time_elapsed_ms = (c_end - c_start) / 1000.0;

	std::cerr << "\n**********************************************************\n";
	std::cout << "CPU time used for GREEDY Baseline " << time_elapsed_ms << " s\n";
	std::cerr << "the final set size:: " << finalSet.size() << std::endl;
	std::cerr << "\n**********************************************************\n";

	hp->~Heap();
	delete[] visited;
	delete[] pruned;

}

*/