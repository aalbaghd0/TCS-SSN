#ifndef INDEX_TRAV
#define INDEX_TRAV
#include "parameterSettings.h"
#include "vertex.h"
#include"heap.h"
/*

	This pice of the code will traverse the index to get the solution

*/


void indexTrav() {

	Heap* hp = new Heap();
	Heap* hp_p = new Heap();
	Heap* temp = new Heap();

	std::unordered_set<int> S;


	int root = 0;

	HeapEntry* he = new HeapEntry();
	
	he->son1 = root;
	he->key = 0;
	hp->insert(he);
	delete he;
	
	// this is an update
	while (hp->used > 0) {
		HeapEntry* he = new HeapEntry();
		hp->remove(he);
		delete he;

		int theNode = he->son1;

		if (tree[theNode].level == 0) { // if a leaf node

			for (std::unordered_set<int>::iterator it = tree[theNode].child.begin(); it != tree[theNode].child.end(); ++it) {

				if (!Lem2(q, *it)) { // if the object cannot be pruned

					S.insert(*it);

				}

			}

		}
		else{ // non-leaf node

			// optain the tree node that contains the query vertex
			int queryNode = get_queryNode(tree, q);
			
			for (std::unordered_set<int>::iterator it = tree[theNode].child.begin(); it != tree[theNode].child.end(); ++it) {
				if (!Lem7(q, *it)) {
					
					HeapEntry* he = new HeapEntry();

					he->son1 = *it;
					he->key = lb_dist_RN(queryNode, *it);

					hp_p->insert(he);
					delete he;
				}
			}
		}

		Heap* temp = new Heap();
		temp = hp;
		hp = hp_p;
		delete hp_p;
		delete temp;

	}

	S = Refine(S);

	//return S;
}

#endif // !INDEX_TRAV

