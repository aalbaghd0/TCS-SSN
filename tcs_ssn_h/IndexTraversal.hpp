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
	int root = 0;

	HeapEntry* he = new HeapEntry();

	he->son1 = root;
	he->key = 0;
	hp->insert(he);
	delete he;

	while (hp->used > 0) {
		HeapEntry* he = new HeapEntry();
		hp->remove(he);
		
		int theNode = he->son1;

		if (tree[theNode].level == 0) {


		}
		else{ // non-leaf node

			int queryNode = get_queryNode(tree, q);
			
			for (std::unordered_set<int>::iterator it = tree[theNode].child.begin(); it != tree[theNode].child.end(); ++it) {

			}
		}


	}


}

#endif // !INDEX_TRAV

