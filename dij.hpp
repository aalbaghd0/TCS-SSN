#ifndef DD
#define DD


#include<iostream>
#include<list>
#include"vertex.h"
#include"heap.h"


void dij(std::list<int> edge[], int V_num, int src, int dst) {
	Heap* hp;
	int* prev = new int[V_num];
	int* dist = new int[V_num];

	for (int v = 0; v < V_num; v++) {
		prev[v] = -1;
	}
	HeapEntry* he = new HeapEntry();
	he->son1 = src;
	he->key = 0;
	dist[src] = 0;
	hp->insert(he);
	delete he;


	while (hp->used > 0) {
		HeapEntry *he = new HeapEntry;
		hp->remove(he);
		for (std::list<int>::iterator it = edge[he->son1].begin(); it != edge[he->son1].end(); ++it) {
			HeapEntry* he2 = new HeapEntry();
			he2->son1 = *it;
			he2->key = dist[he->son1] + 1;
			if (he2->key < dist[*it]) {
				dist[*it] = he2->key;
				prev[*it] = he->son1;
				hp->insert(he2);
				delete he2;
			}

		}
	}
	delete hp;

}
#endif // !1