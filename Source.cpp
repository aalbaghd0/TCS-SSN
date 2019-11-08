#include <iostream>
#include "vertex.h"
#include<fstream>
#include"readFiles.hpp"
#include <list>
#include"heap.h"

const int No_V = 100;
const int No_E = 40000;
char Vfile[10] = "V.txt";
char Efile[10] = "E.txt";
vertex *vrtx = new vertex[No_V];
std::list<int> * e_list = new std::list<int>[No_E];

void dij(std::list<int> edge[], int V_num, int src, int dst);

void main(int args, char arg[]) {
	read_vertices(vrtx, 10, Vfile);
	read_edges(e_list, 10, Efile);
	dij(e_list, No_V, 0, 10);
	
	
};


void dij(std::list<int> edge[], int V_num, int src, int dst) {
	Heap* hp = new Heap();
	
	int* prev = new int[V_num];
	int* dist = new int[V_num];

	for (int v = 0; v < V_num; v++) {
		prev[v] = -1;
	}
	HeapEntry* he = new HeapEntry();
	he->son1 = src;
	he->key = 0;
	he->level = 0;
	dist[src] = 0;
	hp->init(2);
	hp->insert(he);
	delete he;


	while (hp->used > 0) {
		HeapEntry* he = new HeapEntry;
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