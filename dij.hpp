#ifndef DD
#define DD


#include<iostream>
#include<list>
#include"vertex.h"
#include"heap.h"
#include <list>


double rn_Dij(std::list<int> graph[], edges edge[],int V_num, int src, int dst) {
	
	int* dist = new int[V_num];
	int* parent = new int[V_num];
	bool* s = new bool[V_num];
	bool* f = new bool[V_num];

	for (int v = 0; v < V_num; v++) {
		s[v] = false;
		f[v] = false;
		dist[v] = INT_MAX;
		parent[v] = NULL;
	}
	Heap* hp = new Heap();
	HeapEntry* he = new HeapEntry();

	hp->used = 0; //clear the heap

	hp->init(2);
	dist[src] = 0;
	he->key = 0;
	he->son1 = src;
	f[src] = true;
	hp->insert(he);
	parent[src] = NULL;
	double cost = 0;
	delete he;
	while (hp->used > 0) {
		HeapEntry* he = new HeapEntry();
		hp->remove(he);
		int e = he->son1;
		if (dst == e) {
			cost = he->key;
			break;
		}
		s[e] = true;
		f[e] = false;
		delete he;

		for (std::list<int>::iterator it = graph[e].begin(); it != graph[e].end(); ++it) {
			int to = *it;
			if (s[to] == false) {
				double new_dist = dist[e] + edge[to].weight;
				if (new_dist < dist[to]) {
					dist[to] = new_dist;
					HeapEntry* he = new HeapEntry();
					parent[to] = e;
					he->son1 = to;
					he->key = new_dist;
					hp->insert(he);
					delete he;
				}
			}
		}
	}
	hp->~Heap();

	if (cost != 0) {// if cost ==0, no path exsists
		int p = parent[dst];
		std::cerr << std::endl << dst << " " << p << " ";
		while (true) {
			if (p == src) break;
			p = parent[p];
			std::cerr << p << " ";
		}
	}

	return cost;
};



double sn_Dij(std::list<int> graph[], int E_numm, int V_num,int src, int dst) {
	//initialise all vertices as unexplored 
	int * dist = new int[V_num];
	int *parent = new int[V_num];
	bool* s = new bool[V_num];
	bool* f = new bool[V_num];


	for (int v = 0; v < V_num; v++) {
		s[v] = false;
		f[v] = false;
		dist[v] = INT_MAX;
		parent[v] = NULL;
	}
	Heap* hp = new Heap();
	HeapEntry* he = new HeapEntry();

	hp->used = 0; //clear the heap

	hp->init(2);
	dist[src] = 0;
	he->key = 0;
	he->son1 = src;
	f[src] = true;
	hp->insert(he);
	parent[src] = NULL;
	//delete he;
	double cost = 0;
	delete he;


	while (hp->used > 0) {
		HeapEntry* he = new HeapEntry();
		hp->remove(he);
		int e = he->son1;
		if (dst == e) {
			cost = he->key;
			break;
		}
		s[e] = true;
		f[e] = false;
		delete he;
		// explore the OUT set of v 

		//Node_N<the_edge>* temp = graphh[e].beginning;
		//list<the_edge> *temp = graphh[e].front;
		for (std::list<int>::iterator it = graph[e].begin(); it != graph[e].end(); ++it) {
			int to = *it;
			if (!s[to]) {
				double new_dist = dist[e] + 1;
				if (new_dist < dist[to]) {
					dist[to] = new_dist;
					HeapEntry* he = new HeapEntry();
					parent[to] = e;
					he->son1 = to;
					he->key = new_dist;
					hp->insert(he);
					delete he;
				}
			}
		}
	}
	hp->~Heap();

	if (cost != 0) {// if cost ==0, no path exsists
		int p = parent[dst];
		std::cerr << std::endl << dst << " " << p << " ";
		while (true) {
			if (p == src) break;
			p = parent[p];
			std::cerr << p << " ";
		}
	}
	return cost;
}






#endif // !1