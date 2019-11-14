#ifndef DIG_HPP
#define DIG_HPP

#include<iostream>
#include<list>
#include"tcs_ssn_h/vertex.h"
#include"heap.h"
#include <list>
#include "parameterSettings.h"


double rn_Dij(int src, int dst) {
	
	int* dist = new int[No_rn_V];
	int* parent = new int[No_rn_V];
	bool* s = new bool[No_rn_V];
	bool* f = new bool[No_rn_V];

	for (int v = 0; v < No_rn_V; v++) {
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

		for (std::list<int>::iterator it = rnGraph[e].begin(); it != rnGraph[e].end(); ++it) {
			int to = *it;
			if (s[to] == false) {
				double new_dist = dist[e] + rnEdges[to].weight;
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
	delete[] s;
	delete[] f;
	delete[] dist;
	delete[] parent;
	return cost;
};



double sn_Dij(int src, int dst) {
	//initialise all vertices as unexplored 
	double * dist = new double[No_sn_V];
	int *parent = new int[No_sn_V];
	bool* s = new bool[No_sn_V];
	bool* f = new bool[No_sn_V];


	for (int v = 0; v < No_sn_V; v++) {
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
		for (std::list<int>::iterator it = snGraph[e].begin(); it != snGraph[e].end(); ++it) {
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

	delete[] s;
	delete[] f;
	delete[] dist;
	delete[] parent;

	return cost;
}


/*
double inf_score(int src, int dst) {
	//initialise all vertices as unexplored 
	double* dist = new double[No_sn_V];
	int* parent = new int[No_sn_V];
	bool* s = new bool[No_sn_V];
	bool* f = new bool[No_sn_V];
	double sum = 0;

	for (int v = 0; v < No_sn_V; v++) {
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
		for (std::list<int>::iterator it = snGraph[e].begin(); it != snGraph[e].end(); ++it) {
			int to = *it;
			sum = 0.0;
			if (!s[to]) {
				for (int a = 0; a < No_topics; ++a) {
					sum = sum + snEdges[to].topics[a];
				}
				double new_dist = dist[e] * sum ;
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

	delete[] s;
	delete[] f;
	delete[] dist;
	delete[] parent;

	return cost;
}
*/
#endif // !1