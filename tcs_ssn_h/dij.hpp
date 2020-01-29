#ifndef DIG_HPP
#define DIG_HPP

#include<iostream>
#include<list>
#include"tcs_ssn_h/vertex.h"
#include"heap.h"
#include <list>
#include "parameterSettings.h"
#include <fstream>


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

			hash_rn_dist[std::make_pair(src, e)] = dist[e];
			check_hash_rn_dist[std::make_pair(src, e)] = true;

			hash_rn_dist[std::make_pair(e, src)] = dist[e];
			check_hash_rn_dist[std::make_pair(e, src)] = true;

			cost = he->key;
			break;
		}
		s[e] = true;
		f[e] = false;
		delete he;

		for (std::list<int>::iterator it = rnGraph[e].begin(); it != rnGraph[e].end(); ++it) {
			int to = *it;
			if (s[to] == false) {
				double new_dist = dist[e] + rn_edge_info[std::make_pair(e, to)].weight;
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

	/*
	this is a printing function
	*/
	/*
	if (cost != 0) {// if cost ==0, no path exsists
		int p = parent[dst];
		std::cerr << std::endl << dst << " " << p << " ";
		while (true) {
			if (p == src) break;
			p = parent[p];
			std::cerr << p << " ";
		}
	}
	*/
	delete[] s;
	delete[] f;
	delete[] dist;
	delete[] parent;
	return cost;
};



double sn_Dij(int src, int dst) {
	//initialise all vertices as unexplored 
	double* dist = new double[No_sn_V];
	int* parent = new int[No_sn_V];
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


	/*
	this is a printing function
	*/
	/*
	if (cost != 0) {// if cost ==0, no path exsists
		int p = parent[dst];
		std::cerr << std::endl << dst << " " << p << " ";
		while (true) {
			if (p == src) break;
			p = parent[p];
			std::cerr << p << " ";
		}
	}

	*/

	delete[] s;
	delete[] f;
	delete[] dist;
	delete[] parent;

	return cost;
}



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
		int v = he->son1;
		if (dst == v) {
			
			cost = he->key;
			break;
		}
		s[v] = true;
		f[v] = false;
		delete he;
		// explore the OUT set of v 

		//Node_N<the_edge>* temp = graphh[e].beginning;
		//list<the_edge> *temp = graphh[e].front;
		for (std::list<int>::iterator it = snGraph[v].begin(); it != snGraph[v].end(); ++it) {
			int to = *it;
			sum = 0.0;
			if (!s[to]) {
				for (int a = 0; a < No_of_TOPICS; ++a) {
					sum = sum + sn_edge_info[std::make_pair(v, to)].topics[a];
				}
				double new_dist = dist[v] + sum;
				if (new_dist < dist[to]) {
					dist[to] = new_dist;
					HeapEntry* he = new HeapEntry();
					parent[to] = v;
					he->son1 = to;
					he->key = new_dist;
					hp->insert(he);
					delete he;
				}
			}
		}
	}
	hp->~Heap();

	/*
	this is a printing function
	*/
	///*
	if (cost != 0) {// if cost ==0, no path exsists
		int p = parent[dst];
		std::cerr << std::endl << dst << " " << p << " ";
		while (true) {
			if (p == src) break;
			p = parent[p];
			std::cerr << p << " ";
		}
	}
	//*/

	delete[] s;
	delete[] f;
	delete[] dist;
	delete[] parent;

	return cost;
}

/////////////////////////////////////////////////////////////////////////////////////
/*
	Given	::	a social network and a start point vertex
	ENSURES ::  the shortest path distance to all the other social network vertices
*/
////////////////////////////////////////////////////////////////////////////////////
double sn_Dij_to_all_vertices(int src) {
	int dst = INT_MAX;
	//initialise all vertices as unexplored 
	double* dist = new double[No_sn_V];
	int* parent = new int[No_sn_V];
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
		//if (dst == e) {
			//cost = he->key;
			//break;
		//}
		if (dist[e] < INT_MAX) {
			hash_sn_dist[std::make_pair(src, e)] = dist[e];
			check_hash_sn_dist[std::make_pair(src, e)] = true;

			hash_sn_dist[std::make_pair(e, src)] = dist[e];
			check_hash_sn_dist[std::make_pair(e, src)] = true;
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


	/*
	this is a printing function
	*/
	/*
	if (cost != 0) {// if cost ==0, no path exsists
		int p = parent[dst];
		std::cerr << std::endl << dst << " " << p << " ";
		while (true) {
			if (p == src) break;
			p = parent[p];
			std::cerr << p << " ";
		}
	}

	*/

	for (int i = 0; i < No_rn_V; ++i) {
		if (!(check_hash_sn_dist[std::make_pair(src, i)] && check_hash_sn_dist[std::make_pair(i, src)])) {

			check_hash_sn_dist[std::make_pair(src, i)] = true;
			check_hash_sn_dist[std::make_pair(i, src)] = true;
			hash_sn_dist[std::make_pair(src, i)] = +INT_MAX;
			hash_sn_dist[std::make_pair(i, src)] = +INT_MAX;

		}
	}

	delete[] s;
	delete[] f;
	delete[] dist;
	delete[] parent;

	return cost;
}
//////////////////////////////////////////////////////////////////////////
/*
	Given	::	a social network and a start point vertex
	ENSURES ::  the shortest path distance to all the other social network vertices
*/
//////////////////////////////////////////////////////////////////////////
double rn_Dij_to_all_vertices(int src) {

	int dst = INT_MAX;

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
		//if (dst == e) {
			//cost = he->key;
			//break;
		//}
		if (dist[e] < INT_MAX) {
			hash_rn_dist[std::make_pair(src, e)] = dist[e];
			check_hash_rn_dist[std::make_pair(src, e)] = true;

			hash_rn_dist[std::make_pair(e, src)] = dist[e];
			check_hash_rn_dist[std::make_pair(e, src)] = true;

		}
		s[e] = true;
		f[e] = false;
		delete he;

		for (std::list<int>::iterator it = rnGraph[e].begin(); it != rnGraph[e].end(); ++it) {
			int to = *it;
			if (s[to] == false) {
				double new_dist = dist[e] + rn_edge_info[std::make_pair(e, to)].weight;
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

	/*
	this is a printing function
	*/
	/*
	if (cost != 0) {// if cost ==0, no path exsists
		int p = parent[dst];
		std::cerr << std::endl << dst << " " << p << " ";
		while (true) {
			if (p == src) break;
			p = parent[p];
			std::cerr << p << " ";
		}
	}
	*/
	// if we cannot hit the vertex, then we assign the distance to it to maximum
	for (int i = 0; i < No_rn_V; ++i) {
		if (!(check_hash_rn_dist[std::make_pair(src, i)] && check_hash_rn_dist[std::make_pair(i, src)])) {

			check_hash_rn_dist[std::make_pair(src, i)] = true;
			check_hash_rn_dist[std::make_pair(i, src)] = true;
			hash_rn_dist[std::make_pair(src, i)] = +INT_MAX;
			hash_rn_dist[std::make_pair(i, src)] = +INT_MAX;

		}
	}
	delete[] s;
	delete[] f;
 	delete[] dist;
	delete[] parent;
	return cost;
};

double inf_score_to_all_vertices(int src) {
	//initialise all vertices as unexplored 
	int dst = INT_MAX;
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
		int v = he->son1;
		if (dst == v) {

			hash_infScore[std::make_pair(src, v)] = dist[v];
			check_hash_infScore[std::make_pair(src, v)] = true;

			hash_infScore[std::make_pair(v, src)] = dist[v];
			check_hash_infScore[std::make_pair(v, src)] = true;

		}
		s[v] = true;
		f[v] = false;
		delete he;
		// explore the OUT set of v 

		//Node_N<the_edge>* temp = graphh[e].beginning;
		//list<the_edge> *temp = graphh[e].front;
		for (std::list<int>::iterator it = snGraph[v].begin(); it != snGraph[v].end(); ++it) {
			int to = *it;
			sum = 0.0;
			if (!s[to]) {
				for (int a = 0; a < No_of_TOPICS; ++a) {
					sum = sum + sn_edge_info[std::make_pair(v, to)].topics[a];
				}
				double new_dist = dist[v] + sum;
				if (new_dist < dist[to]) {
					dist[to] = new_dist;
					HeapEntry* he = new HeapEntry();
					parent[to] = v;
					he->son1 = to;
					he->key = new_dist;
					hp->insert(he);
					delete he;
				}
			}
		}
	}
	hp->~Heap();

	/*
	this is a printing function
	*/
	///*
	if (cost != 0) {// if cost ==0, no path exsists
		int p = parent[dst];
		std::cerr << std::endl << dst << " " << p << " ";
		while (true) {
			if (p == src) break;
			p = parent[p];
			std::cerr << p << " ";
		}
	}
	//*/

	delete[] s;
	delete[] f;
	delete[] dist;
	delete[] parent;

	return cost;
}



////////////////////////////////////////////////////////////////////////////////////
/*
	Given	::	a social network and a start point vertex TO Print
	ENSURES ::  the shortest path distance to all the other social network vertices
*/
////////////////////////////////////////////////////////////////////////////////////
double To_Print_sn_Dij_to_all_vertices(int src, std::ofstream& f_out) {
	
	
	int dst = INT_MAX;
	//initialise all vertices as unexplored 

	

	double* dist = new double[No_sn_V];
	int* parent = new int[No_sn_V];
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
		//if (dst == e) {
			//cost = he->key;
			//break;
		//}
		
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


	/*
	this is a printing function
	*/
	/*
	if (cost != 0) {// if cost ==0, no path exsists
		int p = parent[dst];
		std::cerr << std::endl << dst << " " << p << " ";
		while (true) {
			if (p == src) break;
			p = parent[p];
			std::cerr << p << " ";
		}
	}

	*/

	for (int i = 0; i < No_rn_V; ++i) {
		if (!(check_hash_sn_dist[std::make_pair(src, i)] && check_hash_sn_dist[std::make_pair(i, src)])) {

			check_hash_sn_dist[std::make_pair(src, i)] = true;
			check_hash_sn_dist[std::make_pair(i, src)] = true;
			hash_sn_dist[std::make_pair(src, i)] = +INT_MAX;
			hash_sn_dist[std::make_pair(i, src)] = +INT_MAX;

		}
	}
	for (int i = 0; i < No_sn_V; ++i)
		f_out << dist[i]<<" ";
	f_out << "\n";

	delete[] s;
	delete[] f;
	delete[] dist;
	delete[] parent;

	//f_out.close();
	return cost;
}

//////////////////////////////////////////////////////////////////////////
/*
	Given	::	a social network and a start point vertex
	ENSURES ::  the shortest path distance to all the other social network vertices
*/
//////////////////////////////////////////////////////////////////////////
double To_Print_rn_Dij_to_all_vertices(int src, std::ofstream& f_out) {

	int dst = INT_MAX;

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
		//if (dst == e) {
			//cost = he->key;
			//break;
		//}
		//if (dist[e] < INT_MAX) {
			//hash_rn_dist[std::make_pair(src, e)] = dist[e];
			//check_hash_rn_dist[std::make_pair(src, e)] = true;

			//hash_rn_dist[std::make_pair(e, src)] = dist[e];
			//check_hash_rn_dist[std::make_pair(e, src)] = true;

		//}
		s[e] = true;
		f[e] = false;
		delete he;

		for (std::list<int>::iterator it = rnGraph[e].begin(); it != rnGraph[e].end(); ++it) {
			int to = *it;
			if (s[to] == false) {
				double new_dist = dist[e] + rn_edge_info[std::make_pair(e, to)].weight;
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

	/*
	this is a printing function
	*/
	/*
	if (cost != 0) {// if cost ==0, no path exsists
		int p = parent[dst];
		std::cerr << std::endl << dst << " " << p << " ";
		while (true) {
			if (p == src) break;
			p = parent[p];
			std::cerr << p << " ";
		}
	}
	*/
	// if we cannot hit the vertex, then we assign the distance to it to maximum
	for (int i = 0; i < No_rn_V; ++i) {
		if (!(check_hash_rn_dist[std::make_pair(src, i)] && check_hash_rn_dist[std::make_pair(i, src)])) {

			check_hash_rn_dist[std::make_pair(src, i)] = true;
			check_hash_rn_dist[std::make_pair(i, src)] = true;
			hash_rn_dist[std::make_pair(src, i)] = +INT_MAX;
			hash_rn_dist[std::make_pair(i, src)] = +INT_MAX;

		}
	}

	for (int i = 0; i < No_rn_V; ++i)
		f_out << dist[i] << " ";
	f_out << "\n";


	delete[] s;
	delete[] f;
	delete[] dist;
	delete[] parent;
	return cost;
};


#endif // !1