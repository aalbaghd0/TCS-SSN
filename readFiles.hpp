#ifndef READ_FILES
#define READ_FILES

#include <iostream>
#include "vertex.h"
#include<fstream>
#include<list>

#pragma warning(disable:4996)

void read_vertices(vertex vrtx[], int size, char Vfile[]) {
	FILE* file = fopen(Vfile, "r");
	int v_id, i, j, r, f = 0;
	j = 0;

	if (file == NULL) {
		std::cerr << "could not open the file to read vertices";
		return;
	}

	while (f < 100) {
		++f;
		fscanf(file, "%d", &v_id);
		j = 0;
		for (i = 1; i <= 12; i++) {
			fscanf(file, "%lf", &vrtx[v_id].ckins[j]);
			j++;
		}
		std::cerr << std::endl;
		for (i = 1; i <= 5; i++) {
			fscanf(file, "%d", &r);
			vrtx[v_id].key[r] = 1;
		}
		for (i = 1; i <= 5; i++) {
			fscanf(file, "%d", &r);
			vrtx[v_id].topic[r] = 1;
		}
	}
	/*
	// print the vrtx array for checking
	for (i = 0; i < 10; i++) {
		for (j = 0; j < 12; j++)
			std::cerr << vrtx[i].ckins[j] << " ";
		for (j = 0; j < 10; j++)
			std::cerr << vrtx[i].key[j] << " ";
		for (j = 0; j < 10; j++)
			std::cerr << vrtx[i].topic[j]<< " ";
		std::cout << std::endl;
	}
	*/
};

void read_edges(std::list<int> edge[], int size, char Efile[]) {
	FILE* file = fopen(Efile, "r");
	
	if (file == NULL) {
		std::cout << "The edge file cannot be open";
		return;
	}
	int from, to;
	int f = 0;
	while (f<100) {
		++f;
		fscanf(file, "%d %d", &from, &to);
		edge[from].push_back(to);
		edge[to].push_back(from);
	}

	/*
	// test the edge list
	for (int i = 0; i < 10; i++) {
		for (std::list<int>::iterator it = edge[i].begin(); it != edge[i].end(); ++it) {
			std::cerr << *it << " ";
		}
		std::cout << "\n";
	}
	*/
}
#endif // !READ_FILES