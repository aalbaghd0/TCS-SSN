#ifndef READ_FILES
#define READ_FILES

#define _CRT_SECURE_NO_DEPRECATE
#pragma warning (disable : 4996)

#include <iostream>
#include "vertex.h"
#include<fstream>
#include<list>
#include"parameterSettings.h"



void read_vertices(vertex vrtx[], int size, char Vfile[]) {
	FILE* file = fopen(Vfile, "r");
	int v_id, i, j, r, f = 0, r1;
	int r2, r3;
	std::pair<double, double> pir;
	j = 0;
	if (file == NULL) {
		std::cout << "could not open the file to read social network vertices";
		return;
	}
	// get ride of the number of lines in the file
	int getRide;
	fscanf(file, "%d", &getRide);

	while (!feof(file)) {
		fscanf(file, "%d", &v_id);
		//std::cerr << v_id << " ";
		for (i = 1; i <= 12; i++) {
			fscanf(file, "%d", &r2);
			vrtx[v_id].ckins[i] = r2;
			std::cerr << vrtx[v_id].ckins[i] << " ";
		}
		std::cerr << "\n";
		for (i = 1; i <= 5; i++) {
			fscanf(file, "%d", &r);
			vrtx[v_id].key[r] = 1;
			//std::cerr << r << " ";
		}
		for (i = 1; i <= 5; i++) {
			fscanf(file, "%d", &r);
			//vrtx[v_id].topic[r] = 1;
			//std::cerr << r << " ";
		}
		//std::cerr << " \n";
	}
	fclose(file);
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

void sn_read_edges(std::list<int> edge[], char Efile[]) {
	FILE* file = fopen(Efile, "r");
	int from, to, e_id = 0;

	if (file == NULL) {
		std::cout << "The edge file cannot be open";
		return;
	}
	// get ride of the number of lines in the file
	int getRide;
	fscanf(file, "%d", &getRide);

	while (!feof(file)) {
		fscanf(file, "%d %d", &from, &to);
		edge[from].push_back(to);
		edge[to].push_back(from);
		
		sn_vrtx[from].nbrs.insert(to);
		sn_vrtx[to].nbrs.insert(from);

		/*
		//get edges associated with each 
		sn_vrtx[from].myedges.insert(e_id);
		sn_vrtx[to].myedges.insert(e_id);
		sn_vrtx[from].myedges.insert((No_sn_E / 2) + e_id);
		sn_vrtx[to].myedges.insert((No_sn_E / 2) + e_id);
		*/

		snEdges[e_id].from = from;
		snEdges[e_id].to = to;

		snEdges[(No_sn_E / 2) + e_id].from = to;
		snEdges[(No_sn_E / 2) + e_id].to = from;

		//hash_edge[std::make_pair(from, to)] = e_id;
		//hash_edge[std::make_pair(to, from)] = (No_sn_E / 2) + e_id;

		e_id++;
	}
	fclose(file);
}

void rn_read_edges(std::list<int> graph[], edges edge[], int E_num,int size, char Efile[]) {
	FILE* file = fopen(Efile, "r");
	int from, to, id;
	int max = 0;
	double w = 0;
	if (file == NULL) {
		std::cout << "The edge file cannot be open";
		return;
	}
	// get ride of the number of lines in the file
	int getRide;
	fscanf(file, "%d", &getRide);


	while (!feof(file)) {
		fscanf(file, "%d %d %d %lf", &id, &from, &to, &w);
		edge[id].weight = w;
		edge[(E_num / 2) + id].weight = w;
		graph[from].push_back(to);
		graph[to].push_back(from);
		if (max < from)
			max = from;
		if (max < to)
			max = to;
	}
	fclose(file);
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