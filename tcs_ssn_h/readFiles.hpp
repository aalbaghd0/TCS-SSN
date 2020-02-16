#ifndef READ_FILES
#define READ_FILES

#define _CRT_SECURE_NO_DEPRECATE
#pragma warning (disable : 4996)

#include <iostream>
#include "vertex.h"
#include<fstream>
#include<list>
#include"parameterSettings.h"
#include"gendef.h"


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
		//std::cerr << v_id << "\n";

		//std::cerr << v_id << " ";
		for (i = 0; i < 4; i++) {
			fscanf(file, "%d", &r2);
			vrtx[v_id].ckins[i] = r2;
			rn_vrtx[r2].myUsers.insert(v_id);
			//std::cerr << vrtx[v_id].ckins[i] << " ";

		}
		//std::cerr << "\n";
		for (i = 1; i <= 5; i++) {
			fscanf(file, "%d", &r);
			vrtx[v_id].key[r] = 1;
			//std::cerr << r << " ";
		}
		for (i = 1; i <= 5; i++) {
			fscanf(file, "%d", &r);
			vrtx[v_id].topic[r] = 1;
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
	int from, to;
	double  pr0, pr1, pr2;
	int max_vrtx = -INT_MAX;
	
	int e_num = 0;


	if (file == NULL) {
		std::cout << "The edge file cannot be open";
		return;
	}
	// get ride of the number of lines in the file
	int getRide;
	fscanf(file, "%d", &getRide);

	while (!feof(file)) {
		fscanf(file, "%d %d %lf %lf %lf", &from, &to, &pr0, &pr1, &pr2);
		edge[from].push_back(to);
		edge[to].push_back(from);

		sn_vrtx[from].nbrs.insert(to);
		sn_vrtx[to].nbrs.insert(from);

		sn_vrtx[from].nbrs_set.insert(to);
		sn_vrtx[to].nbrs_set.insert(from);


		//get edges associated with each 
		sn_vrtx[from].myedges.insert(e_num);
		sn_vrtx[to].myedges.insert(e_num);
		sn_vrtx[from].myedges.insert(e_num + 1);
		sn_vrtx[to].myedges.insert(e_num + 1);

		//get maximum vrtx in the file
		max_vrtx = max(max_vrtx, from);
		max_vrtx = max(max_vrtx, to);


		sn_vrtx[to].in_nbrs.insert(from);
		sn_vrtx[to].out_nbrs.insert(from);

		sn_vrtx[from].in_nbrs.insert(to);
		sn_vrtx[from].out_nbrs.insert(to);



		snEdges[e_num].from = from;
		snEdges[e_num].to = to;

		sn_edge_info[std::make_pair(from, to)].topics[0] = pr0 * 0.8;
		sn_edge_info[std::make_pair(from, to)].topics[1] = pr1 * 0.8;
		sn_edge_info[std::make_pair(from, to)].topics[2] = pr2 * 0.8;

		sn_edge_info[std::make_pair(to, from)].topics[0] = pr0 * 0.8;
		sn_edge_info[std::make_pair(to, from)].topics[1] = pr1 * 0.8;
		sn_edge_info[std::make_pair(to, from)].topics[2] = pr2 * 0.8;
		
		


		snEdges[e_num + 1].from = to;
		snEdges[e_num + 1].to = from;

		

		hash_edge[std::make_pair(from, to)] = e_num;
		hash_edge[std::make_pair(to, from)] = e_num + 1;

		e_num = e_num + 2;

		//std::cerr << "e_id " << e_id << " from: " << from << " to: " << to << std::endl;
		//std::cerr << "e_id " << (No_sn_E / 2) + e_id << " from: " << to << " to: " << from << std::endl;
		
	}
	fclose(file);
	std::cout << "The maximum vrtx is: " << max_vrtx << std::endl;
}

void rn_read_edges(std::list<int> graph[], RN_edges edge[], int E_num, int size, char Efile[]) {
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

		rn_edge_info[std::make_pair(from, to)].weight = w;
		rn_edge_info[std::make_pair(to, from)].weight = w;

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