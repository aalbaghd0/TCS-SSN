#pragma once
#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include<unordered_map>
#include<fstream>


char sn_V_file[100] =			"data/SN_vertices.txt";
char sn_E_file[100] =			"data/SN_edges.txt";
char rn_V_file[100] =			"data/RN_vertices.txt";
char rn_E_file[100] =			"data/RN_edges.txt";


int set_parameters(char filename[]) {
	FILE* file = fopen(filename, "r");
	int read;
	if (file == NULL) {
		std::cout << "could not open the file to read social network vertices";
		return -1;
	}
	else {
		fscanf(file, "%d", &read);
	}
	fclose(file);
	return read;
}



// file sizes
#define	No_sn_V				set_parameters(sn_V_file)
#define No_sn_E				set_parameters(sn_E_file) * 2
#define No_rn_V				set_parameters(rn_V_file)
#define No_rn_E				set_parameters(rn_E_file) * 2
#define No_index_piv			10
#define No_SN_piv				10
#define No_RN_piv				10
#define No_subgraphs			10
#define No_CKINs				6

//set parameters//
#define W1					0.5
#define W2					0.5
#define W3					0.5
#define No_topics			3



//////

int* index_piv = new int[No_index_piv];

// file names



/// some globale structurs 
vertex* sn_vrtx = new vertex[No_sn_V];
vertex* rn_vrtx = new vertex[No_rn_V];

std::list<int>* snGraph = new std::list<int>[No_sn_V];
std::list<int>* rnGraph = new std::list<int>[No_rn_V];

edges* snEdges = new edges[No_sn_E];
edges* rnEdges = new edges[No_rn_E];

node* index = new node[100];

/// some global functions
int uniform(int _min, int _max) {
	//	cout<<_min<<"  "<<_max<<endl;
	int int_r = rand();
	long base = RAND_MAX - 1;
	float f_r = ((float)int_r) / base;
	return (int) (_max - _min) * f_r + _min;
}
//////////////////////////////////////////
typedef std::pair<int, int> pair;
struct pair_hash {
	template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2>& pair) const {
		return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	}
};
//////////////
std::unordered_map<pair, int, pair_hash> hash_edge;



// functions

#endif // !SETTINGS_HPP