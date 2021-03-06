#pragma once
#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include<unordered_map>
#include<fstream>


char sn_V_file[100] = "data/200K_Gow_Cali/New_mapping_By_groups.txt";
char sn_E_file[100] = "data/200K_Gow_Cali/SN_Gwo_edges.txt";
char rn_V_file[100] = "data/200K_Gow_Cali/RN_cali_vertices.txt";
char rn_E_file[100] = "data/200K_Gow_Cali/RN_cali_edges.txt";


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
//////////////////////////////////
/*        parameters        */
#define			SIGMA			0.5
#define			THETA			0.5
#define			No_Hops			3 // the value of d
#define			Ktruss			3 // the truss k 



///////////////////////////////


// file sizes
#define	No_sn_V				set_parameters(sn_V_file)
#define No_sn_E				set_parameters(sn_E_file) * 2
#define No_rn_V				set_parameters(rn_V_file)
#define No_rn_E				set_parameters(rn_E_file) * 2

#define No_index_piv			273
#define No_subgraphs			No_index_piv

#define No_SN_piv				4
#define No_RN_piv				4
#define No_CKINs				3



//set parameters//
#define W1					0.5
#define W2					0.5
#define W3					0.5
#define No_of_TOPICS	    3



//////

int* index_piv = new int[No_index_piv];

// file names



/// some globale structurs 
vertex* sn_vrtx = new vertex[No_sn_V];
rn_vertices* rn_vrtx = new rn_vertices[No_rn_V];



std::list<int>* snGraph = new std::list<int>[No_sn_V];
std::list<int>* snGraph_fnl = new std::list<int>[No_sn_V];
std::list<int>* rnGraph = new std::list<int>[No_rn_V];

SN_edges* snEdges = new SN_edges[No_sn_E];
RN_edges* rnEdges = new RN_edges[No_rn_E];

int* RN_piv_set = new int[No_RN_piv];
int* SN_piv_set = new int[No_SN_piv];


std::set<int>* index = new std::set<int>[1000];
int GlobalIndexIter = 0;
/// some global functions
int uniform(int _min, int _max) {
	//	cout<<_min<<"  "<<_max<<endl;
	int int_r = rand();
	long base = RAND_MAX - 1;
	float f_r = ((float)int_r) / base;
	return (int)(_max - _min) * f_r + _min;
}

/// some global functions
double uniform_dou(int _min, int _max) {
	//	cout<<_min<<"  "<<_max<<endl;
	int int_r = rand();
	long base = RAND_MAX - 1;
	float f_r = ((float)int_r) / base;
	return (int)(_max - _min) * f_r + _min;
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

std::unordered_map<pair, bool, pair_hash> f_edges;

std::unordered_map<pair, SN_edges, pair_hash> sn_edge_info;
std::unordered_map<pair, RN_edges, pair_hash> rn_edge_info;

std::unordered_map<pair, double, pair_hash> hash_sn_dist;
std::unordered_map<pair, bool, pair_hash> check_hash_sn_dist;

std::unordered_map<pair, double, pair_hash> hash_rn_dist;
std::unordered_map<pair, bool, pair_hash> check_hash_rn_dist;

std::unordered_map<pair, double, pair_hash> hash_rnToUser_dist;
std::unordered_map<pair, bool, pair_hash> check_hash_rnToUser_dist;

std::unordered_map<pair, double, pair_hash> hash_infScore;
std::unordered_map<pair, bool, pair_hash> check_hash_infScore;


std::unordered_map<int, int> hash_my_position_in_tree;

///////////////////////////////////////
/*
	the unique function,
	REQUIRES :: an array of integers and an integer
	ENSURES  :: the given integer is not in the array
*/
bool isInTheArray(int arr[], int arr_size, int x) {
	for (int i = 0; i < arr_size; i++) {
		if (arr[i] == x)
			return true;
	}
	return false;
}



/////////////////////////////////
// se the index size and define the tree
int number_nodes(int a, int b, int& c, int PivsInLevelTree[]) {
	int count = 0;
	c = a;
	
	int j = 0;
	PivsInLevelTree[j] = a;

	while (a > 0) {
		j++;

		if (a > b) {
			a = a / b;
			c = c + a;
			count++;
			PivsInLevelTree[j] = a;
		}
		else {
			c++;
			a = a / b;
			count++;
			PivsInLevelTree[j] = 1;
		}
	}

	for (int i = 0; i <= count; i++) {
		std::cerr << PivsInLevelTree[i] << "  ";
	}

	return count;
	std::cerr << "\n" << c;
}
int PivsInLevelTree[1000];

// control the number of pivots in tree layers
#define		NO_INTER_PIVS		8

int INDEXSIZE = 0;
int LengthOfPivLevels = number_nodes(No_index_piv, NO_INTER_PIVS, INDEXSIZE, PivsInLevelTree) + 1;
Gnode* tree = new Gnode[INDEXSIZE];


#endif // !SETTINGS_HPP