#pragma once
#include <bitset>
#include <set>
#include<unordered_set>
#include<unordered_map>
#ifndef VERTEX_H
#define VERTEX_H

//////////////////////////////////////////
typedef std::pair<int, int> pair;
struct pair_hashh {
	template <class T1, class T2>
	std::size_t operator() (const std::pair<T1, T2>& pair) const {
		return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	}
};

const int CKINS = 6 * 2; // 6 ckin locations
const int No_K = 10;
const int No_T = 10;

#define No_of_TOPICS	    3

class rn_vertices {
public:
	std::unordered_set<int> myUsers;
};


class vertex {
public:
	int id;
	int ckins[CKINS];
	std::bitset<No_K> key;
	std::bitset<No_T> topic;
	std::unordered_set<int> nbrs;

	std::unordered_set<int> myedges;
	int truss;
	std::unordered_set<int> out_nbrs;
	std::unordered_set<int> in_nbrs;
	std::unordered_map<pair, double, pair_hashh> rn_distToPiv;
	std::unordered_map<pair, double, pair_hashh> sn_distToPiv;

	/// need to be defined
	double in_inf [No_of_TOPICS];
	double out_inf [No_of_TOPICS];

};

class RN_edges {
public:
	int id;
	int from = 0;
	int to = 0;
	double weight = 0;
};

class SN_edges {
public:
	int id;
	int from = 0;
	int to = 0;
	int weight = 1;
	int sup;
	double topics[No_of_TOPICS];
};

class node {
public:
	std::set<int> chldrn;
	int* s_ptr;
	int* e_ptr;
	int* parent;
};

class Gnode {
public:
	std::unordered_set<int> child;
	std::unordered_set<int> ptr;
	int parent;
	int start;
	int end;
	int n_id;
	int truss;
	int level;

	std::unordered_map<pair, double, pair_hashh> rn_min_dist_to_piv;
	std::unordered_map<pair, double, pair_hashh> rn_max_dist_to_piv;
	int rn_minimum_piv;
	int rn_maximum_piv;

	double sn_min_dist_to_piv;
	double sn_max_dist_to_piv;
	int sn_minimum_piv;
	int sn_maximum_piv;

	double in_topics_prob[No_of_TOPICS];
	double out_topics_prob[No_of_TOPICS];

	std::bitset<No_K> keys;

};
#endif