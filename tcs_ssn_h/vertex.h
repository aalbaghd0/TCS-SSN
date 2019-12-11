#pragma once
#include <bitset>
#include <set>
#include<unordered_set>
#ifndef VERTEX_H
#define VERTEX_H


const int CKINS = 6 * 2; // 6 ckin locations
const int No_K = 10;
const int No_T = 10;

#define No_of_TOPICS	    3


class vertex {
public:
	int id;
	int ckins[CKINS];
	std::bitset<No_K> key;
	std::bitset<No_T> topic;
	std::set<int> nbrs;
	std::set<int> myedges;
	int truss;
	std::unordered_set<int> out_nbrs;
	std::unordered_set<int> in_nbrs;

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
	std::set<int> ptr;
	int parent;
	int start;
	int end;
	int n_id;
	int truss;
	int level;

	std::unordered_map<pair, double, pair_hash> rn_min_dist_to_piv;
	double rn_max_dist_to_piv;
	int rn_minimum_piv;
	int rn_maximum_piv;

	double sn_min_dist_to_piv;
	double sn_max_dist_to_piv;
	int sn_minimum_piv;
	int sn_maximum_piv;

	int in_topics_prob[No_of_TOPICS];
	int out_topics_prob[No_of_TOPICS];

	std::bitset<No_K> keys;

};
#endif