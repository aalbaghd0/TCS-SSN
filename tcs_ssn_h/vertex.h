#pragma once
#include <bitset>
#include <set>
#include<unordered_set>
#ifndef VERTEX_H
#define VERTEX_H


const int CKINS = 6 * 2; // 6 ckin locations
const int No_K= 10;
const int No_T= 10;

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
	std::unordered_set<int> child;
	int start;
	int end;
	int n_id;
};
#endif