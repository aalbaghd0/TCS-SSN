#pragma once
#include <bitset>
#include <set>
#ifndef VERTEX_H
#define VERTEX_H


const int CKINS = 6 * 2; // 6 ckin locations
const int No_K= 10;
const int No_T= 10;
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

class edges {
public:
	int id;
	int from = 0;
	int to = 0;
	int weight = 1;
	int sup;
	double topics[No_T];
};

class node {
public:
	std::set<int> chldrn;
	int* s_ptr;
	int* e_ptr;
	int* parent;
};
#endif