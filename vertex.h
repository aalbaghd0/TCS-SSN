#pragma once
#include <bitset>
#ifndef VERTEX_H
#define VERTEX_H


const int CKINS = 6 * 2; // 6 ckin locations
const int No_K= 10;
const int No_T= 10;
class vertex {
public:
	int id;
	double ckins[CKINS];
	std::bitset<No_K> key;
	std::bitset<No_T> topic;
};

class edges {
public:
	int from;
	int to;
};


#endif