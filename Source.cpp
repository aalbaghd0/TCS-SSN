#include <iostream>
#include "vertex.h"
#include<fstream>
#include"readFiles.hpp"
#include <list>


const int No_V = 100;
const int No_E = 100;
char Vfile[10] = "V.txt";
char Efile[10] = "E.txt";
vertex *vrtx = new vertex[No_V];
std::list<int> * e_list = new std::list<int>[100];

void main(int args, char arg[]) {
	read_vertices(vrtx, 10, Vfile);
	
	
	
	int i, j;
	for (i = 0; i < 10; i++) {
		for (j = 0; j < 12; j++)
			std::cerr << vrtx[i].ckins[j] << " ";
		for (j = 0; j < 10; j++)
			std::cerr << vrtx[i].key[j] << " ";
		for (j = 0; j < 10; j++)
			std::cerr << vrtx[i].topic[j] << " ";
		std::cout << std::endl;
	}
	
	
	read_edges(e_list, 10, Efile);
	for (int i = 0; i < 10; i++) {
		for (std::list<int>::iterator it = e_list[i].begin(); it != e_list[i].end(); ++it) {
			std::cerr << *it << "-->";
		}
		std::cout << "\n";
	}
};
