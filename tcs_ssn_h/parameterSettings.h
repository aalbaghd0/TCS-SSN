#pragma once
#ifndef SETTINGS_HPP
#define SETTINGS_HPP

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
// file names

#endif // !SETTINGS_HPP