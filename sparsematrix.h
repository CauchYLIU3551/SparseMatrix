#pragma once
#ifndef _sparsematrix_h_
#define _sparsematrix_h_

#include<iostream>
#include<cstdio>
#include<stdio.h>
#include<vector>
#include<math.h>

//template <class T1>
class sparsematrix
{
public:

	sparsematrix(); 
	sparsematrix(std::vector<int> n, std::vector<int> j, std::vector<double> v);
	sparsematrix(std::vector<int> n, std::vector<int> j, std::vector<double> v, int row, int col);

	sparsematrix transpose();
	sparsematrix inverse();
	
	sparsematrix diag();
	sparsematrix lowtri();
	sparsematrix uppertri();

	static sparsematrix add(sparsematrix A, sparsematrix B);
	std::vector<double> multiply(std::vector<double> x);
	void showmatrix();
	std::vector<int> show_offset();
	std::vector<int> show_indice();
	std::vector<double> show_value();
	//std::vector<int> rowoffset;
	//std::vector<int> indice;
	//std::vector<double> value;

private:
	std::vector<int> rowoffset;
	std::vector<int> indice;
	std::vector<double> value;
	int col;
	int row;
	int max();
};
//template<class T1>

// get the transpose matrix of the sparse matrix.


#endif