#pragma once
#ifndef _sparsematrix_h_
#define _sparsematrix_h_

#include<sparse/sparsity_pattern.h>

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
	
	// following functions are designed to solve the implement the Jacobi iteration.
	sparsematrix diag();
	sparsematrix diag_inverse();
	sparsematrix lowtri();
	sparsematrix uppertri();
	sparsematrix inverse_num();

	static sparsematrix add(sparsematrix A, sparsematrix B);
	std::vector<double> multiply(std::vector<double> x);
	void showmatrix();
	std::vector<int> show_offset();
	std::vector<int> show_indice();
	std::vector<double> show_value();

	int show_row();
	//std::vector<int> rowoffset;
	//std::vector<int> indice;
	//std::vector<double> value;
	
	int n();// return the number of rows;
	int m();// return the number of cols;
	int get_row_length(int r);
	int n_nonzero_elements();
	//sparsitypattern get_sparsity_pattern();

private:
	SparsityPattern pattern;
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
