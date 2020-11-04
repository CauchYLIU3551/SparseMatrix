#ifndef __sparsity_pattern_h
#define __sparsity_pattern_h

//#include<sparse/sparsematrix.h>
#include<vector>

class SparsityPattern
{
public:
	SparsityPattern();
	SparsityPattern(int m, int n);

private:
	//int max_dim;
	int rows;
	int cols;
	//int max_vec_len;
	//unsigned int max_row_lenth;
	//bool compressed;
	//bool store_diagonal_first_in_row;
};

#endif
