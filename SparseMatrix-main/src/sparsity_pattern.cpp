#include<sparse/sparsity_pattern.h>

#include<iostream>

/*
SparsityPattern::SparsityPattern ()
  	:
  	max_dim(0),
  	max_vec_len(0),
  	compressed(false),
 	store_diagonal_first_in_row(false)
{
  	reinit (0,0,0);
}

SparsityPattern::SparsityPattern(int m,
				 int n,
				 std::vector<unsigned int>&row_lengths)
	:
	max_dim(0),
	max_vec_len(0),
	store_diagonal_fisrt_in_row(m==n)
{
	reinit (m, n, row_lengths);
}
*/

Sparsitypattern::Sparsitypattern():rows(0),cols(0)
{}

Sparsitypattern::Sparsitypattern(int m, int n):rows(m),cols(n)
{}

