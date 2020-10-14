#pragma once
#ifndef _JacobiSolver_h_
#define _JacobiSolver_h_

#include<iostream>
#include<vector>
#include<stdio.h>
#include<sparsematrix.h>

class jacobisolver
{
public:
	typedef sparsematrix matrix;
private:
	jacobisolver();//default constructor;
	jacobisolver(matrix A);
	void solve(std::vector<double> b);
	void GaussSidel();// unfinished function and it will be completed at next step.
};


#endif // !_JacobiSolver_h_
