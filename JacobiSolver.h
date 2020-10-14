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
	jacobisolver();//default constructor;
	jacobisolver(matrix a);
	std::vector<double> solve(std::vector<double> b, double eps);
	void GaussSidel();// unfinished function and it will be completed at next step.
private:
	matrix A;
	//double eps;
};


#endif // !_JacobiSolver_h_
