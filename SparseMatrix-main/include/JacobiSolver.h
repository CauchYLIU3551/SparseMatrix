#pragma once
#ifndef _JacobiSolver_h_
#define _JacobiSolver_h_

#include<iostream>
#include<vector>
#include<stdio.h>
#include<sparse/sparsematrix.h>

class jacobisolver
{
public:
	typedef sparsematrix matrix;
	jacobisolver();//default constructor;
	jacobisolver(matrix a);
	std::vector<double> solve(std::vector<double>x,std::vector<double> b, double eps);
	std::vector<double> GaussSeidel(std::vector<double>x,std::vector<double> b, double eps, int num=10);
	// GaussSidel Solver and the default iteration number=10.
private:
	matrix A;
	//double eps;
};


#endif // !_JacobiSolver_h_
