#ifndef _Solver_h_
#define _Solver_h_

#include<iostream>
#include<vector>
#include<stdio.h>
#include<sparse/sparsematrix.h>


class Solver
{
public:
	Solver();
	Solver(sparsematrix A);
	std::vector<double> solve(std::vector<doube> b, double eps);
	//
	// In fact, in AMGSolver, GaussSidel is used in sparsematrix, and the 
	// object is the projection matrix output from the function projection
	// What's more, can have a try to use the point or & to refer the
	// components. In this way, we do not need to define the return value;
	std::vector<double> GaussSeidel(std::vector<double> b, double eps, int num);
	std::vector<double> GaussSeidel(std::vector<std::vector<double>> IA,std::vector<double> b, double eps=1.0e-08, int num=10)

private:
	sparsematrix A;
};

#endif
