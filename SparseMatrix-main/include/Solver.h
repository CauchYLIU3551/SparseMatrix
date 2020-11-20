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
	std::vector<double> GaussSeidel(std::vector<double> b, double eps, int num);
	std::vector<double> GaussSeidel(std::vector<std::vector<double>> IA,std::vector<double> b, double eps=1.0e-08, int num=10)

private:
	sparsematrix A;
};

#endif
