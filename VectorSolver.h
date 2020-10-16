#pragma once
#ifndef _VectorSolver_h_ 
#define _VectorSolver_h_
#include<iostream>
#include<cstdio>
#include<math.h>
#include<vector>
#include<sparsematrix.h>

class vectorsolver
{
public:
	static std::vector<double> add(std::vector<double> a, std::vector<double>b);
	static std::vector<double> minus(std::vector<double> a, std::vector<double>b);
	double multiply(std::vector<double> a, std::vector<double>b);// This function compute aT times b;
	std::vector<double> multiply(double arg1, std::vector<double>b); //this function compute number arg1 times the matrix b;
	static std::vector<double> transpose();

	double frobenis(std::vector<double> a, std::vector<double>b);
};
#endif // !_VectorSolver_h_ 
