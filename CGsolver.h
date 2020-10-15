#pragma once
#ifndef _CGsolver_h_
#define _CGsolver_h_

#include<iostream>
#include<vector>
#include<stdio.h>
#include<sparsematrix.h>
#include<JacobiSolver.h>


class CGsolver
{
public:
	std::vector<double> solve(std::vector<double> x0,std::vector<double> b,double eps);
	CGsolver();
	CGsolver(sparsematrix a);
private:
	sparsematrix A;
};
#endif