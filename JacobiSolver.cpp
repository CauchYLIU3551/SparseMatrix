#include<iostream>
#include<vector>
#include<stdio.h>
#include<sparsematrix.h>
#include<JacobiSolver.h>

typedef jacobisolver::matrix matrix;

jacobisolver::jacobisolver()
{
	matrix a;
	A = a;
	//eps = 0.0001;
};

jacobisolver::jacobisolver(matrix a) :A(a) {}

void jacobisolver::GaussSidel() {}

double frobenis(std::vector<double> a, std::vector<double>b)
{
	if (a.size() != b.size())
	{
		std::cout << "The two vector are not in the same dimension! Please check it!\n";
	}
	else
	{
		double F = 0;
		int len = a.size();
		for (int i = 0; i < len; i++)
		{
			F += pow(a[i] - b[i], 2);
		}
		F = sqrt(F);
		return F;
	}
}

std::vector<double> vector_add(std::vector<double> a, std::vector<double> b)
{

	if (a.size() != b.size())
	{
		std::cout << "The two vector are not in the same dimension! Please check it!\n";
	}
	else
	{
		int len = a.size();
		for (int i = 0; i < len; i++)
		{
			a[i] = a[i] + b[i];
		}
		return a;
	}
}

void showvector(std::vector<double> A)
{
	for (int i = 0; i < A.size(); i++)
	{
		std::cout << A[i] << " ";
	}
	std::cout << "\n";
}

std::vector<double> jacobisolver::solve(std::vector<double>b, double eps)
{
	sparsematrix DI, L, U;
	DI = A.diag_inverse();
	L = A.lowtri();
	U = A.uppertri();

	sparsematrix M = sparsematrix::add(L, U);
	M = M.inverse_num(); //get the -(L+U) to compute the iteration Matrix;

	//M.showmatrix();

	int dim = A.show_row();
	std::vector<double> x(dim, 1), x0;

	do
	{
		x0 = x;
		x = vector_add(DI.multiply(M.multiply(x0)), DI.multiply(b));
	} while (frobenis(x, x0) > eps);

	std::cout << "The answer is that:::\n";
	for (int i = 0; i < x.size(); i++)
	{
		std::cout << x[i] << " ";
	}
	std::cout << "\n";
	return x;
}