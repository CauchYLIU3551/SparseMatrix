#include<iostream>
#include<vector>
#include<stdio.h>
#include<sparse/sparsematrix.h>
#include<sparse/JacobiSolver.h>

typedef jacobisolver::matrix matrix;

jacobisolver::jacobisolver()
{
	matrix a;
	A = a;
	//eps = 0.0001;
};

jacobisolver::jacobisolver(matrix a) :A(a) {}

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


std::vector<double> jacobisolver::GaussSeidel(std::vector<double>x,std::vector<double> b, double eps, int num) 
{
	std::vector<double> D,value;
	D=A.get_diag();
	std::vector<int>offset,indice;
	offset=A.get_offset();
	indice=A.get_indice();
	value=A.get_value();

        //sparsematrix M = sparsematrix::add(L, U);
        //M = M.inverse_num(); //get the -(L+U) to compute the iteration Matrix;

        //M.showmatrix();

        int row = A.m();
        std::vector<double> x0;
	
	int N=0;
        do
        {
		x0=x;
		for(int i=0;i<row;i++)
		{
			double tmp=0;
			//COMPUTE TMP1 TMP2;
			for (int j=offset[i];j<offset[i+1];j++)
			{
				if(i!=indice[j])
				{
					tmp+=value[j]*x[indice[j]];
				}
			}
			x[i]=(b[i]-tmp)/D[i];
		}
		N++;
                //x0=x;
                //x = vector_add(DI.multiply(M.multiply(x0)), DI.multiply(b));
        } while (N<num&&frobenis(x0,x)>eps);



//      std::cout << "The answer is that:::\n";
//      for (int i = 0; i < x.size(); i++)
//      {
//              std::cout << x[i] << " ";
//      }
//      std::cout << "\n";
        return x;

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

std::vector<double> jacobisolver::solve(std::vector<double>x,std::vector<double>b, double eps)
{
	sparsematrix DI, L, U;
	DI = A.diag_inverse();
	L = A.lowtri();
	U = A.uppertri();

	sparsematrix M = sparsematrix::add(L, U);
	M = M.inverse_num(); //get the -(L+U) to compute the iteration Matrix;

	//M.showmatrix();

	int dim = A.show_row();
	std::vector<double> x0;

	do
	{
		x0=x;
		x = vector_add(DI.multiply(M.multiply(x0)), DI.multiply(b));
	} while (frobenis(x, x0) > eps);

//	std::cout << "The answer is that:::\n";
//	for (int i = 0; i < x.size(); i++)
//	{
//		std::cout << x[i] << " ";
//	}
//	std::cout << "\n";
	return x;
}
