#include<VectorSolver.h>

std::vector<double> vectorsolver::add(std::vector<double> a, std::vector<double>b)
{
	if (a.size() != b.size())
	{
		std::cout << "Vector Add ERROR:: The two matrix are not in the same dimension!\n";
	}
	else
	{
		int len = a.size();
		std::vector<double> c(len);
		for (int i = 0; i < a.size(); i++)
		{
			c[i] = a[i] + b[i];
		}
		return c;
	}
}
std::vector<double> vectorsolver::minus(std::vector<double> a, std::vector<double>b)
{
	if (a.size() != b.size())
	{
		std::cout << "Vector Minus ERROR:: The two matrix are not in the same dimension!\n";
	}
	else
	{
		int len = a.size();
		std::vector<double> c(len);
		for (int i = 0; i < a.size(); i++)
		{
			c[i] = a[i] - b[i];
		}
		return c;
	}
}

double vectorsolver::multiply(std::vector<double> a, std::vector<double>b)
{
	int len = a.size();
	int len2 = b.size();
	if (len != len2)
	{
		std::cout << "Multiply ERROR: The two matrixes are not in the same dimension!\n";
	}
	else
	{
		double sum = 0;
		for (int i = 0; i < len; i++)
		{
			sum += a[i] * b[i];
		}
		return sum;
	}
}

std::vector<double> vectorsolver::multiply(double arg1, std::vector<double>b)
{
	for (int i = 0; i < b.size(); i++)
	{
		b[i] = arg1 * b[i];
	}
	return b;
}

std::vector<double> vectorsolver::transpose()
{}

double vectorsolver::frobenis(std::vector<double> a, std::vector<double>b)
{
	double F = 0;
	int len = a.size();
	for (int i = 0; i < len; i++)
	{
		F += pow(a[i], 2);
	}
	F = sqrt(F);
	return F;

}