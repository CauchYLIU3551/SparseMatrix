#include<CGsolver.h>

CGsolver::CGsolver()
{
	sparsematrix A;
}

CGsolver::CGsolver(sparsematrix a)
{
	A = a;
}

// Following are a series of mathematical computing functions;
// I can try to put all the mathematical functions into One Math .h file!
// This is a function that help us to compute vector a minus vector b;
std::vector<double> vector_minus(std::vector<double> a, std::vector<double>b)
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

std::vector<double> vector_add2(std::vector<double> a, std::vector<double>b)
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

double frobenis(std::vector<double> a)
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

double trans_multiply(std::vector<double> a, std::vector<double>b)
{
	int len = a.size();
	int len2 = b.size();
	if (len != len2)
	{
		std::cout << "Multiply ERROR: The two matrixes are not in the same dimension!\n";
	}
	else
	{
		int sum = 0;
		for (int i = 0; i < len; i++)
		{
			sum += a[i] * b[i];
		}
		return sum;
	}
}
std::vector<double> num_multiply(double arg1, std::vector<double>b)
{
	//std::vector<double> tmp;
	for (int i = 0; i < b.size(); i++)
	{
		b[i] = arg1 * b[i];
	}
	return b;
}

std::vector<double> CGsolver::solve(std::vector<double> x0, std::vector<double> b,double eps)
{
	std::vector<double> r, r0, p;
	r0 = vector_minus(b, A.multiply(x0));
	p = r0;
	int k = 0;
	int Max = b.size();
	double alpha, beta;

	while (k <= Max)
	{
		std::cout << "this is the k times test!" << k << "\n";
		std::cout << "the frobenis error is " << frobenis(r0);
		if (frobenis(r0) <eps)
		{
			return x0;
		}
		else
		{
			//here is a little problem that we can get 1 correct ans at the first step, but lose the control after that.
			alpha=trans_multiply(r0, r0)/trans_multiply(p,A.multiply(p));
			//std::cout << "Check point 1\n";
			x0 = vector_add2(x0, num_multiply(alpha, p));
			//std::cout << "Check point 2\n";
			r = vector_minus(r0, num_multiply(alpha, A.multiply(p)));
			//std::cout << "Check point 3\n";
			beta = trans_multiply(r, r) / trans_multiply(r0, r0);
			p = vector_add2(r,num_multiply(beta, p));
			r0 = r;
			for (int i = 0; i < x0.size(); i++)
			{
				std::cout << x0[i] << " ";
			}
			std::cout << std::endl;
		}
		k++;
	}

	return x0;
}