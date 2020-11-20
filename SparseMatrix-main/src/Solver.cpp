#include<sparse/Solver.h>

Solver::Solver()
{
	sparsematrix a;
	A=a;
};

Solver::Solver(sparsematrix A): A(a) {};



std::vector<double> solve(std::vector<double> b, double eps)
{
	int n=b.size();
	int m=A.m();
	if(n!=m)
	{
		std::cout << "Function solve() ERROR: Please check the dimension of Matrix A and the right side vector b!";
	}
	else
	{
		std::vector<double> x(n,0);
		// Try to define the projection matrix;
		std::vector<std::vector<double>> I12,I21;
		// A2h=I12*A*I21;
		//
	}
}

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


std::vector<double> Solver::GaussSeidel(std::vector<double> b, double eps=1.0e-08, int num=10)
{
        std::vector<double> D,value;
        D=A.get_diag();
        std::vector<int>offset,indice;
        offset=A.get_offset();
        indice=A.get_indice();
        value=A.get_value();

        int row = A.m();
	std::vector<double> x(row,0);
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

std::vector<double> Solver::GaussSeidel(std::vector<std::vector<double>> IA,std::vector<double> b, double eps=1.0e-08, int num=10)
{
        std::vector<double> D,value;
	int n=A[0].size();
	std::vector<double> x(n,0); 
        do
        {
                x0=x;
                for(int i=0;i<n;i++)
                {
                        double tmp=0;
                        //COMPUTE TMP1 TMP2;
                        for (int j=0;i<n;i++)
                        {
				if(j!=i)
				{
					tmp+=A[i][j]*x[j];
				}
                        }
                        x[i]=(b[i]-tmp)/A[i][i];
                }
                num--;
                //x0=x;
                //x = vector_add(DI.multiply(M.multiply(x0)), DI.multiply(b));
        } while (num>0&&frobenis(x0,x)>eps);



//      std::cout << "The answer is that:::\n";
//      for (int i = 0; i < x.size(); i++)
//      {
//              std::cout << x[i] << " ";
//      }
//      std::cout << "\n";
        return x;

}

