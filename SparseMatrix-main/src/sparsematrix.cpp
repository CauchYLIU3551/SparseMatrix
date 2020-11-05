#include<sparse/sparsematrix.h>


// default construction function without input;
sparsematrix::sparsematrix() :rowoffset({ 0 ,0 }), indice({ 0 }), value({ 0 })
{
	col = max();
	row = 1;
};

// construction function with input;
sparsematrix::sparsematrix(std::vector<int> n, std::vector<int> j, std::vector<double> v) :rowoffset(n), indice(j), value(v)
{
	col = max();
	row = rowoffset.size() - 1;
}

// construction function with defined row and col;
sparsematrix::sparsematrix(std::vector<int> n, std::vector<int> j, std::vector<double> v, int ind_row, int ind_col) :rowoffset(n), indice(j), value(v)
{
	row = ind_row;
	col = ind_col;
}

// get the transpose matrix of the sparse matrix.
sparsematrix sparsematrix::transpose()
{

	std::vector<double>v;
	std::vector<int>n, c;
	n.push_back(0);
	int index = 0;

	int num = 0; // compute the number of the values sit in the jth col, inorder to compute the transpose rowoffsets;

	for (int j = 0; j < col; j++)
	{

		for (int i = 0; i < indice.size(); i++)
		{
			if (indice[i] == j) //if the indice start from 1, here will change to indice[k]-1
			{
				v.push_back(value[i]);
				num += 1;
				for (int k = 0; k < col; k++)
				{
					if (rowoffset[k] <= i && i < rowoffset[k + 1])
					{
						c.push_back(k);
						break;
					}
				}
			}

			// find the related index of the value in row and input it into the transpose col vector;

		}
		n.push_back(num);
	}

	sparsematrix A(n, c, v);

	return A;
}


// diag() function is used to get the diagnal element of the matrix and output a Diagnal matrix D
// we can divide the matrix in this way: A=D-L-U; D is diagnal matrix, L is lower-triangular matrix, U is upper-triangular matrix;
sparsematrix sparsematrix::diag()
{
	std::vector<int> r(1, 0), Dindice;
	std::vector<double> Dvalue;
	int flag = 0;

	for (int i = 0; i < row; i++)
	{
		//std::cout << "this is the " << i << "th times test!!!\n";
		for (int j = rowoffset[i]; j < rowoffset[i + 1]; j++)
		{
			//std::cout << "this is the" << indice[j] << "col element!\n";
			if (i == indice[j]) //if the indice start from 1, here will change to indice[k]-1
			{
				Dindice.push_back(indice[j]);
				Dvalue.push_back(value[j]);
				r.push_back(i + 1);
				flag = 1;
				break;
			}

		}
		if (flag != 1)
		{
			std::cout<<"In the "<<i+1<<"th row the element in the diag is zero\n";
			r.push_back(r[i]);
		}

	}

	sparsematrix D(r, Dindice, Dvalue);

	return D;
}

//this is a function to get the inverse matrix of the diagnal matrix D of A 
sparsematrix sparsematrix::diag_inverse()
{
	sparsematrix DI = diag();
	for (int i = 0; i < DI.value.size(); i++)
	{
		DI.value[i] = 1.0 / DI.value[i];
	}
	return DI;
}

// Need to finish!!!!!!
// Compute the lower-triangular L;
sparsematrix sparsematrix::lowtri()
{
	std::vector<int> r(1, 0), Lindice;
	std::vector<double> Lvalue;
	for (int i = 0; i < row; i++)
	{
		int num = 0;
		for (int j = rowoffset[i]; j < rowoffset[i + 1]; j++)
		{
			if (indice[j] < i)
			{
				num++;
				Lindice.push_back(indice[j]);
				Lvalue.push_back(value[j]);
			}
			else
				break;
		}
		r.push_back(r[i] + num);
	}
	sparsematrix L(r, Lindice, Lvalue, row, col); // here I use the index of row and col index to define the sparse matrix otherwise it will lose the col info.
	return L;
}

// Compute the upper-triangular U;
sparsematrix sparsematrix::uppertri()
{
	std::vector<int> r(1, 0), Uindice;
	std::vector<double> Uvalue;
	for (int i = 0; i < row; i++)
	{
		int num = 0;
		for (int j = rowoffset[i]; j < rowoffset[i + 1]; j++)
		{
			if (indice[j] > i)
			{
				num++;
				Uindice.push_back(indice[j]);
				Uvalue.push_back(value[j]);
			}
		}
		r.push_back(r[i] + num);
	}
	sparsematrix U(r, Uindice, Uvalue, row, col); // here I use the index of row and col index to define the sparse matrix otherwise it will lose the col info.
	return U;
}

sparsematrix sparsematrix::inverse_num()
{
	std::vector<double> Mvalue;
	for (int i = 0; i < value.size(); i++)
	{
		Mvalue.push_back(-1 * value[i]);
	}

	sparsematrix M(rowoffset, indice, Mvalue);

	return M;
}

////////////////////////need to finish !
sparsematrix sparsematrix::inverse()
{
	sparsematrix A;
	return A;
}

sparsematrix sparsematrix::add(sparsematrix A, sparsematrix B)
{
	std::vector<double> Cvalue;
	std::vector<int>Cind, Crow;
	Crow.push_back(0);

	if (A.row == B.row && A.col == B.col)
	{
		for (int i = 0; i < A.row; i++)
		{
			//std::cout << "this is No." << i << " test of this !!!!\n";
			int len, len1, len2;
			len1 = A.rowoffset[i + 1] - A.rowoffset[i];
			len2 = B.rowoffset[i + 1] - B.rowoffset[i];

			len = len1 + len2;
			//std::cout << len << std::endl;
			//std::cout << "this is No." << i << " test len1="<<len1<<"!!!!\n";
			//std::cout << "this is No." << i << " test len2="<<len2<<"!!!!\n";
			//std::cout << A.rowoffset[i] << std::endl;


			int j = 0, k = 0;

			while (j < len1 && k < len2)//there is a problem!
			{
				if (A.indice[A.rowoffset[i] + j] < B.indice[B.rowoffset[i] + k])
				{
					Cind.push_back(A.indice[A.rowoffset[i] + j]);
					Cvalue.push_back(A.value[A.rowoffset[i] + j]);
					j++;
				}
				else if (B.indice[B.rowoffset[i] + k] < A.indice[A.rowoffset[i] + j])
				{
					Cind.push_back(B.indice[B.rowoffset[i] + k]);
					Cvalue.push_back(B.value[B.rowoffset[i] + k]);
					k++;
				}
				else
				{
					if (A.value[A.rowoffset[i] + j] + B.value[B.rowoffset[i] + k] != 0)
					{
						Cind.push_back(A.indice[A.rowoffset[i] + j]);
						Cvalue.push_back(A.value[A.rowoffset[i] + j] + B.value[B.rowoffset[i] + k]);
						len = len - 1;
						//std::cout << "ATTENTION!\n";
						//std::cout << A.indice[A.rowoffset[i] + j] << std::endl;
						//std::cout << B.indice[B.rowoffset[i] + k] << std::endl;
						j++;
						k++;
					}
					else
					{
						len = len - 2;
						j++;
						k++;
					}

				}
			}
			while (j < len1)
			{
				Cind.push_back(A.indice[A.rowoffset[i] + j]);
				Cvalue.push_back(A.value[A.rowoffset[i] + j]);
				j++;
			}
			while (k < len2)
			{
				Cind.push_back(B.indice[B.rowoffset[i] + k]);
				Cvalue.push_back(B.value[B.rowoffset[i] + k]);
				k++;
			}
			//std::cout <<"After that"<< len << std::endl;
			Crow.push_back(Crow[i] + len);
		}

		/*
		for (int i = 0; i < Crow.size(); i++)
		{
			std::cout << Crow[i] << " ";
		}
		std::cout << std::endl;
		for (int i = 0; i < Cind.size(); i++)
		{
			std::cout << Cind[i] << " ";
		}
		std::cout << std::endl;
		for (int i = 0; i < Cvalue.size(); i++)
		{
			std::cout << Cvalue[i] << " ";
		}
		std::cout << std::endl;
		*/

		sparsematrix C(Crow, Cind, Cvalue);

		return C;
	}
	else
	{
		std::cout << "Function add() ERROR:The two matrixes are not in same size! Please check it!\n";
	}
}



// compute the value of the sparsematrix times a n-dimensional vector x.
std::vector<double> sparsematrix::multiply(std::vector<double> x)
{
	std::vector<double> b(row, 0);
	//std::cout << rowoffset[row]<<std::endl;
	int num = 0;
	for (int i = 0; i < row; i++)
	{
		for (int j = rowoffset[i]; j < rowoffset[i + 1]; j++)
		{
			//num += 1;
			//std::cout << "this the No." << num << " times \n";
			//std::cout << "the operation is ::::" << value[j] << "times" << x[indice[j]] << std::endl;
			b[i] = b[i] + value[j] * x[indice[j]]; //if the indice start from 1, here will change to indice[k]-1
		}
	}

	return b;
}

// print the value of the whole matrix;
void sparsematrix::showmatrix()
{
	std::cout << "this is the sparse matrix:::" << std::endl;
	std::cout << "the col is :::" << col << std::endl;
	std::cout << "the row is :::" << row << std::endl;
	if (col == 1)
		std::cout << value[0] << std::endl;
	else

		for (int j = 0; j < rowoffset.size() - 1; j++)
		{
			int len = rowoffset[j + 1] - rowoffset[j];
			/*
			for (int i = 0; i < len; i++)
			{
				int num = indice[i + rowoffset[j]];
				while (num > 0)
				{
					num = num - 1;
				}
			}*/
			int k = rowoffset[j];
			int num = 0;
			for (int i = 0; i < col; i++)
			{
				if (k < rowoffset[row])
				{
					if (i == indice[k] && num < len) //if the indice start from 1, here will change to indice[k]-1
					{
						std::cout << value[k] << " ";
						k++;
						num++;
					}
					else
						std::cout << "0 ";
				}
				else
					std::cout << "0 ";
			}
			std::cout << std::endl;
		}

	std::cout << std::endl;
}

// A special function, that helps the class to get the num of col;
int sparsematrix::max()
{
	int num = 0;
	for (int i = 0; i < indice.size(); i++)
	{
		if (indice[i] >= num)
			num = indice[i];
	}
	return num + 1;
};

std::vector<int> sparsematrix::get_offset()
{
	/*
	for (int i = 0; i < rowoffset.size(); i++)
	{
		std::cout << rowoffset[i] << " ";
	}
	std::cout << std::endl;
	*/
	return rowoffset;
}

std::vector<int> sparsematrix::get_indice()
{
	/*
	for (int i = 0; i < indice.size(); i++)
	{
		std::cout << indice[i] << " ";
	}
	std::cout << std::endl;
	*/
	return indice;
}

std::vector<double> sparsematrix::get_value()
{
	/*
	for (int i = 0; i < value.size(); i++)
	{
		std::cout << value[i] << " ";
	}
	std::cout << std::endl;
	*/
	return value;
}

std::vector<double> sparsematrix::get_diag()
{
	if(row==col)
	{
		std::vector<double> D;
		for(int i=0;i<row;i++)
		{
			for(int j=rowoffset[i];j<rowoffset[i+1];j++)
			{
				if(i==indice[j])
				{
					D.push_back(value[j]);
				}
			}
		}
		return D;
	}
	else
	{
		std::cout << "Function get_diag() ERROR:The matrix is not square matrix! Please check it!\n";
	}
}

int sparsematrix::show_row()
{
	return row;
}

int sparsematrix::n()
{
	return row;
}

int sparsematrix::m()
{
	return col;
}

int sparsematrix::get_row_length(int r)
{
	return rowoffset[r]-rowoffset[r-1];
}

int sparsematrix::n_nonzero_elements()
{
	return rowoffset[row];
}
