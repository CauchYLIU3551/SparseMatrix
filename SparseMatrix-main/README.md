# SparseMatrix
This is a basic implement of CRS in C++

This is a simple practice in computing sparse matrix. It is unfinished yet.

In this sparsematrix.h, I use the CRS method to save the matrix and it will provide a series of operations of sparse matrix.
It can finish some simple computation such as :
1.multiply 2. add 3. divide A into L, U, D such that A=L+U+D(it will be helpful in Jacobi iteration) 4. show the values in sparsematrix 5. inverse (unfinised!) 

I have implement a simple version of Jacobi iteration method in the JacobiSolver.cpp;
Now I finish the CG method in CGsolver.cpp;

In this master branch, I add an example using AFEPack. I convert the AFEPack::StiffMatrix to sparsematrix and using the Jacobi and CG methods to solve the equations.

In the future, I hope I can achieve the CG or PCG to solve large scale sparse matrix to solve some 'Ax=b' problem with this .h file!

