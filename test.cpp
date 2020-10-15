// sparsematrix.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <sparsematrix.h>
#include <JacobiSolver.h>
#include <CGsolver.h>

int main()
{
    std::vector<int> n, c;
    std::vector<double> v, x = { 1,2,3,4 };//there is a little problem in multiply!!

    n = { 0,2,4,7,9 };
    c = { 0,1,1,2,0,2,3,1,3 };
    v = { 1,7,2,8,5,3,9,6,4 };
    sparsematrix C(n, c, v);

    n = { 0,1,3,4,6 };
    c = { 0,1,3,2,1,3 };
    v = { 1,2,6,3,5,4 };
    sparsematrix B(n, c, v);

    n = { 0,2,5,8,12,16,19 };
    c = { 0,4,0,1,5,1,2,3,0,2,3,4,1,3,4,5,1,4,5 };
    v = { 10,-2,3,9,3,7,8,7,3,8,7,5,8,9,9,13,4,2,-1 };
    sparsematrix A(n, c, v);

    /*
    A.showmatrix();
    A.show_indice();
    A.show_offset();
    A.show_value();
    */


    //std::cout << "this is the test of the getting diagnal function!\n";
    //A = A.diag();
    //A.showmatrix();
    //A.show_offset();
    //A.show_indice();
    //A.show_value();
    //B.showmatrix();
    //B.show_indice();
    //C.showmatrix();

    //std::cout << "this is the test of the add function!\n";

    //B=sparsematrix::add(B, C);

    //B.show_indice();
    //B.show_offset();
    //B.show_value();
    //B.showmatrix();

    //A.showmatrix();
    //B.showmatrix();


    C.showmatrix();
    C.show_indice();
    C.show_offset();
    C.show_value();

    /*
    std::cout << "After transpose!\n";
    A = C.transpose();
    A.showmatrix();
    A.show_indice();
    A.show_offset();
    A.show_value();
    */


    std::cout << "This is the multiply of the sparse matrix!!\n";

    x = C.multiply(x);
    for (int i = 0; i < x.size(); i++)
    {
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;
    //std::cout << "this is the test of the getting diagnal function!\n";
    //A = C.diag();
    //A.showmatrix();
    //A.show_offset();
    //A.show_indice();
    //A.show_value();
    //n = { 0,3,6,9 };
    //c = { 0,1,2,0,1,2,0,1,2 };
    //v = { 10,-1,-2,-1,10,-2,-1,-1,5 };
    //sparsematrix E(n, c, v);
    //E.showmatrix();


    // here is the test of jacobi iteration;
    //std::vector<double>b = { 72,83,42 };


    //jacobisolver sol(E);
    //sol.solve(b, 0.0001);

    n = { 0,2,5,7 };
    c = { 0,1,0,1,2,1,2 };
    v = { 2,-1,-1,2,-1,-1,2 };
    sparsematrix E(n, c, v);
    E.showmatrix();
    std::vector<double> x0(3, 0), b = {1,0,1.8};
    CGsolver sol2(E);
    x0 = sol2.solve(x0, b, 0.001);
    
    for (int i = 0; i < x0.size(); i++)
    {
        std::cout << x0[i] << " ";
    }
    std::cout << std::endl;


    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
