g++ -c -o possion_equation.o test.cpp -I/usr/local/dealii-8.1.0/base/include -I/usr/local/dealii-8.1.0/lac/include -I/usr/local/dealii-8.1.0/include -I/usr/local/dealii-8.1.0/include/deal.II  -g -O2 -D__SERIALIZATION__ -DMULTITHREAD -pthread -fPIC -std=c++11 -ftemplate-depth-256 -I/usr/local/AFEPack/library/include
g++ -o main  possion_equation.o -L/usr/local/AFEPack/library/lib -lAFEPack -ldeal_II -L/usr/local/dealii-8.1.0/lib  -lm -ldl -L /home/chengyuliu/SparseMatrix/SparseMatrix-main/lib/ -I /home/chengyuliu/SparseMatrix/SparseMatrix-main/include/ -lsparse

Be careful about the path of the CXX lib and include path. Compiling in different environment might need to edit the Makefile;
