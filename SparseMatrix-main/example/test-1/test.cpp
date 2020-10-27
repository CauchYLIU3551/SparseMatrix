////////////////////////////////////////////////////////////////////////////////////////////
// main1.cpp :
//

#include <iostream>
#include <fstream>

#include <AFEPack/AMGSolver.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>
#include <AFEPack/EasyMesh.h>
#include <sparse/sparsematrix.h>
#include <sparse/JacobiSolver.h>

#define PI (4.0*atan(1.0))

double u(const double *);
double f(const double *);

int main(int argc, char * argv[])
{
  EasyMesh mesh;
  mesh.readData(argv[1]);

  TemplateGeometry<2>	triangle_template_geometry;
  triangle_template_geometry.readData("triangle.tmp_geo");
  CoordTransform<2,2>	triangle_coord_transform;
  triangle_coord_transform.readData("triangle.crd_trs");
  TemplateDOF<2>	triangle_template_dof(triangle_template_geometry);
  triangle_template_dof.readData("triangle.1.tmp_dof");
  BasisFunctionAdmin<double,2,2> triangle_basis_function(triangle_template_dof);
  triangle_basis_function.readData("triangle.1.bas_fun");

  std::vector<TemplateElement<double,2,2> > template_element(1);
  template_element[0].reinit(triangle_template_geometry,
			     triangle_template_dof,
			     triangle_coord_transform,
			     triangle_basis_function);

  //using FEMSpace to define a fem_space object with the data from mesh and 
  //template_element which are defined by the above commands;
  FEMSpace<double,2> fem_space(mesh, template_element);
	
  int n_element = mesh.n_geometry(2);
  fem_space.element().resize(n_element);
  for (int i = 0;i < n_element;i ++)
    fem_space.element(i).reinit(fem_space,i,0);

  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  // Maybe I can have a try to replace the sparse matrix in StiffMatrix class
  // instead of in the process of solving the equations.
  // StiffMatrix class inherits BilinearOperator class; while BilinearOperator
  // class inherits SparseMatrix
  
  StiffMatrix<2,double> stiff_matrix(fem_space);
  stiff_matrix.algebricAccuracy() = 4;
  stiff_matrix.build();

  FEMFunction<double,2> solution(fem_space);
  Vector<double> right_hand_side;
  Operator::L2Discretize(&f, fem_space, right_hand_side, 4);

  BoundaryFunction<double,2> boundary(BoundaryConditionInfo::DIRICHLET, 1, &u);
  BoundaryConditionAdmin<double,2> boundary_admin(fem_space);
  boundary_admin.add(boundary);
  boundary_admin.apply(stiff_matrix, solution, right_hand_side);

  std::cout<<"this is the info of the stiff matrix\n";
  std::cout<<stiff_matrix.n()<<"\n";
  std::cout<<stiff_matrix.m()<<"\n";
  SparsityPattern sparse;
  sparse.copy_from(stiff_matrix.getSparsityPattern());
  std::ofstream out2("sparse_matrix");
  stiff_matrix.print(out2);

  SparseMatrixIterators::Iterator<double,true> i=stiff_matrix.begin();
  std::cout<<i->row()<<"\n";
  std::cout<<i->column()<<"\n";
  std::cout<<i->value()<<"\n";
  std::cout<<"Finish the output of the info!\n";

  const SparseMatrixIterators::Iterator<double,true> e=stiff_matrix.end();
  
  std::vector<int> row(stiff_matrix.n()+1),col;
  row[0]=0;
  int num1=stiff_matrix.n_nonzero_elements();//get the number of nonzero element in stiff_matrix;
  //std::vector<double> col(num1),val(num1); 
  std::vector<double> val; 
  // Produce the col_indice and value vector with the same size as the num1;
  // Save the col_indice and value in these two vectors to produce a sparse
  // matrix defined by myself.
  // 
  while (i!=e)
  {
	  col.push_back(i->column());
	  val.push_back(i->value());
	  ++i;
  }
  for(int j=1;j<row.size();j++)
  {
	  row[j]=stiff_matrix.get_row_length(j-1)+row[j-1];
  }
  sparsematrix A(row,col,val);

  //std::cout<<"The last element of row:::"<<row[row.size()-1]<<"\n"; 
  //std::cout<<"This is the rowlen of the 1st row"<<stiff_matrix.get_row_length(0)<<"\n";
  //std::ofstream out ("sparsity_pattern.1");
  //sparse.print_gnuplot (out);
  //std::ofstream out ("sparsity_pattern.2");
  //sparse.print (out);
  
  /////////////
  //following are the process to solve the equations!
  jacobisolver solve_2(A);
  std::vector<double> sol2,b;
  // here are commands to copy the elements from solution and right** to sol2
  // and b;
  
  for(int k=0;k<solution.size();k++)
  {
	  sol2.push_back(solution[k]);
	  b.push_back(right_hand_side[k]);
  }

  //here are some problems that the solution and right_hand_side are 
  //datatype dealii::Vector<double,2>;
  

  //following are commands to print all the elements in solution;
  //std::cout<<"this is the size of the solution:::"<<solution.size()<<"\n"; 
  //std::ofstream out("solution");
  //solution.print(out);


  //std::cout<<"this is the size of the right_hand_side:::"<<right_hand_side.size()<<"\n"; 

  //std::ofstream out("right_hand_side");
  //right_hand_side.print(out);
  sol2=solve_2.solve(sol2,b,1.0e-08);
  //sol2.writeOpenDXData("u2.dx");
  
  AMGSolver solver(stiff_matrix);
  solver.solve(solution, right_hand_side, 1.0e-08, 200);	
  
  std::freopen("sol2","w",stdout);
  for(int k=0;k<sol2.size();k++)
  {
	  std::cout<<sol2[k]<<" ";
  }
  fclose(stdout);
  
  // Notes: solution inherits from dealii::Vector<double>, So it can be computed
  // as a Vector<double> object in AMGSolver::solver; But in error compute it is
  // used as a FEMFunction<double,2> object, So I need to use some commands to 
  // convert the result vector<double>sol2 into AFEPack::FEMFunction<double,2>
  // solution. Then use the solution to complete the error computing.
  //
  // Unfinished target: understanding the commands in FEMFunction.h and 
  // AMGSolver.h. After that replace the classes from deal.ii by some DIY class.


  //std::cout<<"this is the size of the solution:::"<<solution.size()<<"\n";
  //std::ofstream out("solution_after");
  //solution.print(out);
  
  // The class FunctionFunction is defined in Miscellaneous.h

  double err=0,max=0;
  for(int k=0;k<sol2.size();k++)
  {
	  if(err<abs(sol2[k]-solution[k]))
	  {
		  err=abs(sol2[k]-solution[k]);
	  }
	  if(max<solution[k])
	  {
		  max=solution[k];
	  }
  }
  err=err/max;
  std::cout<<"The error between sol2 and the solution::::"<<err<<std::endl;

  solution.writeOpenDXData("u.dx");
  double error = Functional::L2Error(solution, FunctionFunction<double>(&u), 3);
  // L2Error is defined in Functional.h and .template.h; The FunctionFunction is
  // defined in Miscellaneous.h; 
  // L2Error(FEMFunction<value_type, DIM>& f, const Function<value_type>& f1, int algebric_accuracy)
  std::cerr << "\nL2 error = " << error << std::endl;

  //A.showmatrix();

  return 0;
};

double u(const double * p)
{
  return sin(PI*p[0]) * sin(2*PI*p[1]);
};

double f(const double * p)
{
  return 5*PI*PI*u(p);
};

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
