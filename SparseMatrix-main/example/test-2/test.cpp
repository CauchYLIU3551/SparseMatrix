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
#include <sparse/CGsolver.h>


#include <typeinfo>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>

#include <base/exceptions.h>
#include <lac/vector.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>

#include <time.h>


#define PI (4.0*atan(1.0))


double u(const double *);
double f(const double *);

// 
// Try to edit the error function in the Functional.template.h instead of in different cpp file.
//
/*template <class value_type, int DIM>
value_type EditError(FEMFunction<value_type, DIM>& f, const Function<value_type>& f1, int algebric_accuracy)
{
        value_type error = 0;
        FEMSpace<value_type,DIM>& fem_space = f.femSpace();
        typename FEMSpace<value_type,DIM>::ElementIterator the_element = fem_space.beginElement();
        typename FEMSpace<value_type,DIM>::ElementIterator end_element = fem_space.endElement();
        for (;the_element != end_element;the_element ++) {
                double volume = the_element->templateElement().volume();
                const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
                std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
                int n_quadrature_point = quad_info.n_quadraturePoint();
                std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
                std::vector<double> f_value = f.value(q_point, *the_element);
                for (int l = 0;l < n_quadrature_point;l ++) {
                        double df_value = f1.value(q_point[l]) - f_value[l];
                        df_value = fabs(df_value);
                        if (df_value > error) error = df_value;
                }
        }
        return error;
};
*/

clock_t start, end, start2, end2, start3, end3;

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
  //  get the number of elements in FEMSpace and reinit every element in the 
  //  FEMSpace;
  //  Make every element is corresponding to 
  //  the related geometry element images the No.0 template element
  std::cout<< "number of element::"<<n_element<<"\n";

  fem_space.element().resize(n_element);
  for (int i = 0;i < n_element;i ++)
    fem_space.element(i).reinit(fem_space,i,0);

  // build the FEMSpace initial data structure.
  // build Dof and Dof Boundary Mark;
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

  // initialize a function to approximate the original function;
  FEMFunction<double,2> solution(fem_space), testsolution(fem_space);
  // using the testsolution to compute the error and output some image to help
  // understanding the AMGsolve
  
  Vector<double> right_hand_side, right2;
  // right2 is serve for testsolution;

  Operator::L2Discretize(&f, fem_space, right_hand_side, 4);
  Operator::L2Discretize(&f, fem_space, right2, 4);
  
  BoundaryFunction<double,2> boundary(BoundaryConditionInfo::DIRICHLET, 1, &u);
  BoundaryConditionAdmin<double,2> boundary_admin(fem_space);
  boundary_admin.add(boundary);
  boundary_admin.apply(stiff_matrix, solution, right_hand_side);
  // to apply testsolution 
  boundary_admin.apply(stiff_matrix, testsolution, right2);
  
  /*
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
*/

  SparseMatrixIterators::Iterator<double,true> i=stiff_matrix.begin();
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

  // output the sparsity pattern of the stiff_matrix;
  //
  //std::cout<<"The last element of row:::"<<row[row.size()-1]<<"\n"; 
  //std::cout<<"This is the rowlen of the 1st row"<<stiff_matrix.get_row_length(0)<<"\n";
  //std::ofstream out ("sparsity_pattern.1");
  //sparse.print_gnuplot (out);
  //std::ofstream out ("sparsity_pattern.2");
  //sparse.print (out);
  
  /////////////
  //following are the process to solve the equations!
  jacobisolver solve_2(A), solve_GS(A);
  std::vector<double> sol2,sol3,sol4,sol5,b;
  // here are commands to copy the elements from solution and right** to sol2
  // and b;
  
  for(int k=0;k<solution.size();k++)
  {
	  sol2.push_back(solution[k]);
	  sol3.push_back(solution[k]);
	  sol4.push_back(solution[k]);
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
  
  // test my DIY Solver and check the result:::
  //
  //
  sol2=solve_2.solve(sol2,b,1.0e-08);
  //sol2.writeOpenDXData("u2.dx");
  CGsolver solve_3(A);
  
  //get the computing time of CG method;
  start3=clock();
  
  sol3=solve_3.solve(sol3,b,1.0e-08);
  
  end3=clock();
  double endtime3=(double)(end3-start3)/CLOCKS_PER_SEC;


  int numtmp=10;
  std::vector<double> errorn;


  /*
  // to compute the order of convergence of GS order;
  while(numtmp<100)
  {
	sol5=solve_GS.GaussSeidel(sol4,b,1.0e-08,numtmp);
  	for(int k=0;k<solution.size();k++)
  	{
        	  testsolution[k]=sol5[k];
  	}
	errorn.push_back(Functional::L2Error(testsolution, FunctionFunction<double>(&u), 3));
	numtmp+=1;
  }
  double order=0;
  for(int i=20;i<25;i++)
  {
	  order=log(errorn[i+1]/errorn[i])/log(errorn[i]/errorn[i-1]);
	  std::cout<<"The order of GS is ::"<<order<<"\n";
  }
// to compute the order of convergence of CG method;
  while(numtmp<100)
  {
        sol5=solve_3.solve(sol4,b,1.0e-08);
        for(int k=0;k<solution.size();k++)
        {
                  testsolution[k]=sol5[k];
        }
        errorn.push_back(Functional::L2Error(testsolution, FunctionFunction<double>(&u), 3));
        numtmp+=1;
  }
  for(int i=20;i<25;i++)
  {
          order=log(errorn[i+1]/errorn[i])/log(errorn[i]/errorn[i-1]);
          std::cout<<"The order of CG is ::"<<order<<"\n";
  }
*/

  // to compute the rate of convergence of CG method should edit in the 
  // CG function! instead of in this main function. Because CG just give 
  // back the final result without any mid-process result!

  start=clock();  
  
  sol4=solve_GS.GaussSeidel(sol4,b,1.0e-08,200);
  
  end=clock();
  double endtime=(double)(end-start)/CLOCKS_PER_SEC;
//  std::cout<<"The computing time of GS method is :::"<<1000*endtime<<" ms \n";

  for(int k=0;k<solution.size();k++)
  {
          testsolution[k]=sol4[k];
  }

  double err34=0;
  err34=Functional::L2Error(testsolution, FunctionFunction<double>(&u), 3);
  testsolution.writeOpenDXData("GS.dx");


  
  AMGSolver solver(stiff_matrix);
  
  start2=clock();
  
  solver.solve(solution, right_hand_side, 1.0e-08, 200);	
  
  end2=clock();
  double endtime2=(double)(end2-start2)/CLOCKS_PER_SEC;
//  std::cout<<"The computing time of AMG method is :::"<<1000*endtime<<" ms \n";

//  std::cout<<"The computing time of GS method is :::"<<1000*endtime<<" ms \n";




  /*
  std::freopen("sol3","w",stdout);
  for(int k=0;k<sol3.size();k++)
  {
	  std::cout<<sol3[k]<<" ";
  }
  fclose(stdout);
  */

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

  // Following steps are compute the infinite error between the sol2 sol3 and
  // the solution;
/*  double err=0,max=0;
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

  double err2=0,max2=0;
  for(int k=0;k<sol3.size();k++)
  {
          if(err2<abs(sol3[k]-solution[k]))
          {
                  err2=abs(sol3[k]-solution[k]);
          }
          if(max2<solution[k])
          {
                  max2=solution[k];
          }
  }
  err2=err2/max2;

  std::cout<<"The error between sol3 and the solution::::"<<err2<<std::endl;
*/;

  double testError=0;
/*
  FEMSpace<double, 2>& fem = solution.femSpace();
  FEMSpace<double, 2>::ElementIterator the_element = fem.beginElement();
  FEMSpace<double, 2>::ElementIterator end_element = fem.endElement();
  
  int flag=0;
  for (;the_element != end_element;the_element ++) {
        double volume = the_element->templateElement().volume();
	//std::cout<<"volume is ::"<<volume<<"\n";
        const QuadratureInfo<2>& quad_info = the_element->findQuadratureInfo(3);
        std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
        int n_quadrature_point = quad_info.n_quadraturePoint();
        std::vector<double> f_value = solution.value(the_element->local_to_global(quad_info.quadraturePoint()), *the_element);
        if (flag!=1) std::cout<<"n_quadrature_point::"<<n_quadrature_point<<"\n";
	for (int l = 0;l < n_quadrature_point;l ++) {
		if(flag!=1)
		{
			std::cout<<"q_point[l]:"<<the_element->local_to_global(quad_info.quadraturePoint())[0]<<"\n";
			std::cout<<"q_point[l].size:"<<the_element->local_to_global(quad_info.quadraturePoint()).size()<<"\n";
			std::cout<<"u[q_point[l]]:"<<FunctionFunction<double>(&u).value(the_element->local_to_global(quad_info.quadraturePoint())[l])<<"\n";
			std::cout<<"f_value[l]"<<f_value[l]<<"\n";
			std::cout<<"l::"<<l<<"\n";
		}	
		double df_value = FunctionFunction<double>(&u).value(the_element->local_to_global(quad_info.quadraturePoint())[l]) - f_value[l];
                
	       	df_value = fabs(df_value);
        
	      	if (df_value > testError) testError = df_value;
        }
	
	flag=1;
  }

*/

  solution.writeOpenDXData("u.dx");
//  double testerr = EditError(solution, FunctionFunction<double>(&u),3);
//  it is useless to add the Functional namespace function into this 
  double error = Functional::L2Error(solution, FunctionFunction<double>(&u), 3);
  double error2 = Functional::L0Error(solution, FunctionFunction<double>(&u), 3);

  // L2Error is defined in Functional.h and .template.h; The FunctionFunction is
  // defined in Miscellaneous.h; 
  // L2Error(FEMFunction<value_type, DIM>& f, const Function<value_type>& f1, int algebric_accuracy)
  std::cerr << "\nL2 error = " << error << std::endl;

  for(int k=0;k<solution.size();k++)
  {
          solution[k]=sol3[k];
  }
  
  error2 = Functional::L2Error(solution, FunctionFunction<double>(&u), 3);

  std::cout<<"L2 error of GS method is :"<<err34<<"\n";
  std::cerr << "\nL2 error of CGsolver= " << error2 << std::endl;
  
  std::cout<<"The computing time of AMG method is :::"<<1000*endtime2<<" ms \n";
  std::cout<<"The computing time of GS method is :::"<<1000*endtime<<" ms \n";
 // std::cout<<"The computing time of CG method is :::"<<1000*endtime3<<" ms \n";

  // Have a test in FunctionFunction<double>(&u);i
  //
  /*
  double L0err = 0;
  FEMSpace<double,DIM>& fem_space = solution.femSpace();
  FEMSpace<double,DIM>::ElementIterator the_element = fem_space.beginElement();
  FEMSpace<double,DIM>::ElementIterator end_element = fem_space.endElement();
  for (;the_element != end_element;the_element ++) {
      double volume = the_element->templateElement().volume();
      const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(algebric_accuracy);
      std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
      int n_quadrature_point = quad_info.n_quadraturePoint();
      std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
      std::vector<double> f_value = f.value(q_point, *the_element);
      for (int l = 0;l < n_quadrature_point;l ++) {
             double df_value = f1.value(q_point[l]) - f_value[l];
             df_value = fabs(df_value);
             if (df_value > error) error = df_value;
      }
  }
  */


  //A.showmatrix();

  return 0;
};

double u(const double * p)
{
  return sin(10*PI*p[0]) * sin(20*PI*p[1]);
};

double f(const double * p)
{
  return 500*PI*PI*u(p);
};

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
