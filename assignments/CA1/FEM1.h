

//Include files
//Data structures and solvers
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
//Mesh related classes
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
//Finite element implementation classes
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
//Standard C++ libraries
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace dealii;

template <int dim>
class FEM
{
 public:
  //Class functions
  FEM (unsigned int order,unsigned int problem); // Class constructor
  ~FEM(); //Class destructor

  double xi_at_node(unsigned int dealNode);

  double basis_function(unsigned int node, double xi);
  double basis_gradient(unsigned int node, double xi);

  //Solution steps
  void generate_mesh(unsigned int numberOfElements);
  void define_boundary_conds();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results();

  double l2norm_of_error();

  //Class objects
  Triangulation<dim>   triangulation; //mesh
  FESystem<dim>        fe; 	      //FE element
  DoFHandler<dim>      dof_handler;   //Connectivity matrices

  unsigned int	        quadRule;    //quadrature rule, i.e. number of quadrature points
  std::vector<double>	quad_points; //vector of Gauss quadrature points
  std::vector<double>	quad_weight; //vector of the quadrature point weights

  //Data structures
  SparsityPattern       sparsity_pattern; //Sparse matrix pattern
  SparseMatrix<double>  K;		 //Global stiffness (sparse) matrix
  Vector<double>        D, F; 		 //Global vectors - Solution vector (D) and Global force vector (F)
  std::vector<double>   nodeLocation;	 //Vector of the x-coordinate of nodes by global dof number
  std::map<unsigned int,double> boundary_values;	//Map of dirichlet boundary conditions
  double                basisFunctionOrder, prob, L, g1, g2;

  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};

// Class constructor for a vector field
template <int dim>
FEM<dim>::FEM(unsigned int order,unsigned int problem)
:
fe (FE_Q<dim>(order), dim),
  dof_handler (triangulation)
{
  basisFunctionOrder = order;
  if(problem == 1 || problem == 2){
    prob = problem;
  }
  else{
    std::cout << "Error: problem number should be 1 or 2.\n";
    exit(0);
  }

  for (unsigned int i=0; i<dim; ++i){
    nodal_solution_names.push_back("u");
    nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }
}

//Class destructor
template <int dim>
FEM<dim>::~FEM(){
  dof_handler.clear();
}

template <int dim>
double FEM<dim>::xi_at_node(unsigned int dealNode){
  double xi;

  if(dealNode == 0){
    xi = -1.;
  }
  else if(dealNode == 1){
    xi = 1.;
  }
  else if(dealNode <= basisFunctionOrder){
    if (basisFunctionOrder == 3)
      xi = (dealNode == 2 ? -1. : 1.) / sqrt(5.);
    else
      xi = -1. + 2.*(dealNode-1.)/basisFunctionOrder;
  }
  else{
    std::cout << "Error: you input node number "
	      << dealNode << " but there are only "
	      << basisFunctionOrder + 1 << " nodes in an element.\n";
    exit(0);
  }

  return xi;
}

//Define basis functions
template <int dim>
double FEM<dim>::basis_function(unsigned int node, double xi){


  double value = 1.; //Store the value of the basis function in this variable



  const double xi_node = xi_at_node(node);
  for (unsigned int other = 0; other <= basisFunctionOrder; ++other)
    if (other != node)
      value *= (xi - xi_at_node(other)) / (xi_node - xi_at_node(other));

  return value;
}

//Define basis function gradient
template <int dim>
double FEM<dim>::basis_gradient(unsigned int node, double xi){


  double value = 0.; //Store the value of the gradient of the basis function in this variable



  for (unsigned int omitted = 0; omitted <= basisFunctionOrder; ++omitted)
    if (omitted != node) {
      double term = 1. / (xi_at_node(node) - xi_at_node(omitted));
      for (unsigned int other = 0; other <= basisFunctionOrder; ++other)
        if (other != node && other != omitted)
          term *= (xi - xi_at_node(other)) /
                  (xi_at_node(node) - xi_at_node(other));
      value += term;
    }

  return value;
}

template <int dim>
void FEM<dim>::generate_mesh(unsigned int numberOfElements){

  //Define the limits of your domain
  L = 0.1;
  double x_min = 0.;
  double x_max = L;

  Point<dim,double> min(x_min),
    max(x_max);
  std::vector<unsigned int> meshDimensions (dim,numberOfElements);
  GridGenerator::subdivided_hyper_rectangle (triangulation, meshDimensions, min, max);
}

//Specify the Dirichlet boundary conditions
template <int dim>
void FEM<dim>::define_boundary_conds(){
  const unsigned int totalNodes = dof_handler.n_dofs(); //Total number of nodes




  for(unsigned int globalNode=0; globalNode<totalNodes; globalNode++){
    if(nodeLocation[globalNode] == 0){
      boundary_values[globalNode] = g1;
    }
    if(nodeLocation[globalNode] == L){
      if(prob == 1){
	boundary_values[globalNode] = g2;
      }
    }
  }

}

template <int dim>
void FEM<dim>::setup_system(){

  g1 = 0.; g2 = 0.001;

  //Let deal.II organize degrees of freedom
  dof_handler.distribute_dofs (fe);

  MappingQ1<dim,dim> mapping;
  std::vector< Point<dim,double> > dof_coords(dof_handler.n_dofs());
  nodeLocation.resize(dof_handler.n_dofs());
  DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
  for(unsigned int i=0; i<dof_coords.size(); i++){
    nodeLocation[i] = dof_coords[i][0];
  }

  define_boundary_conds();

  sparsity_pattern.reinit (dof_handler.n_dofs(), dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();
  K.reinit (sparsity_pattern);
  F.reinit (dof_handler.n_dofs());
  D.reinit (dof_handler.n_dofs());

  //Define quadrature rule

  quadRule = 4;
  quad_points.resize(quadRule); quad_weight.resize(quadRule);

  quad_points[0] = -0.8611363115940526;
  quad_points[1] = -0.3399810435848563;
  quad_points[2] =  0.3399810435848563;
  quad_points[3] =  0.8611363115940526;
  quad_weight[0] = 0.3478548451374538;
  quad_weight[1] = 0.6521451548625461;
  quad_weight[2] = 0.6521451548625461;
  quad_weight[3] = 0.3478548451374538;

  //Just some notes...
  std::cout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;
}

template <int dim>
void FEM<dim>::assemble_system(){

  K=0; F=0;

  const unsigned int   			dofs_per_elem = fe.dofs_per_cell; //This gives you number of degrees of freedom per element
  FullMatrix<double> 				Klocal (dofs_per_elem, dofs_per_elem);
  Vector<double>      			Flocal (dofs_per_elem);
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);
  double										h_e, x;

  //loop over elements
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){


    elem->get_dof_indices (local_dof_indices);


    h_e = nodeLocation[local_dof_indices[1]] - nodeLocation[local_dof_indices[0]];

    Flocal = 0.;
    for(unsigned int A=0; A<dofs_per_elem; A++){
      for(unsigned int q=0; q<quadRule; q++){
	x = 0;
	for(unsigned int B=0; B<dofs_per_elem; B++){
	  x += nodeLocation[local_dof_indices[B]]*basis_function(B,quad_points[q]);
	}
        const double force_per_length = 1.e11 * x * 1.e-4;
        Flocal[A] += force_per_length * basis_function(A, quad_points[q]) *
                     (h_e / 2.) * quad_weight[q];
      }
    }
    if(prob == 2){
      if(nodeLocation[local_dof_indices[1]] == L){
        for (unsigned int A=0; A<dofs_per_elem; ++A)
          Flocal[A] += 1.e6 * basis_function(A, 1.);
      }
    }

    Klocal = 0;
    for(unsigned int A=0; A<dofs_per_elem; A++){
      for(unsigned int B=0; B<dofs_per_elem; B++){
	for(unsigned int q=0; q<quadRule; q++){
          const double grad_A = basis_gradient(A, quad_points[q]) * 2. / h_e;
          const double grad_B = basis_gradient(B, quad_points[q]) * 2. / h_e;
          Klocal[A][B] += 1.e7 * grad_A * grad_B * (h_e / 2.) * quad_weight[q];
	}
      }
    }

    //Assemble local K and F into global K and F
    //You will need to used local_dof_indices[A]
    for(unsigned int A=0; A<dofs_per_elem; A++){
      F[local_dof_indices[A]] += Flocal[A];

      for(unsigned int B=0; B<dofs_per_elem; B++){
        K.add(local_dof_indices[A], local_dof_indices[B], Klocal[A][B]);

      }
    }

  }

  //Apply Dirichlet boundary conditions

  MatrixTools::apply_boundary_values (boundary_values, K, D, F, false);
}

//Solve for D in KD=F
template <int dim>
void FEM<dim>::solve(){

  //Solve for D
  SparseDirectUMFPACK  A;
  A.initialize(K);
  A.vmult (D, F); //D=K^{-1}*F

}

//Output results
template <int dim>
void FEM<dim>::output_results (){

  //Write results to VTK file
  std::ofstream output1("solution.vtk");
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector(D, nodal_solution_names, DataOut<dim>::type_dof_data,
			   nodal_data_component_interpretation);
  data_out.build_patches();
  data_out.write_vtk(output1);
  output1.close();
}

template <int dim>
double FEM<dim>::l2norm_of_error(){

  double l2norm = 0.;

  const unsigned int   			dofs_per_elem = fe.dofs_per_cell; //This gives you dofs per element
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);
  double u_exact, u_h, x, h_e;

  //loop over elements
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (),
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    elem->get_dof_indices (local_dof_indices);

    //Find the element length
    h_e = nodeLocation[local_dof_indices[1]] - nodeLocation[local_dof_indices[0]];

    for(unsigned int q=0; q<quadRule; q++){
      x = 0.; u_h = 0.;
      for(unsigned int B=0; B<dofs_per_elem; B++){
	x += nodeLocation[local_dof_indices[B]]*basis_function(B,quad_points[q]);
	u_h += D[local_dof_indices[B]]*basis_function(B,quad_points[q]);
      }
      if (prob == 1)
        u_exact = -x*x*x/6. + ((g2-g1)/L + L*L/6.)*x + g1;
      else
        u_exact = -x*x*x/6. + (1.e6/1.e7 + L*L/2.)*x + g1;
      l2norm += (u_exact-u_h)*(u_exact-u_h) * (h_e/2.) * quad_weight[q];

    }
  }

  return sqrt(l2norm);
}
