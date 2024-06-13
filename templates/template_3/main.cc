/*%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-
%                  ma deal                   %
%  Template for Incompressible Hyperelastic  % 
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-*/

//! Include the "large_3field" header for elastic large deformation (nearly-incompressible): 3-field mixed-formulation
#include <large_3field.h>

//! Include RVE generator
#include <rve.h>

//! Use madeal name space
using namespace madeal;
 
int main(){
  
  deallog.depth_console(0);
  
  //! (1) Enter Dimension, Problem dimension and polynomial degree
  const int dim   = 2;
  const int p_dim = 3;
  const int p_deg = 1;

  //! (2) Define FEM object and Entering input file
  Elastic::Large_3field<dim, p_dim> problem1("inputfile.prm", p_deg);

  //! (3) Enforce any input manually
  problem1.output_name = "out"; //The output name
  problem1.s_typ       = 2;     //Solution type (1 plane-stress; 2 plane-strain)

  //! (4) Define the mesh
  // (i) From deal
  problem1.L_m = 1.0;
  problem1.W_m = 1.0;
  Point<dim> p1(0, 0), p2(problem1.L_m, problem1.W_m);
  vector<unsigned int> subd={10,10};
  GridGenerator::subdivided_hyper_rectangle(problem1.mesh, subd, p1, p2);
  // (ii) Import the mesh
  // GridIn<dim> grid_in;
  // grid_in.attach_triangulation(problem1.mesh);
  // std::ifstream input_mesh("mesh.inp");
  // grid_in.read_abaqus(input_mesh);
  // problem1.L_m = 1.0;
  // problem1.W_m = 1.0;
  // (iii) Generate randomly distributed RVE
  // const int inputnr = 7;
  // RVE::CNT<dim> myrve1("inputrve.prm", inputnr);
  // myrve1.generate_rve(false);
  // problem1.mesh.copy_triangulation(myrve1.mesh);
  // problem1.L_m = myrve1.L_m;
  // problem1.W_m = myrve1.L_m;

  //! (5) Specify Loads and BCs
  // (i) dbc=Dirichlet
  LoadAndBCs::BC bc1;
  bc1.type     = LoadAndBCs::dbc;
  bc1.position = LoadAndBCs::left;
  bc1.value    = {0.0, 0.0};
  problem1.BCs.push_back(bc1);
  // (ii) nbc=Neumann
  LoadAndBCs::BC bc2;
  bc2.type     = LoadAndBCs::nbc;
  bc2.position = LoadAndBCs::right;
  bc2.value    = {1.0e4, 1.0e4};
  problem1.BCs.push_back(bc2);
  // (iii) ldb=Linear Displacement
  // LoadAndBCs::BC bc3;
  // bc3.type     = LoadAndBCs::ldb;
  // // F = [F11, F22, F12, F21]
  // bc3.value    = {1.5, 1.0, 0.0, 0.0};
  // problem1.BCs.push_back(bc3);
  // (iv) pbc=Periodic
  // LoadAndBCs::BC bc4;
  // bc4.type     = LoadAndBCs::pbc;
  // // F = [F11, F22, F12, F21]
  // bc4.value    = {1.5, 1.0, 0.0, 0.0};
  // problem1.BCs.push_back(bc4);

  //! (6) Run the FE Analysis
  problem1.run();
  
  return 0;

}