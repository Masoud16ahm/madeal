/*%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%                                                   %
%                      EXAMPLE                      %
%                                                   %
%        < Solve a Linear elastic problem:          %
%               a Beam with Fixed-Free BCs >        %
%                                                   %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%*/

//! Include the "small_def" header
#include <small_def.h>

using namespace madeal;


// //I: Plane-stress
// //
// int main(){
  
//   deallog.depth_console(0);

//   // (1) Enter inputs
//   const int dim   = 2;
//   const int p_deg = 2;
//   Elastic::Small_def<dim> example1("inputfile.prm", p_deg);
//   example1.L_m = 13;
//   example1.W_m = 0.8;
//   example1.H_m = 0.8;
//   Point<dim> p1(0, 0), p2(13.0, 0.8);
//   vector<unsigned int> subd={40,4};
//   GridGenerator::subdivided_hyper_rectangle(example1.mesh, subd, p1, p2);
//   example1.output_name="beam";

//   // (2) Specify Loads and BCs
//   // (i) dbc=Dirichlet
//   LoadAndBCs::BC bc1;
//   bc1.type     = LoadAndBCs::dbc;
//   bc1.position = LoadAndBCs::left;
//   bc1.value    = {0.0, 0.0};
//   example1.BCs.push_back(bc1);
//   // (ii) nbc=Neumann
//   LoadAndBCs::BC bc2;
//   bc2.type     = LoadAndBCs::nbc;
//   bc2.position = LoadAndBCs::top;
//   bc2.value    = {0.0, -2000000};
//   example1.BCs.push_back(bc2);
  
//   // (3) Run the FE Analysis
//   example1.run();

//   // Max deflection comparison
//   //Using solid mechanics: Max u2 = 3.0 m
//   //Our result = 3.0 m

//   return 0;

// }


//II: 3D beam
//  
int main(){
  
  deallog.depth_console(0);


  // (1) Enter inputs
  const int dim   = 3;
  const int p_deg = 2;
  Elastic::Small_def<dim> example1("inputfile3d.prm", p_deg);
  example1.L_m = 13;
  example1.W_m = 0.8;
  example1.H_m = 0.8;
  Point<dim> p1(0, 0, 0), p2(13.0, 0.8, 0.8);
  vector<unsigned int> subd={20,4,4};
  GridGenerator::subdivided_hyper_rectangle(example1.mesh, subd, p1, p2);
  example1.output_name="beam";

  // (2) Specify Loads and BCs
  // (i) dbc=Dirichlet
  LoadAndBCs::BC bc1;
  bc1.type     = LoadAndBCs::dbc;
  bc1.position = LoadAndBCs::left;
  bc1.value    = {0.0, 0.0, 0.0};
  example1.BCs.push_back(bc1);
  // (ii) nbc=Neumann
  LoadAndBCs::BC bc2;
  bc2.type     = LoadAndBCs::nbc;
  bc2.position = LoadAndBCs::top;
  bc2.value    = {0.0, -2000000/0.8, 0.0};
  example1.BCs.push_back(bc2);
  
  // (3) Run the FE Analysis
  example1.run();

  // Max deflection comparison
  //Using solid mechanics: Max u2 = 3.00 m
  //Our result = 2.97 m

  return 0;

}