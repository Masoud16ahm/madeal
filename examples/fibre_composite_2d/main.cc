/*%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%                                                   %
%        EXAMPLE: Fibre-Reinforced composite        %
%                                                   %
%         <  nonLinear elastic problem  >           %
%            Compressible Hyperelastic              %
%                                                   %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%*/

//! Include the "large_classic" header
#include <large_classic.h>
#include <rve.h>

using namespace madeal;


// //I: Linear Displacement BC, Plane-stress & Incompressible matrix
// //
// int main(){
  
//   deallog.depth_console(0);

//   // (1) Enter inputs
//   const int dim   = 2;
//   const int p_dim = 2;
//   const int p_deg = 2;
//   Elastic::Large_classic<dim, p_dim> example("inputfile.prm", p_deg);
//   //Input mesh:
//   // GridIn<dim> grid_in;
//   // grid_in.attach_triangulation(example.mesh);
//   // std::ifstream input_mesh("rve2d.inp");
//   // grid_in.read_abaqus(input_mesh);
//   // example.L_m = 44.72;
//   // example.W_m = 44.72;
//   //Voxel mesh:
//   const int inputnr = 7;
//   RVE::CNT<dim> myrve1("inputrve.prm", inputnr);
//   myrve1.generate_rve(true);
//   example.mesh.copy_triangulation(myrve1.mesh);
//   example.L_m = myrve1.L_m;
//   example.W_m = myrve1.L_m;
//   //
//   example.output_name="rve";

//   // (2) Specify Loads and BCs
//   // (i) ldb=Linear Displacement
//   LoadAndBCs::BC bc;
//   bc.type     = LoadAndBCs::ldb;
//   bc.value    = {1.2, 0.9, 0.2, 0.2};
//   example.BCs.push_back(bc);
  
//   // (3) Run the FE Analysis
//   example.run();

//   return 0;

// }



//II: Periodic BC, Plane-stress & Incompressible matrix
//
int main(){
  
  deallog.depth_console(0);

  // (1) Enter inputs
  const int dim   = 2;
  const int p_dim = 2;
  const int p_deg = 2;
  Elastic::Large_classic<dim, p_dim> example("inputfile.prm", p_deg);
  //Voxel mesh:
  const int inputnr = 7;
  RVE::CNT<dim> myrve1("inputrve.prm", inputnr);
  myrve1.generate_rve(false);
  example.mesh.copy_triangulation(myrve1.mesh);
  example.L_m = myrve1.L_m;
  example.W_m = myrve1.L_m;
  //
  example.output_name="rve";

  // (2) Specify Loads and BCs
  // (i) pbc=Periodic
  LoadAndBCs::BC bc;
  bc.type     = LoadAndBCs::pbc;
  bc.value    = {1.2, 0.9, 0.1, 0.1};
  example.BCs.push_back(bc);
  
  // (3) Run the FE Analysis
  example.run();

  return 0;

}