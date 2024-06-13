/*%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%                                                   %
%               EXAMPLE: Punch problem              %
%                                                   %
%       < Solve a nonLinear elastic problem >       %
%            Incompressible Hyperelastic            %
%                                                   %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%*/

//! Include the "large_3field" header
#include <large_3field.h>

using namespace madeal;


//I: 2D world / Plane-strain
//
int main(){

  deallog.depth_console(0);

  // (1) Enter inputs
  const int dim   = 2;
  const int p_dim = 3; //p-strain(3) or 2D wrold(2)
  const int p_deg = 2;
  Elastic::Large_3field<dim, p_dim> example("inputfile.prm", p_deg); //p-strain(s_typ=2) or 2D wrold(s_typ=0)
  example.L_m = 1.0e-3;
  example.W_m = 1.0e-3;
  Point<dim> p1(0, 0), p2(example.L_m, example.W_m);
  vector<unsigned int> subd={4,4};
  GridGenerator::subdivided_hyper_rectangle(example.mesh, subd, p1, p2);
  example.output_name="punch";

  // (2) Specify Loads and BCs
  // (i) dbc=Dirichlet
  LoadAndBCs::BC bc1;
  bc1.type     = LoadAndBCs::dbc;
  bc1.position = LoadAndBCs::bottom;
  bc1.value    = {0.0, 0.0};
  bc1.x_free   = true;
  example.BCs.push_back(bc1);
  // (ii) nbc=Neumann
  LoadAndBCs::BC bc2;
  bc2.type     = LoadAndBCs::nbc;
  bc2.position = LoadAndBCs::halftop;
  bc2.value    = {0.0, -400000.0};
  example.BCs.push_back(bc2);
  // (iii) dbc=Dirichlet
  LoadAndBCs::BC bc3;
  bc3.type     = LoadAndBCs::dbc;
  bc3.position = LoadAndBCs::left;
  bc3.value    = {0.0, 0.0};
  bc3.y_free   = true;
  example.BCs.push_back(bc3);
  // (iv) dbc=Dirichlet
  LoadAndBCs::BC bc4;
  bc4.type     = LoadAndBCs::dbc;
  bc4.position = LoadAndBCs::top;
  bc4.value    = {0.0, 0.0};
  bc4.y_free   = true;
  example.BCs.push_back(bc4);
  // (v) dbc=Dirichlet
  LoadAndBCs::BC bc5;
  bc5.type     = LoadAndBCs::dbc;
  bc5.position = LoadAndBCs::halftop;
  bc5.value    = {0.0, 0.0};
  bc5.y_free   = true;
  example.BCs.push_back(bc5);
  
  // (3) Run the FE Analysis
  example.run();

  // Vertical displacement (2D world)
  //References = 0.5611 mm
  //Our result = 0.5611 mm; p-strain = 0.5581

}



// //II: 3D punch problem
// //
// int main(){

//   deallog.depth_console(0);

//   // (1) Enter inputs
//   const int dim   = 3;
//   const int p_dim = 3;
//   const int p_deg = 2;
//   Elastic::Large_3field<dim, p_dim> example("inputfile3d.prm", p_deg);
//   example.L_m = 1.0e-3;
//   example.W_m = 1.0e-3;
//   example.H_m = 1.0e-3;
//   Point<dim> p1(0, 0, 0), p2(example.L_m, example.W_m, example.H_m);
//   vector<unsigned int> subd={4,4,4};
//   GridGenerator::subdivided_hyper_rectangle(example.mesh, subd, p1, p2);
//   example.output_name="punch";

//   // (2) Specify Loads and BCs
//   // (i) dbc=Dirichlet
//   LoadAndBCs::BC bc1;
//   bc1.type     = LoadAndBCs::dbc;
//   bc1.position = LoadAndBCs::bottom;
//   bc1.value    = {0.0, 0.0, 0.0};
//   bc1.x_free   = true;
//   bc1.z_free   = true;
//   example.BCs.push_back(bc1);
//   // (ii) nbc=Neumann
//   LoadAndBCs::BC bc2;
//   bc2.type     = LoadAndBCs::nbc;
//   bc2.position = LoadAndBCs::quartertop;
//   bc2.value    = {0.0, -400000000.0, 0.0};
//   example.BCs.push_back(bc2);
//   // (iii) dbc=Dirichlet
//   LoadAndBCs::BC bc3;
//   bc3.type     = LoadAndBCs::dbc;
//   bc3.position = LoadAndBCs::left;
//   bc3.value    = {0.0, 0.0, 0.0};
//   bc3.y_free   = true;
//   bc3.z_free   = true;
//   example.BCs.push_back(bc3);
//   // (iv) dbc=Dirichlet
//   LoadAndBCs::BC bc4;
//   bc4.type     = LoadAndBCs::dbc;
//   bc4.position = LoadAndBCs::top;
//   bc4.value    = {0.0, 0.0, 0.0};
//   bc4.y_free   = true;
//   example.BCs.push_back(bc4);
//   // (v) dbc=Dirichlet
//   LoadAndBCs::BC bc5;
//   bc5.type     = LoadAndBCs::dbc;
//   bc5.position = LoadAndBCs::quartertop;
//   bc5.value    = {0.0, 0.0, 0.0};
//   bc5.y_free   = true;
//   example.BCs.push_back(bc5);
//   // (vi) dbc=Dirichlet
//   LoadAndBCs::BC bc6;
//   bc6.type     = LoadAndBCs::dbc;
//   bc6.position = LoadAndBCs::back;
//   bc6.value    = {0.0, 0.0, 0.0};
//   bc6.x_free   = true;
//   bc6.y_free   = true;
//   example.BCs.push_back(bc6);

//   // (3) Run the FE Analysis
//   example.run();

//   // Vertical displacement
//   //References = 0.7719 mm [4X4X4, p_deg=2]
//   //Our result = 0.7719 mm

// }