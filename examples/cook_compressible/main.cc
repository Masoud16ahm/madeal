/*%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%                                                   %
%              EXAMPLE: Cook's membrane             %
%                                                   %
%       < Solve a nonLinear elastic problem >       %
%             Compressible Hyperelastic             %
%                                                   %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%*/

//! Include the "large_classic" header
#include <large_classic.h>

using namespace madeal;


// A Function for transformation to Cook's Membrane
int dim = 2;
template <int dim>
Point<dim> transform (const Point<dim> &p_i) {
  const double &x = p_i[0];
  const double &y = p_i[1];

  const double y_upper = 44.0e-3 + (16.0e-3/48.0e-3)*x; // Line defining upper edge of beam
  const double y_lower =  0.0 + (44.0e-3/48.0e-3)*x; // Line defining lower edge of beam
  const double theta = y/44.0e-3; // Fraction of height along left side of beam
  const double y_transform = (1-theta)*y_lower + theta*y_upper; // Final transformation

  Point<dim> p_o = p_i;
  p_o[1] = y_transform;

  return p_o;
}


// //I: Plane-strain
// //
// int main(){
  
//   deallog.depth_console(0);

//   // (1) Enter inputs
//   const int dim   = 2;
//   const int p_dim = 3;
//   const int p_deg = 2;
//   Elastic::Large_classic<dim, p_dim> example("inputfile.prm", p_deg);
//   example.L_m = 48.0e-3;
//   example.W_m = 0;
//   const Point<dim> p1 = Point<dim>(0.0,         0.0);
//   const Point<dim> p2 = Point<dim>(example.L_m, 44.0e-3);
//   vector<unsigned int> subd={32, 32};
//   GridGenerator::subdivided_hyper_rectangle(example.mesh, subd, p1, p2);
//   GridTools::transform(&transform<dim>, example.mesh);
//   example.output_name="cook";

//   // (2) Specify Loads and BCs
//   // (i) dbc=Dirichlet
//   LoadAndBCs::BC bc1;
//   bc1.type     = LoadAndBCs::dbc;
//   bc1.position = LoadAndBCs::left;
//   bc1.value    = {0.0, 0.0};
//   example.BCs.push_back(bc1);
//   // (ii) nbc=Neumann
//   LoadAndBCs::BC bc2;
//   bc2.type     = LoadAndBCs::nbc;
//   bc2.position = LoadAndBCs::right;
//   bc2.value    = {0.0, 1.0/(0.016)};
//   example.BCs.push_back(bc2);
  
//   // (3) Run the FE Analysis
//   example.run();
  
//   // Vertical displacement of the tip
//   //References = 13.77 mm [running code, not dealii tutorial][16X16, p_deg=1]
//   //Our result = 13.77 mm

//   return 0;

// }



//II: 3D cook's membrane
//
int main(){
  
  deallog.depth_console(0);

  // (1) Enter inputs
  const int dim   = 3;
  const int p_dim = 3;
  const int p_deg = 2;
  Elastic::Large_classic<dim, p_dim> example("inputfile3d.prm", p_deg);
  example.L_m = 48.0e-3;
  example.W_m = 0;
  example.H_m = 1.0e-3;
  const Point<dim> p1 = Point<dim>(0.0,                 0.0,     0.0);
  const Point<dim> p2 = Point<dim>(example.L_m, 44.0e-3, 1.0e-3);
  vector<unsigned int> subd={16, 16, 1};
  GridGenerator::subdivided_hyper_rectangle(example.mesh, subd, p1, p2);
  GridTools::transform(&transform<dim>, example.mesh);
  example.output_name="cook";

  // (2) Specify Loads and BCs
  // (i) dbc=Dirichlet
  LoadAndBCs::BC bc1;
  bc1.type     = LoadAndBCs::dbc;
  bc1.position = LoadAndBCs::left;
  bc1.value    = {0.0, 0.0, 0.0};
  example.BCs.push_back(bc1);
  // (ii) nbc=Neumann
  LoadAndBCs::BC bc2;
  bc2.type     = LoadAndBCs::nbc;
  bc2.position = LoadAndBCs::right;
  bc2.value    = {0.0, 1000.0/(0.016), 0.0};
  example.BCs.push_back(bc2);
  // (iii) dbc=Dirichlet
  LoadAndBCs::BC bc3;
  bc3.type     = LoadAndBCs::dbc;
  bc3.position = LoadAndBCs::front;
  bc3.value    = {0.0, 0.0, 0.0};
  bc3.x_free   = true;
  bc3.y_free   = true;
  example.BCs.push_back(bc3);
  // (iv) dbc=Dirichlet
  LoadAndBCs::BC bc4;
  bc4.type     = LoadAndBCs::dbc;
  bc4.position = LoadAndBCs::back;
  bc4.value    = {0.0, 0.0, 0.0};
  bc4.x_free   = true;
  bc4.y_free   = true;
  example.BCs.push_back(bc4);
  
  // (3) Run the FE Analysis
  example.run();

  // Vertical displacement
  //References = 13.77 mm
  //Our result = 13.77 mm

  return 0;

}