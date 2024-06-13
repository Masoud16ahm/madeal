/*%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%                    < Test >                     %
%    Large Deformation, Classic (Compressible)    %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%*/

//! Include the headers
#include <large_classic.h>
#include <rve.h>

// Include Google test
#include <gtest/gtest.h>

using namespace madeal;


//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
//
//   the Classic Example: Cook's membrane    
//
//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

// A Function for transformation to Cook's Membrane
int dim = 2;
template <int dim>
Point<dim> transform (const Point<dim> &pt_in) {
  const double &x = pt_in[0];
  const double &y = pt_in[1];

  const double y_upper = 44.0e-3 + (16.0e-3/48.0e-3)*x; // Line defining upper edge of beam
  const double y_lower =  0.0 + (44.0e-3/48.0e-3)*x; // Line defining lower edge of beam
  const double theta = y/44.0e-3; // Fraction of height along left side of beam
  const double y_transform = (1-theta)*y_lower + theta*y_upper; // Final transformation

  Point<dim> pt_out = pt_in;
  pt_out[1] = y_transform;

  return pt_out;
}


/////////////////////
//      TEST1      //
/////////////////////

// cook's membrane: Plane-strain

TEST(ClassicExample, test1)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_dim = 3;
  const int p_deg = 1;
  Elastic::Large_classic<dim, p_dim> exampleclassic1("inputfile1.prm", p_deg);
  exampleclassic1.L_m = 48.0e-3;
  exampleclassic1.W_m = 0;
  const Point<dim> p1 = Point<dim>(0.0,         0.0);
  const Point<dim> p2 = Point<dim>(exampleclassic1.L_m, 44.0e-3);
  vector<unsigned int> subd={16, 16};
  GridGenerator::subdivided_hyper_rectangle(exampleclassic1.mesh, subd, p1, p2);
  GridTools::transform(&transform<dim>, exampleclassic1.mesh);
  exampleclassic1.output_name="cook";

  // (2) Specify Loads and BCs
  // (i) dbc=Dirichlet
  LoadAndBCs::BC bc1;
  bc1.type     = LoadAndBCs::dbc;
  bc1.position = LoadAndBCs::left;
  bc1.value    = {0.0, 0.0};
  exampleclassic1.BCs.push_back(bc1);
  // (ii) nbc=Neumann
  LoadAndBCs::BC bc2;
  bc2.type     = LoadAndBCs::nbc;
  bc2.position = LoadAndBCs::right;
  bc2.value    = {0.0, 1.0/(0.016)};
  exampleclassic1.BCs.push_back(bc2);
  
  // (3) Run the FE Analysis
  exampleclassic1.run();
  
  // TEST sigmaBar
  ASSERT_NEAR(exampleclassic1.sigmaBar[0], -2291.81, 1.0);
  ASSERT_NEAR(exampleclassic1.sigmaBar[1], 34259.5,  1.0);

  // TEST Q     
  ASSERT_NEAR(exampleclassic1.meanQ, 0.00130268, 0.000001);
  ASSERT_NEAR(exampleclassic1.l1Q,   2.34007, 0.001);
  ASSERT_NEAR(exampleclassic1.l2Q,     0.138761, 0.00001);

  // Vertical displacement
  //References = 13.58 mm [running code, not dealii tutorial][16X16, p_deg=1]
  //Our result = 13.58 mm

}


/////////////////////
//      TEST2      //
/////////////////////

// 3D cook's membrane

TEST(ClassicExample, test2)
{ 
  // (1) Enter inputs
  const int dim   = 3;
  const int p_dim = 3;
  const int p_deg = 1;
  Elastic::Large_classic<dim, p_dim> exampleclassic1("inputfile3d.prm", p_deg);
  exampleclassic1.L_m = 48.0e-3;
  exampleclassic1.W_m = 0;
  exampleclassic1.H_m = 1.0e-3;
  const Point<dim> p1 = Point<dim>(0.0,                 0.0,     0.0);
  const Point<dim> p2 = Point<dim>(exampleclassic1.L_m, 44.0e-3, 1.0e-3);
  vector<unsigned int> subd={8, 8, 2};
  GridGenerator::subdivided_hyper_rectangle(exampleclassic1.mesh, subd, p1, p2);
  GridTools::transform(&transform<dim>, exampleclassic1.mesh);
  exampleclassic1.output_name="cook";

  // (2) Specify Loads and BCs
  // (i) dbc=Dirichlet
  LoadAndBCs::BC bc1;
  bc1.type     = LoadAndBCs::dbc;
  bc1.position = LoadAndBCs::left;
  bc1.value    = {0.0, 0.0, 0.0};
  exampleclassic1.BCs.push_back(bc1);
  // (ii) nbc=Neumann
  LoadAndBCs::BC bc2;
  bc2.type     = LoadAndBCs::nbc;
  bc2.position = LoadAndBCs::right;
  bc2.value    = {0.0, 1000.0/(0.016), 0.0};
  exampleclassic1.BCs.push_back(bc2);
  // (iii) dbc=Dirichlet
  LoadAndBCs::BC bc3;
  bc3.type     = LoadAndBCs::dbc;
  bc3.position = LoadAndBCs::front;
  bc3.value    = {0.0, 0.0, 0.0};
  bc3.x_free   = true;
  bc3.y_free   = true;
  exampleclassic1.BCs.push_back(bc3);
  // (iv) dbc=Dirichlet
  LoadAndBCs::BC bc4;
  bc4.type     = LoadAndBCs::dbc;
  bc4.position = LoadAndBCs::back;
  bc4.value    = {0.0, 0.0, 0.0};
  bc4.x_free   = true;
  bc4.y_free   = true;
  exampleclassic1.BCs.push_back(bc4);
  
  // (3) Run the FE Analysis
  exampleclassic1.run();
  
  // TEST sigmaBar
  ASSERT_NEAR(exampleclassic1.sigmaBar[0], -2589.66, 1);
  ASSERT_NEAR(exampleclassic1.sigmaBar[1], 32520.0,  1);

  // TEST Q
  ASSERT_NEAR(exampleclassic1.meanQ, 0.00084857, 0.000001);
  ASSERT_NEAR(exampleclassic1.l1Q,        1.927, 0.01);
  ASSERT_NEAR(exampleclassic1.l2Q,     0.125989, 0.00001);


  // Vertical displacement
  //References = 13.07 mm [running code, not dealii tutorial][8X8X2, p_deg=1]
  //Our result = 13.07 mm

}



//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
//
//   Example 1: A simple block
//
//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

/////////////////////
//      TEST1      //
/////////////////////

// Compressible Plane-stress
// Left:fixed  | Right:Neumann

TEST(Example1, test1)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_dim = 2;
  const int p_deg = 3;
  Elastic::Large_classic<dim, p_dim> example1("inputfile2.prm", p_deg);
  example1.L_m = 10;
  example1.W_m = 10;
  Point<dim> p1(0, 0), p2(example1.L_m, example1.W_m);
  vector<unsigned int> subd={2,2};
  GridGenerator::subdivided_hyper_rectangle(example1.mesh, subd, p1, p2);

  // (2) Specify Loads and BCs
  // (i) dbc=Dirichlet
  LoadAndBCs::BC bc1;
  bc1.type     = LoadAndBCs::dbc;
  bc1.position = LoadAndBCs::left;
  bc1.value    = {0.0, 0.0};
  example1.BCs.push_back(bc1);
  // (ii) nbc=Neumann
  LoadAndBCs::BC bc2;
  bc2.type     = LoadAndBCs::nbc;
  bc2.position = LoadAndBCs::right;
  bc2.value    = {5.0e5, 0.0};
  example1.BCs.push_back(bc2);
  
  // (3) Run the FE Analysis
  example1.run();

  // TEST Q
  ASSERT_NEAR(example1.meanQ, 1.53584, 0.001);
  ASSERT_NEAR(example1.l1Q,   1052.01, 0.1);
  ASSERT_NEAR(example1.l2Q,   66.8367, 0.1);

}



//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
//
//   Example 2: A block with a hole
//              with 2 Macro BCs    
//
//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%


/////////////////////
//      TEST1      //
/////////////////////

// 2D world
// MacroBC: Linear Displacement

TEST(Example2, test1)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_dim = 2;
  const int p_deg = 1;
  Elastic::Large_classic<dim, p_dim> example2("inputfile0.prm", p_deg);
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(example2.mesh);
  std::ifstream input_mesh("hole.inp");
  grid_in.read_abaqus(input_mesh);
  example2.L_m = 44.72;
  example2.W_m = 44.72;

  // (2) Specify Loads and BCs
  // (i) ldb=Linear Displacement
  LoadAndBCs::BC bc;
  bc.type     = LoadAndBCs::ldb;
  bc.value    = {1.4, 1.5, 0.2, 0.2};
  example2.BCs.push_back(bc);
  
  // (3) Run the FE Analysis
  example2.run();

  // TEST Q
  ASSERT_NEAR(example2.meanQ, 13.8259, 0.001);
  ASSERT_NEAR(example2.l1Q,   50630.4, 1.0);
  ASSERT_NEAR(example2.l2Q,   932.719, 0.01);

}


/////////////////////
//      TEST2      //
/////////////////////

// InCompressible Plane-stress
// MacroBC: Periodic BC

TEST(Example2, test2)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_dim = 2;
  const int p_deg = 1;
  Elastic::Large_classic<dim, p_dim> example2("inputfile3.prm", p_deg);
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(example2.mesh);
  std::ifstream input_mesh("hole.inp");
  grid_in.read_abaqus(input_mesh);
  example2.L_m = 44.72;
  example2.W_m = 44.72;

  // (2) Specify Loads and BCs
  // (i) pbc=Periodic
  LoadAndBCs::BC bc;
  bc.type     = LoadAndBCs::pbc;
  bc.value    = {1.0, 1.05, 0.2, 0.2};
  example2.BCs.push_back(bc);
  
  // (3) Run the FE Analysis
  example2.run();

  // TEST Q
  ASSERT_NEAR(example2.meanQ, 3.71977, 0.001);
  ASSERT_NEAR(example2.l1Q,   13627.3, 1.0);
  ASSERT_NEAR(example2.l2Q,   309.785, 0.01);

}



//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
//
//   Example 3: CNT/Matrix periodic voxel mesh
//              with 2 different BCs    
//
//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%


/////////////////////
//      TEST1      //
/////////////////////

// Plane-stress
// MacroBC: Linear Displacement

TEST(Example3, test1)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_dim = 2;
  const int p_deg = 1;
  Elastic::Large_classic<dim, p_dim> example3("inputfile4.prm", p_deg);
  const int inputnr = 7;
  RVE::CNT<dim> myrve1("inputrve.prm", inputnr);
  myrve1.generate_rve(true);
  example3.mesh.copy_triangulation(myrve1.mesh);
  example3.L_m = myrve1.L_m;
  example3.W_m = myrve1.L_m;

  // (2) Specify Loads and BCs
  // (i) ldb=Linear Displacement
  LoadAndBCs::BC bc;
  bc.type     = LoadAndBCs::ldb;
  bc.value    = {1.0, 1.0, 0.3, 0.0};
  example3.BCs.push_back(bc);
  
  // (3) Run the FE Analysis
  example3.run();

  // TEST Q
  ASSERT_NEAR(example3.meanQ, 3.7,     1.0);
  ASSERT_NEAR(example3.l1Q,   51000,  5000);
  ASSERT_NEAR(example3.l2Q,   700.0,   100);

}


/////////////////////
//      TEST2      //
/////////////////////

// Plane-strain
// MacroBC: Periodic BC

TEST(Example3, test2)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_dim = 3;
  const int p_deg = 1;
  Elastic::Large_classic<dim, p_dim> example3("inputfile5.prm", p_deg);
  const int inputnr = 7;
  RVE::CNT<dim> myrve1("inputrve.prm", inputnr);
  myrve1.generate_rve(false);
  example3.mesh.copy_triangulation(myrve1.mesh);
  example3.L_m = myrve1.L_m;
  example3.W_m = myrve1.L_m;

  // (2) Specify Loads and BCs
  // (i) pbc=Periodic
  LoadAndBCs::BC bc;
  bc.type     = LoadAndBCs::pbc;
  bc.value    = {1.7, 1.9, 0.15, 0.15};
  example3.BCs.push_back(bc);
  
  // (3) Run the FE Analysis
  example3.run();

  // TEST Q
  ASSERT_NEAR(example3.meanQ, 23.0,       5);
  ASSERT_NEAR(example3.l1Q,   300000, 50000);
  ASSERT_NEAR(example3.l2Q,   3000,     300); 

}




//////////
// MAIN //
//////////

int main(int argc, char *argv[])
{
  testing::InitGoogleTest(&argc, argv);

  // Fast test option:
  // testing::GTEST_FLAG(filter) = "Fast_Test*:ClassicExample.test1:Example2.test1:Example3.test2";
  
  return RUN_ALL_TESTS();
}
