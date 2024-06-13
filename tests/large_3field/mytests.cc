/*%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%                    < Test >                         %
%  Large Deformation, 3-field (Nearly-Incompressible) %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%*/

//! Include the headers
#include <large_3field.h>
#include <rve.h>

// Include Google test
#include <gtest/gtest.h>

using namespace madeal;


//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
//
//   the Classic Example: Punch problem    
//
//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

// 2D world 

TEST(ClassicExample, test)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_dim = 2;
  const int p_deg = 2;
  Elastic::Large_3field<dim, p_dim> exampleclassic("inputfile0.prm", p_deg);
  exampleclassic.L_m = 1.0e-3;
  exampleclassic.W_m = 1.0e-3;
  Point<dim> p1(0, 0), p2(exampleclassic.L_m, exampleclassic.W_m);
  vector<unsigned int> subd={4,4};
  GridGenerator::subdivided_hyper_rectangle(exampleclassic.mesh, subd, p1, p2);
  exampleclassic.output_name="punch";

  // (2) Specify Loads and BCs
  // (i) dbc=Dirichlet
  LoadAndBCs::BC bc1;
  bc1.type     = LoadAndBCs::dbc;
  bc1.position = LoadAndBCs::bottom;
  bc1.value    = {0.0, 0.0};
  bc1.x_free   = true;
  exampleclassic.BCs.push_back(bc1);
  // (ii) nbc=Neumann
  LoadAndBCs::BC bc2;
  bc2.type     = LoadAndBCs::nbc;
  bc2.position = LoadAndBCs::halftop;
  bc2.value    = {0.0, -400000.0};
  exampleclassic.BCs.push_back(bc2);
  // (iii) dbc=Dirichlet
  LoadAndBCs::BC bc3;
  bc3.type     = LoadAndBCs::dbc;
  bc3.position = LoadAndBCs::left;
  bc3.value    = {0.0, 0.0};
  bc3.y_free   = true;
  exampleclassic.BCs.push_back(bc3);
  // (iv) dbc=Dirichlet
  LoadAndBCs::BC bc4;
  bc4.type     = LoadAndBCs::dbc;
  bc4.position = LoadAndBCs::top;
  bc4.value    = {0.0, 0.0};
  bc4.y_free   = true;
  exampleclassic.BCs.push_back(bc4);
  // (v) dbc=Dirichlet
  LoadAndBCs::BC bc5;
  bc5.type     = LoadAndBCs::dbc;
  bc5.position = LoadAndBCs::halftop;
  bc5.value    = {0.0, 0.0};
  bc5.y_free   = true;
  exampleclassic.BCs.push_back(bc5);
  
  // (3) Run the FE Analysis
  exampleclassic.run();
  
  // TEST
  ASSERT_NEAR(exampleclassic.sigmaBar[0], 7.76127e6, 0.001e6);
  ASSERT_NEAR(exampleclassic.sigmaBar[1], -1.0035e8, 0.001e8);
  ASSERT_NEAR(exampleclassic.sigmaBar[2], 1.59768e7, 0.001e7);

  // Vertical displacement
  //References = 0.5611 mm [3X3, p_deg=2]
  //Our result = 0.5611 mm

}


/////////////////////
//      TEST2      //
/////////////////////

// 3D

TEST(ClassicExample, test2)
{ 
  // (1) Enter inputs
  const int dim   = 3;
  const int p_dim = 3;
  const int p_deg = 1;
  Elastic::Large_3field<dim, p_dim> exampleclassic("inputfile3d.prm", p_deg);
  exampleclassic.L_m = 1.0e-3;
  exampleclassic.W_m = 1.0e-3;
  exampleclassic.H_m = 1.0e-3;
  Point<dim> p1(0, 0, 0), p2(exampleclassic.L_m, exampleclassic.W_m, exampleclassic.H_m);
  vector<unsigned int> subd={4,4,4};
  GridGenerator::subdivided_hyper_rectangle(exampleclassic.mesh, subd, p1, p2);
  exampleclassic.output_name="punch";

  // (2) Specify Loads and BCs
  // (i) dbc=Dirichlet
  LoadAndBCs::BC bc1;
  bc1.type     = LoadAndBCs::dbc;
  bc1.position = LoadAndBCs::bottom;
  bc1.value    = {0.0, 0.0, 0.0};
  bc1.x_free   = true;
  bc1.z_free   = true;
  exampleclassic.BCs.push_back(bc1);
  // (ii) nbc=Neumann
  LoadAndBCs::BC bc2;
  bc2.type     = LoadAndBCs::nbc;
  bc2.position = LoadAndBCs::quartertop;
  bc2.value    = {0.0, -400000000.0, 0.0};
  exampleclassic.BCs.push_back(bc2);
  // (iii) dbc=Dirichlet
  LoadAndBCs::BC bc3;
  bc3.type     = LoadAndBCs::dbc;
  bc3.position = LoadAndBCs::left;
  bc3.value    = {0.0, 0.0, 0.0};
  bc3.y_free   = true;
  bc3.z_free   = true;
  exampleclassic.BCs.push_back(bc3);
  // (iv) dbc=Dirichlet
  LoadAndBCs::BC bc4;
  bc4.type     = LoadAndBCs::dbc;
  bc4.position = LoadAndBCs::top;
  bc4.value    = {0.0, 0.0, 0.0};
  bc4.y_free   = true;
  exampleclassic.BCs.push_back(bc4);
  // (v) dbc=Dirichlet
  LoadAndBCs::BC bc5;
  bc5.type     = LoadAndBCs::dbc;
  bc5.position = LoadAndBCs::quartertop;
  bc5.value    = {0.0, 0.0, 0.0};
  bc5.y_free   = true;
  exampleclassic.BCs.push_back(bc5);
  // (vi) dbc=Dirichlet
  LoadAndBCs::BC bc6;
  bc6.type     = LoadAndBCs::dbc;
  bc6.position = LoadAndBCs::back;
  bc6.value    = {0.0, 0.0, 0.0};
  bc6.x_free   = true;
  bc6.y_free   = true;
  exampleclassic.BCs.push_back(bc6);

  // (3) Run the FE Analysis
  exampleclassic.run();
  
  // TEST
  ASSERT_NEAR(exampleclassic.sigmaBar[0], 9.53534e6,  0.01e6);
  ASSERT_NEAR(exampleclassic.sigmaBar[1], -4.30058e7, 0.01e7);
  ASSERT_NEAR(exampleclassic.sigmaBar[2], 8.269e6,    0.01e6);

  // Vertical displacement
  //References = 0.7719 mm [3X3X3, p_deg=2]
  //Our result = 0.7719 mm [3X3X3, p_deg=2]
  //Our result = 0.7897 mm [3X3X3, p_deg=1]

}



//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
//
//   Example 1: A simple block 
//
//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

/////////////////////
//      TEST1      //
/////////////////////

// Left:fixed  | Right:Neumann

TEST(Example1, test1)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_dim = 3;
  const int p_deg = 3;
  Elastic::Large_3field<dim, p_dim> example1("inputfile2.prm", p_deg);
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
  bc2.value    = {5.0e7, 0.0};
  example1.BCs.push_back(bc2);
  
  // (3) Run the FE Analysis
  example1.run();

  // TEST Q
  ASSERT_NEAR(example1.meanQ, 852272,      10);
  ASSERT_NEAR(example1.l1Q,   1.22175e9, 0.001e9);
  ASSERT_NEAR(example1.l2Q,   2.05142e8, 0.001e8);

}



//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
//
//   Example 2: A block with a hole
//
//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%


/////////////////////
//      TEST1      //
/////////////////////

// MacroBC: Linear Displacement

TEST(Example2, test1)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_dim = 2;
  const int p_deg = 1;
  Elastic::Large_3field<dim, p_dim> example2("inputfile3.prm", p_deg);
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
  bc.value    = {1.7, 1.6, 0.15, 0.15};
  example2.BCs.push_back(bc);
  
  // (3) Run the FE Analysis
  example2.run();

  // TEST Q
  ASSERT_NEAR(example2.meanQ, 1.12288e11, 0.001e11);
  ASSERT_NEAR(example2.l1Q,   7.92977e14, 0.001e14);
  ASSERT_NEAR(example2.l2Q,   1.92325e13, 0.001e13);

}



//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
//
//   Example 3: CNT/Matrix periodic voxel mesh
//
//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

/////////////////////
//      TEST1      //
/////////////////////

// MacroBC: Periodic BC

TEST(Example3, test1)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_dim = 3;
  const int p_deg = 1;
  Elastic::Large_3field<dim, p_dim> example3("inputfile4.prm", p_deg);
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
  bc.value    = {1.0, 1.05, 0.2, 0.2};
  example3.BCs.push_back(bc);
  
  // (3) Run the FE Analysis
  example3.run();

  // TEST Q
  ASSERT_NEAR(example3.meanQ, 2.0e7,  1e7);
  ASSERT_NEAR(example3.l1Q,   5.0e11, 2e11);
  ASSERT_NEAR(example3.l2Q,   7.0e9,  1e9); 

}




//////////
// MAIN //
//////////

int main(int argc, char **argv) 
{
  testing::InitGoogleTest(&argc, argv);
  
  // Fast test option:
  // testing::GTEST_FLAG(filter) = "Fast_Test*:ClassicExample.test1:Example3.test1";
  
  return RUN_ALL_TESTS();
}