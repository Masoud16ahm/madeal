/*%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%                    < Test >                     %
%               Small Deformation                 %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%*/

// Include the headers
#include <small_def.h>
#include <rve.h>

// Include Google test
#include <gtest/gtest.h>

using namespace madeal;


//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
//
//   the Classic Example: a Beam 
//              with Fixed-Free BCs    
//
//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

/////////////////////
//      TEST1      //
/////////////////////

// 3D model

TEST(ClassicExample, test1)
{ 
  // (1) Enter inputs
  const int dim   = 3;
  const int p_deg = 2;
  Elastic::Small_def<dim> exampleclassic("inputfile3d.prm", p_deg);
  exampleclassic.L_m = 13;
  exampleclassic.W_m = 0.8;
  exampleclassic.H_m = 0.8;
  Point<dim> p1(0, 0, 0), p2(13.0, 0.8, 0.8);
  vector<unsigned int> subd={10,4,4};
  GridGenerator::subdivided_hyper_rectangle(exampleclassic.mesh, subd, p1, p2);
  exampleclassic.output_name="beam";

  // (2) Specify Loads and BCs
  // (i) dbc=Dirichlet
  LoadAndBCs::BC bc1;
  bc1.type     = LoadAndBCs::dbc;
  bc1.position = LoadAndBCs::left;
  bc1.value    = {0.0, 0.0, 0.0};
  exampleclassic.BCs.push_back(bc1);
  // (ii) nbc=Neumann
  LoadAndBCs::BC bc2;
  bc2.type     = LoadAndBCs::nbc;
  bc2.position = LoadAndBCs::top;
  bc2.value    = {0.0, -2000000/0.8, 0.0};
  exampleclassic.BCs.push_back(bc2);
  
  // (3) Run the FE Analysis
  exampleclassic.run();
  
  // TEST
  ASSERT_NEAR(exampleclassic.SigmaBar[1], -1.231e6, 0.01e6);
  ASSERT_NEAR(exampleclassic.SigmaBar[2], -2.031e7, 0.01e7);

  // Max deflection comparison
  //Using solid mechanics: Max u2 = 3.00 m
  //Our result = 2.96 m

}


/////////////////////
//      TEST2      //
/////////////////////

// 2D Plane-stress

TEST(ClassicExample, test2)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_deg = 2;
  Elastic::Small_def<dim> exampleclassic("inputfile1.prm", p_deg);
  exampleclassic.L_m = 13;
  exampleclassic.W_m = 0.8;
  Point<dim> p1(0, 0), p2(13, 0.8);
  vector<unsigned int> subd={39,4};
  GridGenerator::subdivided_hyper_rectangle(exampleclassic.mesh, subd, p1, p2);
  exampleclassic.output_name="beam";

  // (2) Specify Loads and BCs
  // (i) dbc=Dirichlet
  LoadAndBCs::BC bc1;
  bc1.type     = LoadAndBCs::dbc;
  bc1.position = LoadAndBCs::left;
  bc1.value    = {0.0, 0.0};
  exampleclassic.BCs.push_back(bc1);
  // (ii) nbc=Neumann
  LoadAndBCs::BC bc2;
  bc2.type     = LoadAndBCs::nbc;
  bc2.position = LoadAndBCs::top;
  bc2.value    = {0.0, -2000000};
  exampleclassic.BCs.push_back(bc2);
  
  // (3) Run the FE Analysis
  exampleclassic.run();
  
  // TEST
  ASSERT_NEAR(exampleclassic.SigmaBar[1], -1.238e6, 0.01e6);
  ASSERT_NEAR(exampleclassic.SigmaBar[2], -2.031e7, 0.01e7);  

  // Max deflection comparison
  //Using solid mechanics: Max u2 = 3.0 m
  //Our result = 3.0 m

}


//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
//
//  Example 1: A simple block
//             with 2 different BCs    
//
//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

/////////////////////
//      TEST1      //
/////////////////////

// Adaptive mesh refinement
// Left:fixed  | Right:Dirichlet

TEST(Example1, test1)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_deg = 1;
  Elastic::Small_def<dim> example1("inputfile2.prm", p_deg);
  example1.L_m = 10;
  example1.W_m = 10;
  Point<dim> p1(0, 0), p2(example1.L_m, example1.W_m);
  vector<unsigned int> subd={10,10};
  GridGenerator::subdivided_hyper_rectangle(example1.mesh, subd, p1, p2);

  // (2) Specify Loads and BCs
  // (i) dbc=Dirichlet
  LoadAndBCs::BC bc1;
  bc1.type     = LoadAndBCs::dbc;
  bc1.position = LoadAndBCs::left;
  bc1.value    = {0.0, 0.0};
  example1.BCs.push_back(bc1);
  // (ii) dbc=Dirichlet
  LoadAndBCs::BC bc2;
  bc2.type     = LoadAndBCs::dbc;
  bc2.position = LoadAndBCs::right;
  bc2.value    = {0.0, 1.5};
  example1.BCs.push_back(bc2);
  
  // (3) Run the FE Analysis
  example1.run();
  
  // TEST Q
  ASSERT_NEAR(example1.meanQ,  0.375, 0.01);
  ASSERT_NEAR(example1.l1Q,  1266.72, 0.1);
  ASSERT_NEAR(example1.l2Q,  35.4877, 0.01);
  ASSERT_NEAR(example1.SigmaBar[2], 2.74792e9, 0.001e9);

}


/////////////////////
//      TEST2      //
/////////////////////

// Left:fixed  | Right:Neumann

TEST(Example1, test2)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_deg = 3;
  Elastic::Small_def<dim> example1("inputfile3.prm", p_deg);
  example1.L_m = 10;
  example1.W_m = 10;
  Point<dim> p1(0, 0), p2(example1.L_m, example1.W_m);
  vector<unsigned int> subd={10,10};
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
  bc2.value    = {0.0, 1.0e7};
  example1.BCs.push_back(bc2);
  
  // (3) Run the FE Analysis
  example1.run();
  
  // TEST Q
  ASSERT_NEAR(example1.meanQ, 0.075192, 0.0001);
  ASSERT_NEAR(example1.l1Q,   191.749,  0.01);
  ASSERT_NEAR(example1.l2Q,   6.08796,  0.001);
  ASSERT_NEAR(example1.SigmaBar[2], 2.5e7, 0.01e7);

}


//%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
//
//   Example 2: CNT/Matrix periodic voxel mesh
//              with 2 different BCs    
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
  const int p_deg = 1;
  Elastic::Small_def<dim> example2("inputfile3.prm", p_deg);
  const int inputnr = 7;
  RVE::CNT<dim> myrve1("inputrve.prm", inputnr);
  myrve1.generate_rve(true);
  example2.mesh.copy_triangulation(myrve1.mesh);
  example2.L_m = myrve1.L_m;
  example2.W_m = myrve1.L_m;

  // (2) Specify Loads and BCs
  // (i) ldb=Linear Displacement
  LoadAndBCs::BC bc;
  bc.type     = LoadAndBCs::ldb;
  bc.value    = {0.1, 0.0, 0.2};
  example2.BCs.push_back(bc);
  
  // (3) Run the FE Analysis
  example2.run();

  // TEST Q
  ASSERT_NEAR(example2.meanQ, 6,     3);
  ASSERT_NEAR(example2.l1Q,   85000, 10000);
  ASSERT_NEAR(example2.l2Q,   850,   100);

}


/////////////////////
//      TEST2      //
/////////////////////

// MacroBC: Periodic BC

TEST(Example2, test2)
{ 
  // (1) Enter inputs
  const int dim   = 2;
  const int p_deg = 1;
  Elastic::Small_def<dim> example2("inputfile3.prm", p_deg);
  const int inputnr = 7;
  RVE::CNT<dim> myrve1("inputrve.prm", inputnr);
  myrve1.generate_rve(false);
  example2.mesh.copy_triangulation(myrve1.mesh);
  example2.L_m = myrve1.L_m;
  example2.W_m = myrve1.L_m;

  // (2) Specify Loads and BCs
  // (i) pbc=Periodic
  LoadAndBCs::BC bc;
  bc.type     = LoadAndBCs::pbc;
  bc.value    = {0.1, 0.3, 0.2};
  example2.BCs.push_back(bc);
  
  // (3) Run the FE Analysis
  example2.run();

  // TEST Q
  ASSERT_NEAR(example2.meanQ, 10,     3);
  ASSERT_NEAR(example2.l1Q,   140000, 50000);
  ASSERT_NEAR(example2.l2Q,   1400,   500); 

}




//////////
// MAIN //
//////////

int main(int argc, char *argv[])
{
  testing::InitGoogleTest(&argc, argv);

  // Fast test option:
  // testing::GTEST_FLAG(filter) = "Fast_Test*:ClassicExample.test2:Example2.test2";

  return RUN_ALL_TESTS();
}
