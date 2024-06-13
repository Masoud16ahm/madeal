/*%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%                                                   %
%                       madeal                      %
%                                                   %
%    < A deal module developed basd on dealii >     %
%                                                   %
%              < by: Masoud Ahmadi >                %
%         < Masoud.Ahmadi@glasgow.ac.uk >           %
%           UNIVERSITY of GLASGOW 2023              %
%                                                   %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%*/


// Include C++ headers
#include <fstream>
#include <iostream>
#include <memory>
#include <chrono>

// Include Deal II headers
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
//
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
//
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
//
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_point_data.h>
//
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
//
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>


// Using functions from std
using std::endl;
using std::vector;
using std::map;
using std::string;
using std::shared_ptr;

// Defining madeal namespace
#define MADEAL_NAMESPACE_OPEN namespace madeal{
#define MADEAL_NAMESPACE_CLOSE }

// Impoert dealii namespace
using namespace dealii;