/*%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%                                                   %
%                << Linear Elastic >>               %
%               Vector-valued approach              %
%                                                   %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%*/

// Include Deal II headers
#include <deal.II/lac/sparse_direct.h>

// Include madeal headers
#include <madeal.h>
#include <store_data.h>
#include <loads_and_bcs.h>

// Using functions from std
using std::cout;


MADEAL_NAMESPACE_OPEN
namespace Elastic{

// Small_def Class
template <int dim>
class Small_def{
public:
  // Constructor
  Small_def(const string &inputfile, int p_deg):
    fe(FE_Q<dim>(p_deg), dim) ,dof_h(mesh), ifname(inputfile)
  {}

  // Public Functions
  void run();

  // Public Variables
  //Mesh
  Triangulation<dim> mesh;
  //Increment Setting
  int nc;     //Number of mesh refining cycles
  //Material properties
  vector <double> E;  //Young's Modulus
  vector <double> v;  //Poisson's Ratio
  //Planner setting
  int s_typ;  //Solution type for planner assumptions (1 plane-stress; 2 plane-strain)
  double z;   //Plate thickness for 2D
  //Dimentions
  double L_m; //Length of Geometry
  double W_m; //Width  of Geometry
  double H_m; //Height of Geometry
  //BCs
  vector <LoadAndBCs::BC> BCs;
  //SigmaBar
  vector<double> SigmaBar;
  //Q values for test
  double meanQ;
  double l1Q;
  double l2Q;
  //Name of output files
  string output_name="output";

private:
  // Private Functions
  //Principal Functions
  void readinput();
  void prepare_structure(int);
  void evolve(int);
  //Other Functions
  void refine_mesh();
  void assemble();
  void mark_bc();
  void apply_ebc();
  void solve();
  void cal_stress();
  void results();

  // Private Variables
  ParameterHandler ph; //Handle the Input Parameters
  //
  CellDataStorage<typename Triangulation<dim>::cell_iterator, StoreData::SD1<dim>> cell_ds; //Cell storage
  //
  FESystem<dim>      fe;    //Finite element class
  DoFHandler<dim>    dof_h; //Degree of freedom handler
  //
  BlockSparsityPattern     sp; //Sparsity pattern
  BlockSparseMatrix<double> K; //Stiffness matrix 
  //
  BlockVector<double> Q;  //Nodal variable list
  BlockVector<double> R;  //the Right hand side matrix
  //
  vector<Tensor<2, dim>> sigma_e;  //Stress tensor of elements
  //
  AffineConstraints<double> cn_hn;  //List of constraints for the hanging nodes
  //
  string ifname;  //Input file name
  //
  //FEM Number Parameters
  int n_dof_e; //Number of Degrees of Freedom per element
  int n_elm_t; //Number of Elements in total
  int n_nod_t; //Number of Nodes in total
  int n_qp_d;  //Number of Quadrature Points in every direction
  int n_dof_t; //Number of Degrees of Freedom in total
  int n_qp_e;  //Number of Quadrature Points per element
  int n_qp_f;  //Number of Quadrature Points per face
  int save_n;  //Save results number

};


// Public Functions
template <int dim>
void Small_def<dim>::run(){

  // Start measuring run time
  auto start_t = std::chrono::high_resolution_clock::now();

  // Read inputs and assign them
  readinput();
  
  cout << "......................................." <<endl;
  cout << " Linear Elastic Analysis in " << dim << "D." << endl;
  cout << " Isotropic Material Model. " << endl;
  cout << "......................................." <<endl;
  
  //Mesh refining cycle
  save_n=0; //Set results number to 0
  //
  for (int cyc=0; cyc<nc; cyc++) {

    cout << endl << " | Cycle: " << cyc+1 << "/" << nc <<" | "<<endl;

    prepare_structure(cyc);
    evolve(cyc);
  
  }

  // Set values for test
  meanQ = Q.mean_value();
  l1Q   = Q.l1_norm();
  l2Q   = Q.l2_norm();

  // Done!
  auto end_t = std::chrono::high_resolution_clock::now();
  auto run_t = std::chrono::duration_cast<std::chrono::nanoseconds>(end_t - start_t);
  cout<<endl<< "......................................." <<endl;
  cout << "Run time: " << run_t.count() * 1e-9;
  //
  cout<<endl<< "......................................." <<endl;
  cout      << ".......... <   D O N E !   > .........." <<endl;
  cout      << "......................................." <<endl<<endl;

}

/*---------------------------------------------------------------
-------      P R I N C I P A L     F U N C T I O N S      -------
---------------------------------------------------------------*/

/*  READ INPUTS  */
template <int dim>
void Small_def<dim>::readinput(){

  //Declear entries
  ph.declare_entry("Number of meshing refine cycle", "1", Patterns::Integer());
  ph.declare_entry("Young's Modulus1", "210e9", Patterns::Double());
  ph.declare_entry("Young's Modulus2", "210e9", Patterns::Double());
  ph.declare_entry("Poisson's Ratio1", "0.3", Patterns::Double());
  ph.declare_entry("Poisson's Ratio2", "0.3", Patterns::Double());
  ph.declare_entry("Solution type (1 plane-stress; 2 plane-strain)", "1", Patterns::Integer());
  ph.declare_entry("Plate thickness", "1.0", Patterns::Double());

  //Open input file
  ph.parse_input(ifname);

  //Assign entries to variables
  nc = ph.get_integer("Number of meshing refine cycle");
  //
  E.push_back(ph.get_double("Young's Modulus1"));
  E.push_back(ph.get_double("Young's Modulus2"));
  v.push_back(ph.get_double("Poisson's Ratio1"));
  v.push_back(ph.get_double("Poisson's Ratio2"));
  //
  s_typ = ph.get_integer("Solution type (1 plane-stress; 2 plane-strain)");
  z     = ph.get_double("Plate thickness");
  //
  if (dim==3)
    AssertThrow(z==1.0, ExcMessage("The thickness cannot be used for 3D problems!"));

}

/*  PREPARE STRUCTURE  */
template <int dim>
void Small_def<dim>::prepare_structure(int cyc){

  // Refine the mesh if it's not the first cycle
  if (cyc!=0)
    refine_mesh();
  
  //Mark the boundaries
  mark_bc();

  // Enumerate DoF on mesh
  dof_h.distribute_dofs(fe);
  // Renumber degrees of freedom based on their vector component
  vector<unsigned int> block_component(dim, 0);
  DoFRenumbering::component_wise(dof_h, block_component);
  
  // Compute hanging node constraints
  cn_hn.clear(); //clear the current set of constraints from the last system
  DoFTools::make_hanging_node_constraints(dof_h, cn_hn);
  //
  //Apply Macro-Strain Periodic BC
  for (auto &bc : BCs)
    if (bc.type==LoadAndBCs::pbc){
      bc.value.push_back(bc.value.back()); // eps21=eps12
      LoadAndBCs::MacroStrainPBC(bc.value, L_m, W_m, dof_h, fe, cn_hn, true);}
  //
  cn_hn.close();

  // Initialize cell storage and number of the points
  for (const auto &cell : dof_h.active_cell_iterators())
    cell_ds.initialize(cell, 1);

  // Assign number of FE parameters
  n_dof_e = fe.dofs_per_cell; 
  n_qp_d  = fe.degree+1;
  n_dof_t = dof_h.n_dofs();
  n_elm_t = mesh.n_active_cells();
  n_nod_t = n_dof_t/fe.dofs_per_vertex;
  
  cout << " | Number of nodes:              " << n_nod_t <<" | "<<endl;
  cout << " | Number of degrees of freedom: " << n_dof_t <<" | "<<endl;
  cout << " | Number of Elements:           " << n_elm_t <<" | "<<endl;

  // Evaluate the number of non-zero shape functions for a particular vector component
  vector<types::global_dof_index> n_dof_vc = DoFTools::count_dofs_per_fe_block(dof_h, block_component);
  int n_u_dof = n_dof_vc[0];   //Number of u degree of freedoms

  // Sparsity
  BlockDynamicSparsityPattern dsp(1, 1);
  dsp.block(0, 0).reinit(n_u_dof, n_u_dof);
  dsp.collect_sizes();
  //
  DoFTools::make_sparsity_pattern(dof_h, dsp, cn_hn, true);
  sp.copy_from(dsp);

  // Resize KQR
  K.reinit(sp);
  //
  Q.reinit(1);
  Q.block(0).reinit(n_u_dof);
  Q.collect_sizes();
  //
  R.reinit(1);
  R.block(0).reinit(n_u_dof);
  R.collect_sizes();
  //
  sigma_e.resize(n_elm_t);

}

/*  EVOLVE  */
template <int dim>
void Small_def<dim>::evolve(int cyc){

  //Assemble the global Stiffness matrix and Right hand side
  assemble();

  //Modify K and R for essential boundary conditions
  apply_ebc();
  
  //Solve: (the resulting linear equations: R=KQ )
  solve();
  
  //Compute the values of constrained nodes from the values of the unconstrained ones
  cn_hn.distribute(Q);

  //Calculate stress
  cal_stress();

  //Print results
  results();
  cout<<endl<<" SigmaBar: ["<<SigmaBar[0]<<", "<<SigmaBar[1]<<", "<<SigmaBar[2]<<"]"<<endl;

  cout << " Completed: % " << (100*((cyc+1.0)/nc)) <<endl;
  save_n++;

}


/*---------------------------------------------------------------
-------           O T H E R    F U N C T I O N S          -------
---------------------------------------------------------------*/


/*  REFINE MESH  */
template <int dim>
void Small_def<dim>::refine_mesh(){
  Vector<float> eepc(n_elm_t); //Estimated error per cell
  KellyErrorEstimator<dim>::estimate(dof_h, QGauss<dim - 1>(n_qp_d), {}, Q, eepc);

  GridRefinement::refine_and_coarsen_fixed_number(mesh, eepc, 0.3, 0.03);
  
  mesh.execute_coarsening_and_refinement();

}


/*  ASSEMBLING STIFNESS MATRIX  AND R */
template <int dim>
void Small_def<dim>::assemble(){

  // Determine the gaussian points and weight factors for elements and faces
  QGauss<dim>   egauss(n_qp_d);
  QGauss<dim-1> fgauss(n_qp_d);
  n_qp_e =      egauss.size();
  n_qp_f =      fgauss.size();

  // Calculate the Cartesian derivatives of shape functions and area for all elements and faces in the solution
  FEValues<dim>     fe_values_e(fe, egauss, update_values | update_gradients      | update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_values_f(fe, fgauss, update_values | update_normal_vectors | update_quadrature_points | update_JxW_values);

  // Define element stiffness and rhs
  FullMatrix<double> k_el(n_dof_e, n_dof_e);
  Vector<double>     r_el(n_dof_e);

  // Create a list for storing DOFs of the element
  vector<types::global_dof_index> dofs_el(n_dof_e);

  // Create lists for storing the b, t and P
  vector<Tensor<1, dim>>  b(n_qp_e);
  vector<Tensor<1, dim>>  t(n_qp_f);
  // Tensor<1, dim>          P;

  // Finding P point for applying point load
  // bool PpointFound=false;
  
  // Extractors that get at the displacement components of vector-valued shape functions
  FEValuesExtractors::Vector displacement(0);

  for (const auto &cell : dof_h.active_cell_iterators()){

    // Initialize FE values of the element
    fe_values_e.reinit(cell);

    // Store element data
    const vector<shared_ptr<StoreData::SD1<dim>>> csgd = cell_ds.get_data(cell);
    int mat;
    cell->material_id()==0 ? mat=0 : mat=cell->material_id()-1;
    csgd[0]->setup_ed(E,v,mat,s_typ,z);

    // Get the list of DOFs of the element
    cell->get_dof_indices(dofs_el);

    // Initialize element stiffness and r
    k_el = 0;
    r_el = 0;

    // Get the values of b at the quadrature points of the element
    LoadAndBCs::b_calculate(fe_values_e.get_quadrature_points(),b);

    // Getting the values of intended data for the element
    double la = csgd[0]->la_el;
    double mu = csgd[0]->mu_el;

    // Calculate element stiffness and R(b contribution) over each QP of the element
    for (int q=0; q<n_qp_e; q++){

      for (int I=0; I<n_dof_e; I++){
        Tensor<1, dim> N_Iu    = fe_values_e[displacement].value(I, q);
        Tensor<2, dim> grdN_Iu = fe_values_e[displacement].gradient(I, q);
        double         divN_Iu = fe_values_e[displacement].divergence(I, q);
        
        for (int J=0; J<n_dof_e; J++){
          Tensor<2, dim> grdN_Ju = fe_values_e[displacement].gradient(J, q);
          double         divN_Ju = fe_values_e[displacement].divergence(J, q);

          //K of element
          k_el(I,J) += ( la * divN_Iu * divN_Ju + mu * double_contract<0,0,1,1>(grdN_Iu, grdN_Ju) + mu * double_contract<0,1,1,0>(grdN_Iu, grdN_Ju) )  * z * fe_values_e.JxW(q);

        }

        //R of element
        r_el(I) += b[q] * N_Iu * z * fe_values_e.JxW(q);

      }

    }
    

    // BCs for t and P
    for (auto &bc : BCs)
    if (bc.type==LoadAndBCs::nbc){
      for (const auto &face : cell->face_iterators()){
        if (face->boundary_id()==bc.position){

          fe_values_f.reinit(cell, face);

          // Get the values of t at the quadrature points of the face
          for (int iq=0; iq<n_qp_f; iq++)
            for (int iiq=0; iiq<dim; iiq++)
              t[iq][iiq] = bc.value[iiq];
          
          // Calculate R(t contribution) over each QP of the face
          for (int qf = 0; qf<n_qp_f; qf++){
            for (int I = 0; I<n_dof_e; I++){
              Tensor<1, dim> N_Iu = fe_values_f[displacement].value(I, qf);            

              r_el(I) += t[qf] * N_Iu * fe_values_f.JxW(qf);
            
            }
          }

          // Calculate R(P contribution) over each vertices of the face
          /*LoadAndBCs::P_calculate(P);
          Point<dim> p1;
          p1[0]=10.0; p1[1]=0.4; if(dim==3){p1[2]=0.2;}
          for (const auto v : face->vertex_indices()){
            if ( ((face->vertex(v) - p1).norm_square() < 1e-3) & (!PpointFound) ){
              for (int i=0; i<dim; i++)
                R(face->vertex_dof_index(v, i)) += P[i];
              PpointFound=true;
            }
          }*/

        }
      }
    }
    

    // Form global stiffness and R
    cn_hn.distribute_local_to_global(k_el, r_el, dofs_el, K, R);

  }

}


/*  MARK THE BOUNDARY  */
template <int dim>
void Small_def<dim>::mark_bc(){
  
  //If applying linear displacement BC, skip (all 0 by default)
  bool skipall = false;
  for (auto &bc : BCs)
    if (bc.type==LoadAndBCs::ldb)
      skipall=true;

  if (skipall==false){
  
  //Get the size of the geometry and its corners 
  const double x_max = L_m, y_max = W_m, z_max = H_m, x_min = 0.0, y_min = 0.0, z_min = 0.0;
  Point<2>  c_dl(x_min,y_min), c_dr(x_max,y_min), c_ul(x_min,y_max), c_ur(x_max,y_max); //Corner Points (l,r,u,d: left,right,up,down)
  bool skipcell, checkcorners;
  checkcorners=false;
  skipcell = false;
  for (auto &bc : BCs)
    if (bc.type==LoadAndBCs::pbc)
      checkcorners = true;

  //Mark the boundaries
  for (const auto &cell : mesh.active_cell_iterators())
    if (cell->at_boundary()){

      //Check if cell is at corner [for PBC only]
      if constexpr (dim==2){
      if (checkcorners==true){
        for (const auto v : cell->vertex_indices()){
          if (cell->vertex(v).distance(c_dl)<1e-9)
            skipcell=true;
          else if (cell->vertex(v).distance(c_dr)<1e-9)
            skipcell=true;
          else if (cell->vertex(v).distance(c_ul)<1e-9)
            skipcell=true;
          else if (cell->vertex(v).distance(c_ur)<1e-9)
            skipcell=true;
        }
      }
      }

      //Mark the requested boundaries
      if (skipcell==false)
        for (const auto &face : cell->face_iterators())
        if (face->at_boundary()) {
          if (std::fabs(face->center()(0) - (x_min)) < 1e-9)
            face->set_boundary_id(1);
          else if (std::fabs(face->center()(0) - (x_max)) < 1e-9)
            face->set_boundary_id(2);
          else if (std::fabs(face->center()(1) - (y_min)) < 1e-9)
            face->set_boundary_id(3);
          else if (std::fabs(face->center()(1) - (y_max)) < 1e-9)
            face->set_boundary_id(4);

          if (dim==3){
            if (std::fabs(face->center()(2) - (z_min)) < 1e-9)
              face->set_boundary_id(5);
            else if (std::fabs(face->center()(2) - (z_max)) < 1e-9)
              face->set_boundary_id(6);
          }

        }

      //check half bcs
      if (skipcell==false)
        for (auto &bc : BCs)
          if (bc.position==7)
            for (const auto &face : cell->face_iterators())
              if (face->at_boundary())
                if (std::fabs(face->center()(1) - (y_max)) < 1e-9 && face->center()(0) < (x_max/2.0) )
                  face->set_boundary_id(7);

    }
  
  }
  
}


/*  APPLY Essential BCs  */
template <int dim>
void Small_def<dim>::apply_ebc(){  

  // Extractors
  const FEValuesExtractors::Vector displacement(0);

  // Apply Dirichlet BC
  for (auto &bc : BCs)
    if (bc.type==LoadAndBCs::dbc){
      map<types::global_dof_index, double> cn_bc1;
      VectorTools::interpolate_boundary_values(dof_h, bc.position, Functions::ConstantFunction<dim>(bc.value), cn_bc1);
      MatrixTools::apply_boundary_values(cn_bc1, K, Q, R);
    }

  // Macro-Strain Linear BC
  for (auto &bc : BCs)
    if (bc.type==LoadAndBCs::ldb){
      map<types::global_dof_index, double> cn_bc1;
      bc.value.push_back(bc.value.back()); // eps21=eps12
      VectorTools::interpolate_boundary_values(dof_h, 0, LoadAndBCs::MacroStrainLBC<dim>(bc.value), cn_bc1, fe.component_mask(displacement));
      MatrixTools::apply_boundary_values(cn_bc1, K, Q, R);
    }

}


/*  SOLVE THE LINEAR SYSTEMS  */
template <int dim>
void Small_def<dim>::solve(){

  // (I) Direct Solver
  // SparseDirectUMFPACK Direct_Solver;
  // Direct_Solver.initialize(K);
  // Direct_Solver.vmult(Q, R);
  
  // (II) Iterative CG Solver
  auto &M = K.block(0, 0);
  auto &F = R.block(0);
  auto &U = Q.block(0);
  auto op_M = linear_operator(M);
  //
  ReductionControl reduction_control_M(10000, 1.0e-18, 1.0e-10);
  SolverCG<>       solver_M(reduction_control_M);
  PreconditionJacobi<SparseMatrix<double>> preconditioner_M;
  preconditioner_M.initialize(M);
  auto op_M_inv = inverse_operator(op_M, solver_M, preconditioner_M);
  //
  U = op_M_inv * F;

}

/*  CALCULATE STRESS  */ 
template <int dim>
void Small_def<dim>::cal_stress(){

  // Define SigmaBarVol
  Tensor<2, dim> SigmaBarVol;
  SigmaBar.resize((dim-1)*3);

  // Determine the gaussian points and weight factors
  QGauss<dim> gauss(n_qp_d);

  // Calculate the Cartesian derivatives of shape functions and area for all elements in the solution
  FEValues<dim> fe_values_e(fe, gauss, update_gradients | update_JxW_values);

  // Create a list for storing DOFs of the element
  vector<types::global_dof_index> dofs_el(n_dof_e);

  // Reset volume and sigma of elements
  for (const auto &cell : dof_h.active_cell_iterators()) ///THIS NEEDS UPDATE
    sigma_e[cell->active_cell_index()] = 0;
    
  for (const auto &cell : dof_h.active_cell_iterators()){
    
    // Initialize FE values of the element
    fe_values_e.reinit(cell);

    // Element data
    const vector<shared_ptr<StoreData::SD1<dim>>> csgd = cell_ds.get_data(cell);

    // Get the list of DOFs of the element
    cell->get_dof_indices(dofs_el);

    // Getting the C of the element
    Tensor<4, dim> C = csgd[0]->C_el;

    // Get the u value of nodes of the element [This one also can be improved]
    vector<double> Q_nodes(n_dof_e);
    for (int I=0; I<n_dof_e; I++){
      Q_nodes[I] = Q(dofs_el[I]);
    }

    // Calculate stress of every IP of the element
    for (int q=0; q<n_qp_e; q++){

      // Extractors that get at the displacement
      FEValuesExtractors::Vector displacement(0);

      Tensor<2, dim> grd_u_qp;

      //Loop over each dof
      for (int I=0; I<n_dof_e; I++){
        
        Tensor<2, dim> grd_N_Iu  = fe_values_e[displacement].gradient(I, q);
        
        grd_u_qp += grd_N_Iu * Q_nodes[I];
        
      }

      //Sigma += C*grad_u [This is sum of sigma which will be divided by n_qp_e when outputing results]
      sigma_e[cell->active_cell_index()] += double_contract<0,0,1,1>(C, grd_u_qp);

      //SigmaBarVol
      SigmaBarVol  += double_contract<0,0,1,1>(C, grd_u_qp) * fe_values_e.JxW(q);

    }

  }

  //SigmaBar
  SigmaBar[0] = SigmaBarVol[0][0]/GridTools::volume(mesh);
  SigmaBar[1] = SigmaBarVol[1][1]/GridTools::volume(mesh);
  SigmaBar[2] = SigmaBarVol[0][1]/GridTools::volume(mesh);
  
}

/*  PRINT RESULTS  */
template <int dim>
void Small_def<dim>::results(){
  
  DataOut<dim> res;

  // Vector u
  vector<string> rtitles(dim, "u");
  vector<DataComponentInterpretation::DataComponentInterpretation> interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
  // Write higher orders
  DataOutBase::VtkFlags res_flags;
  res_flags.write_higher_order_cells = true;
  res.set_flags(res_flags);
  
  res.add_data_vector(dof_h, Q, rtitles, interpretation);

  // Stress sigma
  Vector<double> sigmaxx(n_elm_t), sigmayy(n_elm_t), sigmaxy(n_elm_t), sigmazz(n_elm_t), sigmayz(n_elm_t), sigmazx(n_elm_t);
  switch (dim){
    case 1:
      Assert(dim >= 2, ExcNotImplemented());
      break;
    case 2:
      for (int i=0; i<n_elm_t; i++){
        sigmaxx(i)=sigma_e[i][0][0]/((double) n_qp_e);
        sigmayy(i)=sigma_e[i][1][1]/((double) n_qp_e);
        sigmaxy(i)=sigma_e[i][0][1]/((double) n_qp_e);
      }
      res.add_data_vector(sigmaxx, "sigma_xx");
      res.add_data_vector(sigmayy, "sigma_yy");
      res.add_data_vector(sigmaxy, "sigma_xy");
      break;
    case 3:
      for (int i=0; i<n_elm_t; i++){
        sigmaxx(i)=sigma_e[i][0][0]/((double) n_qp_e);
        sigmayy(i)=sigma_e[i][1][1]/((double) n_qp_e);
        sigmazz(i)=sigma_e[i][2][2]/((double) n_qp_e);
        sigmaxy(i)=sigma_e[i][0][1]/((double) n_qp_e);
        sigmayz(i)=sigma_e[i][1][2]/((double) n_qp_e);
        sigmazx(i)=sigma_e[i][2][0]/((double) n_qp_e);
      }
      res.add_data_vector(sigmaxx, "sigma_xx");
      res.add_data_vector(sigmayy, "sigma_yy");
      res.add_data_vector(sigmazz, "sigma_zz");
      res.add_data_vector(sigmaxy, "sigma_xy");
      res.add_data_vector(sigmayz, "sigma_yz");
      res.add_data_vector(sigmazx, "sigma_zx");
      break;
  }

  res.build_patches(fe.degree);
  
  string i=std::to_string(save_n);
  string k=std::to_string(dim);
  std::ofstream output(output_name+"_"+k+"d_"+i+".vtk");
  res.write_vtk(output);

}


}
MADEAL_NAMESPACE_CLOSE