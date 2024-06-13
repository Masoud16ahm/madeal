/*%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
%                                                   %
%       << Large Deformation Hyperelastic >>        %
%                Classic formulation                %
%                                                   %
%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%*/

// Include Deal II headers
#include <deal.II/lac/sparse_direct.h>
//
#include <deal.II/physics/elasticity/standard_tensors.h>
#include <deal.II/physics/elasticity/kinematics.h>
//
#include <deal.II/base/table_handler.h>

// Include madeal headers
#include <madeal.h>
#include <store_data.h>
#include <loads_and_bcs.h>

// Using functions from std
using std::cout;


MADEAL_NAMESPACE_OPEN
namespace Elastic{

// Large_classic Class
template <int dim, int p_dim>
class Large_classic{
public:
  // Constructor
  Large_classic(const string &inputfile, int p_deg):
    fe(FE_Q<dim>(p_deg), dim) ,dof_h(mesh), ifname(inputfile)
  {}

  // Public Functions
  void run();

  // Public Variables
  //Mesh
  Triangulation<dim> mesh;
  //Increment Setting
  int nc;      //Number of mesh refining cycles
  int ns;      //Number of time increment steps
  int save_r;  //Print frequency to output the results to file (Save result)
  double tol;  //Newton iteration tolerance
  int max_iter;//Maximum Newton iteration
  int lsolver; //Linear Solver for KQ=R (1 Direct; 2 Iterative)

  //Material properties
  vector <double> G; //Shear Modulus (mu)
  vector <double> k; //Bulk Modulus
  vector <double> v; //Poisson's Ratio
  //2D setting
  int s_typ;  //Solution type (0 flatworld; 1 plane-stress; 2 plane-strain)
  double z;   //Plate thickness
  //Dimentions
  double L_m; //Length of Geometry
  double W_m; //Width  of Geometry
  double H_m; //Height of Geometry
  //BCs
  vector <LoadAndBCs::BC> BCs;
  //Volume
  double Volume, volume;
  //sigmaBar
  vector<double> sigmaBar;
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
  void assemble(double, double);
  void mark_bc();
  void apply_ebc(int);
  void solve();
  bool convergance_check(int, double&);
  void results(int);

  // Private Variables
  ParameterHandler ph; //Handle the Input Parameters
  //
  CellDataStorage<typename Triangulation<dim>::cell_iterator, StoreData::SD2<dim, p_dim>> Element_ds; //Element data storage
  CellDataStorage<typename Triangulation<dim>::cell_iterator, StoreData::SD2<dim, p_dim>> QP_ds;      //Quadrature Point data storage
  //
  FESystem<dim>   fe;    //Finite element class
  DoFHandler<dim> dof_h; //Degree of freedom handler
  //
  BlockSparsityPattern     sp; //Sparsity pattern
  BlockSparseMatrix<double> K; //Stiffness matrix 
  //
  BlockVector<double> Q;     //Nodal variable list
  BlockVector<double> Q_stp; //Nodal variable list for time step
  BlockVector<double> R;     //the Right hand side matrix
  //
  vector<Tensor<2, p_dim>> sigma_e; //Stress tensor of elements
  TableHandler results_table;       //Table of results
  //
  AffineConstraints<double> cn_hn;   //List of constraints for the hanging nodes
  AffineConstraints<double> cn_hn2;  //List of constraints for the hanging nodes 2
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
template <int dim, int p_dim>
void Large_classic<dim, p_dim>::run(){

  // Start measuring run time
  auto start_t = std::chrono::high_resolution_clock::now();
  
  // Read inputs and assign them
  readinput();
  
  cout << "......................................." <<endl;
  cout << " Nonlinear Hyperelastic Analysis in ";
  //
  if (dim==3)
    cout << "3D" << endl;
  else if (dim==2){
    if (p_dim==2 && s_typ==1)
      cout << "Plane-Stress" << endl;
    else if (p_dim==3 && s_typ==2)
      cout << "Plane-Strain" << endl;
    else if (p_dim==2 && s_typ==0)
      cout << "2D flatland" << endl;
  }
  //
  cout << " usin classic Single-field Formulation. " << endl;
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
template <int dim, int p_dim>
void Large_classic<dim, p_dim>::readinput(){

  //Declear entries
  ph.declare_entry("Number of meshing refine cycle", "1", Patterns::Integer());
  ph.declare_entry("Shear Modulus 1", "0", Patterns::Double());
  ph.declare_entry("Shear Modulus 2", "0", Patterns::Double());
  ph.declare_entry("Bulk Modulus 1", "0", Patterns::Double());
  ph.declare_entry("Bulk Modulus 2", "0", Patterns::Double());
  ph.declare_entry("Poisson's Ratio 1", "0", Patterns::Double());
  ph.declare_entry("Poisson's Ratio 2", "0", Patterns::Double());
  ph.declare_entry("Solution type (1 plane-stress; 2 plane-strain)", "2", Patterns::Integer());
  ph.declare_entry("Plate thickness", "1.0", Patterns::Double());
  ph.declare_entry("Number of time increment steps", "10", Patterns::Double());
  ph.declare_entry("Save result", "1", Patterns::Integer());
  ph.declare_entry("Newton iteration tolerance", "1.0e-12", Patterns::Double());
  ph.declare_entry("Maximum Newton iteration", "10", Patterns::Integer());
  ph.declare_entry("Linear Solver (1 Direct; 2 Iterative)", "1", Patterns::Selection("1|2"));

  //Open input file 
  ph.parse_input(ifname);

  //Assign entries to variables
  nc      = ph.get_integer("Number of meshing refine cycle");
  ns      = ph.get_integer("Number of time increment steps");
  tol     = ph.get_double ("Newton iteration tolerance");
  max_iter= ph.get_integer("Maximum Newton iteration");
  save_r  = ph.get_integer("Save result");
  lsolver = ph.get_integer("Linear Solver (1 Direct; 2 Iterative)");
  //
  G.push_back(ph.get_double("Shear Modulus 1"));
  k.push_back(ph.get_double("Bulk Modulus 1"));
  v.push_back(ph.get_double("Poisson's Ratio 1"));
  G.push_back(ph.get_double("Shear Modulus 2"));
  k.push_back(ph.get_double("Bulk Modulus 2"));
  v.push_back(ph.get_double("Poisson's Ratio 2"));
  //
  s_typ = ph.get_integer("Solution type (1 plane-stress; 2 plane-strain)");
  z     = ph.get_double("Plate thickness");
  //
  AssertThrow(dim==3 || (dim==2 && p_dim==2 && s_typ==0) || (dim==2 && p_dim==2 && s_typ==1) || (dim==2 && p_dim==3 && s_typ==2), ExcMessage("The solution type is not valid!"));
  if (dim==3)
    AssertThrow(z==1.0, ExcMessage("The thickness cannot be used for 3D problems!"));
  
}

/*  PREPARE STRUCTURE  */
template <int dim, int p_dim>
void Large_classic<dim, p_dim>::prepare_structure(int cyc){

  // Create geometry and mesh it, refine the mesh if it's not the first cycle
  if (cyc!=0)
    refine_mesh();

  //Initial Volume
  Volume = GridTools::volume(mesh) * z;
  
  //Mark the boundaries
  mark_bc();

  // Enumerate DoF on mesh
  dof_h.distribute_dofs(fe);
  // Renumber degrees of freedom based on their vector component
  vector<unsigned int> block_component(dim, 0);
  DoFRenumbering::component_wise(dof_h, block_component);
  
  // Compute hanging node constraints
  cn_hn.clear(); //clear the current set of constraints from the last system
  cn_hn2.clear();
  DoFTools::make_hanging_node_constraints(dof_h, cn_hn);
  DoFTools::make_hanging_node_constraints(dof_h, cn_hn2);
  //
  //Apply Macro-Strain Periodic BC
  for (auto &bc : BCs)
    if (bc.type==LoadAndBCs::pbc){
      vector <double> valuepertime(bc.value.size());
      
      // F - I 
      for (long unsigned int i=0; i<bc.value.size(); i++){
        if (i<dim)
          valuepertime[i] = (bc.value[i]-1.0)/ns;
        else
          valuepertime[i] = bc.value[i]/ns;}

      LoadAndBCs::MacroStrainPBC(valuepertime, L_m, W_m, dof_h, fe, cn_hn, true);
      LoadAndBCs::MacroStrainPBC(valuepertime, L_m, W_m, dof_h, fe, cn_hn2, false);
      }
  //
  cn_hn.close();
  cn_hn2.close();

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

  // Sparsity
  BlockDynamicSparsityPattern dsp(n_dof_vc, n_dof_vc);
  dsp.collect_sizes();
  //
  DoFTools::make_sparsity_pattern(dof_h, dsp, cn_hn, true);
  sp.copy_from(dsp);

  // Resize KQR
  K.reinit(sp);
  //
  Q.reinit(n_dof_vc);
  Q.collect_sizes();
  //
  Q_stp.reinit(n_dof_vc);
  Q_stp.collect_sizes();
  //
  R.reinit(n_dof_vc);
  R.collect_sizes();
  //
  sigma_e.resize(n_elm_t);

}


/*  EVOLVE  */
template <int dim, int p_dim>
void Large_classic<dim, p_dim>::evolve(int cyc){
  
  //Output time=0 result
  n_qp_e=1;
  results(-1);
  results_table.add_value("step", 0.0);
  results_table.add_value("sigmaBar_xx", 0.0);
  results_table.add_value("sigmaBar_yy", 0.0);
  results_table.add_value("sigmaBar_xy", 0.0);
  if (p_dim==3){
    results_table.add_value("sigmaBar_zz", 0.0);
    results_table.add_value("sigmaBar_yz", 0.0);
    results_table.add_value("sigmaBar_zx", 0.0);
    results_table.add_value("sigmaBar_vm", 0.0);
  }
  save_n++;

  //Time evolution
  for (int i=0; i<ns; i++){

    //Newton Iteration
    cout << endl << " Newton Iteration: ";
    double time = (1.0+i)/ns;
    double l2norm_0;

    for (int iter=0; iter<max_iter; iter++){

      //Assemble the global Stiffness matrix and Right hand side
      R=0; K=0;
      assemble(time, iter);

      //Modify K and R for essential boundary conditions
      Q_stp = 0;
      apply_ebc(iter);
      
      //Solve: (the resulting linear equations: R=KQ )
      solve();

      //Compute the values of constrained nodes from the values of the unconstrained ones
      if (iter==0)
        cn_hn.distribute(Q_stp);
      else
        cn_hn2.distribute(Q_stp);
      
      //Update the results
      Q += Q_stp;

      //Convergance check
      if (convergance_check(iter, l2norm_0) == true)
        break;
      AssertThrow(iter+1 < max_iter, ExcMessage("Newton Iteration was not converged for maximum iteration!"));

    }

    
    //Print Status
    cout << " Completed: % " << (100*((cyc+1.0)/nc)*((i+1.0)/ns)) <<endl;

    //Print results
    if ((i+1) % save_r == 0 || i==0){
      
      results(i);

      //Print v/V
      cout<<" v / V : "<<volume/Volume<<endl;

      save_n++;

    }
  
  }

}

/*---------------------------------------------------------------
-------           O T H E R    F U N C T I O N S          -------
---------------------------------------------------------------*/

/*  REFINE MESH  */
template <int dim, int p_dim>
void Large_classic<dim, p_dim>::refine_mesh(){
  Vector<float> eepc(n_elm_t); //Estimated error per cell
  KellyErrorEstimator<dim>::estimate(dof_h, QGauss<dim - 1>(n_qp_d), {}, Q, eepc);

  GridRefinement::refine_and_coarsen_fixed_number(mesh, eepc, 0.3, 0.03);
  
  mesh.execute_coarsening_and_refinement();

}


/*  ASSEMBLING STIFNESS MATRIX  AND R */
template <int dim, int p_dim>
void Large_classic<dim, p_dim>::assemble(double time, double iter){

  // Determine the gaussian points and weight factors for elements and faces
  QGauss<dim>   egauss(n_qp_d);
  QGauss<dim-1> fgauss(n_qp_d);
  n_qp_e =      egauss.size();
  n_qp_f =      fgauss.size();

  // Calculate the Cartesian derivatives of shape functions and area for all elements and faces in the solution
  FEValues<dim>     fe_values_e(fe, egauss, update_values | update_gradients         | update_quadrature_points | update_JxW_values);
  FEFaceValues<dim> fe_values_f(fe, fgauss, update_values | update_quadrature_points | update_JxW_values);

  // Define element stiffness and rhs
  FullMatrix<double> k_el(n_dof_e, n_dof_e);
  Vector<double>     r_el(n_dof_e);

  // Define constants and variables
  //For element: 
  double mu, kappa;
  //For QPs:
  double detF;
  Tensor<2, p_dim> tau, F_inv, sigma, S, F_t;
  Tensor<4, p_dim> Jc, C;

  // Define sigmaBarVol
  Tensor<2, p_dim> sigmaBarVol;
  sigmaBar.resize((p_dim-1)*3);

  // Create a list for storing DOFs of the element
  vector<types::global_dof_index> dofs_el(n_dof_e);

  // Initialize data storage for element and QP
  for (const auto &cell : dof_h.active_cell_iterators()){
    Element_ds.initialize(cell, 1);
    QP_ds.initialize(cell, n_qp_e);
  }

  // Create lists for storing the b, t and P [Loads]
  vector<Tensor<1, dim>>  b(n_qp_e);
  vector<Tensor<1, dim>>  t(n_qp_f);
  // Tensor<1, dim>       P;

  // Finding P point for applying point load
  // bool PpointFound=false;
  
  // Extractors that get at the displacement components of vector-valued shape functions
  FEValuesExtractors::Vector displacement(0);

  // Reset volume and sigma
  volume = 0;
  sigmaBarVol = 0;
  for (const auto &cell : dof_h.active_cell_iterators()) ///THIS NEEDS UPDATE
    sigma_e[cell->active_cell_index()] = 0;

  for (const auto &cell : dof_h.active_cell_iterators()){

    // Initialize FE values of the element
    fe_values_e.reinit(cell);
    
    // Initialize Data storage of the element
    const vector<shared_ptr<StoreData::SD2<dim, p_dim>>> esgd  = Element_ds.get_data(cell);
    const vector<shared_ptr<StoreData::SD2<dim, p_dim>>> qpsgd = QP_ds.get_data(cell);
    
    // Store element data
    int mat;
    cell->material_id()==0 ? mat=0 : mat=cell->material_id()-1;
    esgd[0]->setup_ed(G, k, v, mat, s_typ);

    // Getting the values of intended data for the element
    mu   = esgd[0]->mu_el;
    kappa= esgd[0]->kappa_el;

    // Get the list of DOFs of the element
    cell->get_dof_indices(dofs_el);

    // Get the u value of nodes of the element [This one also can be improved]
    vector<double> Q_nodes(n_dof_e);
    for (int I=0; I<n_dof_e; I++)
      Q_nodes[I] = Q(dofs_el[I]);
    
    // Initialize element stiffness and r
    k_el = 0;
    r_el = 0;

    // Get the values of b at the quadrature points of the element
    LoadAndBCs::b_calculate(fe_values_e.get_quadrature_points(), b, time);

    // Calculate element stiffness and R over each QP of the element
    for (int q=0; q<n_qp_e; q++){
      
      // Store QP data
      qpsgd[q]->setup_qpd(n_dof_e, fe_values_e, q, Q_nodes, s_typ, kappa, mu);

      // Getting the values of intended data for the QP
      detF  = qpsgd[q]->detF_qp;
      Jc    = qpsgd[q]->Jc_iso_qp + detF * qpsgd[q]->c_vol_qp;
      C     = qpsgd[q]->C;
      //
      tau   = qpsgd[q]->tau_qp;
      sigma = qpsgd[q]->sigma_qp;
      S     = qpsgd[q]->S_qp;
      F_inv = qpsgd[q]->F_inv_qp;
      F_t   = qpsgd[q]->F_t_qp;

      // Create Shape functions and derivatives
      Tensor<1, p_dim> N_Iu;
      Tensor<2, p_dim> GrdN_Iu, grdN_Iu, FtGrdN_Iu, FtGrdN_Ju;
      Tensor<2, p_dim> GrdN_Ju, grdN_Ju;

      // Loop over dofs
      for (int I=0; I<n_dof_e; I++){
        

        // Get the shape functions and derivatives
        if constexpr (dim==2){
          Tensor<1, 2> N_Iu_2d    = fe_values_e[displacement].value(I, q);
          Tensor<2, 2> GrdN_Iu_2d = fe_values_e[displacement].gradient(I, q);
          for (int i=0; i<2; i++){
            N_Iu[i] = N_Iu_2d[i];
            for (int j=0; j<2; j++)
              GrdN_Iu[i][j] = GrdN_Iu_2d[i][j];
          }
          if constexpr (p_dim==3)
            GrdN_Iu[2][2] = 0.0;
        }
        else if constexpr (dim==3){
          N_Iu    = fe_values_e[displacement].value(I, q);
          GrdN_Iu = fe_values_e[displacement].gradient(I, q);
        }
        grdN_Iu = GrdN_Iu * F_inv;
        FtGrdN_Iu = F_t * GrdN_Iu;

        
        for (int J=0; J<n_dof_e; J++){
          
          
          // Get the shape functions and derivatives
          if constexpr (dim==2){
            Tensor<2, 2> GrdN_Ju_2d = fe_values_e[displacement].gradient(J, q);
            for (int i=0; i<2; i++)
              for (int j=0; j<2; j++)
                GrdN_Ju[i][j] = GrdN_Ju_2d[i][j];
            
            if constexpr (p_dim==3)
              GrdN_Ju[2][2] = 0.0;
          }
          else if constexpr (dim==3){
            GrdN_Ju = fe_values_e[displacement].gradient(J, q);
          }
          grdN_Ju = GrdN_Ju * F_inv;
          FtGrdN_Ju = F_t * GrdN_Ju;

          //-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
          //
          // 3 forms of L are considered:
          //
          // L1 = Intg ( Grad du . Grad Du : S     +  Ft Grad du : C : Ft Grad Du ) dV
          //
          // L2 = Intg ( grad du . grad Du : sigma +     grad du : c : grad Du ) dv
          //
          // L3 = Intg ( grad du . grad Du : tau   +     grad du :Jc : grad Du ) dV
          //
          //-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

          //K of element
          //L1:
          if (dim==2 and p_dim==2 and s_typ==1)
            k_el(I,J) += ( double_contract<0,0,1,1>(FtGrdN_Iu, double_contract<2,0,3,1>(C, FtGrdN_Ju))  + (  double_contract<0,0,1,1>(contract<0,0>(GrdN_Iu, GrdN_Ju), S) ) ) * z * fe_values_e.JxW(q);
          //L2 & L3
          if (dim!=2 or p_dim!=2 or s_typ!=1){
            k_el(I,J) += ( double_contract<0,0,1,1>(grdN_Iu, double_contract<2,0,3,1>(Jc/detF, grdN_Ju))  + (  double_contract<0,0,1,1>(contract<0,0>(grdN_Iu, grdN_Ju), sigma) ) ) * z * detF*fe_values_e.JxW(q);
            // k_el(I,J) += ( double_contract<0,0,1,1>(grdN_Iu, double_contract<2,0,3,1>(Jc, grdN_Ju))  + (  double_contract<0,0,1,1>(contract<0,0>(grdN_Iu, grdN_Ju), tau) ) ) * z * fe_values_e.JxW(q);
          }
          
        }


        Tensor<1, dim> N_Iu_dim = fe_values_e[displacement].value(I, q);
        
        //R of element
        //L1:
        if (dim==2 and p_dim==2 and s_typ==1)
          r_el(I) -= (double_contract<0,0,1,1>(FtGrdN_Iu, S) + b[q] * N_Iu_dim ) * z * fe_values_e.JxW(q);
        //L2 & L3
        if (dim!=2 or p_dim!=2 or s_typ!=1){
          r_el(I) -= (double_contract<0,0,1,1>(grdN_Iu, sigma) + b[q] * N_Iu_dim ) * z * detF*fe_values_e.JxW(q);
          // r_el(I) -= (double_contract<0,0,1,1>(grdN_Iu, tau) + b[q] * N_Iu_dim ) * z * fe_values_e.JxW(q);
        }

      }

      // Calculate Stress [This is sum of sigma which will be divided by n_qp_e when outputing results]
      sigma_e[cell->active_cell_index()] += sigma;

      //Calculate volume
      volume += detF * z * fe_values_e.JxW(q);

      //sigmaBarVol
      sigmaBarVol += sigma * fe_values_e.JxW(q);

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
              t[iq][iiq] = time * bc.value[iiq];
          
          // Calculate R(t contribution) over each QP of the face
          for (int qf = 0; qf<n_qp_f; qf++){
            for (int I = 0; I<n_dof_e; I++){
              Tensor<1, dim> N_Iu_dim = fe_values_f[displacement].value(I, qf);            

              r_el(I) += t[qf] * N_Iu_dim * fe_values_f.JxW(qf);
            
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
    if (iter==0)
      cn_hn.distribute_local_to_global(k_el, r_el, dofs_el, K, R);
    else
      cn_hn2.distribute_local_to_global(k_el, r_el, dofs_el, K, R);


  }

  //sigmaBar
  sigmaBar[0] = sigmaBarVol[0][0]/GridTools::volume(mesh);
  sigmaBar[1] = sigmaBarVol[1][1]/GridTools::volume(mesh);
  sigmaBar[2] = sigmaBarVol[0][1]/GridTools::volume(mesh);
  if (p_dim==3){
    sigmaBar[3] = sigmaBarVol[2][2]/GridTools::volume(mesh);
    sigmaBar[4] = sigmaBarVol[0][2]/GridTools::volume(mesh);
    sigmaBar[5] = sigmaBarVol[1][2]/GridTools::volume(mesh);
    sigmaBar[6] = pow(0.5*( pow(sigmaBar[0]-sigmaBar[1],2) + pow(sigmaBar[1]-sigmaBar[3],2) + pow(sigmaBar[3]-sigmaBar[0],2) ) 
                    + 3.0*( pow(sigmaBar[2],2) + pow(sigmaBar[4],2) + pow(sigmaBar[5],2) ), 0.5);
  }

}


/*  MARK THE BOUNDARY  */
template <int dim, int p_dim>
void Large_classic<dim, p_dim>::mark_bc(){
  
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
  for (auto &bc : BCs)
    if (bc.type==LoadAndBCs::pbc)
      checkcorners = true;

  //Mark the boundaries
  for (const auto &cell : mesh.active_cell_iterators())
    if (cell->at_boundary()){

      //Check if cell is at corner [for PBC only]
      skipcell = false;
      if (checkcorners==true){
      if constexpr (dim==2){
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
      
      //check for extra bcs
      if (skipcell==false)
        for (auto &bc : BCs){
          if (bc.position==7)
            for (const auto &face : cell->face_iterators())
              if (face->at_boundary())
                if (std::fabs(face->center()(1) - (y_max)) < 1e-9 && face->center()(0) < (x_max/2.0) )
                  face->set_boundary_id(7);

          if (dim==3 && bc.position==8)
            for (const auto &face : cell->face_iterators())
              if (face->at_boundary())
                if (std::fabs(face->center()(1) - (y_max)) < 1e-9 && face->center()(0) < (x_max/2.0) && face->center()(2) < (z_max/2.0) )
                  face->set_boundary_id(8);
          }

    }
  
  }
  
}


/*  APPLY Essential BCs  */
template <int dim, int p_dim>
void Large_classic<dim, p_dim>::apply_ebc(int iter){

  // Extractors
  const FEValuesExtractors::Vector displacement(0);
  const FEValuesExtractors::Scalar x_displacement(0);
  const FEValuesExtractors::Scalar y_displacement(1);
  const FEValuesExtractors::Scalar z_displacement(2);

  // Apply Dirichlet BC
  for (auto &bc : BCs)
    if (bc.type==LoadAndBCs::dbc){
      if (bc.value == vector<double> {0.0,0.0} || bc.value == vector<double> {0.0,0.0,0.0}){        
        map<types::global_dof_index, double> cn_bc1, cn_bc2, cn_bc3;
        if (bc.x_free==false){
          VectorTools::interpolate_boundary_values(dof_h, bc.position, Functions::ZeroFunction<dim>(dim), cn_bc1, fe.component_mask(x_displacement));
          MatrixTools::apply_boundary_values(cn_bc1, K, Q_stp, R);}
        if (bc.y_free==false){
          VectorTools::interpolate_boundary_values(dof_h, bc.position, Functions::ZeroFunction<dim>(dim), cn_bc2, fe.component_mask(y_displacement));
          MatrixTools::apply_boundary_values(cn_bc2, K, Q_stp, R);}
        if (bc.z_free==false && dim==3){
          VectorTools::interpolate_boundary_values(dof_h, bc.position, Functions::ZeroFunction<dim>(dim), cn_bc3, fe.component_mask(z_displacement));
          MatrixTools::apply_boundary_values(cn_bc3, K, Q_stp, R);}
      }
      else{
        map<types::global_dof_index, double> cn_bc1, cn_bc2, cn_bc3;
        if (iter==0){
          vector <double> valuepertime(bc.value.size());
          for (long unsigned int i=0; i<bc.value.size(); i++)
            valuepertime[i] = bc.value[i]/ns;
          
          if (bc.x_free==false){
            VectorTools::interpolate_boundary_values(dof_h, bc.position, Functions::ConstantFunction<dim>(valuepertime), cn_bc1, fe.component_mask(x_displacement));
            MatrixTools::apply_boundary_values(cn_bc1, K, Q_stp, R);}
          if (bc.y_free==false){
            VectorTools::interpolate_boundary_values(dof_h, bc.position, Functions::ConstantFunction<dim>(valuepertime), cn_bc2, fe.component_mask(y_displacement));
            MatrixTools::apply_boundary_values(cn_bc2, K, Q_stp, R);}
          if (bc.z_free==false && dim==3){
            VectorTools::interpolate_boundary_values(dof_h, bc.position, Functions::ConstantFunction<dim>(valuepertime), cn_bc3, fe.component_mask(z_displacement));
            MatrixTools::apply_boundary_values(cn_bc3, K, Q_stp, R);}
        }
        else{
          if (bc.x_free==false){
            VectorTools::interpolate_boundary_values(dof_h, bc.position, Functions::ZeroFunction<dim>(dim), cn_bc1, fe.component_mask(x_displacement));
            MatrixTools::apply_boundary_values(cn_bc1, K, Q_stp, R);}
          if (bc.y_free==false){
            VectorTools::interpolate_boundary_values(dof_h, bc.position, Functions::ZeroFunction<dim>(dim), cn_bc2, fe.component_mask(y_displacement));
            MatrixTools::apply_boundary_values(cn_bc2, K, Q_stp, R);}
          if (bc.z_free==false && dim==3){
            VectorTools::interpolate_boundary_values(dof_h, bc.position, Functions::ZeroFunction<dim>(dim), cn_bc3, fe.component_mask(z_displacement));
            MatrixTools::apply_boundary_values(cn_bc3, K, Q_stp, R);}
        }
        
      }
    }


  // Macro-Strain Linear BC
  for (auto &bc : BCs)
    if (bc.type==LoadAndBCs::ldb){
      if (iter==0){
        map<types::global_dof_index, double> cn_bc1;
        vector <double> valuepertime(bc.value.size());
        
        // F - I 
        for (long unsigned int i=0; i<bc.value.size(); i++){
          if (i<dim)
            valuepertime[i] = (bc.value[i]-1.0)/ns;
          else
            valuepertime[i] = bc.value[i]/ns;}

        // Enforcing J=1 by changing F_(dim dim)
        // if (v>0.47){
        // double g1 = ( 1.0 + ((time-(1.0/ns))*bc.value[2])*((time-(1.0/ns))*bc.value[3]) ) / ( (time-(1.0/ns))*(bc.value[0]-1.0) + 1.0 );
        // double g2 = ( 1.0 + (time*bc.value[2])*(time*bc.value[3]) ) / ( time*(bc.value[0]-1.0) + 1.0 );
        // valuepertime[dim-1] = g2 - g1;}
        
        VectorTools::interpolate_boundary_values(dof_h, 0, LoadAndBCs::MacroStrainLBC<dim>(valuepertime), cn_bc1, fe.component_mask(displacement));
        MatrixTools::apply_boundary_values(cn_bc1, K, Q_stp, R);
      }
      else{
        map<types::global_dof_index, double> cn_bc1;
        VectorTools::interpolate_boundary_values(dof_h, 0, Functions::ZeroFunction<dim>(dim), cn_bc1, fe.component_mask(displacement));
        MatrixTools::apply_boundary_values(cn_bc1, K, Q_stp, R);
      }
    }

}


/*  SOLVE THE LINEAR SYSTEMS  */
template <int dim, int p_dim>
void Large_classic<dim, p_dim>::solve(){
  
  // [I] Direct Solver
  if (lsolver==1){
    SparseDirectUMFPACK Direct_Solver;
    Direct_Solver.initialize(K);
    Direct_Solver.vmult(Q_stp, R);
  }

  // [II] Iterative CG Solver
  else if (lsolver==2){
    auto &M = K.block(0, 0);
    auto &F = R.block(0);
    auto &U = Q_stp.block(0);
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

}


/*  CHECK THE CONVERGANCE  */
template <int dim, int p_dim>
bool Large_classic<dim, p_dim>::convergance_check(int iter, double &l2norm_0){
  
  //Print Status 
  cout << " " << iter+1 <<":";
  
  //Initialize
  bool converged = false;
  vector<unsigned int> block_component(dim, 0);
  vector<types::global_dof_index> n_dof_vc = DoFTools::count_dofs_per_fe_block(dof_h, block_component);
  int n_u_dof = n_dof_vc[0];  //Number of u degree of freedoms
  BlockVector<double> error_ud;
  error_ud.reinit(n_dof_vc);
  error_ud.collect_sizes();
  
  //Get
  for (int i = 0; i < n_u_dof; ++i)
    if (!cn_hn.is_constrained(i))
        error_ud(i) = Q_stp.block(0)[i];
  if (iter==0)
    l2norm_0 = error_ud.l2_norm();

  //Print Normalized Error
  std::cout<<std::abs(error_ud.l2_norm()/l2norm_0)<<" ..."; 
  
  //Check!
  if (std::abs(error_ud.l2_norm()/l2norm_0)<tol){
    cout<<" Converged!"<<endl;
    converged = true;
  }
  
  return converged;

}


/*  PRINT RESULTS  */
template <int dim, int p_dim>
void Large_classic<dim, p_dim>::results(int i){
  
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
  Vector<double> sigmaxx(n_elm_t), sigmayy(n_elm_t), sigmaxy(n_elm_t), sigmazz(n_elm_t), sigmayz(n_elm_t), sigmazx(n_elm_t), sigmavm(n_elm_t);
  switch (p_dim){
    case 1:
      Assert(dim >= 2, ExcNotImplemented());
      break;
    case 2:
      for (int j=0; j<n_elm_t; j++){
        sigmaxx(j)=sigma_e[j][0][0]/((double) n_qp_e);
        sigmayy(j)=sigma_e[j][1][1]/((double) n_qp_e);
        sigmaxy(j)=sigma_e[j][0][1]/((double) n_qp_e);
        if (s_typ==0)
          sigmavm(j)=pow(0.75*( pow(sigmaxx(j)-sigmayy(j),2)) +  3.0*pow(sigmaxy(j),2), 0.5);
        else if (s_typ==1)
          sigmavm(j)=pow(0.5*( pow(sigmaxx(j)-sigmayy(j),2) + pow(sigmayy(j),2) + pow(sigmaxx(j),2) ) 
                     + 3.0*( pow(sigmaxy(j),2) ), 0.5);
      }
      res.add_data_vector(sigmaxx, "sigma_xx");
      res.add_data_vector(sigmayy, "sigma_yy");
      res.add_data_vector(sigmaxy, "sigma_xy");
      res.add_data_vector(sigmavm, "sigma_vonMises");
      break;
    case 3:
      for (int j=0; j<n_elm_t; j++){
        sigmaxx(j)=sigma_e[j][0][0]/((double) n_qp_e);
        sigmayy(j)=sigma_e[j][1][1]/((double) n_qp_e);
        sigmazz(j)=sigma_e[j][2][2]/((double) n_qp_e);
        sigmaxy(j)=sigma_e[j][0][1]/((double) n_qp_e);
        sigmayz(j)=sigma_e[j][1][2]/((double) n_qp_e);
        sigmazx(j)=sigma_e[j][2][0]/((double) n_qp_e);
        sigmavm(j)=pow(0.5*( pow(sigmaxx(j)-sigmayy(j),2) + pow(sigmayy(j)-sigmazz(j),2) + pow(sigmazz(j)-sigmaxx(j),2) ) 
                     + 3.0*( pow(sigmaxy(j),2) + pow(sigmayz(j),2) + pow(sigmazx(j),2) ), 0.5);
      }
      res.add_data_vector(sigmaxx, "sigma_xx");
      res.add_data_vector(sigmayy, "sigma_yy");
      res.add_data_vector(sigmazz, "sigma_zz");
      res.add_data_vector(sigmaxy, "sigma_xy");
      res.add_data_vector(sigmayz, "sigma_yz");
      res.add_data_vector(sigmazx, "sigma_zx");
      res.add_data_vector(sigmavm, "sigma_vonMises");
      break;
  }

  res.build_patches(fe.degree);
  
  string k;
  if (dim==3)
    k = "3d";
  else if (dim==2){
    if (p_dim==2 && s_typ==1)
      k = "pstress";
    else if (p_dim==3 && s_typ==2)
      k = "pstrain";
    else if (p_dim==2 && s_typ==0)
      k = "2d";
  }
  //Output results to a vtk file
  string a=std::to_string(save_n);
  std::ofstream output(output_name+"_"+k+"_"+a+".vtk");
  res.write_vtk(output);
  //
  //Print and Save sigmaBar
  if (i>-1){
    cout<<" sigmaBar: [sigma_xx: "<<sigmaBar[0]<<", sigma_yy: "<<sigmaBar[1]<<", sigma_xy: "<<sigmaBar[2]<<"]"<<endl;
    results_table.add_value("step", ((i+1.0)/ns) );
    results_table.add_value("sigmaBar_xx", sigmaBar[0]);
    results_table.add_value("sigmaBar_yy", sigmaBar[1]);
    results_table.add_value("sigmaBar_xy", sigmaBar[2]);
    if (p_dim==3){
      results_table.add_value("sigmaBar_zz", sigmaBar[3]);
      results_table.add_value("sigmaBar_yz", sigmaBar[4]);
      results_table.add_value("sigmaBar_zx", sigmaBar[5]);
      results_table.add_value("sigmaBar_vm", sigmaBar[6]);
    }
  }
  //Output sigmaBar to a dat file
  if (i+1==ns){
    std::ofstream out_res(output_name+"_"+k+"_"+"sigmabar.dat");
    results_table.write_text(out_res, TableHandler::TextOutputFormat::simple_table_with_separate_column_description);
    out_res.close();
  }

}


}
MADEAL_NAMESPACE_CLOSE