/*---------------------------------------------------------------
-------                L O A D S and B C s                -------
---------------------------------------------------------------*/

// Include Deal II headers
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q1.h>

// NameSpaces
using std::vector;


MADEAL_NAMESPACE_OPEN
namespace LoadAndBCs{


// BC structure
enum BCtype {free, dbc, nbc, ldb, pbc};
enum BCposition {left=1, right=2, bottom=3, top=4, back=5, front=6, //Mains
                 halftop=7, quartertop=8}; //Extras
//
struct BC{
  public:
  BCtype          type;
  BCposition      position;
  vector<double>  value;
  bool            x_free=false;
  bool            y_free=false;
  bool            z_free=false;
};


/*  b = BODY FORCE */
// b1: Base
template <int dim>
void b_calculate(const vector<Point<dim>> &qps, vector<Tensor<1, dim>> &b){
  
  int n_qp_e = qps.size();

  for (int i=0; i<n_qp_e; i++)
    for (int ii=0; ii<dim; ii++)
      b[i][ii]=0;

}

// b2: With time
template <int dim>
void b_calculate(const vector<Point<dim>> &qps, vector<Tensor<1, dim>> &b, double &time){
  
  int n_qp_e = qps.size();

  for (int i=0; i<n_qp_e; i++)
    for (int ii=0; ii<dim; ii++)
      b[i][ii]=0.0*time;

}


/* P = POINT LOAD */
// P1: Base
// template <int dim>
// void P_calculate(Tensor<1, dim> &P){

//   //for (int i=0; i<dim; i++)
//     //P[i] = 1000.0;
  
//   P[0]=0.0;
//   P[1]=1000000.0;

// }

// // P2: With time
// template <int dim>
// void P_calculate(Tensor<1, dim> &P, double &time){
  
//   P[0]=0.0*time;
//   P[1]=0.0*time;

// }


/* MACRO-STRAIN LINEAR BC */
template <int dim>
class MacroStrainLBC:public Function<dim>{
public:
  MacroStrainLBC(vector<double> MacStrain)
    :Function<dim>(dim), MStrain(MacStrain)
  {}

  // For (n_ev)-field formulation (dim+n_ev)
  MacroStrainLBC(vector<double> MacStrain, int n_ev)
    :Function<dim>(dim+n_ev), MStrain(MacStrain)
  {}

  vector<double> MStrain;
  
  virtual void vector_value(const Point<dim> &p, Vector<double> &value) const override;
};
//
template <int dim>
void MacroStrainLBC<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const{
  values(0) = MStrain[0]*p(0) + MStrain[2]*p(1);
  values(1) = MStrain[3]*p(0) + MStrain[1]*p(1);
}


/* MACRO-STRAIN PERIODIC BC */
template <int dim>
void MacroStrainPBC(vector<double> &MacStrain, double L_m, double W_m, DoFHandler<dim> &dof_h, FESystem<dim> &fe, AffineConstraints<double> &cn_hn, bool full){

  const FEValuesExtractors::Vector displacement(0);
  const FEValuesExtractors::Scalar x_displacement(0);
  const FEValuesExtractors::Scalar y_displacement(1);

  // Periodicity on Boundaries (except corners)
  DoFTools::make_periodicity_constraints<dim, dim>(dof_h, 1, 2, 0, cn_hn, fe.component_mask(displacement));
  DoFTools::make_periodicity_constraints<dim, dim>(dof_h, 3, 4, 1, cn_hn, fe.component_mask(displacement));
  
  // Inhomogeneity on Boundaries (except corners)
  if (full==true){
  IndexSet right_x_dofs, right_y_dofs, left_x_dofs, left_y_dofs;
  IndexSet up_x_dofs,    up_y_dofs,    down_x_dofs, down_y_dofs;
  std::set<types::boundary_id> right_boundary, left_boundary, up_boundary, down_boundary;
  right_boundary.insert(2);
  left_boundary.insert(1);
  up_boundary.insert(4);
  down_boundary.insert(3);
  DoFTools::extract_boundary_dofs(dof_h, fe.component_mask(x_displacement), right_x_dofs, right_boundary);
  DoFTools::extract_boundary_dofs(dof_h, fe.component_mask(y_displacement), right_y_dofs, right_boundary);
  DoFTools::extract_boundary_dofs(dof_h, fe.component_mask(x_displacement), left_x_dofs,  left_boundary);
  DoFTools::extract_boundary_dofs(dof_h, fe.component_mask(y_displacement), left_y_dofs,  left_boundary);
  DoFTools::extract_boundary_dofs(dof_h, fe.component_mask(x_displacement), up_x_dofs,    up_boundary);
  DoFTools::extract_boundary_dofs(dof_h, fe.component_mask(y_displacement), up_y_dofs,    up_boundary);
  DoFTools::extract_boundary_dofs(dof_h, fe.component_mask(x_displacement), down_x_dofs,  down_boundary);
  DoFTools::extract_boundary_dofs(dof_h, fe.component_mask(y_displacement), down_y_dofs,  down_boundary);
  //
  for(auto i=right_x_dofs.begin(); i!=right_x_dofs.end(); i++)
    if (cn_hn.is_constrained(*i))
      cn_hn.set_inhomogeneity(*i, MacStrain[0]*L_m);
  for(auto i=right_y_dofs.begin(); i!=right_y_dofs.end(); i++)
    if (cn_hn.is_constrained(*i))
      cn_hn.set_inhomogeneity(*i, MacStrain[3]*W_m);
  for(auto i=up_x_dofs.begin(); i!=up_x_dofs.end(); i++)
    if (cn_hn.is_constrained(*i))
      cn_hn.set_inhomogeneity(*i, MacStrain[2]*L_m);
  for(auto i=up_y_dofs.begin(); i!=up_y_dofs.end(); i++)
    if (cn_hn.is_constrained(*i)) 
      cn_hn.set_inhomogeneity(*i, MacStrain[1]*W_m);

  for(auto i=left_x_dofs.begin(); i!=left_x_dofs.end(); i++)
    if (cn_hn.is_constrained(*i))
      cn_hn.set_inhomogeneity(*i, -MacStrain[0]*L_m);
  for(auto i=left_y_dofs.begin(); i!=left_y_dofs.end(); i++)
    if (cn_hn.is_constrained(*i))
      cn_hn.set_inhomogeneity(*i, -MacStrain[3]*W_m);
  for(auto i=down_x_dofs.begin(); i!=down_x_dofs.end(); i++)
    if (cn_hn.is_constrained(*i))
      cn_hn.set_inhomogeneity(*i, -MacStrain[2]*L_m);
  for(auto i=down_y_dofs.begin(); i!=down_y_dofs.end(); i++)
    if (cn_hn.is_constrained(*i))
      cn_hn.set_inhomogeneity(*i, -MacStrain[1]*W_m);
  
  }


  // Periodicity of Corners
  Point<dim>  c_dl(0,0), c_dr(L_m,0), c_ul(0,W_m), c_ur(L_m,W_m); //Corner Points (l,r,u,d: left,right,up,down)
  vector<int> cID_dl, cID_dr, cID_ul, cID_ur; //Corners IDs
  //
  std::map<types::global_dof_index, Point<dim>> support_points_map;
  for (const auto &cell : dof_h.active_cell_iterators())
    for (const auto v : cell->vertex_indices())
      for (int d=0; d<dim; d++)
        support_points_map[cell->vertex_dof_index(v,d)]=cell->vertex(v);
  

  for (auto const &[key, val] : support_points_map){

    if (val.distance(c_dl)<1e-12)
      cID_dl.emplace_back(key);
    else if (val.distance(c_dr)<1e-12)
      cID_dr.emplace_back(key);
    else if (val.distance(c_ul)<1e-12)
      cID_ul.emplace_back(key);
    else if (val.distance(c_ur)<1e-12)
      cID_ur.emplace_back(key);
    
  }
  //A
  cn_hn.add_line(cID_dl[0]);
  cn_hn.add_line(cID_dl[1]);
  if (full==true){
    cn_hn.set_inhomogeneity(cID_dl[0], 0);
    cn_hn.set_inhomogeneity(cID_dl[1], 0);}
  //D
  cn_hn.add_line(cID_dr[0]);
  cn_hn.add_line(cID_dr[1]);
  if (full==true){
  cn_hn.add_entry(cID_dr[0], cID_dl[0], 1);
  cn_hn.add_entry(cID_dr[1], cID_dl[1], 1);
  cn_hn.set_inhomogeneity(cID_dr[0], MacStrain[0]*L_m);
  cn_hn.set_inhomogeneity(cID_dr[1], MacStrain[3]*W_m);}
  //B
  cn_hn.add_line(cID_ul[0]);
  cn_hn.add_line(cID_ul[1]);
  if (full==true){
  cn_hn.add_entry(cID_ul[0], cID_dl[0], 1);
  cn_hn.add_entry(cID_ul[1], cID_dl[1], 1);
  cn_hn.set_inhomogeneity(cID_ul[0], MacStrain[2]*L_m);
  cn_hn.set_inhomogeneity(cID_ul[1], MacStrain[1]*W_m);}
  //C
  cn_hn.add_line(cID_ur[0]);
  cn_hn.add_line(cID_ur[1]);
  if (full==true){
  cn_hn.add_entry(cID_ur[0], cID_ul[0], 1);
  cn_hn.add_entry(cID_ur[1], cID_ul[1], 1);
  cn_hn.set_inhomogeneity(cID_ur[0], (MacStrain[0])*L_m);
  cn_hn.set_inhomogeneity(cID_ur[1], (MacStrain[3])*W_m);}
  
}


}

MADEAL_NAMESPACE_CLOSE