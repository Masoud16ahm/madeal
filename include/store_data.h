/*---------------------------------------------------------------
-------                S T O R E   D A T A                -------
---------------------------------------------------------------*/

// Include Deal II headers
#include <deal.II/physics/elasticity/standard_tensors.h>
#include <deal.II/physics/elasticity/kinematics.h>

// NameSpaces
using std::vector;


MADEAL_NAMESPACE_OPEN
namespace StoreData{


//Newton-Raphson method for Plane-stress large deformation
class Newton_pstress{
public:
  //Constructor
  Newton_pstress(double kappa, double mu, double F11, double F12, double F21, double F22):
    kappa(kappa), mu(mu)
  {
    //Get C components
    C11 = F11*F11 + F21*F21;
    C12 = F11*F12 + F21*F22;
    C22 = F12*F12 + F22*F22;
  }

  //Main function: Calculate C33 and dC33/dC
  double cal_F33_dC33dC(SymmetricTensor<2, 2> &dC33_dCbar){
    
    //I:calculate C33
    double C33_old   = 1.0; //Initial guess
    double tolerance = 1e-5;
    int    max_iter  = 27;
    int    iteration = 1;

    //Newton iteration
    while (iteration < max_iter){
      
      double f  = get_f(C33_old); 
      double fp = get_fp(C33_old);
      //
      C33 = C33_old - f / fp;
      if (C33<0)
        C33=0.01;
      
      if (std::abs(C33-C33_old) < tolerance && std::abs(get_f(C33)) < 0.001){
        // std::cout << "Converged after " << iteration << " iterations." << std::endl;
        break;
      }
      
      C33_old = C33;
      iteration++;

    }
    
    AssertThrow(iteration < max_iter, ExcMessage("Newton-Raphson for plane-stress condition did not converge!"));

    //II:Calculate dC33/dC components
    double dS33_dC33 = get_fp(C33);
    get_beta();
    //
    dC33_dCbar = - beta / dS33_dC33;

    return pow(C33,0.5);

  }


private:
  //Parameters
  double kappa, mu;
  double C11, C12, C22, C33;
  SymmetricTensor<2, 2> I2 = unit_symmetric_tensor<2>();
  SymmetricTensor<2, 2> beta;

  //Functions:
  
  // MATERIAL model:
  //
  // (I,i)  << neo-Hookean >> 
  // << Psi = (mu/2) [tr(b_hat)-3] + kappa/4 [J^2 - 1 - 2lnJ] >>
  //
  // f=S33 = mu/J^(2/3) + [ k/2 [J^2-1] - [mu trC]/[3 J^(2/3)] ]/C33
  // J^2 = det(C)
  //
  // NOTE! << if changing the material model, do not forget to change the Store Data 2 for Large Def Classic 1-field: SD2 >>
  double get_f(double X){
    // X:C33
    double J2 =X*(C11*C22-C12*C12);
    return mu/pow(J2,1.0/3.0) + (0.5*kappa * (J2-1.0) - (mu*(C11+C22+X))/(3.0 * pow(J2,1.0/3.0)))/X;
  }
  //
  double get_fp(double X){
    // X:C33
    double J2=X*(C11*C22-C12*C12);
    return (9.0 * kappa * pow(J2, 1.0/3.0) + 4.0 * mu * (2.0 * C11 + 2.0 * C22 - X)) / (18.0 * X*X * pow(J2, 1.0/3.0));
  }

  //calculate beta
  void get_beta(){

    double detCbar = C11*C22-C12*C12;
    double J = sqrt(C33*detCbar);
    SymmetricTensor<2, 2> cbar; 
    cbar[0][0] = C11;
    cbar[0][1] = C12;
    cbar[1][0] = C12;
    cbar[1][1] = C22;
    auto Cbarinv   = invert(cbar);
    double trC = C11 + C22 + C33;
    auto dJ23dCbar = - (2.0/3.0) * pow(J,-5.0/3.0) * Cbarinv * detCbar; 
    auto dgdCbar   = (kappa/2.0) * C33 * Cbarinv * detCbar - (mu/3.0) * (pow(J,-2.0/3.0)*I2 + trC * dJ23dCbar);
    beta = mu * dJ23dCbar + (1.0/C33) * dgdCbar;

  }

};




// Stored Data classes [What to be stored in each element]

// Store Data 1: Small Def
template <int dim>
class SD1{
public:

  //Store mu, lambda and C
  double mu_el;
  double la_el;
  SymmetricTensor<4, dim> C_el;

  //Setup element data
  void setup_ed(vector<double> E, vector<double> v, int mat, int s_typ, double &z){

    //Calculate mu and lambda based on E and v
    mu_el = (E[mat])/(2.0*(1.0+v[mat]));
    la_el = (E[mat]*v[mat])/((1.0+v[mat])*(1.0-2.0*v[mat]));

    //Plane stress and strain 
    if ((dim==2) & (s_typ==1)){la_el=la_el*(mu_el)/(mu_el+0.5*la_el);}
    
    //Thickness for 3D model
    if (dim==3){z=1;}

    //Calculate C
    C_el = la_el * outer_product(I2, I2) + 2 * mu_el * I4;

  }

private:
  // Define the Identity tensor
  SymmetricTensor<2, dim> I2 = unit_symmetric_tensor<dim>();
  SymmetricTensor<4, dim> I4 = identity_tensor<dim>();

};



// Store Data 2: Large Def Classic 1-field
template <int dim, int p_dim>
class SD2{
public:

  //[1]Store Data in ELEMENT:
  double mu_el, kappa_el;

  //[2]Store Data in QUADRATURE POINTs:
  double detF_qp;
  Tensor<2, p_dim> F_inv_qp, F_t_qp;
  SymmetricTensor<4, p_dim> Jc_iso_qp, c_vol_qp, C_qp, C;
  SymmetricTensor<2, p_dim> S_qp, sigma_qp, tau_qp;

  //[1] Setup ELEMENT data
  void setup_ed(vector<double> G, vector<double> k, vector<double> v, int mat, int s_typ){

    //Calculate mu & kappa based on input parameters
    mu_el    = G[mat];
    kappa_el = k[mat];
    if (k[mat]==0)
      kappa_el = ( 2.0*G[mat] * (1.0+v[mat]) )/ ( 3.0*(1.0-2.0*v[mat]) );
    else
      AssertThrow(v[mat]==0 || k[mat]==0, ExcInternalError());
    
    //Error for nearly-incompressible non plane stress
    if ( (v[mat]>0.47 || (3.0*kappa_el-2.0*mu_el)/(2.0*(3.0*kappa_el+mu_el))>0.47) & (s_typ == 2 || p_dim==3 || dim==3) )
      AssertThrow(v[mat]<0.47, ExcMessage("Nearly-incompressible materials in a non Plane-stress condition don't vonverge to correct answer with Classic-formulation!! please try 3field-mixed formulation."));

  }


  //[2] Setup QUADRATURE POINTs data
  void setup_qpd(int n_dof_e, FEValues<dim> &fe_values_e, int q, vector<double> Q_nodes, int s_typ, double kappa, double mu){
    
    // Extractors that get at the displacement
    FEValuesExtractors::Vector displacement(0);

    //Reset values at qps
    Grd_u_qp    = 0.0;
    Grd_u_qp_2d = 0.0;

    //Loop over each dof
    for (int I=0; I<n_dof_e; I++){
      
      Tensor<2, dim> Grd_N_Iu = fe_values_e[displacement].gradient(I, q);
      
      if constexpr (dim==2)
        Grd_u_qp_2d += Grd_N_Iu * Q_nodes[I];
      else if constexpr (dim==3)
        Grd_u_qp += Grd_N_Iu * Q_nodes[I];
      
    }

    // Convert Grd_u_qp_2d to Grd_u_qp
    if (dim==2)
      for (int i=0; i<2; i++)
        for (int j=0; j<2; j++)
          Grd_u_qp[i][j] = Grd_u_qp_2d[i][j];

    //F33 for Planar assumptions
    if constexpr (dim==2){
      
      //Plane Stress
      if (s_typ==1){
        double F11 = Grd_u_qp[0][0] + 1.0;
        double F12 = Grd_u_qp[0][1];
        double F21 = Grd_u_qp[1][0]; 
        double F22 = Grd_u_qp[1][1] + 1.0;
        
        //Get F33 and dC33_dCbar 
        Newton_pstress n_pstress(kappa, mu, F11, F12, F21, F22);
        F33 = n_pstress.cal_F33_dC33dC(dC33_dCbar);
      }

      //Plane Strain
      else if (s_typ==2)
        Grd_u_qp[2][2] = 0.0;
    
    }

    //Deformation gradient F
    F_qp     = Physics::Elasticity::Kinematics::F(Grd_u_qp);
    detF_qp  = determinant(F_qp);
    F_inv_qp = invert(F_qp);
    F_t_qp   = transpose(F_qp);
    AssertThrow(detF_qp > 0.0, ExcMessage("The J=det(F) should be positive!"));


    // MATERIAL models:
    //
    // (I)  << neo-Hookean >> 
    // << Psi_iso = (mu/2) [tr(b_hat)-3] >>
    //
    // i)   nH1: << Psi_vol = kappa/4 [J^2 - 1 - 2lnJ] >>
    // ii)  nH2: << Psi_vol = kappa/2 [J - 1]^2 >>
    //
    // NOTE! << If changing the material model, do not forget to change the plane-stress calss: functions_pstress >>
    //

    //
    // (A) plane-stress case: L1
    //

    if constexpr (dim==2 and p_dim==2) 
    if (s_typ==1){

      // Material model: I, i
      //Get the needed values!
      SymmetricTensor<2, dim> Cbar_IJ = Physics::Elasticity::Kinematics::C(F_qp);
      SymmetricTensor<2, dim> Cbar_IJ_inv = invert(Cbar_IJ);
      SymmetricTensor<4, dim> dCbinv_dCbinv = Physics::Elasticity::StandardTensors<dim>::dC_inv_dC(F_qp);
      detF_qp = detF_qp*F33; // Correct J for plane-stress
      //
      double C33 = F33*F33;
      double detCbar = determinant(Cbar_IJ);
      double dJdC33  = 1.0/(2.0*detF_qp) * detCbar;
      //
      auto dJdCbar = 1.0/(2.0*detF_qp) * C33*Cbar_IJ_inv*detCbar;
      auto term    = (2.0/3.0)*mu*pow(detF_qp,-5.0/3.0) * (I2 - C33*Cbar_IJ_inv);
      auto dSdCbar = C33*mu*pow(detF_qp,-2.0/3.0)*(-dCbinv_dCbinv) - outer_product(term, dJdCbar);
      auto dSdC33  = -mu*pow(detF_qp,-2.0/3.0)*Cbar_IJ_inv - term*dJdC33;
      
      //Calculate S,sigma & C 
      S_qp = mu * pow(detF_qp,-2.0/3.0) * (I2 - C33*Cbar_IJ_inv);
      sigma_qp = symmetrize ((1.0/detF_qp) * F_qp * S_qp * transpose(F_qp));
      // 
      C = 2.0*(dSdCbar + outer_product(dSdC33, dC33_dCbar));

    }

    
    //
    // (B) Non plane-sress cases: L2 & L3
    //
    
    if (dim!=2 or p_dim!=2 or s_typ!=1){

      //Get the needed values!
      F_hat_qp = Physics::Elasticity::Kinematics::F_iso(F_qp);
      b_hat_qp = Physics::Elasticity::Kinematics::b(F_hat_qp);

      // Material model: I, i
      //p = d Psi_vol / d J
      p_qp  = (kappa/2.0) * (detF_qp - (1.0/detF_qp));
      //d p / d J
      d_p_J = (kappa/2.0) * (1.0 + 1.0/(detF_qp*detF_qp));
      //d Psi_iso / d b_hat
      d_Psi_iso_b_hat  = mu/2.0;
      //d2 Psi_iso / d b_hat d b_hat
      d2_Psi_iso_b_hat = 0;

      //Calculate tau, sigma & c 
      tau_hat_qp = 2.0 * b_hat_qp * d_Psi_iso_b_hat;
      c_hat      = 4.0*b_hat_qp * d2_Psi_iso_b_hat * b_hat_qp;
      tau_iso    = Dev * tau_hat_qp;
      tau_qp     = tau_iso + p_qp*detF_qp*I2;   
      sigma_qp   = (1.0/detF_qp) * tau_qp;
      c_vol_qp   = (p_qp + detF_qp*d_p_J)*outer_product(I2, I2) - 2.0*p_qp*I4;
      Jc_iso_qp  = Dev * c_hat * Dev +  (2.0/p_dim) * trace(tau_hat_qp) * Dev 
                   - (2.0/p_dim) * (outer_product(tau_iso, I2) + outer_product(I2, tau_iso));
    
    }

  }


private:
  //Define the Identity tensor
  SymmetricTensor<2, p_dim> I2 = unit_symmetric_tensor<p_dim>();
  SymmetricTensor<4, p_dim> I4 = identity_tensor<p_dim>();
  SymmetricTensor<4, p_dim> Dev= Physics::Elasticity::StandardTensors<p_dim>::dev_P;

  //Private variables
  double p_qp, d_p_J, d_Psi_iso_b_hat, d2_Psi_iso_b_hat, F33;
  Tensor<2, 2> Grd_u_qp_2d;
  Tensor<2, p_dim> Grd_u_qp, F_qp, F_hat_qp;
  SymmetricTensor<2, dim> dC33_dCbar;
  SymmetricTensor<2, p_dim> tau_hat_qp, tau_iso, b_hat_qp;
  SymmetricTensor<4, p_dim> c_hat;

};



// Store Data 3: Large Def mixed 3-fields
template <int dim, int p_dim>
class SD3{
public:

  //[1]Store Data in ELEMENT:
  double mu_el, kappa_el;

  //[2]Store Data in QUADRATURE POINTs:
  double detF_qp, pt_qp, Jt_qp, p_qp, d_p_Jt;
  Tensor<2, p_dim> F_inv_qp;
  SymmetricTensor<4, p_dim> Jc_iso_qp, c_vol_qp;
  SymmetricTensor<2, p_dim> tau_qp, sigma_qp;


  //[1] Setup ELEMENT data
  void setup_ed(vector<double> G, vector<double> k, vector<double> v, int mat){

    //Calculate mu & kappa based on input parameters
    mu_el    = G[mat];
    kappa_el = k[mat];
    if (k[mat]==0)
      kappa_el = ( 2.0*G[mat] * (1.0+v[mat]) ) / ( 3.0*(1.0-2.0*v[mat]) );
    else
      AssertThrow(v[mat]==0 || k[mat]==0, ExcInternalError());

  }


  //[2] Setup QUADRATURE POINTs data
  void setup_qpd(int n_dof_e, FEValues<dim> &fe_values_e, int q, vector<double> Q_nodes, int s_typ, double kappa, double mu){
    
    //Extractors that get at the displacement
    FEValuesExtractors::Vector displacement(0);
    FEValuesExtractors::Scalar ptilde(dim);
    FEValuesExtractors::Scalar Jtilde(dim+1);

    //Reset values at qps
    Jt_qp = 0.0;
    pt_qp = 0.0;
    Grd_u_qp    = 0.0;
    Grd_u_qp_2d = 0.0;
    
    //Loop over each dof
    for (int I=0; I<n_dof_e; I++){
      
      double N_Ip              = fe_values_e[ptilde].value(I, q);
      double N_Ij              = fe_values_e[Jtilde].value(I, q);
      Tensor<2, dim> Grd_N_Iu  = fe_values_e[displacement].gradient(I, q);
      
      pt_qp += N_Ip * Q_nodes[I];
      Jt_qp += N_Ij * Q_nodes[I];

      if constexpr (dim==2)
        Grd_u_qp_2d += Grd_N_Iu * Q_nodes[I];
      else if constexpr (dim==3)
        Grd_u_qp += Grd_N_Iu * Q_nodes[I];
      
    }

    // Convert Grd_u_qp_2d to Grd_u_qp
    if (dim==2)
      for (int i=0; i<2; i++)
        for (int j=0; j<2; j++)
          Grd_u_qp[i][j] = Grd_u_qp_2d[i][j];

    //Plane Strain or Stress
    if (dim==2 and p_dim==3){
      //Plane Stress
      AssertThrow(s_typ == 2, ExcMessage("Please use the Classic single-field model for the plane-stress condition!"));
      
      //Plane Strain
      Grd_u_qp[2][2] = 0.0;
    }
    
    // Get the needed values!
    F_qp     = Physics::Elasticity::Kinematics::F(Grd_u_qp);
    F_hat_qp = Physics::Elasticity::Kinematics::F_iso(F_qp);
    b_hat_qp = Physics::Elasticity::Kinematics::b(F_hat_qp);
    detF_qp  = determinant(F_qp);
    F_inv_qp = invert(F_qp);
    AssertThrow(detF_qp > 0.0, ExcInternalError());
    
    // MATERIAL model:
    //
    // (I)  << neo-Hookean >> 
    // << Psi_iso = (mu/2) [tr(b_hat)-3] >>
    //
    // i)   nH1: << Psi_vol = (kappa/4) [J^2 - 1 - 2lnJ] >>
    // ii)  nH2: << Psi_vol = kappa/2 [J - 1]^2 >>
    //

    // Material I, i:
    //p = d Psi_vol / d J_tilde
    p_qp   = (kappa/2.0) * (Jt_qp - (1.0/Jt_qp));
    //d p / d J_tilde
    d_p_Jt = (kappa/2.0) * (1.0 + 1.0/(Jt_qp*Jt_qp));
    //d Psi_iso / d b_hat
    d_Psi_iso_b_hat  = mu/2.0;
    //d2 Psi_iso / d b_hat d b_hat
    d2_Psi_iso_b_hat = 0;


    //Get the needed values!
    tau_hat_qp = 2.0 * b_hat_qp * d_Psi_iso_b_hat;
    c_hat      = 4.0*b_hat_qp * d2_Psi_iso_b_hat * b_hat_qp;
    tau_iso    = Dev * tau_hat_qp;
    tau_qp     = p_qp*detF_qp*I2 + tau_iso;
    sigma_qp   = (1.0/detF_qp) * tau_qp;
    c_vol_qp   = p_qp*(outer_product(I2, I2) - 2.0*I4);
    Jc_iso_qp  = Dev * c_hat * Dev + (2.0/p_dim) * trace(tau_hat_qp) * Dev - (2.0/p_dim) * (outer_product(tau_iso, I2) + outer_product(I2, tau_iso));

  }


private:
  //Define the Identity tensor
  SymmetricTensor<2, p_dim> I2 = unit_symmetric_tensor<p_dim>();
  SymmetricTensor<4, p_dim> I4 = identity_tensor<p_dim>();
  SymmetricTensor<4, p_dim> Dev= Physics::Elasticity::StandardTensors<p_dim>::dev_P;

  //Private variables
  Tensor<2, 2> Grd_u_qp_2d;
  double d_Psi_iso_b_hat, d2_Psi_iso_b_hat;
  Tensor<2, p_dim> Grd_u_qp, F_qp, F_hat_qp;
  SymmetricTensor<2, p_dim> tau_hat_qp, tau_iso, b_hat_qp;
  SymmetricTensor<4, p_dim> c_hat;

};




}
MADEAL_NAMESPACE_CLOSE