/*---------------------------------------------------------------
-------                     R   V   E                     -------
---------------------------------------------------------------*/


// Include C++ headers
#include <cmath>
#include <random>
#include <algorithm>


MADEAL_NAMESPACE_OPEN
namespace RVE{


// Rounding function
double roundit(double value, int decimal){
    
    if (decimal==1)
        value = std::floor(value * 10.0) / 10.0;
    else if (decimal==2)
        value = std::floor(value * 100.0) / 100.0;
    else if (decimal==3)
        value = std::floor(value * 1000.0) / 1000.0;

    return value;
}



// CNT-Reinforced RVE Class
template <int dim>
class CNT{
public:
  // Constructor
  CNT(const string &inputrve, const int &inputnr):
    ifname(inputrve), Nref(inputnr)
  {}

  // Public Functions
  void generate_rve(bool);

  // Public Variables
  Triangulation<dim> mesh;
  double L_m;  //Length of Matrix

private:
  // Private Functions
  void readinput();
  void distribute_cnts();
  void rve_linenodes();
  void mark_boundaries(bool);
  void refine_mesh();
  void assign_mid();

  // Private Variables
  ParameterHandler ph; //Handle the Input Parameters
  string ifname;       //Input file name
  int    NofF;         //Number of Fibres
  int    Nref;         //Number of refinement steps
  double AR;           //Aspect Ratio of Fibres
  double VF;           //Volume Fraction of Fibres
  double D_f;          //Diameter of Fibres
  double L_f;          //Length of Fibres
  //
  vector<Point<dim>> cnt_centres;
  vector<double> cnt_angles;
  vector<Point<dim>> cntlines, cntsolids;

};


template <int dim>
void CNT<dim>::generate_rve(bool skipmb){
  
  readinput();
  distribute_cnts();
  rve_linenodes();
  
  
  vector<unsigned int> A;
  Point<dim> P1, P2; 
  if (dim==2){
    A={1,1};
    P1={0,0};
    P2={L_m,L_m};}
  else if (dim==3){
    A={1,1,1};
    P1={0,0,0};
    P2={L_m,L_m,L_m};}
  //
  GridGenerator::subdivided_hyper_rectangle(mesh, A, P1, P2);
  

  mark_boundaries(skipmb);
  refine_mesh();
  assign_mid();
  
}


template <int dim>
void CNT<dim>::readinput(){

  //Declear entries
  ph.declare_entry("Number of Fibres", "20", Patterns::Integer());
  ph.declare_entry("Aspect Ratio of Fibres", "10.0", Patterns::Double());
  ph.declare_entry("Volume Fraction of Fibres", "0.10", Patterns::Double());
  ph.declare_entry("Diameter of Fibres", "1.0", Patterns::Double());

  //Open input file
  ph.parse_input(ifname);

  //Assign entries to variables
  NofF = ph.get_integer("Number of Fibres");
  AR   = ph.get_double("Aspect Ratio of Fibres");
  VF   = ph.get_double("Volume Fraction of Fibres");
  D_f  = ph.get_double("Diameter of Fibres");

}


// Random CNT distributer
template <int dim>
void CNT<dim>::distribute_cnts(){

  // Calculate Matrix and Fibre Lengths
  L_f = AR*D_f;
  L_m = RVE::roundit(pow(((NofF*L_f*D_f)/VF), 0.5), 2);

  // Set random engine
  std::random_device rnd_device;
  std::uniform_real_distribution<double> rnd_d1(0,M_PI);
  std::uniform_real_distribution<double> rnd_d2(0,L_m);

  // Create some vectors and other objects
  vector<Point<dim>> thePoints, newPoints;
  //
  double theta, x_f, y_f;
  double x_fX, y_fX, x_fY, y_fY, x_fC, y_fC;
  //
  int xdir_status, ydir_status;
  
  // Distribute the CNTs:

  for (int iter=0; iter<NofF; iter++){

    // Generate proper rotation angle and position
    bool Overlaping=true;
    while(Overlaping==true){
      Overlaping=false;
      theta = RVE::roundit(rnd_d1(rnd_device),1);
      x_f   = RVE::roundit(rnd_d2(rnd_device),1);
      y_f   = RVE::roundit(rnd_d2(rnd_device),1);

      // Set position of out of box Fibres as the original one
      x_fX  = x_f; y_fX  = y_f;
      x_fY  = x_f; y_fY  = y_f;
      x_fC  = x_f; y_fC  = y_f;

      //Calculate the normal vectors
      double NX = std::cos(theta);
      double NY = std::sin(theta);

      //Check if Fibre is out of box
      //Reset Status (0:in, 1:+out, -1:-out)
      xdir_status = 0;
      ydir_status = 0;
      //
      //X direction
      if ( x_f > L_m-(std::abs(NX)*(L_f/2)+D_f) ){
        xdir_status = 1;
        x_fX = x_f-L_m;}
      else if ( x_f < (std::abs(NX)*(L_f/2)+D_f) ){
        xdir_status = -1;
        x_fX = x_f+L_m;}
      //
      //Y direction
      if ( y_f > L_m-(NY*(L_f/2)+D_f) ){
        ydir_status = 1;
        y_fY = y_f-L_m;}
      else if ( y_f < (NY*(L_f/2)+D_f) ){
        ydir_status = -1;
        y_fY = y_f+L_m;}
      //
      //Corners
      if (xdir_status == 1){
        if (ydir_status == 1){
          x_fC = x_f-L_m;
          y_fC = y_f-L_m;}
        else if (ydir_status == -1){
          x_fC = x_f-L_m;
          y_fC = y_f+L_m;}
      }
      else if (xdir_status == -1){
        if (ydir_status == 1){
          x_fC = x_f+L_m;
          y_fC = y_f-L_m;}
        else if (ydir_status == -1){
          x_fC = x_f+L_m;
          y_fC = y_f+L_m;}
      }

      //Reset the list of new Points
      newPoints.clear();

      //Centre point
      newPoints.push_back(Point<dim> (x_f,y_f));
      double CPx, CPy;

      //Other points
      int alpha = int(std::floor(L_f/D_f)+1);

      for (int fib=0; fib<4; fib++){

        bool search=false;
        //
        //Main Fibre
        if (fib==0){
          CPx = x_f;
          CPy = y_f;
          search=true;
        }
        //Out of Box xdir
        else if (fib==1 && xdir_status != 0){
          CPx = x_fX;
          CPy = y_fX;
          search=true;
        }
        //Out of Box ydir
        else if (fib==2 && ydir_status != 0){
          CPx = x_fY;
          CPy = y_fY;
          search=true;
        }
        //Out of Box corner
        else if (fib==3 && (xdir_status != 0 && ydir_status != 0)){
          CPx = x_fC;
          CPy = y_fC;
          search=true;
        }
        if (search==true){
          for (int i=1; i<alpha+1; i++){
            float i_float=i*1.0;
            //one side
            double lx_dummy = CPx + (L_f/2)*(i_float/alpha)*NX;
            double ly_dummy = CPy + (L_f/2)*(i_float/alpha)*NY;
            newPoints.push_back(Point<dim> (lx_dummy,ly_dummy));
            //other side
            double rx_dummy = CPx - (L_f/2)*(i_float/alpha)*NX;
            double ry_dummy = CPy - (L_f/2)*(i_float/alpha)*NY;
            newPoints.push_back(Point<dim> (rx_dummy,ry_dummy));
          }
        }
      }

      //Check overlaping
      if (iter!=0){

        for (auto &newpoint : newPoints){

          Point<dim> TestP = newpoint;

          for (auto &thepoint : thePoints){
            if ( thepoint.distance(TestP) < 1.2*D_f){
              Overlaping=true;
              break;
            }
          }
        }
      }

    }

    //Add new Points of Fibres to Point list
    for (auto &newpoint : newPoints)
      thePoints.push_back(newpoint);

    //Save centre points and rotation angle of the fibres
    cnt_centres.push_back(Point<dim> (x_f,  y_f));
    cnt_angles.push_back(theta);

    if (xdir_status!=0){
      cnt_centres.push_back(Point<dim> (x_fX, y_fX));
      cnt_angles.push_back(theta);}
    if (ydir_status!=0){
       cnt_centres.push_back(Point<dim> (x_fY, y_fY));
       cnt_angles.push_back(theta);}
    if (xdir_status != 0 && ydir_status != 0){
       cnt_centres.push_back(Point<dim> (x_fC, y_fC));
       cnt_angles.push_back(theta);}

  }

}


//RVE lines
template <int dim>
void CNT<dim>::rve_linenodes(){

  Point<dim> pCNT, pointA, pointB, pointC, pointD, movingP;
  double theta, NX, NY;
  double Lsc= L_m/(std::pow(2.0,Nref));
  int alpha = int(1.2*(std::floor(L_f/(2*Lsc))+2));
  int beta  = int(1.2*(std::floor(D_f/(2*Lsc))+2));

  for (long unsigned int i=0; i<cnt_centres.size(); i++){

    pCNT  = cnt_centres[i];
    theta = cnt_angles[i];
    NX = std::cos(theta);
    NY = std::sin(theta);
    
    // I) CNT lines
    pointA(0) = pCNT(0)-(D_f/2.0*NY);
    pointA(1) = pCNT(1)+(D_f/2.0*NX);
    pointB(0) = pCNT(0)+(D_f/2.0*NY);
    pointB(1) = pCNT(1)-(D_f/2.0*NX);
    pointC(0) = pCNT(0)+(L_f/2.0*NX);
    pointC(1) = pCNT(1)+(L_f/2.0*NY);
    pointD(0) = pCNT(0)-(L_f/2.0*NX);
    pointD(1) = pCNT(1)-(L_f/2.0*NY);

    //Line A
    cntlines.push_back(pointA);
    for (int i=1; i<alpha+1; i++){
      float i_float=i*1.0;
      //one side
      double lx_dummy = pointA(0) + (L_f/2)*(i_float/alpha)*NX;
      double ly_dummy = pointA(1) + (L_f/2)*(i_float/alpha)*NY;
      if (lx_dummy<=L_m && lx_dummy>=0 && ly_dummy<=L_m && ly_dummy>=0)
        cntlines.push_back(Point<dim> (lx_dummy,ly_dummy));
      //other side
      double rx_dummy = pointA(0) - (L_f/2)*(i_float/alpha)*NX;
      double ry_dummy = pointA(1) - (L_f/2)*(i_float/alpha)*NY;
      if (rx_dummy<=L_m && rx_dummy>=0 && ry_dummy<=L_m && ry_dummy>=0)  
        cntlines.push_back(Point<dim> (rx_dummy,ry_dummy));
    }

    //Line B
    cntlines.push_back(pointB);
    for (int i=1; i<alpha+1; i++){
      float i_float=i*1.0;
      //one side
      double lx_dummy = pointB(0) + (L_f/2)*(i_float/alpha)*NX;
      double ly_dummy = pointB(1) + (L_f/2)*(i_float/alpha)*NY;
      if (lx_dummy<=L_m && lx_dummy>=0 && ly_dummy<=L_m && ly_dummy>=0)
        cntlines.push_back(Point<dim> (lx_dummy,ly_dummy));
      //other side
      double rx_dummy = pointB(0) - (L_f/2)*(i_float/alpha)*NX;
      double ry_dummy = pointB(1) - (L_f/2)*(i_float/alpha)*NY;
      if (rx_dummy<=L_m && rx_dummy>=0 && ry_dummy<=L_m && ry_dummy>=0)
        cntlines.push_back(Point<dim> (rx_dummy,ry_dummy));
    }

    //Line C
    cntlines.push_back(pointC);
    for (int i=1; i<beta+1; i++){
      float i_float=i*1.0;
      //one side
      double lx_dummy = pointC(0) - (D_f/2)*(i_float/beta)*NY;
      double ly_dummy = pointC(1) + (D_f/2)*(i_float/beta)*NX;
      if (lx_dummy<=L_m && lx_dummy>=0 && ly_dummy<=L_m && ly_dummy>=0)
        cntlines.push_back(Point<dim> (lx_dummy,ly_dummy));
      //other side
      double rx_dummy = pointC(0) + (D_f/2)*(i_float/beta)*NY;
      double ry_dummy = pointC(1) - (D_f/2)*(i_float/beta)*NX;
      if (rx_dummy<=L_m && rx_dummy>=0 && ry_dummy<=L_m && ry_dummy>=0)
        cntlines.push_back(Point<dim> (rx_dummy,ry_dummy));
    }

    //Line D
    cntlines.push_back(pointD);
    for (int i=1; i<beta+1; i++){
      float i_float=i*1.0;
      //one side
      double lx_dummy = pointD(0) - (D_f/2)*(i_float/beta)*NY;
      double ly_dummy = pointD(1) + (D_f/2)*(i_float/beta)*NX;
      if (lx_dummy<=L_m && lx_dummy>=0 && ly_dummy<=L_m && ly_dummy>=0)
        cntlines.push_back(Point<dim> (lx_dummy,ly_dummy));
      //other side
      double rx_dummy = pointD(0) + (D_f/2)*(i_float/beta)*NY;
      double ry_dummy = pointD(1) - (D_f/2)*(i_float/beta)*NX;
      if (rx_dummy<=L_m && rx_dummy>=0 && ry_dummy<=L_m && ry_dummy>=0)
        cntlines.push_back(Point<dim> (rx_dummy,ry_dummy));
    }


    // II) CNT solids
    
    for (int i=0; i<2*beta+1; i++){
      
      float i_float=i*1.0;
      movingP(0) = pointA(0) + (D_f/2)*(i_float/beta)*NY;
      movingP(1) = pointA(1) - (D_f/2)*(i_float/beta)*NX;

      cntsolids.push_back(movingP);
      for (int ii=1; ii<alpha+1; ii++){
        float ii_float=ii*1.0;
        //one side
        double lx_dummy = movingP(0) + (L_f/2)*(ii_float/alpha)*NX;
        double ly_dummy = movingP(1) + (L_f/2)*(ii_float/alpha)*NY;
        if (lx_dummy<=L_m && lx_dummy>=0 && ly_dummy<=L_m && ly_dummy>=0)
          cntsolids.push_back(Point<dim> (lx_dummy,ly_dummy));
        //other side
        double rx_dummy = movingP(0) - (L_f/2)*(ii_float/alpha)*NX;
        double ry_dummy = movingP(1) - (L_f/2)*(ii_float/alpha)*NY;
        if (rx_dummy<=L_m && rx_dummy>=0 && ry_dummy<=L_m && ry_dummy>=0)
          cntsolids.push_back(Point<dim> (rx_dummy,ry_dummy));
      }
    
    }

  }

}

// Mark the Boundaries
template <int dim>
void CNT<dim>::mark_boundaries(bool skipmb){
  
  const double x_max =L_m, y_max = L_m, x_min = 0.0, y_min = 0.0;
  Point<dim>  c_dl(x_min,y_min), c_dr(x_max,y_min), c_ul(x_min,y_max), c_ur(x_max,y_max); //Corner Points (l,r,u,d: left,right,up,down)

  for (const auto &cell : mesh.active_cell_iterators())
    if (cell->at_boundary()){

      if (skipmb==false)
        for (const auto &face : cell->face_iterators())
        if (face->at_boundary()) {
          if (std::fabs(face->center()(0) - (x_min)) < 1e-12)
            face->set_boundary_id(1);
          else if (std::fabs(face->center()(0) - (x_max)) < 1e-12)
            face->set_boundary_id(2);
          else if (std::fabs(face->center()(1) - (y_min)) < 1e-12)
            face->set_boundary_id(3);
          else if (std::fabs(face->center()(1) - (y_max)) < 1e-12)
            face->set_boundary_id(4);
        }

    }  
  
}


// Refinement
template <int dim>
void CNT<dim>::refine_mesh(){

  for (int step=0; step<Nref; step++){
    for (auto &cell : mesh.active_cell_iterators()){
      for (auto &cntline : cntlines){

          if (cell->at_boundary()){
            cell->set_refine_flag();
            break;}
          else if (cell->point_inside(cntline)){
            cell->set_refine_flag();
            break;
          }

      }
    }
  mesh.execute_coarsening_and_refinement();
  }
  
}


// Material IDs
template <int dim>
void CNT<dim>::assign_mid(){

  for (auto &cell : mesh.active_cell_iterators()){
    cell->set_material_id(1);
    for (auto &cntsolid : cntsolids){

      if (cell->point_inside(cntsolid)){
        cell->set_material_id(2);
        break;
      }

    }
  }
  
}


}

MADEAL_NAMESPACE_CLOSE