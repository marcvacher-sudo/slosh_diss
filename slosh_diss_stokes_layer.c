#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "navier-stokes/perfs.h"
#include "reduced.h"

double h;
double h_1;
double R;
double mu_b;
double rho_b;
double grav;

double omega_m;
scalar delta_m[];
double phi_0;
scalar phi_m[];
scalar U_0_m[];

#define reduce_var 1

FILE * fpmax; //

int main() {

  L0=1.; //size of the box
  rho2 = 1.;
  rho1=1000;
  mu1 = 0.25;
  mu2 = 0.01*mu1;
  h=L0/2; 
  grav=9.81;
  
  h_1=0.02*L0; //size of the perturbation

  TOLERANCE = 1e-3 [*];
  
  /**
  These boundary conditions are overriden by the event 'viscous_term'.
  */

  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

  u.n[top] = dirichlet(0.);
  u.t[top] = dirichlet(0.);

  u.n[left] = dirichlet(0.);
  u.n[right] = dirichlet(0.);
  u.t[left] = dirichlet(0.);
  u.t[right] = dirichlet(0.);

  N=128;
  origin (-L0/2, -h);
  init_grid(N);
  
  #if reduce_var 
    G.y=-grav;
  #endif

  fpmax =  fopen("log_wall_model.dat", "w");

  run();
}

/**
The initial state is a sloshing mode $m=3$, of amplitude $h_1$.
*/

event init (t = 0) {
  mu_b=mu2/mu1;
  rho_b=rho2/rho1;
  R=h/L0;

  double k_m=3.*pi/L0;

  fraction (f, h_1*sin(k_m*x)-y);
  omega_m=sqrt(grav*k_m*tanh(k_m*h));

  #if !reduce_var 
    const face vector G[] = {0,-grav};
    a=G;
  #endif

  event("properties");

  phi_0=grav*h_1/omega_m;

  foreach() {
    double mu_c=0.;
    double fm_c=0.;
    foreach_dimension() {
      mu_c+=mu.x[]+mu.x[-1];
      fm_c+=fm.x[]+fm.x[-1];
    }
    delta_m[]=sqrt(2.*mu_c/fm_c/(omega_m*rho[]));
    phi_m[]=phi_0*cosh(k_m*(y+h))/cosh(k_m*h)*cos(k_m*x);
    U_0_m[]=k_m*phi_m[];
    u.x[]=0.;//f[]*phi_0*cosh(k_m*(y+h))/cosh(k_m*h)*k_m*cos(k_m*x);
    u.y[]=0.;//f[]*phi_0*sinh(k_m*(y+h))/cosh(k_m*h)*k_m*sin(k_m*x);
  }

}

/** 
Function to return the distance to the closest wall.
*/

double distance_to_wall(double x,double y) {
  return min(L0/2.-x,min(x+L0/2.,y));
}

/**
Condition for the velocity on the solid. 
*/

double wall_condition_velocity(Point point,double x2, double y2,double t0) {
  if (delta_m[]==0.) {
    return 0.;
  }
  //fprintf(stderr,"%g %g %g %g %g %g %g %g\n",x,y,t,u.x[],u.y[],-U_0_m[]*(cos(omega_m*t0)+sin(omega_m*t0))/delta_m[],U_0_m[],delta_m[]);
  return -U_0_m[]*(cos(omega_m*t0)+sin(omega_m*t0))/delta_m[];
}

/** Wall model */

event viscous_term(i++) {

    if (i>10) {
    u.n[left]=dirichlet(0.);//neumann(wall_condition_velocity(x,y,t));
    u.t[left]=neumann(wall_condition_velocity(point,x,y,t));

    u.n[right]=dirichlet(0.);//neumann(wall_condition_velocity(x,y,t));
    u.t[right]=neumann(wall_condition_velocity(point,x,y,t));

    u.n[bottom]=dirichlet(0.);//neumann(wall_condition_velocity(x,y,t));
    u.t[bottom]=neumann(wall_condition_velocity(point,x,y,t));

    u.n[top]=dirichlet(0.);//neumann(wall_condition_velocity(x,y,t));
    u.t[top]=neumann(wall_condition_velocity(point,x,y,t));
  }
}

/**Computing total energy with potential energy having a shift to be of the same order as the kinetic energy */

event logfile (i++) { 

    double ke=0.;
  double gpe=0.;
  foreach() {
    double norm2=0.;
    foreach_dimension() {
      norm2+=sq(u.x[]);
    }
    ke+=norm2*f[]*dv();
    gpe+=y*f[]*dv();
  }
  double Ec=ke/2.;
  double Ep=grav*(gpe+L0*sq(h)/2.);

  static double E0=0.;
  if (i==10) {
    E0=Ec+Ep;
  }

  if (i>10) {
    fprintf (stderr, "%d %g %g\n", i, t,(Ec+Ep)/E0);
    fprintf (fpmax, "%d %g %g\n", i, t,(Ec+Ep)/E0);
  }
}

/**
We generate a video of the fraction to visualize the decreasing sloshing, as well as videos using "harmonic.h".
*/

event ppm_output (t = 0; t += 0.01; t <= 15) {

  char name[80];
  sprintf (name, "f.mp4");
  output_ppm (f, file = name, n = 512, min = 0, max = 1, linear = true);

  //u.y
  char name_u[80];
  sprintf (name_u, "u_y.mp4");
  output_ppm (u.y, file = name_u, n = 512,min = -0.3, max = 0.3, linear = true);
  
  //u.x
  char name_u_2[80];
  sprintf (name_u_2, "u_x.mp4");
  output_ppm (u.x, file = name_u_2, n = 512,min = -0.3, max = 0.3, linear = true);
}

event profile (t = end) {
  printf ("-----END-----\n");
}

