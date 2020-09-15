#pragma once

namespace HM
{

struct SimParams
{
   //General params.
   int    nt = 0; //Num of timesteps.
   int    nx = 0; //Num of x-steps.
   int    ny = 0; //Num of y-steps.
   double lx = 0; //Box size in x-dim.
   double ly = 0; //Box size in y-dim.
   double dt = 0; //Size of timestep.
   
   //Params specific to model (Haotian's)
   double tau      = 0; //(T_e/T_i)
   double eta      = 0; //(r_n/r_t) - see below for defs.
   double mRat     = 0; //(m_e/m_i)
   double rnByRhoI = 0; //(r_n/rho_i), rho == gyroradius.
   //(grad_x(n_e)/n_e) = -1/r_n
   //(grad_x(T_e)/T_e) = -1/r_t

   bool loadGENE = false; //Possibility for loading data from GENE ASCII output.
};

} //Close namespace.