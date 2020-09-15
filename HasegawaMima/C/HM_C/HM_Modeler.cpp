#include "HM_Modeler.hpp"

#include <iostream>
#include <math.h>
namespace HM
{

Modeler::Modeler(const SimParams& simParams)
   : simParams_(simParams)
{
   CreateGrids();
}

Modeler::~Modeler()
{
   if (grid_)         {delete grid_;}
   if (spectralGrid_) {delete spectralGrid_;}
   if (aliasedGrid_)  {delete aliasedGrid_;}
}

void Modeler::CreateGrids()
{
   gridSize_ = simParams_.nx * simParams_.ny;
   if (gridSize_ > 0)
   {
      grid_         = new double [gridSize_];
      spectralGrid_ = new double [gridSize_];
      aliasedGrid_  = new double [gridSize_];

      double x[simParams_.nx]   = {0};
      double y[simParams_.ny]   = {0};
      double kx[simParams_.nx]  = {0};
      double ky[simParams_.ny]  = {0};
      double kxd[simParams_.nx] = {0};
      double kyd[simParams_.ny] = {0};
      double dx = simParams_.lx/simParams_.nx;
      double dy = simParams_.ly/simParams_.ny;

      //Load x vec's.
      for (int i = 0; i < simParams_.nx; ++i)
      {
         x[i]  = i*dx;
      }
      
      //Load y vec's.
      for (int i = 0; i < simParams_.ny; ++i)
      {
         y[i] = i*dy;
      }
   }
   else
   {
      std::cout << "Grid size must be greather than 0!!!" << std::endl;
   }
}

} //Close namespace