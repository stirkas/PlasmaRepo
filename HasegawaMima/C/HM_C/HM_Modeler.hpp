#pragma once

#include "HM_Definitions.hpp"

namespace HM
{

class Modeler
{
public:
   Modeler(const SimParams& simParams);
   ~Modeler();

private:
   void CreateGrids();
   SimParams simParams_;
   int gridSize_ = 0;
   double* grid_         = nullptr;
   double* spectralGrid_ = nullptr;
   double* aliasedGrid_  = nullptr; //Also spectral, for resolving nonlinearity issues.
};

} //Close namespace.