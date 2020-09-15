#include <iostream>
#include <blas/cblas.h>
#include <lapacke.h>
#include "HM_Definitions.hpp"
#include "HM_Modeler.hpp"

int main()
{
   HM::SimParams params {.nt = 1, .nx = 4, .ny = 4, .lx = 4, .ly = 4};
   HM::Modeler modeler(params);

   return 0;
}