#include "HM_Definitions.hpp"
#include "HM_Modeler.hpp"

int main()
{
   HM::SimParams params {.nt = 4};
   HM::Modeler modeler(params);

   return 0;
}