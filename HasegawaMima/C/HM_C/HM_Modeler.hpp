#pragma once

#include "HM_Definitions.hpp"

namespace HM
{

class Modeler
{
public:
   Modeler(const SimParams& simParams);

private:
   SimParams simParams_;
};

} //Close namespace.