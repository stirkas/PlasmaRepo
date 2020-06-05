#pragma once

#include "HM_Definitions.hpp"

namespace HM
{

class Initializer
{
public:
   Initializer(const SimParams& simParams);

private:
   SimParams simParams_;
   int caseNum_;
};

} //Close namespace.
