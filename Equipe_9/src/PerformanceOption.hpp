#pragma once
#include "Option.hpp"

class PerformanceOption : public Option
{
  public:
     PnlVect* payoffCoefficientsVector;
    double payoff(const PnlMat* path);
};