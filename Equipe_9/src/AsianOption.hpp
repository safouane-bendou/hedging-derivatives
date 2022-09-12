#pragma once
#include "Option.hpp"

class AsianOption : public Option
{
  public:
    double strike;
    PnlVect* payoffCoefficientsVector;

    double payoff(const PnlMat* path);
};