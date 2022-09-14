#pragma once
#include "Option.hpp"

class BasketOption : public Option
{
  public:
    double strike_;
    PnlVect* payoffCoefficientsVector_;
    double payoff(const PnlMat* path);
};