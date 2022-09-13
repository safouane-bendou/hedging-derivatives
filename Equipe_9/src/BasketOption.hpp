#pragma once
#include "Option.hpp"

class BasketOption : public Option
{
  public:
    double strike;
    PnlVect* payoffCoefficientsVector;
    double payoff(const PnlMat* path);
};