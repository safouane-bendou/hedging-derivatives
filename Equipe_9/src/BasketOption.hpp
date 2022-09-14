#pragma once
#include "Option.hpp"

class BasketOption : public Option
{
  public:
    double strike_;
    PnlVect* payoffCoefficientsVector_;

    BasketOption(double T, int nbTimeSteps, int size, double strike, PnlVect * payoffCoefficientsVector);
    double payoff(const PnlMat* path);
};