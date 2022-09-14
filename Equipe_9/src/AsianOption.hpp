#pragma once
#include "Option.hpp"

class AsianOption : public Option
{
  public:
    double strike_;
    PnlVect* payoffCoefficientsVector_;

    AsianOption(double T, int nbTimeSteps, int size, double strike, PnlVect* payoffCoefficientsVector);
    double payoff(const PnlMat* path);
};