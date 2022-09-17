#pragma once
#include "Option.hpp"

class PerformanceOption : public Option
{
  public:
    PnlVect* payoffCoefficientsVector_;


    PerformanceOption(double T, int nbTimeSteps, int size, PnlVect* payoffCoefficientsVector);
    double payoff(const PnlMat* path);
};