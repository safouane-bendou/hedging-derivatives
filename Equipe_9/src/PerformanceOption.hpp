#pragma once
#include "Option.hpp"

class PerformanceOption : public Option
{
  public:
    PnlVect* payoffCoefficientsVector_;// vecteur contenant les coefficients du payoff de l'option

    // Constructeur de la classe Basket Option prenant quatre arguments
    PerformanceOption(double T, int nbTimeSteps, int size, PnlVect* payoffCoefficientsVector);

    // Fonction calculant le payoff de l'option 
    double payoff(const PnlMat* path);

    // Destructeur de la classe Basket Option
    ~PerformanceOption();
};