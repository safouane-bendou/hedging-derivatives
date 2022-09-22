#pragma once
#include "Option.hpp"

class AsianOption : public Option
{
  public:
    double strike_; // Prix d'exercice de l'option
    PnlVect* payoffCoefficientsVector_; // vecteur contenant les coefficients du payoff de l'option


    // Constructeur de la classe Asian Option prenant 5 arguments
    AsianOption(double T, int nbTimeSteps, int size, double strike, PnlVect* payoffCoefficientsVector);

    // Fonction calculant le payoff de l'option
    double payoff(const PnlMat* path);

     // Destructeur de la classe Asian Option
    ~AsianOption();
};