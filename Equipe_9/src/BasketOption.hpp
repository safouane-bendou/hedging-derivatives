#pragma once
#include "Option.hpp"

class BasketOption : public Option
{
  public:
    double strike_; // Prix d'exercice de l'option
    PnlVect* payoffCoefficientsVector_; // vecteur contenant les coefficients du payoff de l'option


    // Constructeur de la classe Basket Option prenant 5 arguments
    BasketOption(double T, int nbTimeSteps, int size, double strike, PnlVect * payoffCoefficientsVector);


    // Fonction calculant le payoff de l'option 
    double payoff(const PnlMat* path);
    

  // Destructeur de la classe Basket Option
  ~BasketOption();
};