#include <iostream>
#include <string>
#include "../BlackScholesModel.hpp"
#include "../BasketOption.hpp"
#include "../PerformanceOption.hpp"
#include "../AsianOption.hpp"
#include "../MonteCarlo.hpp"
#include "jlparser/parser.hpp"
#include "../PricingResults.hpp"

using namespace std;

int main(int argc, char** argv)
{
    double T, r, strike, rho;
    PnlVect *spot, *sigma, *payoffCoefficientsVector;
    string type;
    int size, nbTimeSteps ;
    size_t nbSamples;

    char* infile = argv[1];
    Param* P = new Parser(infile);

    P->extract("option type", type);
    P->extract("maturity", T);
    P->extract("option size", size);
    P->extract("spot", spot, size);
    P->extract("volatility", sigma, size);
    P->extract("payoff coefficients", payoffCoefficientsVector, size);
    P->extract("interest rate", r);
    P->extract("correlation", rho);
    P->extract("sample number", nbSamples);
    P->extract("timestep number", nbTimeSteps);
     double fdStep = 0.1;
     double prix;
     double prix_std_dev;
    Option * NewOption;

    if(type == "basket"){

      P->extract("strike", strike);
      NewOption = new BasketOption(T, nbTimeSteps, size, strike, payoffCoefficientsVector);
    
    }
    else if (type == "asian"){

      P->extract("strike", strike);
      NewOption = new AsianOption(T, nbTimeSteps, size, strike, payoffCoefficientsVector);
    }
    else if (type == "performance"){
      NewOption = new PerformanceOption(T, nbTimeSteps, size, payoffCoefficientsVector);
    }

    PnlRng* rng = pnl_rng_create(0);
    BlackScholesModel * NewModel = new BlackScholesModel(size, r, rho, sigma, spot);
    MonteCarlo * NewMonteCarlo = new MonteCarlo(NewModel, NewOption, rng, fdStep, nbSamples);
    NewMonteCarlo->price(prix, prix_std_dev);
    
    //Deltas:
    PnlVect* delta = pnl_vect_create(size);
    PnlVect* delta_std_dev = pnl_vect_create(size);
    NewMonteCarlo->delta(delta, delta_std_dev);
    PricingResults res(prix, prix_std_dev, delta, delta_std_dev);
    std::cout << res << std::endl;

}