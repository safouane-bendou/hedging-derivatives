#include <iostream>
#include <string>
#include "../BlackScholesModel.hpp"
#include "../PerformanceOption.hpp"
#include "../AsianOption.hpp"
#include "../BasketOption.hpp"
#include "../MonteCarlo.hpp"
#include "jlparser/parser.hpp"
#include "../HedgingResults.hpp"

using namespace std;

int main(int argc, char** argv)
{
    double T, r, strike, rho;
    PnlVect *spot, *sigma, *payoffCoefficientsVector;
    string type;
    int size, nbTimeSteps, hedging ;
    size_t nbSamples;

    char* infile = argv[2];
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
    P->extract("hedging dates number", hedging);
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
    double erreur_couverture;
    PnlMat * marketData= pnl_mat_create_from_file(argv[1]);
    NewMonteCarlo->pAndL(marketData, prix, erreur_couverture);
    HedgingResults res(prix, prix_std_dev, erreur_couverture);
    std::cout << res << std::endl;

}