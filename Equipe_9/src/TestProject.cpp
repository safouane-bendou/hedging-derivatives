#include <iostream>
#include <string>
#include "BasketOption.hpp"
#include "BlackScholesModel.hpp"
#include "MonteCarlo.hpp"

using namespace std;

int
main(int argc, char** argv)
{
    double T = 1;
    int nbTimeSteps = 1;
    int size=40;
    double r=0.04879;
    double rho=0.7;
    double fdStep=0;
    long nbSamples=50000;
    double strike = 100;
    double prix;
    double std_dev;

    PnlVect* sigma = pnl_vect_create_from_scalar(size, 0.2);
    PnlVect* spot = pnl_vect_create_from_scalar(size, 100);
    PnlVect* payoffCoefficientsVector = pnl_vect_create_from_scalar(size, 0.025);
    PnlRng* rng = pnl_rng_create(0);

    BasketOption * NewOption = new BasketOption(T, nbTimeSteps, size, strike, payoffCoefficientsVector);
    BlackScholesModel * NewModel = new BlackScholesModel(size, r, rho, sigma, spot);
    MonteCarlo NewMonteCarlo = new MonteCarlo(NewModel, NewOption, rng, fdStep, nbSamples);
    NewMonteCarlo.price(prix, std_dev);
    cout << "Le prix est: " << prix << "et l'écart-type est:" << std_dev;


}