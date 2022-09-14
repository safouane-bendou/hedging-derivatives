#include <iostream>
#include <string>
#include "BlackScholesModel.hpp"
#include "AsianOption.hpp"
#include "BasketOption.hpp"
#include "MonteCarlo.hpp"

using namespace std;

int
main(int argc, char** argv)
{
    double T = 1.5;
    int nbTimeSteps = 150;
    int size = 2;
    double r = 0.02;
    double rho = 0;
    double fdStep = 0;
    long nbSamples = 50000;
    double strike = 100;
    double price;
    double std_dev;

    PnlVect* sigma = pnl_vect_create_from_scalar(size, 0.2);
    PnlVect* spot = pnl_vect_create_from_scalar(size, 100);
    PnlVect* payoffCoefficientsVector = pnl_vect_create_from_scalar(size, 0.025);
    PnlRng* rng = pnl_rng_create(0);

    AsianOption * NewOption = new AsianOption(T, nbTimeSteps, size, strike, payoffCoefficientsVector);
    BlackScholesModel * NewModel = new BlackScholesModel(size, r, rho, sigma, spot);
    MonteCarlo * NewMonteCarlo = new MonteCarlo(NewModel, NewOption, rng, fdStep, nbSamples);
    NewMonteCarlo->price(price, std_dev);
    cout << "Le prix est: " << price << " et l'écart-type est: " << std_dev << "    ";


}