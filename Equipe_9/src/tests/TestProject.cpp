#include <iostream>
#include <string>
#include "../BlackScholesModel.hpp"
#include "../AsianOption.hpp"
#include "../BasketOption.hpp"
#include "../PerformanceOption.hpp"
#include "../MonteCarlo.hpp"

using namespace std;

int main(int argc, char** argv)
{
    int size = 5;
    double strike = 100;
    PnlVect* spot = pnl_vect_create_from_scalar(size, 100);
    double T = 2;
    PnlVect* sigma = pnl_vect_create_from_scalar(size, 0.2);
    double r = 0.03;
    double rho = 0.5;
    PnlVect* payoffCoefficientsVector = pnl_vect_create_from_scalar(size, 0.2);
    int nbTimeSteps = 12;
    long nbSamples = 50000;

    double fdStep = 0;
    double prix;
    double std_dev;

    PnlRng* rng = pnl_rng_create(0);
    BasketOption * NewOption = new BasketOption(T, nbTimeSteps, size, strike, payoffCoefficientsVector);
    //PerformanceOption * NewOption = new PerformanceOption(T, nbTimeSteps, size, payoffCoefficientsVector);
    BlackScholesModel * NewModel = new BlackScholesModel(size, r, rho, sigma, spot);
    //cout << pnl_vect_get(NewModel->modifiedVolatility, 0);
    MonteCarlo * NewMonteCarlo = new MonteCarlo(NewModel, NewOption, rng, fdStep, nbSamples);
    NewMonteCarlo->price(prix, std_dev);
    cout << "Le prix est: " << prix << " et l'écart-type est: " << std_dev << "    ";


}