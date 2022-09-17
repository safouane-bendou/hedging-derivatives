#include <iostream>
#include <string>
#include "../BlackScholesModel.hpp"
#include "../BasketOption.hpp"
#include "../MonteCarlo.hpp"

using namespace std;

int main(int argc, char** argv)
{
    int size = 40;
    double strike = 100;
    PnlVect* spot = pnl_vect_create_from_scalar(size, 100);
    double T = 3;
    PnlVect* sigma = pnl_vect_create_from_scalar(size, 0.2);
    double r = 0.04879;
    double rho = 0;
    PnlVect* payoffCoefficientsVector = pnl_vect_create_from_scalar(size, 0.025);
    int nbTimeSteps = 1;
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