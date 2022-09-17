#include <iostream>
#include <string>
#include "../BlackScholesModel.hpp"
#include "../AsianOption.hpp"
#include "../MonteCarlo.hpp"

using namespace std;

int main(int argc, char** argv)
{
    int size = 2;
    double strike = 100;
    PnlVect* spot = pnl_vect_create_from_scalar(size, 100);
    double T = 1.5;
    PnlVect* sigma = pnl_vect_create_from_scalar(size, 0.2);
    double r = 0.02;
    double rho = 0;
    PnlVect* payoffCoefficientsVector = pnl_vect_create_from_scalar(size, 0.5);
    int nbTimeSteps = 150;
    long nbSamples = 50000;

    double fdStep = 0.1;
    double prix;
    double std_dev;

    PnlRng* rng = pnl_rng_create(0);
    AsianOption * NewOption = new AsianOption(T, nbTimeSteps, size, strike, payoffCoefficientsVector);
    BlackScholesModel * NewModel = new BlackScholesModel(size, r, rho, sigma, spot);
    //cout << pnl_vect_get(NewModel->modifiedVolatility, 0);
    MonteCarlo * NewMonteCarlo = new MonteCarlo(NewModel, NewOption, rng, fdStep, nbSamples);
    NewMonteCarlo->price(prix, std_dev);
    cout << "Le prix est: " << prix << " et l'écart-type est: " << std_dev << endl;
    PnlVect* delta = pnl_vect_create(size);
    PnlVect* std_dev_delta = pnl_vect_create(size);
    NewMonteCarlo->delta(delta, std_dev_delta);
    cout << "Composition Delta : ";
    for(int i = 0; i < size; i++)
    {
        cout << pnl_vect_get(delta, i) << ", ";
    }
    cout << endl;
    cout << "Ecart-type Deltas : ";
    for(int i = 0; i < size; i++)
    {
        cout << pnl_vect_get(std_dev_delta, i) << ", ";
    }





}