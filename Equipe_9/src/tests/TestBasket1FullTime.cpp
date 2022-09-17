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
    PnlVect* spot = pnl_vect_create_from_scalar(size, 100.0);
    double T = 3;
    double t = 1; //une date inf à T 
    PnlVect* sigma = pnl_vect_create_from_scalar(size, 0.2);
    double r = 0.04879;
    double rho = 0;
    PnlVect* payoffCoefficientsVector = pnl_vect_create_from_scalar(size, 0.025);
    int nbTimeSteps = 1;
    long nbSamples = 50000;
    double fdStep = 0.1;
    double prix1, prix2, std_dev1, std_dev2;
    PnlMat * past = pnl_mat_create_from_scalar(1, size, 100.0); // the past matrix has only one row, this row is exactly the spot values since there is no discretisation

    PnlRng* rng = pnl_rng_create(0);
    BasketOption* newOption1 = new BasketOption(T, nbTimeSteps, size, strike, payoffCoefficientsVector);
    BasketOption* newOption2 = new BasketOption(T - t, nbTimeSteps, size, strike, payoffCoefficientsVector);
    BlackScholesModel* NewModel = new BlackScholesModel(size, r, rho, sigma, spot);
    MonteCarlo* NewMonteCarlo1 = new MonteCarlo(NewModel, newOption1, rng, fdStep, nbSamples);
    MonteCarlo* NewMonteCarlo2 = new MonteCarlo(NewModel, newOption2, rng, fdStep, nbSamples);
    NewMonteCarlo1->price(past, t, prix1, std_dev1);
    NewMonteCarlo2->price(prix2, std_dev2);
    cout << "Le prix attendu est: " << prix2 << " et l'écart-type attendu est: " << std_dev2 << endl ;
    cout << "Le prix en t est: " << prix1 << " et l'écart-type est: " << std_dev1 << endl;
    



    //Deltas :
    PnlVect* delta2 = pnl_vect_create(size);
    PnlVect* std_dev_delta2 = pnl_vect_create(size);
    NewMonteCarlo2->delta(delta2, std_dev_delta2);
    cout << "Composition Delta attendue: ";
    for(int i = 0; i < size; i++)
    {
        cout << pnl_vect_get(delta2, i) << ", ";
    }
    cout << endl;
    cout << "Ecart-type Deltas attendus: ";
    for(int i = 0; i < size; i++)
    {
        cout << pnl_vect_get(std_dev_delta2, i) << ", ";
    }
    cout << endl << endl;

    PnlVect* delta1 = pnl_vect_create(size);
    PnlVect* std_dev_delta1 = pnl_vect_create(size);
    NewMonteCarlo1->delta(past, t, delta1, std_dev_delta1);
    cout << "Composition Delta en t: ";
    for(int i = 0; i < size; i++)
    {
        cout << pnl_vect_get(delta1, i) << ", ";
    }
    cout << endl;
    cout << "Ecart-type Deltas en t ";
    for(int i = 0; i < size; i++)
    {
        cout << pnl_vect_get(std_dev_delta1, i) << ", ";
    }





}