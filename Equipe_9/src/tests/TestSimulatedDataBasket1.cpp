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
    PnlVect* trend = pnl_vect_create_from_scalar(size, 0.04);
    int H = 1000;
    double rho = 0;
    PnlVect* payoffCoefficientsVector = pnl_vect_create_from_scalar(size, 0.025);
    int nbTimeSteps = 1;
    long nbSamples = 50000;
    double fdStep = 0.1;
    double prix;
    double std_dev;
    PnlMat* simulatedData = pnl_mat_create(H+1, size);  

    PnlRng* rng = pnl_rng_create(0);
    BasketOption * NewOption = new BasketOption(T, nbTimeSteps, size, strike, payoffCoefficientsVector);
    BlackScholesModel * NewModel = new BlackScholesModel(size, r, rho, sigma, spot, trend);
    NewModel->simul_market(simulatedData, H, T, rng);

    PnlVect* testVector = pnl_vect_create(size);
    for(int i = 0; i < H + 1; i++){  
        pnl_mat_get_row (testVector, simulatedData, i);
        cout << i << " : " ; 
            for(int d = 0; d < size; d++)
        {
            cout << pnl_vect_get(testVector, d) << ", ";     
        }
        cout << "\n";
        cout << "\n";
    }
    
}
    
    