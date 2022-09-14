#include "PerformanceOption.hpp"

using namespace std;



PerformanceOption::PerformanceOption(double T, int nbTimeSteps, int size, PnlVect* payoffCoefficientsVector)
{
    T_ = T;
    nbTimeSteps_ = nbTimeSteps;
    size_ = size;
    payoffCoefficientsVector_ = payoffCoefficientsVector;
}

double PerformanceOption::payoff(const PnlMat* path)
{
    
  
    double overallSum= 1;
    double sum = 0;

    for(int i=0; i<nbTimeSteps_; i++) { 

        PnlVect* ActionPricesDateII = pnl_vect_create(size_);
        PnlVect* ActionPricesDateI = pnl_vect_create(size_);
        pnl_mat_get_row(ActionPricesDateII, path, i+1);
        pnl_mat_get_row(ActionPricesDateI, path, i);

        sum = (pnl_vect_scalar_prod(payoffCoefficientsVector, ActionPricesDateII)) / pnl_vect_scalar_prod(payoffCoefficientsVector, ActionPricesDateI) - 1;
        if (sum > 0) {
            overallSum += sum;
        }

    }
 
    return overallSum;
    
}

PerformanceOption::~PerformanceOption(){};