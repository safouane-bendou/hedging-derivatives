#include "PerformanceOption.hpp"

using namespace std;

double
BasketOption::payoff(const PnlMat* path)
{
    
  
    double overallSum= 0;

    for(int i=0; i<nbTimesSteps_; i++) { // or nb+1

        PnlVect* ActionPricesInT = pnl_vect_create(size_);
    
        pnl_mat_get_col(ActionPricesInT, i);

        double firstSum = 0;
        double secondSum= 0;
        double result=0;
        for(int d = 0; d < size_+; d++) {
            firstSum += pnl_vect_get(ActionPricesInT, d) * pnl_vect_get(payoffCoefficientsVector, d);
        }

        for(int d =0; i<size_+; i++) {
            secondSum += pnl_vect_get(ActionPricesInt, d-1) * pnl_vect_get(payoffCoefficientsVector, d);
        }
        
        result=firstSum/ secondSum -1;
        if (result < 0) {
            result= 0;
        }

        overallsum += result;

    }
    double interimPayoff = overallSum + 1;
    return interimPayoff;
    
}