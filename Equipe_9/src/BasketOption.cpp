#include "BasketOption.hpp"

using namespace std;

double
BasketOption::payoff(const PnlMat* path)
{
    
    PnlVect* ActionPricesInT = pnl_vect_create(size_);
    
    pnl_mat_get_row(ActionPricesInT, path, nbTimeSteps_);
    double sum = 0;
    for(int i = 0; i < size_; i++) {
        sum += pnl_vect_get(ActionPricesInT, i) * pnl_vect_get(payoffCoefficientsVector, i);
    }
    double interimPayoff = (sum - strike);
    if (interimPayoff >= 0) {
        return interimPayoff;
    } 
    else {
        return 0;    
    }   
}