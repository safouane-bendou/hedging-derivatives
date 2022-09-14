#include "BasketOption.hpp"

BasketOption::BasketOption(double T, int nbTimeSteps, int size, double strike, PnlVect * payoffCoefficientsVector) {
    T_ = T;
    nbTimeSteps_ = nbTimeSteps;
    size_ = size;
    strike_ = strike;
    payoffCoefficientsVector_ = payoffCoefficientsVector;
}

double BasketOption::payoff(const PnlMat* path)
{
    
    PnlVect* ActionPricesInT = pnl_vect_create(size_);
    
    pnl_mat_get_row(ActionPricesInT, path, nbTimeSteps_);
    double sum = pnl_vect_scalar_prod(payoffCoefficientsVector_, ActionPricesInT);
    
    double Payoff = sum - strike_;
    if (Payoff >= 0)
    {
        return Payoff;
    } 
    else
    {
        return 0;    
    }   
}
//BasketOption::~BasketOption(){};