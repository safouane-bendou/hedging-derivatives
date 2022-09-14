#include"AsianOption.hpp"

AsianOption::AsianOption(double T, int nbTimeSteps, int size, double strike, PnlVect* payoffCoefficientsVector)
{
    T_ = T;
    nbTimeSteps_ = nbTimeSteps;
    size_ = size;
    strike_ = strike;
    payoffCoefficientsVector_ = payoffCoefficientsVector;
}
double AsianOption::payoff(const PnlMat* path)
{
    int n = nbTimeSteps_ + 1;
    double sum = 0;
    for (int i = 0; i < size_; i++) {

        PnlVect* ActionPrices = pnl_vect_create(n);
        pnl_mat_get_col(ActionPrices, path, i);
        double SumofPrices = pnl_vect_sum(ActionPrices);
        sum += SumofPrices* pnl_vect_get(payoffCoefficientsVector_, i);

    }
    double ourPayoff = (sum / n) - strike_;
    if (ourPayoff >= 0) {
        return ourPayoff;
    } else {
        return 0;
    }

}
//AsianOption::~AsianOption(){};