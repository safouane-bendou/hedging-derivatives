#include"AsianOption.hpp"
using namespace std;

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
    double sum = 0;
    for (int d = 0; d < size_; d++)
    {
        PnlVect* sharePrices = pnl_vect_create(nbTimeSteps_ + 1);
        pnl_mat_get_col(sharePrices, path, d);
        double sumofPrices = pnl_vect_sum(sharePrices);
        sum += sumofPrices * pnl_vect_get(payoffCoefficientsVector_, d);
    }
    double ourPayoff = (sum / (double) (nbTimeSteps_ + 1)) - strike_;
    if (ourPayoff > 0)
    {
        return ourPayoff;
    } 
    else
    {
        return 0;
    }

}
//AsianOption::~AsianOption(){};