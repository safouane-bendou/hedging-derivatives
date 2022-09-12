#include"AsianOption.hpp"

using namespace std;

double AsianOption::payoff(const PnlMat* path){
    int n = nbTimeSteps_ + 1;
    double sum = 0;
    for (int i = 0; i < size_; i++) {

        PnlVect* ActionPrices = pnl_vect_create(n);
        pnl_mat_get_row(ActionPrices, path, i);
        double SumofPrices = pnl_vect_sum(ActionPrices);
        sum += SumofPrices* pnl_vect_get(payoffCoefficientsVector, i);

    }
    double ourPayoff = (sum / n) - strike;
    if (ourPayoff >= 0) {
        return ourPayoff;
    } else {
        return 0;
    }

}