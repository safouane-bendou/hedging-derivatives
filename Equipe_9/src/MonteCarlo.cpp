#include"MonteCarlo.hpp"
#include <math.h>
using namespace std;

void MonteCarlo::price(double& prix, double& std_dev){
    // Calcul du prix
    int n = (opt_->nbTimeSteps_)+1;
    int nombreActifs = (mod_->size_);
     PnlMat* path = pnl_mat_create (n, nombreActifs);   
    PnlVect* ActionPrices = pnl_vect_create(n);
    double sum = 0;
    double sumDesCarres =0
    for(int i =0; i < nbSamples_; i++){
        mod_->asset( path, (opt_->T_), (n-1), (rng_->rng));
        sum += opt_->payoff(path);
        sumDesCarres += pow(opt_->payoff(path),2);
    }
        double ourPrice = (sum/nbSamples_)*exp(-(mod_->r_)*(opt_->T_));
        prix = ourPrice; // Valeur de l'option 

        // Calcul de l'écart-type:
        double s = pow((sum/nbSamples_),2);
        double volatility = (sumDesCarres/nbSamples_ - s)*exp(-2*(mod_->r_)*(opt_->T_));
        std_dev = pow(volatility,0.5); //valeur de l'écart-type
}


