#include "MonteCarlo.hpp"
using namespace std;

MonteCarlo::MonteCarlo(BlackScholesModel* mod, Option* opt, PnlRng* rng, double fdStep, long nbSamples)
{
    mod_ = mod;
    opt_ = opt;
    rng_ = rng;
    fdStep_ = fdStep;
    nbSamples_ = nbSamples;
}


/**
 * Calcule le prix de l'option à la date 0
 *
 * @param[out] prix valeur de l'estimateur Monte Carlo
 * @param[out] ic écart type de l'estimateur
 */
void MonteCarlo::price(double& prix, double& std_dev){
    // Calcul du prix
    int n = (opt_->nbTimeSteps_)+1;
    int nombreActifs = (mod_->size_);
    PnlMat* path = pnl_mat_create (n, nombreActifs);   
    
    double sum = 0;
    double sumDesCarres = 0;
    for(int round = 0; round < nbSamples_; round++){
        mod_->asset(path, opt_->T_, n-1, rng_);
        //cout << MGET(path, 130, 5);
        sum += opt_->payoff(path);
        sumDesCarres += pow(opt_->payoff(path), 2);  
    }
    double ourPrice = (sum/(double)nbSamples_)*exp(-(mod_->r_)*(opt_->T_));
    prix = ourPrice; // Valeur de l'option 

    // Calcul de l'écart-type:1
    double s = pow(sum/nbSamples_, 2);
    double volatility = (sumDesCarres/nbSamples_ - s)*exp(-2*(mod_->r_)*(opt_->T_));
    std_dev = sqrt(volatility/nbSamples_); //valeur de l'écart-type
       
}

/**
 * Calcule le prix de l'option à la date t
 *
 * @param[in]  past contient la trajectoire du sous-jacent
 * jusqu'à l'instant t
 * @param[in] t date à laquelle le calcul est fait
 * @param[out] prix contient le prix
 * @param[out] std_dev contient l'écart type de l'estimateur
 */
void MonteCarlo::price(const PnlMat* past, double t, double& prix, double& std_dev)
{
    int n = (opt_->nbTimeSteps_)+1;
    int nombreActifs = (mod_->size_);
    PnlMat* path = pnl_mat_create (n, nombreActifs);  
    double sum = 0;
    double sumDesCarres = 0;
    for(int round = 0; round < nbSamples_; round++)
    {
        mod_->asset(path, t, opt_->T_, n-1, rng_, past);
        sum += opt_->payoff(path);
        sumDesCarres += pow(opt_->payoff(path), 2);  
    }
    double ourPrice = (sum/(double)nbSamples_)*exp(-(mod_->r_)*(opt_->T_-t));
    prix = ourPrice; // Valeur de l'option 

    // Calcul de l'écart-type
    double s = pow(sum/nbSamples_,2);
    double volatility = (sumDesCarres/nbSamples_ - s)*exp(-2*(mod_->r_)*(opt_->T_-t));
    std_dev = sqrt(volatility/nbSamples_); //valeur de l'écart-type

}




/**
 * Calcule le delta de l'option à la date 0
 *
 * @param[in] t date à laquelle le calcul est fait
 * @param[out] delta contient le vecteur de delta
 * @param[out] std_dev contient l'écart type de l'estimateur
 */
void MonteCarlo::delta(PnlVect* delta, PnlVect* std_dev)
{
    // Calcul du delta en t = 0
    PnlVect* ourSpots = mod_->spot_;
    int nombreActifs = (mod_->size_);
    double timestep = (opt_->T_)/(opt_->nbTimeSteps_);
    PnlMat* path = pnl_mat_create (opt_->nbTimeSteps_ + 1, nombreActifs);
    PnlMat* shift_path =  pnl_mat_create (opt_->nbTimeSteps_ + 1, nombreActifs);

    for(int round = 0; round < nbSamples_; round++)
    {
        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
        for(int d = 0; d < nombreActifs; d++)
        {
            mod_->shiftAsset(shift_path, path, d, fdStep_, 0, timestep);
            double payoffUp = opt_->payoff(shift_path);
            mod_->shiftAsset(shift_path, path, d, -fdStep_, 0, timestep);
            double payoffDown = opt_->payoff(shift_path);
            double shareDelta = pnl_vect_get(delta, d);
            shareDelta += payoffUp - payoffDown;
            pnl_vect_set(delta, d, shareDelta);
            double deviation = pnl_vect_get(std_dev, d);
            deviation += pow(payoffUp - payoffDown, 2);
            pnl_vect_set(std_dev, d, deviation);
        }

    }

    for(int d = 0; d < nombreActifs; d++)
    {
        double initialValuePerAction = pnl_vect_get(ourSpots, d);
        double shareDelta = pnl_vect_get(delta, d);
        shareDelta = exp(-(mod_->r_) * (opt_->T_)) * shareDelta / (2 * nbSamples_ * fdStep_ * initialValuePerAction);
        pnl_vect_set(delta, d, shareDelta); //valeurs des deltas
        double deviation = pnl_vect_get(std_dev, d);
        deviation = exp(-2 * (mod_->r_) * (opt_->T_)) * deviation / (2 * nbSamples_ * fdStep_ * initialValuePerAction);
        deviation -= pow(shareDelta, 2);
        deviation = sqrt(deviation / nbSamples_);
        pnl_vect_set(std_dev, d, deviation); // valeurs des ecart-type 
    }   
}



/**
 * Calcule le delta de l'option à la date t
 *
 * @param[in] past contient la trajectoire du sous-jacent
 * jusqu'à l'instant t
 * @param[in] t date à laquelle le calcul est fait
 * @param[out] delta contient le vecteur de delta
 * @param[out] std_dev contient l'écart type de l'estimateur
 */
void MonteCarlo::delta(const PnlMat* past, double t, PnlVect* delta, PnlVect* std_dev)
{
    // calcul du delta en t
    PnlVect* ourSpots = mod_->spot_;
    int nombreActifs = (mod_->size_);
    double timestep = (opt_->T_)/(opt_->nbTimeSteps_);
    PnlMat* path = pnl_mat_create (opt_->nbTimeSteps_ + 1, nombreActifs);
    PnlMat* shift_path =  pnl_mat_create (opt_->nbTimeSteps_ + 1, nombreActifs);

    for(int round = 0; round < nbSamples_; round++)
    {
        mod_->asset(path, opt_->T_, t, opt_->nbTimeSteps_, rng_, past);
        for(int d = 0; d < nombreActifs; d++)
        {
            mod_->shiftAsset(shift_path, path, d, fdStep_, 0, timestep);
            double payoffUp = opt_->payoff(shift_path);
            mod_->shiftAsset(shift_path, path, d, -fdStep_, 0, timestep);
            double payoffDown = opt_->payoff(shift_path);
            double shareDelta = pnl_vect_get(delta, d);
            shareDelta += payoffUp - payoffDown;
            pnl_vect_set(delta, d, shareDelta);
            double deviation = pnl_vect_get(std_dev, d);
            deviation += pow(payoffUp - payoffDown, 2);
            pnl_vect_set(std_dev, d, deviation);
        }

    }

    for(int d = 0; d < nombreActifs; d++)
    {
        double initialValuePerAction = pnl_vect_get(ourSpots, d);
        double shareDelta = pnl_vect_get(delta, d);
        shareDelta = exp(-(mod_->r_) * (opt_->T_ - t)) * shareDelta / (2 * nbSamples_ * fdStep_ * initialValuePerAction);
        pnl_vect_set(delta, d, shareDelta); //valeurs des deltas
        double deviation = pnl_vect_get(std_dev, d);
        deviation = exp(-2 * (mod_->r_) * (opt_->T_ - t)) * deviation / (2 * nbSamples_ * fdStep_ * initialValuePerAction);
        deviation -= pow(shareDelta, 2);
        deviation = sqrt(deviation / nbSamples_);
        pnl_vect_set(std_dev, d, deviation); // valeurs des ecart-type 
    }   
}











