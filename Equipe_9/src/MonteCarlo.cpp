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
        sum += opt_->payoff(path);
        sumDesCarres += pow(opt_->payoff(path), 2);  
    }
    double ourPrice = (sum/nbSamples_)*exp(-(mod_->r_)*(opt_->T_));
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

    double timestep = (opt_->T_)/(opt_->nbTimeSteps_);
    PnlMat* path = pnl_mat_create (opt_->nbTimeSteps_ + 1, mod_->size_);
    PnlMat* shift_path_up =  pnl_mat_create (opt_->nbTimeSteps_ + 1, mod_->size_);
    PnlMat* shift_path_down =  pnl_mat_create (opt_->nbTimeSteps_ + 1, mod_->size_);

    for(int round = 0; round < nbSamples_; round++)
    {
        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
        for(int d = 0; d < mod_->size_; d++)
        {
            double shareDelta = pnl_vect_get(delta, d);
            double deviation = pnl_vect_get(std_dev, d);
            mod_->shiftAsset(shift_path_up, path, d, fdStep_, 0, timestep);
            mod_->shiftAsset(shift_path_down, path, d, -fdStep_, 0, timestep);
            double payoffUp = opt_->payoff(shift_path_up);
            double payoffDown = opt_->payoff(shift_path_down);
            shareDelta += payoffUp - payoffDown;
            deviation += pow(payoffUp - payoffDown, 2);
            pnl_vect_set(delta, d, shareDelta);
            pnl_vect_set(std_dev, d, deviation);
        }

    }
    PnlVect* ourSpots = mod_->spot_;
    for(int d = 0; d < mod_->size_; d++)
    {
        double initialValuePerAction = pnl_vect_get(ourSpots, d);
        double shareDelta = pnl_vect_get(delta, d);
        shareDelta = exp(-(mod_->r_) * (opt_->T_)) * shareDelta / (2 * nbSamples_ * fdStep_ * initialValuePerAction);
        pnl_vect_set(delta, d, shareDelta); //valeurs des deltas
        double deviation = pnl_vect_get(std_dev, d);
        deviation = exp(-2 * (mod_->r_) * (opt_->T_)) * deviation / (nbSamples_ * pow(2 * fdStep_ * initialValuePerAction, 2));
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
    double timestep = (opt_->T_)/(opt_->nbTimeSteps_);
    PnlMat* path = pnl_mat_create (opt_->nbTimeSteps_ + 1, mod_->size_);
    PnlMat* shift_path_up =  pnl_mat_create (opt_->nbTimeSteps_ + 1, mod_->size_);
    PnlMat* shift_path_down =  pnl_mat_create (opt_->nbTimeSteps_ + 1, mod_->size_);
    for(int round = 0; round < nbSamples_; round++)
    {
        mod_->asset(path, t, opt_->T_, opt_->nbTimeSteps_, rng_, past);
        for(int d = 0; d < mod_->size_; d++)
        {
            double shareDelta = pnl_vect_get(delta, d);
            double deviation = pnl_vect_get(std_dev, d);

            mod_->shiftAsset(shift_path_up, path, d, fdStep_, t, timestep);
            mod_->shiftAsset(shift_path_down, path, d, -fdStep_, t, timestep);
            double payoffUp = opt_->payoff(shift_path_up);
            double payoffDown = opt_->payoff(shift_path_down);
            shareDelta += payoffUp - payoffDown;
            deviation += pow(payoffUp - payoffDown, 2);
            pnl_vect_set(delta, d, shareDelta);
            pnl_vect_set(std_dev, d, deviation);
        }
    }
    PnlVect* ourShares = pnl_vect_create(mod_->size_);
    pnl_mat_get_row(ourShares, past, past->m - 1);
    for(int d = 0; d < mod_->size_; d++)
    {
        double saved = pnl_vect_get(ourShares, d);
        double shareDelta = pnl_vect_get(delta, d);
        shareDelta = exp(-(mod_->r_) * (opt_->T_ - t)) * shareDelta / (2 * nbSamples_ * fdStep_ * saved);
        pnl_vect_set(delta, d, shareDelta); //valeurs des deltas
        double deviation = pnl_vect_get(std_dev, d);
        deviation = exp(-2 * (mod_->r_) * (opt_->T_ - t)) * deviation / (nbSamples_ * pow(2 * fdStep_ * saved, 2));
        deviation -= pow(shareDelta, 2);
        deviation = sqrt(deviation / nbSamples_);
        pnl_vect_set(std_dev, d, deviation); // valeurs des ecart-type 
    }   
}




void MonteCarlo::makeReguralizedPast(PnlMat * past, PnlVect * shares, int i, double h)
{
    int reguralizedIndex = (i - 1) * opt_->nbTimeSteps_ / h;
    ((i - 1) * opt_->nbTimeSteps_ / h == reguralizedIndex) ? pnl_mat_add_row(past, past->m, shares) : pnl_mat_set_row(past, shares, past->m - 1);
}





void MonteCarlo::pAndL(PnlMat * marketData, double &premium, double &pnlError)
{
    PnlMat * past = pnl_mat_create(1, marketData->n);
    PnlVect * initialValue = pnl_vect_create(marketData->n);
    PnlMat * deltaMatrix = pnl_mat_create(marketData->m, marketData->n);
    PnlVect * spots = pnl_vect_create(marketData->n);
    PnlVect * initDelta = pnl_vect_create(marketData->n);
    PnlVect * init_delta_std_dev = pnl_vect_create(marketData->n);
    PnlVect * differenceDeltas = pnl_vect_create(marketData->n);
    PnlVect * currentShares = pnl_vect_create(marketData->n);
    pnl_mat_get_row(spots, marketData, 0);
    pnl_mat_set_row(past, spots, 0);

    MonteCarlo::delta(initDelta, init_delta_std_dev);
    cout << pnl_vect_get(initDelta, 0);
    pnl_mat_set_row(deltaMatrix, initDelta, 0);
    double value = premium - pnl_vect_scalar_prod(initDelta, spots);   
    for(int i = 1; i < marketData->m; i++)
    {
        PnlVect * currentDeltas = pnl_vect_create(marketData->n);
        PnlVect * delta_std_dev = pnl_vect_create(marketData->n);
        PnlVect * previousDeltas = pnl_vect_create(marketData->n);
        pnl_mat_get_row(currentShares, marketData, i);
        makeReguralizedPast(past, currentShares, i, marketData->m);
        MonteCarlo::delta(past, i  * opt_->T_ / marketData->m, currentDeltas, delta_std_dev);
        pnl_mat_get_row(previousDeltas, deltaMatrix, i - 1);
        pnl_mat_set_row(deltaMatrix, currentDeltas, i);
        pnl_vect_clone(differenceDeltas, currentDeltas);
        pnl_vect_minus_vect(differenceDeltas, previousDeltas);
        value = value * exp((mod_->r_ * opt_->T_)/marketData->m);
        value -= pnl_vect_scalar_prod(differenceDeltas, currentShares);
        //pnl_vect_free(&delta_std_dev);
    }
    PnlVect * lastDeltas = pnl_vect_create(marketData->n);
    pnl_mat_get_row(lastDeltas, deltaMatrix, marketData->m - 1);
    PnlVect * lastShares = pnl_vect_create(marketData->n);
    pnl_mat_get_row(lastShares, marketData, marketData->m - 1);
    pnlError = value + pnl_vect_scalar_prod(lastDeltas, lastShares) - opt_->payoff(past);

}






