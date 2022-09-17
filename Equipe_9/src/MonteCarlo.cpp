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
    double s = pow(sum/nbSamples_,2);
    double volatility = (sumDesCarres/nbSamples_ - s)*exp(-2*(mod_->r_)*(opt_->T_));
    std_dev = sqrt(volatility/nbSamples_); //valeur de l'écart-type
       
}

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

void MonteCarlo::delta(PnlVect* delta, PnlVect* std_dev)
{
    // Calcul du delta en t = 0
    PnlVect* ourSpots = mod_->spot_;
    int n = (opt_->nbTimeSteps_)+1;
    int nombreActifs = (mod_->size_);
    double timestep = (opt_->T_)/(n-1);
    PnlMat* path = pnl_mat_create (n, nombreActifs);
    PnlMat* shift_path =  pnl_mat_create (n, nombreActifs);
    PnlMat* shift_path2 =  pnl_mat_create (n, nombreActifs);
    double sum = 0;
    double sumDesCarres = 0;
    for(int roundOnAction = 0; roundOnAction < nombreActifs; roundOnAction++)
    {
        double initialValuePerAction = pnl_vect_get(ourSpots, roundOnAction);
        for(int round = 0; round < nbSamples_; round++)
        {
            mod_->asset(path, opt_->T_, n-1, rng_);
            mod_->shiftAsset(shift_path, path, roundOnAction, fdStep_, 0, timestep);
            mod_->shiftAsset(shift_path2, path, roundOnAction, -fdStep_, 0, timestep);
            sum += opt_->payoff(shift_path)- opt_->payoff(shift_path2);
            sumDesCarres+= pow(opt_->payoff(shift_path)- opt_->payoff(shift_path2), 2);
        } 
        double deltaPerAction = (sum*exp(-2*(mod_->r_)*(opt_->T_)))/(2*nbSamples_*fdStep_*initialValuePerAction);
        pnl_vect_set(delta, roundOnAction, deltaPerAction); //valeurs des deltas
        double s = pow((sum/(2*nbSamples_*initialValuePerAction*fdStep_)),2);
        double volatilityPerAction = exp(-2*(mod_->r_)*(opt_->T_))*(sumDesCarres/(2*nbSamples_*initialValuePerAction*fdStep_)-s);
        double std_devElement = sqrt(volatilityPerAction/nbSamples_); 
        pnl_vect_set(std_dev, roundOnAction, std_devElement); // valeurs des ecart-type 
    }   
}

void MonteCarlo::delta(const PnlMat* past, double t, PnlVect* delta, PnlVect* std_dev)
{
    // calcul du delta en t
    int n = (opt_->nbTimeSteps_)+1;
    int nombreActifs = (mod_->size_);
    double timestep = (opt_->T_)/(n-1);
    PnlVect* actionPricesInt = pnl_vect_create(nombreActifs);
    pnl_mat_get_row(actionPricesInt, past, past->m -1);
    PnlMat* path = pnl_mat_create (n, nombreActifs);
    PnlMat* shift_path_up = pnl_mat_create (n, nombreActifs);
    PnlMat* shift_path_down = pnl_mat_create (n, nombreActifs);
    double sum = 0;
    double sumDesCarres = 0;
    for(int roundOnAction = 0; roundOnAction < nombreActifs; roundOnAction++)
    {
        double valueIntPerAction = pnl_vect_get(actionPricesInt, roundOnAction);
        for(int round = 0; round < nbSamples_; round++)
        {
            mod_->asset(path, t, opt_->T_, n-1, rng_, past);
            mod_->shiftAsset(shift_path_up, path, roundOnAction, fdStep_, t, timestep);
            mod_->shiftAsset(shift_path_down, path, roundOnAction, -fdStep_, t, timestep);
            sum += opt_->payoff(shift_path_up) - opt_->payoff(shift_path_down);
            sumDesCarres+= pow(opt_->payoff(shift_path_up)- opt_->payoff(shift_path_down), 2);
        }
        double deltaPerAction = (sum*exp(-2*(mod_->r_)*(opt_->T_-t)))/(2*nbSamples_*fdStep_*valueIntPerAction);
        pnl_vect_set(delta, roundOnAction, deltaPerAction); //valeurs des deltas
        double s = pow((sum/(2*nbSamples_*valueIntPerAction*fdStep_)),2);
        double volatilityPerAction = exp(-2*(mod_->r_)*(opt_->T_-t))*(sumDesCarres/(2*nbSamples_*valueIntPerAction*fdStep_)-s);
        double std_devElement = sqrt(volatilityPerAction/nbSamples_); 
        pnl_vect_set(std_dev, roundOnAction, std_devElement); // valeurs des ecart-type 
    }   
}

