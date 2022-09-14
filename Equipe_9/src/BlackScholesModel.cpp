#include <iostream>
#include "BlackScholesModel.hpp"




BlackScholesModel::BlackScholesModel(int size, double r, double rho, PnlVect* sigma, PnlVect* spot)
{
    size_ = size;
    r_ = r;
    rho_ = rho;
    sigma_ = sigma;
    spot_ = spot; 
}


/**
 * Génère une trajectoire du modèle et la stocke dans path
 *
 * @param[out] path contient une trajectoire du modèle.
 * C'est une matrice de taille (nbTimeSteps+1) x d
 * @param[in] T  maturité
 * @param[in] nbTimeSteps nombre de dates de constatation
 */
void BlackScholesModel::asset(PnlMat* path, double T, int nbTimeSteps, PnlRng* rng)
{
    //Create Cholesky
    double volatility;
    double scaleCholeskyGaussian;
    double computedSpot;
    double timeStep = T / nbTimeSteps;
    PnlMat *cholesky = pnl_mat_create_from_scalar(size_, size_, rho_);
    for(int i = 0; i < size_; i++)
    {
        pnl_mat_set_diag(cholesky, rho_, 1);
    }
    pnl_mat_chol(cholesky);
    PnlMat * gaussian = pnl_mat_create(nbTimeSteps + 1, size_);
    pnl_mat_rng_normal(gaussian, nbTimeSteps + 1, size_, rng);
    pnl_mat_set_row(path, spot_, 0);
    //initiating
    PnlVect * choleskyComponent = pnl_vect_create(size_);
    PnlVect * gaussianVector = pnl_vect_create(nbTimeSteps + 1);
    PnlVect * nextSpots = pnl_vect_create(size_);
    PnlVect * currentSpots = pnl_vect_create(size_);
    pnl_mat_get_row(currentSpots, path, 0);
    for(int i = 1; i < nbTimeSteps + 1; i ++)
    {
        pnl_mat_get_row(gaussianVector, gaussian, i); //or i - 1, check that later
        //compute components of nextspots ( vector ; each element is an underlying share)
        for(int d = 0; d < size_; d++)
        {
            pnl_mat_get_row(choleskyComponent, cholesky, d);
            volatility = pnl_vect_get(sigma_, d);
            scaleCholeskyGaussian = pnl_vect_scalar_prod(choleskyComponent, gaussianVector);
            computedSpot = pnl_vect_get(currentSpots, d) * exp((r_ - volatility * volatility / 2) * timeStep + volatility * sqrt(timeStep) * scaleCholeskyGaussian);
            pnl_vect_set(nextSpots, d, computedSpot);
        }
        pnl_mat_set_row(path, nextSpots, i);
        * currentSpots = * nextSpots;
    }
    
}

/**
 * Calcule une trajectoire du modèle connaissant le
 * passé jusqu' à la date t
 *
 * @param[out] path  contient une trajectoire du sous-jacent
 * donnée jusqu'à l'instant t par la matrice past
 * @param[in] t date jusqu'à laquelle on connait la trajectoire.
 * t n'est pas forcément une date de discrétisation
 * @param[in] nbTimeSteps nombre de pas de constatation
 * @param[in] T date jusqu'à laquelle on simule la trajectoire
 * @param[in] past trajectoire réalisée jusqu'a la date t
 */
void BlackScholesModel::asset(PnlMat* path, double t, double T, int nbTimeSteps, PnlRng* rng, const PnlMat* past)
{   
    double volatility;
    double scaleCholeskyGaussian;
    double computedSpot;
    double timeStep = T / nbTimeSteps;
    PnlMat *cholesky = pnl_mat_create_from_scalar(size_, size_, rho_);
    for(int d = 0; d < size_; d++)
    {
        pnl_mat_set_diag(cholesky, d, 1);
    }
    pnl_mat_chol(cholesky);
    pnl_mat_clone(path, past);
    int startingStep = past->m;
    PnlMat * gaussian = pnl_mat_create(nbTimeSteps + 1 - startingStep, size_);
    pnl_mat_rng_normal(gaussian, nbTimeSteps + 1, size_, rng);
    pnl_mat_set_row(path, spot_, 0);
    //initiating
    PnlVect * choleskyComponent = pnl_vect_create(size_);
    PnlVect * gaussianVector = pnl_vect_create(nbTimeSteps + 1 - startingStep);
    PnlVect * nextSpots = pnl_vect_create(size_);
    PnlVect * currentSpots = pnl_vect_create(size_);
    pnl_mat_get_row(currentSpots, path, startingStep);
    for(int i = startingStep + 1; i < nbTimeSteps + 1; i++)
    {
        pnl_mat_get_row(gaussianVector, gaussian, i - startingStep); //or i - 1, check that later
        //compute components of nextspots ( vector ; each element is an underlying share)
        for(int d = 0; d < size_; d++)
        {
            pnl_mat_get_row(choleskyComponent, cholesky, d);
            volatility = pnl_vect_get(sigma_, d);
            scaleCholeskyGaussian = pnl_vect_scalar_prod(choleskyComponent, gaussianVector);
            computedSpot = pnl_vect_get(currentSpots, d) * exp((r_ - volatility * volatility / 2) * (i * T / nbTimeSteps - t) + volatility * sqrt(i * T / nbTimeSteps - t) * scaleCholeskyGaussian);
            pnl_vect_set(nextSpots, d, computedSpot);
        }
        pnl_mat_set_row(path, nextSpots, i);
        * currentSpots = * nextSpots;
    }
}