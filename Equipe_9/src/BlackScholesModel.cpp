#include <iostream>
#include "BlackScholesModel.hpp"

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
    pnlMat *cholesky = pnl_mat_create_from_scalar(size_, size_, rho_);
    for(int i = 0; i < size_; i++)
    {
        pnl_mat_set_diag(cholesky, rho_, 1);
    }
    int choleskyDone = pnl_mat_chol(cholesky);
    pnlMat * gaussian = pnl_mat_create(nbTimeSteps + 1, size_);
    pnl_mat_rng_normal(gaussian, nbTimeSteps + 1, size_, rng);
    pnl_mat_set_row(path, spot_, 0);
    //initiating
    pnlVect * choleskyComponent = pnl_vect_create(size_);
    pnlVect * gaussianVector = pnl_vect_create(nbTimeSteps + 1);
    pnlVect * nextSpots = pnl_vect_create(size_);
    pnlVect * currentSpots = pnl_vect_create(size_);
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
            computedSpot = pnl_vect_get(currentSpots, d) * exp((r - volatility * volatility / 2) * timeStep + volatility * sqrt(timeStep) * scaleCholeskyGaussian);
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
void asset(PnlMat* path, double t, double T, int nbTimeSteps, PnlRng* rng, const PnlMat* past)
{
    
}