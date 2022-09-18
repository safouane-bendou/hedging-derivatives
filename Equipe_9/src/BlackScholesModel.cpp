#include <iostream>
#include "BlackScholesModel.hpp"
using namespace std;




void BlackScholesModel::choleskyComposition(PnlMat* cholesky)
{
    for(int d = 0; d < size_; d++)
    {
        pnl_mat_set_diag(cholesky, 1, d);
    }
    pnl_mat_chol(cholesky);
}

BlackScholesModel::BlackScholesModel(int size, double r, double rho, PnlVect* sigma, PnlVect* spot)
{
    size_ = size;
    r_ = r;
    rho_ = rho;
    sigma_ = pnl_vect_create(size_);
    pnl_vect_clone(sigma_, sigma);
    spot_ = pnl_vect_create(size_);
    pnl_vect_clone(spot_, spot);
    cholesky = pnl_mat_create_from_scalar(size_, size_, rho_);
    choleskyComposition(cholesky);
    choleskyComponent = pnl_vect_create(size_);
    gaussianVector = pnl_vect_create(size_);
    currentShares = pnl_vect_create(size_);
    nextShares = pnl_vect_create(size_);
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
    double timeStep = T / nbTimeSteps;
    PnlMat * gaussian = pnl_mat_create(nbTimeSteps, size_);
    pnl_mat_rng_normal(gaussian, nbTimeSteps + 1, size_, rng);
    pnl_mat_set_row(path, spot_, 0);
    pnl_mat_get_row(currentShares, path, 0);
    for(int i = 1; i < nbTimeSteps + 1; i++)
    {
        pnl_mat_get_row(gaussianVector, gaussian, i);
        for(int d = 0; d < size_; d++)
        {
            pnl_mat_get_row(choleskyComponent, cholesky, d);
            double deviation = pnl_vect_get(sigma_, d);
            double scaleCholeskyGaussian = pnl_vect_scalar_prod(choleskyComponent, gaussianVector);
            double computedShare = pnl_vect_get(currentShares, d) * exp((r_ - deviation * deviation / 2) * timeStep + deviation * sqrt(timeStep) * scaleCholeskyGaussian);
            pnl_vect_set(nextShares, d, computedShare);
        }
        pnl_mat_set_row(path, nextShares, i);
        currentShares = nextShares;
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
    double timeStep = T / nbTimeSteps;
    int startingStep;
    PnlVect * savedSpot = pnl_vect_create(size_);
    pnl_mat_get_row(savedSpot, past, past->m - 1);   // get last row of past
    PnlMat * uniformPast = pnl_mat_create(past->m - 1, size_);
    pnl_mat_extract_subblock(uniformPast, past, 0, past->m - 1, 0, size_);
    if(t / timeStep == (int)(t / timeStep) or t / timeStep < 1)
    {
        pnl_mat_set_subblock(path, past, 0, 0);
        startingStep = past->m;
    }
    else
    {
        pnl_mat_set_subblock(path, uniformPast, 0, 0);
        startingStep = uniformPast->m;
    }

    PnlMat * gaussian = pnl_mat_create(nbTimeSteps + 1 - startingStep, size_);
    pnl_mat_rng_normal(gaussian, nbTimeSteps + 1, size_, rng);
    for(int i = startingStep; i < nbTimeSteps + 1; i++)
    {
        pnl_mat_get_row(gaussianVector, gaussian, i);
        for(int d = 0; d < size_; d++)
        {
            pnl_mat_get_row(choleskyComponent, cholesky, d);
            double deviation = pnl_vect_get(sigma_, d);
            double scaleCholeskyGaussian = pnl_vect_scalar_prod(choleskyComponent, gaussianVector);
            double computedShare = pnl_vect_get(savedSpot, d) * exp((r_ - deviation * deviation / 2) * (i * timeStep - t) + deviation * sqrt(i * timeStep - t) * scaleCholeskyGaussian);
            pnl_vect_set(nextShares, d, computedShare);
        }
        pnl_mat_set_row(path, nextShares, i);
    }
}


/**
 * Shift d'une trajectoire du sous-jacent
 *
 * @param[in]  path contient en input la trajectoire
 * du sous-jacent
 * @param[out] shift_path contient la trajectoire path
 * dont la composante d a été shiftée par (1+h)
 * à partir de la date t.
 * @param[in] t date à partir de laquelle on shift
 * @param[in] h pas de différences finies
 * @param[in] d indice du sous-jacent à shifter
 * @param[in] timestep pas de constatation du sous-jacent
 */
void BlackScholesModel::shiftAsset(PnlMat* shift_path, const PnlMat* path, int d, double h, double t, double timestep)
{
    int startingStep;
    pnl_mat_clone(shift_path, path);
    if(t / timestep == (int)(t / timestep))
    {
        startingStep = (int)(t / timestep);
    }
    else
    {
        startingStep = (int)(t / timestep) + 1;
    }

    PnlVect* shiftComponent = pnl_vect_create(size_);
    pnl_mat_get_col(shiftComponent, shift_path, d);
    for(int i = startingStep; i < shiftComponent->size; i++)
    {
        double shifted = (1 + h) * pnl_vect_get(shiftComponent, i);
        pnl_vect_set(shiftComponent, i, shifted);
    }
    pnl_mat_set_col(shift_path, shiftComponent, d);

}
















