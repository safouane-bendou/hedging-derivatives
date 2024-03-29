#pragma once
#include <iostream>

#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

/// \brief Modèle de Black Scholes
class BlackScholesModel
{
  public:
    int size_;       /// nombre d'actifs du modèle
    double r_;       /// taux d'intérêt
    double rho_;     /// paramètre de corrélation
    PnlVect* sigma_; /// vecteur de volatilités
    PnlVect* spot_;  /// valeurs initiales des sous-jacents
    PnlVect* trend_;
    PnlMat* cholesky;
    PnlVect* choleskyComponent;
    PnlVect* gaussianVector;
    PnlVect* currentShares;
    PnlVect* nextShares;

    /**
     Constructeur à quatre paramètres définissant le modèle BlackScholes 
     */

    BlackScholesModel(int size, double r, double rho, PnlVect* sigma, PnlVect* spot);



    /**
     Constructeur à cinq paramètres prenant en argument le trend afin de simuler les données de marché
     */

    BlackScholesModel(int size, double r, double rho, PnlVect* sigma, PnlVect* spot, PnlVect* trend);



    /*
    calcule la composition de la matrice de cholesky
    */
    void choleskyComposition(PnlMat* cholesky);


    

    /**
     * Génère une trajectoire du modèle et la stocke dans path
     *
     * @param[out] path contient une trajectoire du modèle.
     * C'est une matrice de taille (nbTimeSteps+1) x d
     * @param[in] T  maturité
     * @param[in] nbTimeSteps nombre de dates de constatation
     */
    void asset(PnlMat* path, double T, int nbTimeSteps, PnlRng* rng);

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
    void asset(PnlMat* path, double t, double T, int nbTimeSteps, PnlRng* rng, const PnlMat* past);

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
    void shiftAsset(PnlMat* shift_path, const PnlMat* path, int d, double h, double t, double timestep);



    /**
     * renvoie une simulation du marché
     * @param[out] simulatedData  contient la trajectoire du marché simulé
     * @param[in] H nombre de dates
     * @param[in] T date jusqu'à laquelle on simule la trajectoire
     */
    void simul_market(PnlMat* simulatedData, double H, double T, PnlRng* rng);

    //Destructeur de la classe
    ~BlackScholesModel();
};
