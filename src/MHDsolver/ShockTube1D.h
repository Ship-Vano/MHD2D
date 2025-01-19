//
// Created by Иван on 4/30/2024.
//

#ifndef MAGNETTOPRJCT_SHOCKTUBE1D_H
#define MAGNETTOPRJCT_SHOCKTUBE1D_H

#include"MHDSolver1D.h"

class ShockTube1D {
public:
    double gam_hcr;

    double rhoL;
    double uL;
    double vL;
    double wL;
    double pL;
    double BxL;
    double ByL;
    double BzL;

    double x_mid;

    double rhoR;
    double uR;
    double vR;
    double wR;
    double pR;
    double BxR;
    double ByR;
    double BzR;

    std::function<std::vector<double>(double)> leftBound;
    bool leftBound_is_set = false;

    std::function<std::vector<double>(double)> rightBound;
    bool rightBound_is_set = false;

    std::function<std::vector<double>(double)> initDistrib;
    bool initDistrib_is_set = false;

    ShockTube1D(const double &gam_hcr_init, const std::vector<double> &L_state, const std::vector<double> &R_state, const double &mid_point);
    ShockTube1D(std::string file_loc);

private:
    void file_init(std::string file_loc); // функция инициализации через файл
};


#endif //MAGNETTOPRJCT_SHOCKTUBE1D_H
