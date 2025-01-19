//
// Created by Иван on 4/29/2024.
//

#ifndef MAGNETTOPRJCT_MHDPROBLEM1D_H
#define MAGNETTOPRJCT_MHDPROBLEM1D_H

#include <cmath>
#include <iostream>
#include <functional>

class MHDProblem1D {
public:
    bool periodic_boundaries = false;
    double gam_hcr; // heat capacity rate
    double x0;  // начало отсчёта по пространству (координата левого конца)
    double L;   // характерный пространственный рамзер (длина струны)
    double X;   // координата правого конца
    double t0;  // время отсчёта
    double T;   // время окончания
    double tau; // шаг по времени
    double h;   // шаг по пространству
    double gam_courant; // число Куранта
    int num_time_steps; // количество шагов по времени
    int num_space_steps;// количество шагов по пространству
    MHDProblem1D(double gam_hcr_init, double x0_init, double L_init, double t0_init, double T_init, double h_init, double tau_init, double gam_courant_init, bool what_is_L_init);

    // Начальное состояние системы
    std::function<std::vector<double>(double)> initStateFunc;
    bool initStateFunc_is_set = false;

    // Левые граничные условия
    std::function<std::vector<double>(double)> leftBoundaryFunction;
    bool leftBoundaryFunction_is_set = false;

    // Правые граничные условия
    std::function<std::vector<double>(double)> rightBoundaryFunction;
    bool rightBoundaryFunction_is_set = false;
};

#endif //MAGNETTOPRJCT_MHDPROBLEM1D_H

