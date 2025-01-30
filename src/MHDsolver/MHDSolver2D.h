//
// Created by Иван on 10/21/2024.
//

#ifndef MAGNETTOPRJCT_MHDSOLVER2D_H
#define MAGNETTOPRJCT_MHDSOLVER2D_H

#include "NetGeometry.h"
#include "MHDSolver1D.h"
#include "omp.h"

class MHDSolver2D {
public:
    // mesh
    World geometryWorld;

    // physical quantities
    double gam_hcr = 2.0;
    double startTime = 0.0;  // время отсчёта
    double finalTime = 0.1;   // время окончания
    double tau = 0.0001; // шаг по времени
    double min_tau = 1e-7;
    double cflNum = 0.4; // число Куранта
    int MAX_ITERATIONS = 1000000;
    int iterationsPerFrame = 10;
    int task_type = 1;
    bool periodicBoundaries = false;
    bool debugDivergence = false;

    //states
    std::vector<std::vector<double>> nodeUs; // state U in nodes
    std::vector<std::vector<double>> elemUs; // state U in elements
    std::vector<std::vector<double>> edgeUs; // state U in edges
    std::vector<std::vector<double>> initElemUs;
    std::vector<std::vector<double>> initEdgeUs;
    std::vector<double> initBns;
    std::vector<double> bNs; //Bns at edges

    std::vector<double> rotateStateFromNormalToAxisX(std::vector<double>& U, const std::vector<double>& n);
    std::vector<double> rotateStateFromAxisToNormal(std::vector<double>& U, const std::vector<double>& n);
    double tau_from_cfl2D(const double& sigma, const double& hx, const double& hy, const std::vector<std::vector<double>>& states, const double& gam_hcr);
    double tau_from_cfl2D(const double& sigma, const double& min_h, std::vector<std::vector<double>>& edgeStates, const double& gam_hcr,
                          const EdgePool& ep);
    void setInitElemUs();
    void runSolver();

    // Начальное состояние системы
    std::function<std::vector<double>(double)> initStateFunc;
    bool initStateFunc_is_set = false;

    // Левые граничные условия
    std::function<std::vector<double>(double)> leftBoundaryFunction;
    bool leftBoundaryFunction_is_set = false;

    // Правые граничные условия
    std::function<std::vector<double>(double)> rightBoundaryFunction;
    bool rightBoundaryFunction_is_set = false;
    std::vector<double> applyLimiter(const std::vector<double>& U_left, const std::vector<double>& U_center, const std::vector<double>& U_right);
    MHDSolver2D(const World& world);
};

void writeVTU(const std::string& filename, const World& geometryWorld, const std::vector<std::vector<double>>& elemUs);

double computeDivergence(const std::vector<std::vector<double>>& elemUs, const EdgePool& edgePool);

/*                       0      1      2      3    4   5   6   7
 * U (general state):  rho,  rho*u, rho*v, rho*w,  e,  Bx, Bz, By
 * gasU                rho,  rho*u, rho*v, rho*w,  e
 * magU                 Bx,    Bz,    By
 * */

//void solverHLL2D(const World& world);

#endif //MAGNETTOPRJCT_MHDSOLVER2D_H
