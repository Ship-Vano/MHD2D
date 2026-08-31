//
// Created by Иван on 10/21/2024.
//

#ifndef MAGNETTOPRJCT_MHDSOLVER2D_H
#define MAGNETTOPRJCT_MHDSOLVER2D_H

#include "NetGeometry.h"
#include "MHDSolver1D.h"

#ifdef USE_CUDA
#include "MHDgpu.cuh"
#endif

#include "omp.h"

#include <numbers>

struct DivergenceMetrics {
    double maxAbs = 0.0;
    double maxFlux = 0.0;
    // |Phi_K| / (B_ref * perimeter_K), with B_ref fixed at initialization.
    double maxScaled = 0.0;
    double referenceField = 0.0;
    // Secondary diagnostic; omitted from acceptance gates in near-zero B.
    double maxLocalScaled = 0.0;
};

class MHDSolver2D {
public:
    // mesh
    World geometryWorld;

    // physical quantities
    double gam_hcr = 2.0; // показатель адиабаты
    double startTime = 0.0;  // время отсчёта
    double finalTime = 0.1;   // время окончания
    double tau = 0.0001; // шаг по времени
    double min_tau = 1e-7; // минимальный шаг по времени
    double cflNum = 0.4; // число Куранта
    double configuredCfl = -1.0; // negative preserves a task default
    double fieldLoopU = 2.0;
    double fieldLoopV = 1.0;
    double fieldLoopRadius = 0.3;
    double fieldLoopAmplitude = 1.0e-3;
    int MAX_ITERATIONS = 1000000; // максимальное число итераций
    int iterationsPerFrame = 10; // число итераций для записи в файл
    int task_type = 1; // тип задачи
    bool periodicBoundaries = false; // флаг: периодические г.у.
    bool freeFLowBoundaries = false; // флаг: г.у. на свободный поток на границе
    bool freeFLowBoundaries2 = false; // флаг: г.у. на свободный поток на границе 2-го типа (граничные ячейки - четырёхугольники)
    bool debugDivergence = false;
    int innerElemCount = 0;
    int innerEdgeCount = 0;
    bool ghostOutput = false;

    //states
    std::vector<std::vector<double>> nodeUs; // state U in nodes
    std::vector<std::vector<double>> elemUs; // state U in elements
    std::vector<std::vector<double>> ghostElemUs; // state U in ghost elements (X_x)
    std::vector<std::vector<double>> edgeUs; // state U in edges
    std::vector<std::vector<double>> initElemUs;
    std::vector<std::vector<double>> initGhostElemUs; // state U in ghost elements (X_x)
    std::vector<std::vector<double>> initEdgeUs;
    std::vector<double> initBns;
    std::vector<double> bNs; //Bns at edges
    std::vector<double> initGhostBNs; // state U in ghost elements (X_x)
    std::vector<double> ghostBNs; //Bns at ghost edges

    std::vector<double> rotateStateFromNormalToAxisX(vector<double> &U, const Vec2 &n);
    std::vector<double> rotateStateFromAxisToNormal(std::vector<double>& U, const Vec2 &n);
    double tau_from_cfl2D(const double& sigma, const double& hx, const std::vector<std::vector<double>>& states, const double& gam_hcr);
    double tau_from_cfl2D(const double& sigma, const double& min_h, std::vector<std::vector<double>>& edgeStates, const double& gam_hcr,
                          const EdgePool& ep);
    // Conservative first-order finite-volume CFL bound for arbitrary triangles.
    double stableTimeStepUnstructured(const std::vector<std::vector<double>>& states,
                                      const std::vector<double>& faceNormalB,
                                      const double& gam_hcr);
    void setInitElemUs();
    void runSolver();

    void setInitCylindricElemUs();
    void runCylindricSolver();
    void applyBoundaryConditions(NeighbourService& ns);
    void applyGhostCells2Type();
    void applyZeroRConditions(const ElementPool& elPool, const EdgePool& edgePool, const NodePool& nodePool, const std::vector<std::vector<double>>& elemUs_prev);

    double computeDivergence();
    DivergenceMetrics computeDivergenceMetrics(double referenceField = -1.0);
    void checkNan(bool& foundNan);

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
    void writeVTU(const std::string& filename, const bool& ghost);

    void runGPUSolver();
};

//void writeVTU(const std::string& filename, const World& geometryWorld, const std::vector<std::vector<double>>& elemUs);

double computeDivergence(const std::vector<std::vector<double>>& elemUs, const EdgePool& edgePool);
// Sign relating stored nodeInd1 -> nodeInd2 to the CT tangent
// t_CT=(-n_y,n_x).  It is +1 for the customary CCW physical-element
// winding, but must be applied explicitly for a valid clockwise input mesh.
double edgeCtOrientation(const Edge& edge, const NodePool& nodePool);
double advanceFaceNormalBFromCt(const Edge& edge, const NodePool& nodePool,
                                double previousNormalB, double dt,
                                double firstNodeG, double secondNodeG);
// RT0 reconstruction at the element centroid from edge-normal magnetic
// degrees of freedom.  faceNormalB is oriented with Edge::neighbourInd1;
// this routine applies the opposite sign for neighbourInd2 explicitly.
// Physical elements must reference physical edges only.
Vec2 reconstructCellMagneticFieldRT0(const Element& element,
                                     const NodePool& nodePool,
                                     const EdgePool& edgePool,
                                     const std::vector<double>& faceNormalB,
                                     int physicalEdgeCount);
std::vector<int> findCommonElements(const std::vector<int>& v1, const std::vector<int>& v2);
/*                       0      1      2      3    4   5   6   7
 * U (general state):  rho,  rho*u, rho*v, rho*w,  e,  Bx, By, Bz
 * gasU                rho,  rho*u, rho*v, rho*w,  e
 * magU                 Bx,    By,    Bz
 * */

//void solverHLL2D(const World& world);

#endif //MAGNETTOPRJCT_MHDSOLVER2D_H
