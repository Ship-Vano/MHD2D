//
// Created by Иван on 10/21/2024.
//

#include "MHDSolver2D.h"

double gas_energy(const double& gam_hcr, const double& gas_p, const double& rho, const double& u, const double& v, const double& w){
    return gas_p / (gam_hcr - 1.0) + rho * (u * u + v * v + w * w) / 2.0;
}

double gas_p(const double& gam_hcr, const double& e, const double& rho, const double& u, const double& v, const double& w, const double &Bx, const double &By, const double &Bz){
    return (gam_hcr - 1.0) * (e -  0.5*(Bx * Bx + By * By + Bz * Bz) - (rho/2.0) * (u*u + v*v + w*w));
}

// Вектор состояния из параметров
std::vector<double>
state_from_primitive_vars2D(const double &rho, const double &u, const double &v, const double &w, const double &p,
                            const double &Bx, const double &By, const double &Bz, const double &gam_hcr) {
    std::vector<double> U(8,0.0);

    double mx = rho * u;
    double my = rho * v;
    double mz = rho * w;
    double e = energy(gam_hcr, p, rho, u, v, w, Bx, By, Bz);
    //double e = gas_energy(gam_hcr, p, rho, u, v, w);
    U[0] = rho;
    U[1] = mx;
    U[2] = my;
    U[3] = mz;
    U[4] = e;
    U[5] = Bx;
    U[6] = By;
    U[7] = Bz;

    return U;
}

std::vector<double> state_from_primitive_vars2D(const std::vector<double>& primitiveVars){
    /*rho  u   v   w   p   Bx   By   Bz  gam_hcr*/
    return state_from_primitive_vars2D(primitiveVars[0], primitiveVars[1], primitiveVars[2],
                                       primitiveVars[3], primitiveVars[4], primitiveVars[5],
                                       primitiveVars[6], primitiveVars[7], primitiveVars[8]);
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// конструктор решателя
MHDSolver2D::MHDSolver2D(const World &world): geometryWorld(world), nodeUs(world.getNodePool().nodeCount),
    elemUs(world.getElementPool().elCount - world.ghostElemCount), edgeUs(world.getEdgePool().edgeCount - world.ghostElemCount * 2), initElemUs(world.getElementPool().elCount - world.ghostElemCount),
    initBns(world.getEdgePool().edgeCount - world.ghostElemCount*2), bNs(world.getEdgePool().edgeCount),
    ghostElemUs(world.ghostElemCount), ghostBNs(2*world.ghostElemCount){}

// Minmod-лимитер (пригодится позже)
double minmod(double a, double b) {
    if (a * b > 0.0) {
        return std::copysign(std::min(std::abs(a), std::abs(b)), a);
    }
    return 0.0;
}

// возмжное применение minmod-лимитера (не тестировал, набросок)ы
std::vector<double> MHDSolver2D::applyLimiter(const std::vector<double>& U_left, const std::vector<double>& U_center, const std::vector<double>& U_right) {
    std::vector<double> limited(U_center.size());
    for (size_t i = 0; i < U_center.size(); ++i) {
        double slope_left = U_center[i] - U_left[i];
        double slope_right = U_right[i] - U_center[i];
        double slope_center = (U_right[i] - U_left[i]) / 2.0;
        double limited_slope = minmod(slope_left, slope_center);
        limited[i] = U_center[i] - 0.5 * limited_slope;
    }
    return limited;
}

/*вращения*/
// n = {cos(j), sin(j)} = {n.x, n.y}

// ({cos(j), sin(j), 0}, {-sin(j), cos(j), 0}, {0,0,1}) --- around OZ clockwise
// rotate 1: from normal to OX:      {v.x * n.x + vy * n.y, - v.x * n.y + v.y * n.x, v.z}
// вращаем {u,v,w} и {Bx, By, Bz}, остальные остаются на месте
std::vector<double> MHDSolver2D::rotateStateFromAxisToNormal(std::vector<double> &U, const std::vector<double>& n) {
    std::vector<double> res(U);
    res[1] =  U[1]*n[0] + U[2]*n[1];
    res[2] = -U[1]*n[1] + U[2]*n[0];
    res[5] =  U[5]*n[0] + U[6]*n[1];
    res[6] = -U[5]*n[1] + U[6]*n[0];

    /*double vel0 = std::sqrt(U[1] * U[1] + U[2] * U[2] + U[3] * U[3]);
    double vel1 = std::sqrt(res[1] * res[1] + res[2] * res[2] + res[3] * res[3]);
    if(std::abs(vel0 - vel1) > 1e-12){
        //std::cout << "Unequal velocity modules: " <<std::scientific << vel0 << " vs " << vel1 << std::endl;
        //std::cin.get();
    }*/
    //double p = pressure(res, gam_hcr);
    //double energy_rotated = energy(gam_hcr, p, res[0], res[1]/res[0], res[2]/res[0], res[3]/res[0], res[5], res[6], res[7]);
    //res[4] = energy_rotated;
    /*if(std::abs(energy_rotated - U[4]) > 1e-12) {
        std::cout << "Unequal energy modules: " << std::scientific << U[4] << " vs " << energy_rotated << std::endl;
        std::cin.get();
    }*/
    /*double B0 = std::sqrt(U[5] * U[5] + U[6] * U[6] + U[7] * U[7]);
    double B1 = std::sqrt(res[5] * res[5] + res[6] * res[6] + res[7] * res[7]);
    if(std::abs(B0 - B1) > 1e-12){
        //std::cout << "Unequal B modules: " << std::scientific << B0 << " vs " << B1 << std::endl;
        //std::cin.get();
    }*/
    return res;
}

// ({cos(j), -sin(j), 0}, {sin(j), cos(j), 0}, {0,0,1}) --- around OZ counterclockwise
// rotate 2: from OX to normal:      {v.x * n.x - vy * n.y,  v.x * n.y + v.y * n.x, v.z}
std::vector<double> MHDSolver2D::rotateStateFromNormalToAxisX(vector<double> &U, const vector<double>& n) {
    std::vector<double> res(U);
    res[1] =  U[1]*n[0] - U[2]*n[1];
    res[2] =  U[1]*n[1] + U[2]*n[0];
    res[5] =  U[5]*n[0] - U[6]*n[1];
    res[6] =  U[5]*n[1] + U[6]*n[0];
    //double p = pressure(res, gam_hcr);
    //double energy_rotated = energy(gam_hcr, p, res[0], res[1]/res[0], res[2]/res[0], res[3]/res[0], res[5], res[6], res[7]);
    //res[4] = energy_rotated;
    return res;
}

// условие Куранта (имплементация №1) (элементы)
double MHDSolver2D::tau_from_cfl2D(const double& sigma, const double& hx, const std::vector<std::vector<double>>& states, const double& gam_hcr) {
    double max_speed = 0.0;
    for (const auto& state : states) {
        // для cfast брать модули в-в v и B
        double u = std::fabs(state[1] / state[0]);
        double v = std::fabs(state[2] / state[0]);
        double w = std::fabs(state[3]/state[0]);
        double cf = cfast(state, gam_hcr);
        //double local_speed = (u + cf) / hx + (v + cf) / hy;
        //max_speed = std::max(max_speed, local_speed);
        max_speed = std::max(max_speed, std::sqrt(u*u + v*v + w*w) + cf);
    }
    if(max_speed < 1e-16){
        max_speed = 1e-14;
    }
    double tau = sigma * hx / max_speed;
    const double max_tau = 1e-1; // Define maximum allowable time step
    return std::min(tau, max_tau);
}

// условие Куранта (имплементация №2) (Рёбра)
double MHDSolver2D::tau_from_cfl2D(const double& sigma, const double& min_h, std::vector<std::vector<double>>& edgeStates, const double& gam_hcr,
                      const EdgePool& ep) {
    double max_speed = 0.0;
    for (const auto& edge : ep.edges) {
        // не вращать, а передавать моули в-в
        const vector<double> rotatedState = rotateStateFromAxisToNormal(edgeStates[edge.ind], edge.normalVector);
        double u = std::fabs(rotatedState[1] / rotatedState[0]);
        double cf = cfast(rotatedState, gam_hcr);
        double local_speed = (u + cf);
        max_speed = std::max(max_speed, local_speed);
    }
    max_speed = std::max(max_speed, 1e-14); // Prevent division by zero
    double optimal_tau = sigma * min_h / max_speed;
    return std::min(optimal_tau, 1e-3);
}

// задание начальных условий
void MHDSolver2D::setInitElemUs() {

    ElementPool ep = geometryWorld.getElementPool();
    innerElemCount = ep.elCount - geometryWorld.ghostElemCount;
    std::cout << "Element count = " << ep.elCount << " ; innerElCount = " << innerElemCount << " ; ghostElCOunt = " << geometryWorld.ghostElemCount << std::endl;
    EdgePool edgp = geometryWorld.getEdgePool();
    innerEdgeCount = edgp.edgeCount - geometryWorld.ghostElemCount * 2;

    // Brio-Wu problem
    if(task_type == 1) {
        cflNum = 0.9;
        freeFLowBoundaries = true;
        periodicBoundaries = false;
        std::cout << "SOLVING TASKTYPE 1 (BRIO-WU TEST)" << std::endl;
                                     /*rho  u   v   w   p   Bx   By   Bz gam_hcr*/
        std::vector<double> BrioWu_L1{1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0, gam_hcr};
        std::vector<double> BrioWu_R1{0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0, gam_hcr};

        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = ep.elements[i];
            std::vector<double> centroid = elem.centroid2D;
            int elInd = elem.ind;
            if (centroid[0] < 0.5) {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_L1);
            } else {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_R1);
            }
        }

        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for(int i = innerElemCount; i < ep.elCount; ++i ){
            Element elem = ep.elements[i];
            std::vector<double> centroid = elem.centroid2D;
            int elInd = elem.ind;
            int ghostInd = elInd - innerElemCount;
            if (centroid[0] < 0.5) {
                initGhostElemUs[ghostInd] = state_from_primitive_vars(BrioWu_L1);
            } else {
                initGhostElemUs[ghostInd] = state_from_primitive_vars(BrioWu_R1);
            }
        }

        initBns.resize(edgp.edgeCount, 0.0);
        initGhostBNs.resize(geometryWorld.ghostElemCount*2, 0.0);
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < edgp.edgeCount; ++i) {
            Edge edge = edgp.edges[i];
            std::vector<double> midP = edge.midPoint;
            int edgeInd = edge.ind;
            std::vector<double> normal = edge.normalVector;
            std::vector<double> state;
            if (midP[0] < 0.5) {
                state = state_from_primitive_vars(BrioWu_L1);
            } else {
                state = state_from_primitive_vars(BrioWu_R1);
            }
            double Bn = state[5] * normal[0] + state[6] * normal[1];
            initBns[edgeInd] = Bn;
            if(i < innerEdgeCount) {
                initBns[edgeInd] = Bn;
            }else{
                int ghostInd = edgeInd - innerEdgeCount;
                initGhostBNs[ghostInd] = Bn;
            }
            initEdgeUs[edgeInd] = state;
        }
        NeighbourService ns = geometryWorld.getNeighbourService();
        for(const auto& [boundary, ghost]: ns.boundaryToGhostElements){
            int ghostInd = ghost - innerElemCount;
            initGhostElemUs[ghostInd] =  initElemUs[boundary];
        }
    }

    // alfven wave test
    else if(task_type == 2){
        std::cout << "SOLVING TASKTYPE 2 (ALFVEN WAVE TEST)" << std::endl;
        double r = 6.0;
        double nx = 0.0/*1.0 / std::sqrt(r * r + 1)*/;
        double ny = 1.0/*r / std::sqrt(r * r + 1)*/;
        double rho0 = 1.0;
        double p0 = 1.0;
        double v0 = 0.0;
        double B0 = 1.0/ (std::sqrt(4.0 * M_PI));
        double xi = 0.2;
        gam_hcr = 5.0/3.0;
        cflNum = 0.4;
        periodicBoundaries = true;
        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        initBns.resize(edgp.edgeCount, 0.0);
        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < ep.elCount; ++i) {
            Element elem = ep.elements[i];
            std::vector<double> centroid = elem.centroid2D;
            double x = centroid[0];
            double y = centroid[1];
            double phase = 2.0 * M_PI / ny * (nx * x + ny * y);
            double u_0 = v0 * nx - xi * ny * std::cos(phase);
            double v_0 = v0 * ny + xi * nx * std::cos(phase);
            double w_0 = xi * std::sin(phase);
            double Bx_0 = B0 * nx + xi * ny * std::sqrt(4 * M_PI * rho0) * std::cos(phase);
            double By_0 = B0 * ny - xi * nx * std::sqrt(4 * M_PI * rho0) * std::cos(phase);
            double Bz_0 = - xi * std::sqrt(4 * M_PI * rho0) * std::sin(phase);
            if(i < innerElemCount) {
                initElemUs[elem.ind] = state_from_primitive_vars2D(rho0, u_0, v_0, w_0, p0, Bx_0, By_0, Bz_0, gam_hcr);
            }
            else{
                initGhostElemUs[elem.ind-innerElemCount] = state_from_primitive_vars2D(rho0, u_0, v_0, w_0, p0, Bx_0, By_0, Bz_0, gam_hcr);
            }
        }
        NeighbourService ns = geometryWorld.getNeighbourService();
        for(const auto& [boundary, ghost]: ns.boundaryToGhostElements){
            int ghostInd = ghost - innerElemCount;
            initGhostElemUs[ghostInd] =  initElemUs[boundary];
        }

        initGhostBNs.resize(geometryWorld.ghostElemCount*2, 0.0);
        for(int i = 0; i < edgp.edgeCount; ++i){
            Edge edge = edgp.edges[i];
            double x = edge.midPoint[0];
            double y = edge.midPoint[1];
            double phase = 2.0 * M_PI / ny * (nx * x + ny * y);
            double u_0 = v0 * nx - xi * ny * std::cos(phase);
            double v_0 = v0 * ny + xi * nx * std::cos(phase);
            double w_0 = xi * std::sin(phase);
            double Bx_0 = B0 * nx + xi * ny * std::sqrt(4 * M_PI * rho0) * std::cos(phase);
            double By_0 = B0 * ny - xi * nx * std::sqrt(4 * M_PI * rho0) * std::cos(phase);
            double Bz_0 = - xi * std::sqrt(4 * M_PI * rho0) * std::sin(phase);
            double Bn = Bx_0 * edge.normalVector[0] + By_0 * edge.normalVector[1];
            initBns[edge.ind] = Bn;
            if(i < innerEdgeCount) {
                initBns[edge.ind] = Bn;
            }
            else{
                initGhostBNs[edge.ind - innerEdgeCount] = Bn;
            }
            initEdgeUs[edge.ind] = state_from_primitive_vars2D(rho0, u_0, v_0, w_0, p0, Bx_0, By_0, Bz_0, gam_hcr);
        }
    }
    else if(task_type == 3){
        std::cout << "SOLVING TASKTYPE 3 (boundary and ghost test)" << std::endl;
        double rho0 = 1.0;
        double p0 = 1.0;
        double v0 = 0.0;
        double B0 = 1.0/ (std::sqrt(4.0 * M_PI));
        double xi = 0.2;
        gam_hcr = 5.0/3.0;
        cflNum = 0.4;
        periodicBoundaries = true;
        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        initBns.resize(edgp.edgeCount, 0.0);
        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < ep.elCount; ++i) {
            Element elem = ep.elements[i];
            std::vector<double> centroid = elem.centroid2D;
            double x = centroid[0];
            double y = centroid[1];
            double u_0 = v0;
            double v_0 = v0;
            double w_0 = 0.0;
            double Bx_0 = B0;
            double By_0 = B0 ;
            double Bz_0 = B0;
            if(i < innerElemCount) {
                if(elem.is_boundary){
                    initElemUs[elem.ind] = state_from_primitive_vars(1000.0, u_0, v_0, w_0, p0, Bx_0, By_0, Bz_0, gam_hcr);
                }else {
                    initElemUs[elem.ind] = state_from_primitive_vars(rho0, u_0, v_0, w_0, p0, Bx_0, By_0, Bz_0, gam_hcr);
                }
            }
            else{
                initGhostElemUs[elem.ind-innerElemCount] = state_from_primitive_vars(1000.0, u_0, v_0, w_0, p0, Bx_0, By_0, Bz_0, gam_hcr);
            }
        }
        NeighbourService ns = geometryWorld.getNeighbourService();
        for(const auto& [top, bot]: ns.boundaryElemTopToBottom){
            initElemUs[top][0] = 5000.0;
            initElemUs[bot][0] = -5000.0;
            int ghostIndTop = ns.boundaryToGhostElements[top] - innerElemCount;
            int ghostIndBot = ns.boundaryToGhostElements[bot] - innerElemCount;
            initGhostElemUs[ghostIndBot] =  initElemUs[bot];
            initGhostElemUs[ghostIndTop] = initElemUs[top];
        }
        for(const auto& [left, right]: ns.boundaryElemLeftToRight){
            initElemUs[left][0] = 2000.0;
            initElemUs[right][0] = -2000.0;
            int ghostIndLeft = ns.boundaryToGhostElements[left] - innerElemCount;
            int ghostIndRight = ns.boundaryToGhostElements[right] - innerElemCount;
            initGhostElemUs[ghostIndRight] =  initElemUs[right];
            initGhostElemUs[ghostIndLeft] = initElemUs[left];
        }

        for(const auto& [top, bot]: ns.boundaryElemTopToBottom){
            int ghostIndTop = ns.boundaryToGhostElements[top] - innerElemCount;
            int ghostIndBot = ns.boundaryToGhostElements[bot] - innerElemCount;
            initGhostElemUs[ghostIndTop] =  initElemUs[bot];
            initGhostElemUs[ghostIndBot] = initElemUs[top];
        }
        for(const auto& [left, right]: ns.boundaryElemLeftToRight){
            int ghostIndLeft = ns.boundaryToGhostElements[left] - innerElemCount;
            int ghostIndRight = ns.boundaryToGhostElements[right] - innerElemCount;
            initGhostElemUs[ghostIndLeft] =  initElemUs[right];
            initGhostElemUs[ghostIndRight] = initElemUs[left];
        }

        initGhostBNs.resize(geometryWorld.ghostElemCount*2, 0.0);
        for(int i = 0; i < edgp.edgeCount; ++i){
            Edge edge = edgp.edges[i];
            double u_0 = v0;
            double v_0 = v0;
            double w_0 = 0.0;
            double Bx_0 = B0;
            double By_0 = B0 ;
            double Bz_0 = B0;
            double Bn = Bx_0 * edge.normalVector[0] + By_0 * edge.normalVector[1];
            initBns[edge.ind] = Bn;
            if(i < innerEdgeCount) {
                initBns[edge.ind] = Bn;
            }
            else{
                initGhostBNs[edge.ind - innerEdgeCount] = Bn;
            }
            initEdgeUs[edge.ind] = state_from_primitive_vars(rho0, u_0, v_0, w_0, p0, Bx_0, By_0, Bz_0, gam_hcr);
        }
    }
    else if(task_type == 4){
        std::cout << "SOLVING TASKTYPE 4: rotating cylinder" << std::endl;
        double r0 = 0.1;
        double r1 = 0.115;
        double v0 = 2.0;
        double p = 1.0;
        double Bx = 5.0/std::sqrt(4*M_PI);
        gam_hcr = 1.4;
        finalTime = 0.15;
        cflNum = 0.5;
        periodicBoundaries = false;
        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        initBns.resize(edgp.edgeCount, 0.0);
        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = ep.elements[i];
            std::vector<double> centroid = elem.centroid2D;
            double x = centroid[0];
            double y = centroid[1];
            double r = std::sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5));
            if(r < r0){
                initElemUs[elem.ind] = state_from_primitive_vars2D(10.0, -v0*(y-0.5)/r0, v0*(x-0.5)/r0, 0.0, p, Bx, 0.0, 0.0, gam_hcr);
            }else if(r < r1){
                double f = (r1 - r) / (r1 - r0);
                initElemUs[elem.ind] = state_from_primitive_vars2D(1.0 + 9 * f, -f*v0*(y-0.5)/r0, f*v0*(x-0.5)/r0, 0.0, p, Bx, 0.0, 0.0, gam_hcr);
            }
            else{
                initElemUs[elem.ind] = state_from_primitive_vars2D(1.0, 0.0, 0.0, 0.0, p, Bx, 0.0, 0.0, gam_hcr);
            }
        }
        NeighbourService ns = geometryWorld.getNeighbourService();
        for(const auto& [boundary, ghost]: ns.boundaryToGhostElements){
            int ghostInd = ghost - innerElemCount;
            initGhostElemUs[ghostInd] =  initElemUs[boundary];
        }
        initGhostBNs.resize(geometryWorld.ghostElemCount*2, 0.0);
        for(int i = 0; i < edgp.edgeCount; ++i){
            Edge edge = edgp.edges[i];
            double Bn = Bx * edge.normalVector[0] +  0.0 * edge.normalVector[1];
            initBns[edge.ind] = Bn;
            if(i < innerEdgeCount) {
                initBns[edge.ind] = Bn;
            }
            else{
                initGhostBNs[edge.ind - innerEdgeCount] = Bn;
            }
        }
    }
    else if(task_type == 5){
        std::cout << "SOLVING TASKTYPE 5: Orszang-Tang vortex" << std::endl;
        double rho = 25.0 / (36.0 * M_PI);
        double p = 5.0 / (12.0 * M_PI);
        gam_hcr = 5.0/3.0;
        finalTime = 0.5;
        cflNum = 0.5;
        periodicBoundaries = true;
        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        initBns.resize(edgp.edgeCount, 0.0);
        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = ep.elements[i];
            std::vector<double> centroid = elem.centroid2D;
            double x = centroid[0];
            double y = centroid[1];
            double u = (-1) * std::sin(2.0 * M_PI * y) ;
            double v = std::sin(2.0 * M_PI * x);
            double w = 0.0;
            double Bx = (-1) * std::sin(2.0 * M_PI  * y)/ (std::sqrt(4.0 * M_PI));
            double By = std::sin(4.0 * M_PI * x)/ (std::sqrt(4.0 * M_PI));
            double Bz = 0.0;
            initElemUs[elem.ind] = state_from_primitive_vars2D(rho, u, v, w, p, Bx, By, Bz, gam_hcr);
        }
        NeighbourService ns = geometryWorld.getNeighbourService();
        for(const auto& [boundary, ghost]: ns.boundaryToGhostElements){
            int ghostInd = ghost - innerElemCount;
            initGhostElemUs[ghostInd] =  initElemUs[boundary];
        }
        initGhostBNs.resize(geometryWorld.ghostElemCount*2, 0.0);
        for(int i = 0; i < edgp.edgeCount; ++i){
            Edge edge = edgp.edges[i];
            double x = edge.midPoint[0];
            double y = edge.midPoint[1];
            double Bx = (-1) * std::sin(2.0 * M_PI  * y)/ (std::sqrt(4.0 * M_PI));
            double By = std::sin(4.0 * M_PI * x)/ (std::sqrt(4.0 * M_PI));
            double Bn = Bx * edge.normalVector[0] +  By * edge.normalVector[1];
            initBns[edge.ind] = Bn;
            if(i < innerEdgeCount) {
                initBns[edge.ind] = Bn;
            }
            else{
                initGhostBNs[edge.ind - innerEdgeCount] = Bn;
            }
        }
    } else if(task_type == 6) {
        cflNum = 0.4;
        freeFLowBoundaries = true;
        periodicBoundaries = false;
        std::cout << "SOLVING TASKTYPE 1 (BRIO-WU TEST)" << std::endl;
                                    /*rho   u   v    w     p    Bx   By   Bz gam_hcr*/
        std::vector<double> BrioWu_L1{1.0, 0.0, 0.0, 0.0, 1.0,  1.0, 0.75, 0.0, gam_hcr};
        std::vector<double> BrioWu_R1{0.125, 0.0, 0.0, 0.0, 0.1, -1.0, 0.75, 0.0, gam_hcr};

        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = ep.elements[i];
            std::vector<double> centroid = elem.centroid2D;
            int elInd = elem.ind;
            if (centroid[1] < 0.5) {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_L1);
            } else {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_R1);
            }
        }

        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for(int i = innerElemCount; i < ep.elCount; ++i ){
            Element elem = ep.elements[i];
            std::vector<double> centroid = elem.centroid2D;
            int elInd = elem.ind;
            int ghostInd = elInd - innerElemCount;
            if (centroid[1] < 0.5) {
                initGhostElemUs[ghostInd] = state_from_primitive_vars(BrioWu_L1);
            } else {
                initGhostElemUs[ghostInd] = state_from_primitive_vars(BrioWu_R1);
            }
        }

        initBns.resize(edgp.edgeCount, 0.0);
        initGhostBNs.resize(geometryWorld.ghostElemCount*2, 0.0);
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < edgp.edgeCount; ++i) {
            Edge edge = edgp.edges[i];
            std::vector<double> midP = edge.midPoint;
            int edgeInd = edge.ind;
            std::vector<double> normal = edge.normalVector;
            std::vector<double> state;
            if (midP[1] < 0.5) {
                state = state_from_primitive_vars(BrioWu_L1);
            } else {
                state = state_from_primitive_vars(BrioWu_R1);
            }
            double Bn = state[5] * normal[0] + state[6] * normal[1];
            initBns[edgeInd] = Bn;
            if(i < innerEdgeCount) {
                initBns[edgeInd] = Bn;
            }else{
                int ghostInd = edgeInd - innerEdgeCount;
                initGhostBNs[ghostInd] = Bn;
            }
            initEdgeUs[edgeInd] = state;
        }
    }
}

void MHDSolver2D::runSolver() {
    omp_set_num_threads(1);
    // service
    EdgePool edgePool = geometryWorld.getEdgePool();
    ElementPool elPool = geometryWorld.getElementPool();
    NodePool nodePool = geometryWorld.getNodePool();
    NeighbourService ns = geometryWorld.getNeighbourService();

    // инициализируем состояния
    setInitElemUs();
    elemUs = initElemUs;
    edgeUs = initEdgeUs;
    bNs = initBns;
    ghostElemUs = initGhostElemUs;
    ghostBNs = initGhostBNs;

    // сделать старые дубликаты состояний (предыдущие состояния) чтобы в новые записывать расчёты
    std::vector<std::vector<double>> elemUs_prev(initElemUs.size(), std::vector<double>(8, 0.0));
    elemUs_prev = initElemUs;
    std::vector<std::vector<double>> edgeUs_prev(initEdgeUs.size(), std::vector<double>(8, 0.0));
    edgeUs_prev = initEdgeUs;
    std::vector<std::vector<double>> ghostElemUs_prev(initGhostElemUs.size(), std::vector<double>(8, 0.0));
    ghostElemUs_prev = initGhostElemUs;

    // проверка на нормировку нормалей
    for(const auto edge: edgePool.edges){
        double norm = std::sqrt(edge.normalVector[0]*edge.normalVector[0] + edge.normalVector[1]*edge.normalVector[1]);
        if(std::abs(1.0-norm) > 1e-15){
            std::cout << "bad normal! |1 - norm| = " << std::abs(1.0-norm) << " edge.is_ghost = " << edge.is_ghost <<std::endl;
        }
    }

    // дивергенция магнитного поля
    double divergence = computeDivergence();
    std::cout << "Init divergence = " << divergence << std::endl;

    double h = edgePool.minEdgeLen;
    std::cout << "Min h = " << h << std::endl;

    double currentTime = startTime;
    bool foundNan = false; // флаг для поиска NaN-значений
    int iterations = 0; // число текущих итераций

    // соновной цикл по времени
    while(currentTime < finalTime) {
        if(iterations >= MAX_ITERATIONS){
            std::cout << "iterations limit!" << std::endl;
            break;
        }

        elemUs_prev = elemUs;

        //ghostElemUs_prev.swap(ghostElemUs);

        tau = std::max(min_tau, tau_from_cfl2D(cflNum, h,  elemUs, gam_hcr));

        currentTime += tau;

        // подбор времени (если "перепрыгнули" за финальное время)
        if (currentTime > finalTime) {
            tau -= (currentTime - finalTime);
            currentTime = finalTime;
        }

        if(debugDivergence){
            divergence = computeDivergence();
            std::cout << "Divergence = " << divergence << std::endl;
        }

        // запись в файл временного результата
        if (iterations % iterationsPerFrame == 0) {
            std::cout << std::setprecision(10) << "t = " << currentTime << std::endl;
            writeVTU("OutputData/tmpres_" + std::to_string(iterations) + ".vtu", ghostOutput);
        }


        // вычисляем потоки, проходящие через каждое ребро
        // инициализируем вектор потоков через рёбра // MHD (HLLD) fluxes (from one element to another "<| -> |>")
        std::vector<std::vector<double>> fluxes(edgePool.edges.size(), std::vector<double>(8, 0.0));
        std::vector<std::vector<double>> unrotated_fluxes(edgePool.edges.size(), std::vector<double>(8, 0.0));
//#pragma parallel for
        for (const auto &edge: edgePool.edges) {
            Element neighbour1 = elPool.elements[edge.neighbourInd1];
            if(neighbour1.is_boundary && edge.neighbourInd2 == -1){
                std::vector<double> U1 = rotateStateFromAxisToNormal(elemUs_prev[neighbour1.ind], edge.normalVector);
                int ghostInd = ns.boundaryToGhostElements[neighbour1.ind] - innerElemCount;
                std::vector<double> U2 = rotateStateFromAxisToNormal(ghostElemUs[ghostInd], edge.normalVector);
                fluxes[edge.ind] = HLLD_flux(U1, U2, gam_hcr);
                unrotated_fluxes[edge.ind] = fluxes[edge.ind];
                fluxes[edge.ind] = rotateStateFromNormalToAxisX(fluxes[edge.ind], edge.normalVector);
            }
            else if(neighbour1.is_ghost && edge.neighbourInd2 == -1){
                //continue; //TODO
                std::vector<double> U1 = rotateStateFromAxisToNormal(ghostElemUs[neighbour1.ind - innerElemCount], edge.normalVector);
                fluxes[edge.ind] = HLLD_flux(U1, U1, gam_hcr);
                unrotated_fluxes[edge.ind] = HLLD_flux(U1, U1, gam_hcr);
                fluxes[edge.ind] = rotateStateFromNormalToAxisX(fluxes[edge.ind], edge.normalVector);
            }
            else{
                std::vector<double> U1 = rotateStateFromAxisToNormal(elemUs_prev[neighbour1.ind], edge.normalVector);
                Element neighbour2 = elPool.elements[edge.neighbourInd2];
                std::vector<double> U2 = rotateStateFromAxisToNormal(elemUs_prev[neighbour2.ind], edge.normalVector);
                fluxes[edge.ind] = HLLD_flux(U1, U2, gam_hcr);
                unrotated_fluxes[edge.ind] = fluxes[edge.ind];
                fluxes[edge.ind] = rotateStateFromNormalToAxisX(fluxes[edge.ind], edge.normalVector);
            }
        }

        if(freeFLowBoundaries) {
            for (const auto &[top, bot]: ns.boundaryElemTopToBottom) {
                Element elGhostTop = elPool.elements[ns.boundaryToGhostElements[top]];
                for (const auto &edgeInd: elGhostTop.edgeIndexes) {
                    Edge edge = edgePool.edges[edgeInd];
                    if (edge.is_ghost) {
                        int corspEdge = ns.edgeToGhostEdges[edge.ind];
                        unrotated_fluxes[edge.ind] = (-1.0) * unrotated_fluxes[corspEdge];
                        fluxes[edge.ind] = rotateStateFromAxisToNormal(unrotated_fluxes[edge.ind], edge.normalVector);
                    }
                }
                Element elGhostBot = elPool.elements[ns.boundaryToGhostElements[bot]];
                for (const auto &edgeInd: elGhostBot.edgeIndexes) {
                    Edge edge = edgePool.edges[edgeInd];
                    if (edge.is_ghost) {
                        int corspEdge = ns.edgeToGhostEdges[edge.ind];
                        unrotated_fluxes[edge.ind] = (-1.0) * unrotated_fluxes[corspEdge];
                        fluxes[edge.ind] = rotateStateFromAxisToNormal(unrotated_fluxes[edge.ind], edge.normalVector);
                    }
                }
            }
        }
        // по явной схеме обновляем газовые величины
        //#pragma omp parallel for
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = elPool.elements[i];
            std::vector<double> fluxSum(8, 0.0);
            if(elem.edgeIndexes.size() != 3){
                std::cout << "Bad edge vector size != 3" << std::endl;
            }
            if(!elem.is_boundary || freeFLowBoundaries == false) {
                for (int edgeIndex: elem.edgeIndexes) {
                    Edge edge_j = edgePool.edges[edgeIndex];
                    if (std::abs(edge_j.length) < 1e-16) {
                        std::cout << "BAD EDGE LENGTH = 0!! " << edge_j.length << std::endl;
                        Node node1 = nodePool.getNode(edge_j.nodeInd1);
                        Node node2 = nodePool.getNode(edge_j.nodeInd2);
                        std::cout << "edge: " << edge_j.ind << std::endl;
                        std::cout << "node1 # " << edge_j.nodeInd1 << " { " << node1.x << ", " << node1.y << " } \n";
                        std::cout << "node2 # " << edge_j.nodeInd2 << " { " << node2.x << ", " << node2.y << " } \n";
                        Element elem1 = elPool.elements[edge_j.neighbourInd1];
                        std::cout << "neigel1 # " << elem1.ind << " isBoundary = " << elem1.is_boundary
                                  << " , isGhost = " << elem1.is_ghost << std::endl;
                        std::cout << "neigel2 # " << edge_j.neighbourInd2 << std::endl;
                        std::cin.get();

                    }
                    if (edge_j.neighbourInd1 == i) {
                        fluxSum = fluxSum + edge_j.length * fluxes[edgeIndex];
                    } else if (edge_j.neighbourInd2 == i) {
                        fluxSum = fluxSum - edge_j.length * fluxes[edgeIndex];
                    } else {
                        //std::cerr << "No matching edge..."  << std::endl;
                    }
                }
                elemUs[i] = elemUs_prev[i] - tau / elem.area * fluxSum;
            }
            else{
                int ghostInd = ns.boundaryToGhostElements[i];
                int ghostIndForState = ns.boundaryToGhostElements[i] - innerElemCount;
                Element ghostEl = elPool.elements[ghostInd];
                std::vector<int> commonEdges = findCommonElements(ghostEl.edgeIndexes, elem.edgeIndexes);
                if(commonEdges.size() != 1){
                    std::cout << "NO COMMON EDGES!" << std::endl;
                    std::cout << elem.edgeIndexes[0] << ",  "<< elem.edgeIndexes[1] << ",  "<< elem.edgeIndexes[2] << "; vs/s ghost: " << ghostEl.edgeIndexes[0] << ",  " << ghostEl.edgeIndexes[1] << ",  " << ghostEl.edgeIndexes[2] << std::endl;
                }
                int commonEdgeInd = commonEdges[0];
                double area = elem.area * 2.0;

                for(const auto& edgeInd: elem.edgeIndexes){
                    if(edgeInd == commonEdgeInd){
                        continue;
                    }
                    Edge edge_j = edgePool.edges[edgeInd];
                    if (std::abs(edge_j.length) < 1e-16) {
                        std::cout << "BAD EDGE LENGTH = 0!! " << edge_j.length << std::endl;
                        Node node1 = nodePool.getNode(edge_j.nodeInd1);
                        Node node2 = nodePool.getNode(edge_j.nodeInd2);
                        std::cout << "edge: " << edge_j.ind << std::endl;
                        std::cout << "node1 # " << edge_j.nodeInd1 << " { " << node1.x << ", " << node1.y << " } \n";
                        std::cout << "node2 # " << edge_j.nodeInd2 << " { " << node2.x << ", " << node2.y << " } \n";
                        Element elem1 = elPool.elements[edge_j.neighbourInd1];
                        std::cout << "neigel1 # " << elem1.ind << " isBoundary = " << elem1.is_boundary
                                  << " , isGhost = " << elem1.is_ghost << std::endl;
                        std::cout << "neigel2 # " << edge_j.neighbourInd2 << std::endl;
                        std::cin.get();
                    }
                    if (edge_j.neighbourInd1 == i) {
                        fluxSum = fluxSum + edge_j.length * fluxes[edgeInd];
                    } else if (edge_j.neighbourInd2 == i) {
                        fluxSum = fluxSum - edge_j.length * fluxes[edgeInd];
                    } else {
                        std::cerr << "No matching edge..."  << std::endl;
                    }
                }
                for(const auto& ghostEdgeInd: ghostEl.edgeIndexes){
                    if(ghostEdgeInd == commonEdgeInd){
                        continue;
                    }
                    Edge edge_j = edgePool.edges[ghostEdgeInd];
                    if (std::abs(edge_j.length) < 1e-16) {
                        std::cout << "BAD EDGE LENGTH = 0!! " << edge_j.length << std::endl;
                        Node node1 = nodePool.getNode(edge_j.nodeInd1);
                        Node node2 = nodePool.getNode(edge_j.nodeInd2);
                        std::cout << "edge: " << edge_j.ind << std::endl;
                        std::cout << "node1 # " << edge_j.nodeInd1 << " { " << node1.x << ", " << node1.y << " } \n";
                        std::cout << "node2 # " << edge_j.nodeInd2 << " { " << node2.x << ", " << node2.y << " } \n";
                        Element elem1 = elPool.elements[edge_j.neighbourInd1];
                        std::cout << "neigel1 # " << elem1.ind << " isBoundary = " << elem1.is_boundary
                                  << " , isGhost = " << elem1.is_ghost << std::endl;
                        std::cout << "neigel2 # " << edge_j.neighbourInd2 << std::endl;
                        std::cin.get();
                    }
                    if (edge_j.neighbourInd1 == ghostEl.ind) {
                        fluxSum = fluxSum + edge_j.length * fluxes[ghostEdgeInd];
                    } else if (edge_j.neighbourInd2 == ghostEl.ind) {
                        fluxSum = fluxSum - edge_j.length * fluxes[ghostEdgeInd];
                    } else {
                        std::cerr << "No matching edge..."  << std::endl;
                    }
                }
                elemUs[i] = elemUs_prev[i] - tau / area * fluxSum;
            }
        }


        // г.у (фикт ячейки) поставить перед выч-м потоков
        if(periodicBoundaries){
            // периодическое г.у.
            for(const auto& [top, bot]: ns.boundaryElemTopToBottom){
                int ghostIndTop = ns.boundaryToGhostElements[top] - innerElemCount;
                int ghostIndBot = ns.boundaryToGhostElements[bot] - innerElemCount;
                ghostElemUs[ghostIndTop] = elemUs[bot];
                ghostElemUs[ghostIndBot] = elemUs[top];
            }
            for(const auto& [left, right]: ns.boundaryElemLeftToRight){
                int ghostIndLeft = ns.boundaryToGhostElements[left] - innerElemCount;
                int ghostIndRight = ns.boundaryToGhostElements[right] - innerElemCount;
                ghostElemUs[ghostIndLeft] =  elemUs[right];
                ghostElemUs[ghostIndRight] = elemUs[left];
            }
        }
        else{
            // копируем значения в соответствующие фантомные ячейки (условие free flow)
            for(const auto& [boundary, ghost] : ns.boundaryToGhostElements){
                int ghostInd = ghost - innerElemCount;
                ghostElemUs[ghostInd] = elemUs[boundary];
            }
        }

        for(auto& edge: edgePool.edges){
            Element neighbour1 = elPool.elements[edge.neighbourInd1];
            if(neighbour1.is_ghost){
                edgeUs[edge.ind] = ghostElemUs[neighbour1.ind - innerElemCount];
            }else if(neighbour1.is_boundary){
                std::vector<double> U1 = elemUs[neighbour1.ind];
                int ghostInd = ns.boundaryToGhostElements[neighbour1.ind] - innerElemCount;
                std::vector<double> U2 = ghostElemUs[ghostInd];
                edgeUs[edge.ind] = (1.0/(std::sqrt(U1[0]) + std::sqrt(U2[0]))) * (std::sqrt(U1[0]) * U1 + std::sqrt(U2[0]) * U2);
            }else{
                Element neighbour2 = elPool.elements[edge.neighbourInd2];
                std::vector<double> U1 = elemUs[neighbour1.ind];
                std::vector<double> U2 = elemUs[neighbour2.ind];
                edgeUs[edge.ind] = (1.0/(std::sqrt(U1[0]) + std::sqrt(U2[0]))) * (std::sqrt(U1[0]) * U1 + std::sqrt(U2[0]) * U2);
            }
        }

        //!!! добавить г.у.
        //корректируем магнитные величины
        //находим узловые значения нужных магнитных разностей //(v x B)z в узлах
        std::vector<double> nodeMagDiffs(nodePool.nodeCount, 0.0); //(v x B)z в узлах
//#pragma omp parallel for
        for (const auto &node: nodePool.nodes) {
            int tmp_count = 0;
            for (const auto &neighbourEdgeInd: ns.getEdgeNeighborsOfNode(node.ind)) {
                // добавить проверку на совпадения направлений в-ра скорости и нарпавляющего ребра
                Edge neighbourEdge = edgePool.edges[neighbourEdgeInd];
                if(elPool.elements[neighbourEdge.neighbourInd1].is_ghost && neighbourEdge.neighbourInd2 == -1){
                    //std::cout << "SEE A GHOST!" << std::endl;
                    //continue;
                }
                Node node1 = nodePool.getNode(neighbourEdge.nodeInd1);
                Node node2 = nodePool.getNode(neighbourEdge.nodeInd2);
                int inflowOrientation = 1;
                double u = edgeUs[neighbourEdgeInd][1]/edgeUs[neighbourEdgeInd][0];
                double v = edgeUs[neighbourEdgeInd][2]/edgeUs[neighbourEdgeInd][0];
                if(neighbourEdge.nodeInd1 == node.ind){
                    //inflowOrientation = sgn((u*(node1.x - node2.x) + v*(node1.y - node2.y)));
                    //std::cout << "scalar mult = " << (u*(node1.x - node2.x) + v*(node1.y - node2.y)) <<" , orient = " << inflowOrientation <<std::endl;
                    //std::cin.get();
                }
                else if(neighbourEdge.nodeInd2 == node.ind){
                    //inflowOrientation =  sgn((u*(node2.x - node1.x) + v*(node2.y - node1.y)));
                    //std::cout << "scalar mult = " <<  (u*(node2.x - node1.x) + v*(node2.y - node1.y))  <<" , orient = " << inflowOrientation <<std::endl;
                    //std::cin.get();
                }
                if(inflowOrientation > 0) {
                    nodeMagDiffs[node.ind] += unrotated_fluxes[neighbourEdgeInd][6];
                    ++tmp_count;
                }
            }
            if (tmp_count) {
                nodeMagDiffs[node.ind] /= tmp_count;
            }
        }

        //находим новое значение Bn в ребре
        std::vector<double> bNs_prev(bNs);
        for (int i = 0; i < edgePool.edgeCount; ++i) {
            Edge edge = edgePool.edges[i];
            bNs[edge.ind] = bNs_prev[edge.ind] + (tau / edge.length) * (nodeMagDiffs[edge.nodeInd2] -
                                                        nodeMagDiffs[edge.nodeInd1]);
        }

        //сносим Bn в центр элемента
//#pragma omp parallel for
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = elPool.elements[i];
            std::vector<double> centroid = getElementCentroid2D(elem, nodePool);
            if(elem.edgeIndexes.empty()){
                std::cerr << "Empty element! (no edges)" << std::endl;
            }
            double temp_sum_Bx = 0.0;
            double temp_sum_By = 0.0;
            for (const auto &edgeInd: elem.edgeIndexes) {
                Edge edge = edgePool.edges[edgeInd];
                std::vector<double> cTe = edge.midPoint - centroid;
                double scmult = cTe[0]*edge.normalVector[0] + cTe[1]*edge.normalVector[1];
                if (scmult > 0.0/*edge.neighbourInd1 == elem.ind*/) {
                    // у первого соседа в эдже заданы ноды в порядке полодительного обхода и нормаль тоже
                    const auto nodeInElemInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(),
                                                         edge.nodeInd1);
                    int node_before_ind =
                            nodeInElemInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim - 1] : *(
                                    nodeInElemInd - 1);
                    //std::cout << "edge of nodes = {"<< edge.nodeInd1 << " , " << edge.nodeInd2 << " } current ind = " << edge.nodeInd1 <<" Elem's node indexes: " << elem.nodeIndexes[0] << " "<< elem.nodeIndexes[1] << " "<< elem.nodeIndexes[2] << ", node_beforeInd = " << node_before_ind << std::endl;
                    //std::cin.get();
                    Node node_before = nodePool.getNode(node_before_ind);
                    if(edgeInd < innerEdgeCount) {
                        temp_sum_Bx += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[0] - node_before.x);
                        temp_sum_By += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[1] - node_before.y);
                    }
                    else{
                        continue;
                        temp_sum_Bx += ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid[0] - node_before.x);
                        temp_sum_By += ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid[1] - node_before.y);
                    }
                } else if(scmult < 0.0/*edge.neighbourInd2 != -1*/){
                    // а вот для второго нужно умножать на -1 и в обратном порядке
                    const auto nodeInElemInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(),
                                                         edge.nodeInd2);
                    int node_before_ind =
                            nodeInElemInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim - 1] : *(
                                    nodeInElemInd - 1);
                    Node node_before = nodePool.getNode(node_before_ind);
                    if(edgeInd < innerEdgeCount) {
                        temp_sum_Bx -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[0] - node_before.x);
                        temp_sum_By -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[1] - node_before.y);
                    }
                    else{
                        continue;
                        temp_sum_Bx -= ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid[0] - node_before.x);
                        temp_sum_By -= ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid[1] - node_before.y);
                    }
                }
            }
            if(elem.ind < innerElemCount) {
                elemUs[elem.ind][5] = temp_sum_Bx;
                elemUs[elem.ind][6] = temp_sum_By;
            }
            else{
               // ghostElemUs_prev[elem.ind - innerElemCount][5] = temp_sum_Bx;
               // ghostElemUs_prev[elem.ind - innerElemCount][6] = temp_sum_By;
            }
        }

        // г.у (фикт ячейки) поставить перед выч-м потоков
        if(periodicBoundaries){
            // периодическое г.у.
            for(const auto& [top, bot]: ns.boundaryElemTopToBottom){
                int ghostIndTop = ns.boundaryToGhostElements[top] - innerElemCount;
                int ghostIndBot = ns.boundaryToGhostElements[bot] - innerElemCount;
                ghostElemUs[ghostIndTop] =  elemUs[bot];
                ghostElemUs[ghostIndBot] = elemUs[top];
            }
            for(const auto& [left, right]: ns.boundaryElemLeftToRight){
                int ghostIndLeft = ns.boundaryToGhostElements[left] - innerElemCount;
                int ghostIndRight = ns.boundaryToGhostElements[right] - innerElemCount;
                ghostElemUs[ghostIndLeft] =  elemUs[right];
                ghostElemUs[ghostIndRight] = elemUs[left];
            }
        }
        else{
            // копируем значения в соответствующие фантомные ячейки (условие free flow)
            for(const auto& [boundary, ghost] : ns.boundaryToGhostElements){
                int ghostInd = ghost - innerElemCount;
                ghostElemUs[ghostInd] = elemUs[boundary];
            }
        }

//#pragma parallel for
        for(const auto& u: elemUs){
            for(const auto& val: u){
                if(std::isnan(val)){
                    foundNan = true;
                    break;
                }
            }
        }
        if(foundNan){
            std::cout << "Found a nan VALUE!!! Exiting the solver..." << std::endl;
            break;
        }
        ++iterations;
    }
    writeVTU("OutputData/tmpres_" + std::to_string(iterations) + ".vtu", ghostOutput);
    std::cout << "Final time = " << currentTime << "; iterations = "<< iterations<<  std::endl;
    divergence = computeDivergence();
    std::cout << "Final divergence = " << divergence << std::endl;
}


void MHDSolver2D::writeVTU(const std::string& filename, const bool& ghost) {
    const NodePool& np = geometryWorld.getNodePool();
    const ElementPool& ep = geometryWorld.getElementPool();

    const int ghostElemCount = geometryWorld.ghostElemCount;
    const int ghostNodeCount = geometryWorld.ghostNodeCount;

    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << np.nodeCount - ghostNodeCount * (1 - ghost)
         << "\" NumberOfCells=\"" << ep.elCount - ghostElemCount * (1 - ghost) << "\">\n";

    // Write points (nodes)
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& node : np.nodes) {
        if (!node.is_ghost || ghost) {
            file << "          " << node.x << " " << node.y << " " << node.z << "\n";
        }
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // Write cells (elements)
    file << "      <Cells>\n";

    // Connectivity (fix 1-based index issue)
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& element : ep.elements) {
        if (!element.is_ghost || ghost) {
            for (const auto &nodeIndex : element.nodeIndexes) {
                file << "          " << (nodeIndex) << " ";  // Convert to 0-based indexing
            }
            file << "\n";
        }
    }
    file << "        </DataArray>\n";

    // Offsets
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (const auto& element : ep.elements) {
        if (!element.is_ghost || ghost) {
            offset += element.dim;
            file << "          " << offset << "\n";
        }
    }
    file << "        </DataArray>\n";

    // Types
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (const auto& element : ep.elements) {
        if (!element.is_ghost || ghost) {
            file << "          " << ((element.dim == 3) ? 5 : 9) << "\n";  // 5 for Triangle, 9 for Quadrilateral
        }
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // Write solution data (elemUs) - Fix NumberOfComponents to 8
    file << "      <CellData Scalars=\"elemUs\">\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"9\" Name=\"elemUs\" format=\"ascii\">\n";
    for (const auto& U : elemUs) {
        for (const auto& value : U) {
            file << "          " << value << " ";
        }
        file << "          " << gas_p(gam_hcr, U[4], U[0], U[1]/U[0], U[2]/U[0], U[3]/U[0], U[5], U[6], U[7]) << " ";
        file << "\n";
    }
    if (ghost) {
        for (const auto& U : ghostElemUs) {
            for (const auto& value : U) {
                file << "          " << value << " ";
            }
            file << "          " << gas_p(gam_hcr, U[4], U[0], U[1]/U[0], U[2]/U[0], U[3]/U[0], U[5], U[6], U[7]) << " ";
            file << "\n";
        }
    }
    file << "        </DataArray>\n";
    file << "      </CellData>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
}


double MHDSolver2D::computeDivergence() {
    double max_divergence = 0.0;
    ElementPool elPool = geometryWorld.getElementPool();
    NodePool nodePool = geometryWorld.getNodePool();
    EdgePool edgePool = geometryWorld.getEdgePool();
    // 1/S Sum(l*Bn)
    for (int i = 0; i < innerElemCount; ++i) {
        Element elem = elPool.elements[i];
        std::vector<double> centroid = getElementCentroid2D(elem, nodePool);
        if (elem.edgeIndexes.empty()) {
            std::cerr << "Empty element! (no edges)" << std::endl;
        }
        double divergence = 0.0;
        for (const auto &edgeInd: elem.edgeIndexes) {
            Edge edge = edgePool.edges[edgeInd];
            std::vector<double> cTe = edge.midPoint - centroid;
            double scmult = cTe[0]*edge.normalVector[0] + cTe[1]*edge.normalVector[1];
            if (scmult > 0) {
                divergence += edge.length * bNs[edgeInd];
            } else {
                divergence -= edge.length * bNs[edgeInd];
            }
        }
        divergence /= elem.area;
        if(divergence > max_divergence){
            max_divergence = divergence;
        }
    }

    return max_divergence;
}

std::vector<int> findCommonElements(const std::vector<int>& v1, const std::vector<int>& v2) {
    std::unordered_set<int> elements(v1.begin(), v1.end());
    std::vector<int> common;

    for (int num : v2) {
        if (elements.count(num)) {
            common.push_back(num);
        }
    }

    return common;
}


