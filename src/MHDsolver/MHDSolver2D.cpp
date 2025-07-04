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
// n = {cos(j), sin(j)} = {n.pos.x, n.pos.y}

// ({cos(j), sin(j), 0}, {-sin(j), cos(j), 0}, {0,0,1}) --- around OZ clockwise
// rotate 1: from normal to OX:      {v.pos.x * n.pos.x + vy * n.pos.y, - v.pos.x * n.pos.y + v.pos.y * n.pos.x, v.z}
// вращаем {u,v,w} и {Bx, By, Bz}, остальные остаются на месте
std::vector<double> MHDSolver2D::rotateStateFromAxisToNormal(std::vector<double> &U, const Vec2 &n) {
    std::vector<double> res(U);
    res[1] =  U[1]*n.x + U[2]*n.y;
    res[2] = -U[1]*n.y + U[2]*n.x;
    res[5] =  U[5]*n.x + U[6]*n.y;
    res[6] = -U[5]*n.y + U[6]*n.x;

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
// rotate 2: from OX to normal:      {v.pos.x * n.pos.x - vy * n.pos.y,  v.pos.x * n.pos.y + v.pos.y * n.pos.x, v.z}
std::vector<double> MHDSolver2D::rotateStateFromNormalToAxisX(vector<double> &U, const Vec2 &n) {
    std::vector<double> res(U);
    res[1] =  U[1]*n.x - U[2]*n.y;
    res[2] =  U[1]*n.y + U[2]*n.x;
    res[5] =  U[5]*n.x - U[6]*n.y;
    res[6] =  U[5]*n.y + U[6]*n.x;
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
    NodePool np = geometryWorld.getNodePool();

    // Brio-Wu problem
    if(task_type == 1) {
        cflNum = 0.9;
        freeFLowBoundaries = true;
        periodicBoundaries = false;
        std::cout << "SOLVING TASKTYPE 1 (BRIO-WU TEST)" << std::endl;
                                     /*rho  u   v   w   p      Bx   By   Bz gam_hcr*/
        std::vector<double> BrioWu_L1{1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0, gam_hcr};
        std::vector<double> BrioWu_R1{0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0, gam_hcr};

        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = ep.elements[i];
            Vec2 centroid = elem.centroid2D;
            int elInd = elem.ind;
            if (centroid.x < 0.5) {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_L1);
            } else {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_R1);
            }
        }

        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for(int i = innerElemCount; i < ep.elCount; ++i ){
            Element elem = ep.elements[i];
            Vec2 centroid = elem.centroid2D;
            int elInd = elem.ind;
            int ghostInd = elInd - innerElemCount;
            if (centroid.x < 0.5) {
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
            Vec2 midP = edge.midPoint;
            int edgeInd = edge.ind;
            Vec2 normal = edge.normalVector;
            std::vector<double> state;
            if (midP.x < 0.5) {
                state = state_from_primitive_vars(BrioWu_L1);
            } else {
                state = state_from_primitive_vars(BrioWu_R1);
            }
            double Bn = state[5] * normal.x + state[6] * normal.y;
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
        double B0 = 1.0/ (std::sqrt(4.0 * std::numbers::pi_v<double>));
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
            Vec2 centroid = elem.centroid2D;
            double x = centroid.x;
            double y = centroid.y;
            double phase = 2.0 * std::numbers::pi_v<double> / ny * (nx * x + ny * y);
            double u_0 = v0 * nx - xi * ny * std::cos(phase);
            double v_0 = v0 * ny + xi * nx * std::cos(phase);
            double w_0 = xi * std::sin(phase);
            double Bx_0 = B0 * nx + xi * ny * std::sqrt(4 * std::numbers::pi_v<double> * rho0) * std::cos(phase);
            double By_0 = B0 * ny - xi * nx * std::sqrt(4 * std::numbers::pi_v<double> * rho0) * std::cos(phase);
            double Bz_0 = - xi * std::sqrt(4 * std::numbers::pi_v<double> * rho0) * std::sin(phase);
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
            double x = edge.midPoint.x;
            double y = edge.midPoint.y;
            double phase = 2.0 * std::numbers::pi_v<double> / ny * (nx * x + ny * y);
            double u_0 = v0 * nx - xi * ny * std::cos(phase);
            double v_0 = v0 * ny + xi * nx * std::cos(phase);
            double w_0 = xi * std::sin(phase);
            double Bx_0 = B0 * nx + xi * ny * std::sqrt(4 * std::numbers::pi_v<double> * rho0) * std::cos(phase);
            double By_0 = B0 * ny - xi * nx * std::sqrt(4 * std::numbers::pi_v<double> * rho0) * std::cos(phase);
            double Bz_0 = - xi * std::sqrt(4 * std::numbers::pi_v<double> * rho0) * std::sin(phase);
            double Bn = Bx_0 * edge.normalVector.x + By_0 * edge.normalVector.y;
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
        double B0 = 1.0/ (std::sqrt(4.0 * std::numbers::pi_v<double>));
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
            Vec2 centroid = elem.centroid2D;
            double x = centroid.x;
            double y = centroid.y;
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
            double Bn = Bx_0 * edge.normalVector.x + By_0 * edge.normalVector.y;
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
        double Bx = 5.0/std::sqrt(4*std::numbers::pi_v<double>);
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
            Vec2 centroid = elem.centroid2D;
            double x = centroid.x;
            double y = centroid.y;
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
            double Bn = Bx * edge.normalVector.x +  0.0 * edge.normalVector.y;
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
        double rho = 25.0 / (36.0 * std::numbers::pi_v<double>);
        double p = 5.0 / (12.0 * std::numbers::pi_v<double>);
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
            Vec2 centroid = elem.centroid2D;
            double x = centroid.x;
            double y = centroid.y;
            double u = (-1) * std::sin(2.0 * std::numbers::pi_v<double> * y) ;
            double v = std::sin(2.0 * std::numbers::pi_v<double> * x);
            double w = 0.0;
            double Bx = (-1) * std::sin(2.0 * std::numbers::pi_v<double>  * y)/ (std::sqrt(4.0 * std::numbers::pi_v<double>));
            double By = std::sin(4.0 * std::numbers::pi_v<double> * x)/ (std::sqrt(4.0 * std::numbers::pi_v<double>));
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
            double x = edge.midPoint.x;
            double y = edge.midPoint.y;
            double Bx = (-1) * std::sin(2.0 * std::numbers::pi_v<double>  * y)/ (std::sqrt(4.0 * std::numbers::pi_v<double>));
            double By = std::sin(4.0 * std::numbers::pi_v<double> * x)/ (std::sqrt(4.0 * std::numbers::pi_v<double>));
            double Bn = Bx * edge.normalVector.x +  By * edge.normalVector.y;
            initBns[edge.ind] = Bn;
            if(i < innerEdgeCount) {
                initBns[edge.ind] = Bn;
            }
            else{
                initGhostBNs[edge.ind - innerEdgeCount] = Bn;
            }
        }
    }
    else if(task_type == 6) {
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
            Vec2 centroid = elem.centroid2D;
            int elInd = elem.ind;
            if (centroid.y < 0.5) {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_L1);
            } else {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_R1);
            }
        }

        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for(int i = innerElemCount; i < ep.elCount; ++i ){
            Element elem = ep.elements[i];
            Vec2 centroid = elem.centroid2D;
            int elInd = elem.ind;
            int ghostInd = elInd - innerElemCount;
            if (centroid.y < 0.5) {
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
            Vec2 midP = edge.midPoint;
            int edgeInd = edge.ind;
            Vec2 normal = edge.normalVector;
            std::vector<double> state;
            if (midP.y < 0.5) {
                state = state_from_primitive_vars(BrioWu_L1);
            } else {
                state = state_from_primitive_vars(BrioWu_R1);
            }
            double Bn = state[5] * normal.x + state[6] * normal.y;
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
    else if(task_type == 7) {
        cflNum = 0.4;
        freeFLowBoundaries = true;
        periodicBoundaries = false;
        finalTime = 0.2;
        std::cout << "SOLVING TASKTYPE 7 (SOD'S PROBLEM)" << std::endl;
                                     /*rho   u   v    w     p    Bx   By   Bz gam_hcr*/
        std::vector<double> BrioWu_L1{1.0, 0.0, 0.0, 0.0, 1.0,  0.0, 0.0, 0.0, gam_hcr};
        std::vector<double> BrioWu_R1{0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, gam_hcr};

        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = ep.elements[i];
            Vec2 centroid = elem.centroid2D;
            int elInd = elem.ind;
            if (centroid.x < 0.5) {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_L1);
            } else {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_R1);
            }
        }

        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for(int i = innerElemCount; i < ep.elCount; ++i ){
            Element elem = ep.elements[i];
            Vec2 centroid = elem.centroid2D;
            int elInd = elem.ind;
            int ghostInd = elInd - innerElemCount;
            if (centroid.x < 0.5) {
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
            Vec2 midP = edge.midPoint;
            int edgeInd = edge.ind;
            Vec2  normal = edge.normalVector;
            std::vector<double> state;
            if (midP.x < 0.5) {
                state = state_from_primitive_vars(BrioWu_L1);
            } else {
                state = state_from_primitive_vars(BrioWu_R1);
            }
            double Bn = state[5] * normal.x + state[6] * normal.y;
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
    /*тесты Athena*/
    else if(task_type == 9){
        std::cout << "SOLVING TASKTYPE 9: The Magnetic Field Loop Test" << std::endl;
        double rho = 1.0;
        double p = 1.0;
        gam_hcr = 5.0/3.0;
        finalTime = 2.0;
        cflNum = 0.1;
        periodicBoundaries = true;
        freeFLowBoundaries = false;
        freeFLowBoundaries2 = false;
        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        initBns.resize(edgp.edgeCount, 0.0);
        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        initGhostBNs.resize(geometryWorld.ghostElemCount*2, 0.0);
        double R0 = 0.3;  // Радиус петли
        double A = 1.0e-3;   // Амплитуда
        double xc = 0.0, yc = 0.0;  // Центр петли
        std::vector<double> AzField(np.nodeCount, 0.0);
        for (auto& node : np.nodes) {
            double dx = node.pos.x - xc;
            double dy = node.pos.y - yc;
            double r = std::sqrt(dx*dx + dy*dy);
            AzField[node.ind] = (r <= R0) ? A * (R0 - r) : 0.0;
        }
        // сначала рёбра:
        for(int i = 0; i < edgp.edgeCount; ++i){
            Edge edge = edgp.edges[i];
            Node& node1 = np.nodes[edge.nodeInd1];
            Node& node2 = np.nodes[edge.nodeInd2];
            double Az1 = AzField[node1.ind];
            double Az2 = AzField[node2.ind];
            double Bn = (Az2-Az1)/edge.length;
            initBns[edge.ind] = Bn;
            if(i < innerEdgeCount) {
                initBns[edge.ind] = Bn;
            }
            else{
                initGhostBNs[edge.ind - innerEdgeCount] = Bn;
            }
        }
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = ep.elements[i];
            Vec2 centroid = elem.centroid2D;
            double x = centroid.x;
            double y = centroid.y;
            double u = 2.0;//std::sin(std::numbers::pi_v<double>/3.0);
            double v = 1.0;//std::cos(std::numbers::pi_v<double>/3.0);
            double w = 0.0;
            double Bz = 0.0;
            double Bx = 0.0;
            double By = 0.0;
            double temp_sum_Bx = 0.0;
            double temp_sum_By = 0.0;
            for (const auto &edgeInd: elem.edgeIndexes) {
                Edge edge = edgp.edges[edgeInd];
                if (edge.neighbourInd1 == elem.ind) {
                    // у первого соседа в эдже заданы ноды в порядке положительного обхода и нормаль тоже
                    const auto nodeInElemInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(),
                                                         edge.nodeInd1);
                    int node_before_ind =
                            nodeInElemInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim - 1] : *(
                                    nodeInElemInd - 1);
                    //std::cout << "edge of nodes = {"<< edge.nodeInd1 << " , " << edge.nodeInd2 << " } current ind = " << edge.nodeInd1 <<" Elem's node indexes: " << elem.nodeIndexes[0] << " "<< elem.nodeIndexes[1] << " "<< elem.nodeIndexes[2] << ", node_beforeInd = " << node_before_ind << std::endl;
                    //std::cin.get();
                    Node node_before = np.getNode(node_before_ind);
                    if(edgeInd < innerEdgeCount) {
                        Bx += initBns[edgeInd] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        By += initBns[edgeInd] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
                    }
                    else{
                        continue;
                        Bx += initGhostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        By += initGhostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
                    }
                } else if(edge.neighbourInd2 != -1){
                    // а вот для второго нужно умножать на -1 и в обратном порядке
                    const auto nodeInElemInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(),
                                                         edge.nodeInd2);
                    int node_before_ind =
                            nodeInElemInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim - 1] : *(
                                    nodeInElemInd - 1);
                    Node node_before = np.getNode(node_before_ind);
                    if(edgeInd < innerEdgeCount) {
                        Bx -= initBns[edgeInd] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        By -= initBns[edgeInd] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
                    }
                    else{
                        continue;
                        Bx -= initGhostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        By -= initGhostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
                    }
                }
            }
            initElemUs[elem.ind] = state_from_primitive_vars2D(rho, u, v, w, p, Bx, By, Bz, gam_hcr);
        }
        NeighbourService ns = geometryWorld.getNeighbourService();
        for(const auto& [boundary, ghost]: ns.boundaryToGhostElements){
            int ghostInd = ghost - innerElemCount;
            initGhostElemUs[ghostInd] =  initElemUs[boundary];
        }
    }
    else if(task_type == 8){
        std::cout << "SOLVING TASKTYPE 8: The Circular Polarized Alfven Wave Test " << std::endl;
        double rho = 1.0;
        double p = 0.1;
        gam_hcr = 5.0/3.0;
        finalTime = 1.0;
        cflNum = 0.1;
        double alpha = std::numbers::pi_v<double>/6.0; //The wave propagates along the diagonal of the grid, at an angle θ = tan-1(0.5) ≈ 26.6 degrees with respect to the x-axis
        periodicBoundaries = true;
        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        initBns.resize(edgp.edgeCount, 0.0);
        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = ep.elements[i];
            Vec2 centroid = elem.centroid2D;
            double x = centroid.x;
            double y = centroid.y;
            double x_par = x * std::cos(alpha) + y * std::sin(alpha);
            double v_perp = 0.1 * std::sin(2.0 * std::numbers::pi_v<double> * x_par);
            double v_par = 0.0;
            double B_perp = 0.1 * std::sin(2.0 * std::numbers::pi_v<double> * x_par);
            double B_par = 1.0;
            double u = v_par * std::cos(alpha) - v_perp * std::sin(alpha);
            double v = v_perp * std::cos(alpha) + v_par * std::sin(alpha);
            double w = 0.1 * std::cos(2.0 * std::numbers::pi_v<double> * x_par);
            double Bz = 0.1 * std::cos(2.0 * std::numbers::pi_v<double> * x_par);
            double Bx = B_par * std::cos(alpha) - B_perp * std::sin(alpha);
            double By = B_perp * std::cos(alpha) + B_par * std::sin(alpha);
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
            double x = edge.midPoint.x;
            double y = edge.midPoint.y;
            double x_par = x * std::cos(alpha) + y * std::sin(alpha);
            double B_perp = 0.1 * std::sin(2.0 * std::numbers::pi_v<double> * x_par);
            double B_par = 1.0;
            double Bx = B_par * std::cos(alpha) - B_perp * std::sin(alpha);
            double By = B_perp * std::cos(alpha) + B_par * std::sin(alpha);
            double Bn = Bx * edge.normalVector.x +  By * edge.normalVector.y;
            initBns[edge.ind] = Bn;
            if(i < innerEdgeCount) {
                initBns[edge.ind] = Bn;
            }
            else{
                initGhostBNs[edge.ind - innerEdgeCount] = Bn;
            }
        }
    }
}

void MHDSolver2D::runSolver() {
    omp_set_num_threads(4);
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
        double sqr = edge.normalVector*edge.normalVector;
        double norm = std::sqrt(sqr);
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

        elemUs_prev = elemUs; //TODO: implement swap and test it

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
        #pragma parallel for
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
                        std::cout << "node1 # " << edge_j.nodeInd1 << " { " << node1.pos.x << ", " << node1.pos.y << " } \n";
                        std::cout << "node2 # " << edge_j.nodeInd2 << " { " << node2.pos.x << ", " << node2.pos.y << " } \n";
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
                        std::cout << "node1 # " << edge_j.nodeInd1 << " { " << node1.pos.x << ", " << node1.pos.y << " } \n";
                        std::cout << "node2 # " << edge_j.nodeInd2 << " { " << node2.pos.x << ", " << node2.pos.y << " } \n";
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
                        std::cout << "node1 # " << edge_j.nodeInd1 << " { " << node1.pos.x << ", " << node1.pos.y << " } \n";
                        std::cout << "node2 # " << edge_j.nodeInd2 << " { " << node2.pos.x << ", " << node2.pos.y << " } \n";
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
                    //inflowOrientation = sgn((u*(node1.pos.x - node2.pos.x) + v*(node1.pos.y - node2.pos.y)));
                    //std::cout << "scalar mult = " << (u*(node1.pos.x - node2.pos.x) + v*(node1.pos.y - node2.pos.y)) <<" , orient = " << inflowOrientation <<std::endl;
                    //std::cin.get();
                }
                else if(neighbourEdge.nodeInd2 == node.ind){
                    //inflowOrientation =  sgn((u*(node2.pos.x - node1.pos.x) + v*(node2.pos.y - node1.pos.y)));
                    //std::cout << "scalar mult = " <<  (u*(node2.pos.x - node1.pos.x) + v*(node2.pos.y - node1.pos.y))  <<" , orient = " << inflowOrientation <<std::endl;
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
            Vec2 centroid = getElementCentroid2D(elem, nodePool);
            if(elem.edgeIndexes.empty()){
                std::cerr << "Empty element! (no edges)" << std::endl;
            }
            double temp_sum_Bx = 0.0;
            double temp_sum_By = 0.0;
            for (const auto &edgeInd: elem.edgeIndexes) {
                Edge edge = edgePool.edges[edgeInd];
                Vec2 cTe = edge.midPoint - centroid;
                double scmult = cTe*edge.normalVector;
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
                        temp_sum_Bx += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        temp_sum_By += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
                    }
                    else{
                        continue;
                        temp_sum_Bx += ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        temp_sum_By += ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
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
                        temp_sum_Bx -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        temp_sum_By -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
                    }
                    else{
                        continue;
                        temp_sum_Bx -= ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        temp_sum_By -= ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
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
            file << "          " << node.pos.x << " " << node.pos.y << " " << 0.0 << "\n";
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
        Vec2 centroid = getElementCentroid2D(elem, nodePool);
        if (elem.edgeIndexes.empty()) {
            std::cerr << "Empty element! (no edges)" << std::endl;
        }
        double divergence = 0.0;
        for (const auto &edgeInd: elem.edgeIndexes) {
            Edge edge = edgePool.edges[edgeInd];
            Vec2 cTe = edge.midPoint - centroid;
            double scmult = cTe*edge.normalVector;
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

void MHDSolver2D::setInitCylindricElemUs() {
    //TODO: implement initial distribution with a respect to cylindrical geometry
    ElementPool ep = geometryWorld.getElementPool();
    innerElemCount = ep.elCount - geometryWorld.ghostElemCount;
    std::cout << "Element count = " << ep.elCount << " ; innerElCount = " << innerElemCount << " ; ghostElCOunt = " << geometryWorld.ghostElemCount << std::endl;
    EdgePool edgp = geometryWorld.getEdgePool();
    NodePool np = geometryWorld.getNodePool();
    innerEdgeCount = edgp.edgeCount - geometryWorld.ghostElemCount * 2;
    if(task_type == 0) {
        cflNum = 0.1;
        gam_hcr=1.4;
        freeFLowBoundaries = true;
        periodicBoundaries = false;
        freeFLowBoundaries2 = false;
        finalTime = 0.1;
        std::cout << "SOLVING TASKTYPE 0 (SOD'S PROBLEM)" << std::endl;
                                    /*rho   v_z   v_r    v_phi     p    Bz   Br   Bphi gam_hcr*/
        std::vector<double> BrioWu_L1{1.0, 0.0, 0.0, 0.0, 1.0,  0.0, 0.0, 0.0, gam_hcr};
        std::vector<double> BrioWu_R1{0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, gam_hcr};

        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = ep.elements[i];
            Vec2 centroid = elem.centroid2D;
            int elInd = elem.ind;
            if (centroid.y < 0.2) {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_L1);
            } else {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_R1);
            }
        }

        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for(const auto& [boundary, ghost] : geometryWorld.getNeighbourService().boundaryToGhostElements){
            int ghostInd = ghost - innerElemCount;
            initGhostElemUs[ghostInd] = initElemUs[boundary];
        }

        initBns.resize(edgp.edgeCount, 0.0);
        initGhostBNs.resize(geometryWorld.ghostElemCount*2, 0.0);
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < edgp.edgeCount; ++i) {
            Edge edge = edgp.edges[i];
            Vec2 midP = edge.midPoint;
            int edgeInd = edge.ind;
            Vec2 normal = edge.normalVector;
            std::vector<double> state;
            if (midP.y < 0.2) {
                state = state_from_primitive_vars(BrioWu_L1);
            } else {
                state = state_from_primitive_vars(BrioWu_R1);
            }
            double Bn = state[5] * normal.x + state[6] * normal.y;
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
    else if(task_type == 1) {
        cflNum = 0.1;
        freeFLowBoundaries = true;
        periodicBoundaries = false;
        freeFLowBoundaries2 = false;
        finalTime = 0.1;
        gam_hcr = 5.0/3.0;
        std::cout << "SOLVING TASKTYPE 1 (BRIO-WU PROBLEM)" << std::endl;
                                    /*rho   v_z  v_r  v_phi p    Bz   Br   Bphi gam_hcr*/
        std::vector<double> BrioWu_L1{1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.75, 0.0, gam_hcr};
        std::vector<double> BrioWu_R1{0.125, 0.0, 0.0, 0.0, 0.1,-1.0, 0.75, 0.0, gam_hcr};

        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = ep.elements[i];
            Vec2 centroid = elem.centroid2D;
            int elInd = elem.ind;
            if (centroid.y < 0.5) {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_L1);
            } else {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_R1);
            }
        }

        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for(int i = innerElemCount; i < ep.elCount; ++i ){
            Element elem = ep.elements[i];
            Vec2 centroid = elem.centroid2D;
            int elInd = elem.ind;
            int ghostInd = elInd - innerElemCount;
            if (centroid.y < 0.5) {
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
            Vec2 midP = edge.midPoint;
            int edgeInd = edge.ind;
            Vec2 normal = edge.normalVector;
            std::vector<double> state;
            if (midP.y < 0.5) {
                state = state_from_primitive_vars(BrioWu_L1);
            } else {
                state = state_from_primitive_vars(BrioWu_R1);
            }
            double Bn = state[5] * normal.x + state[6] * normal.y;
            initBns[edgeInd] = Bn;
            if(i < innerEdgeCount) {
                initBns[edgeInd] = Bn;
            }else{
                int ghostInd = edgeInd - innerEdgeCount;
                initGhostBNs[ghostInd] = Bn;
            }
            initEdgeUs[edgeInd] = state;
        }
    } else if(task_type == 2) {
        cflNum = 0.5;
        gam_hcr= 5.0/3.0;
        freeFLowBoundaries = true;
        periodicBoundaries = false;
        freeFLowBoundaries2 = false;
        finalTime = 1.0;
        std::cout << "SOLVING TASKTYPE 2 (NOH'S PROBLEM)" << std::endl;
                                    /*rho   v_z   v_r    v_phi     p    Bz   Br   Bphi gam_hcr*/
        //std::vector<double> BrioWu_L1{1.0, 0.0, 0.0, 0.0, 1.0,  0.0, 0.0, 0.0, gam_hcr};
        //std::vector<double> BrioWu_R1{0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, gam_hcr};
         // Параметры ИЗ ТАБЛИЦЫ B1
        double rho_0 = 3.1831e-3;//0.3;   
        double u_0 = -3.24;//-1.0;    
        double p0 = 0.0001;     
        double B_0 = 0.6;
        
        //double rho_0 = 0.0046;//0.3;   
        //double u_0 = -0.2;//-1.0;    
        //double p0 = 1.52945e-3;     
        //double B_0 = 0.598479;
        // double rho_0 = 1.0;
        // double u_0 = -2.02271;
        // double B_0 = 1.83267;
        // double p0 = 0.267274;
        std::cout << "r=1: rho=" << rho_0 
          << " v_r=" << u_0 
          << " B_phi=" << B_0 
          << " p=" << p0 << std::endl;

        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = ep.elements[i];
            Vec2 centroid = elem.centroid2D;
            int elInd = elem.ind;
                                /*0rho      1v_z  2v_r    3v_phi   4p   5Bz   6Br   7Bphi 8gam_hcr*/
            std::vector<double> ini{0.0,    0.0,  0.0,    0.0,    0.0,  0.0,  0.0,  0.0, gam_hcr};
           
            // Радиальная координата (предполагаем y >= 0)
            double r = centroid.y;
            
            // ПРАВИЛЬНЫЕ РАСПРЕДЕЛЕНИЯ:
            ini[0] = rho_0 * r * r;    // ρ = ρ₀*(r/Rₘ)²
            ini[2] = u_0;              // v_r = u_0 (радиальная скорость)
            ini[7] = B_0 * r;          // B_phi = B₀*(r/Rₘ)
            ini[4] = p0 *  r * r;       // p = p₀*(r/Rₘ)²

            initElemUs[elInd] = state_from_primitive_vars(ini);
        }

        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for(const auto& [boundary, ghost] : geometryWorld.getNeighbourService().boundaryToGhostElements){
            int ghostInd = ghost - innerElemCount;
            initGhostElemUs[ghostInd] = initElemUs[boundary];
        }

        initBns.resize(edgp.edgeCount, 0.0);
        initGhostBNs.resize(geometryWorld.ghostElemCount*2, 0.0);
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < edgp.edgeCount; ++i) {
            Edge edge = edgp.edges[i];
            Vec2 midP = edge.midPoint;
            int edgeInd = edge.ind;
            Vec2 normal = edge.normalVector;
            std::vector<double> ini{0.0,    0.0,  0.0,    0.0,    0.0,  0.0,  0.0,  0.0, gam_hcr};
             
            // Радиальная координата (предполагаем y >= 0)
            double r = midP.y;
            
            // ПРАВИЛЬНЫЕ РАСПРЕДЕЛЕНИЯ:
            ini[0] = rho_0 * r * r;    // ρ = ρ₀*(r/Rₘ)²
            ini[2] = u_0;              // v_r = u_0 (радиальная скорость)
            ini[7] = B_0 * r;          // B_phi = B₀*(r/Rₘ)
            ini[4] = p0 *  r * r;       // p = p₀*(r/Rₘ)²

            if(np.getNode(edge.nodeInd1).pos.y < 1e-15 && np.getNode(edge.nodeInd2).pos.y < 1e-15){
                ini[0] = 0.0;
                ini[7] = 0.0;
                ini[4] = 0.0;
            }

            initEdgeUs[edgeInd] = state_from_primitive_vars(ini);
            
            double Bn = initEdgeUs[edgeInd][5] * normal.x +initEdgeUs[edgeInd][6] * normal.y;
            if(i < innerEdgeCount) {
                initBns[edgeInd] = Bn;
            }else{
                int ghostInd = edgeInd - innerEdgeCount;
                initGhostBNs[ghostInd] = Bn;
            }
        }
    } else if(task_type == 3){
        std::cout << "SOLVING TASKTYPE 3 (MHD BLAST)" << std::endl;
        cflNum = 0.2;
        gam_hcr= 1.4;
        freeFLowBoundaries = true;
        periodicBoundaries = false;
        freeFLowBoundaries2 = false;
        finalTime = 0.01;
        double Bz = 0.0;//50.0/std::sqrt(4*std::numbers::pi_v<double>);
        std::cout << "Bz = " << Bz << std::endl;
        double R = 0.1;
        double R_sqr = R*R;
        double rho = 1.0;
        double p_ambient = 0.1;
        double p_inner = 1000.0;

        initElemUs.resize(innerElemCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = ep.elements[i];
            Vec2 centroid = elem.centroid2D;
            int elInd = elem.ind;
                                /*0rho      1v_z  2v_r    3v_phi   4p   5Bz   6Br   7Bphi 8gam_hcr*/
            std::vector<double> ini{rho,    0.0,  0.0,    0.0,    0.0,  Bz,  0.0,  0.0,  gam_hcr};
           
            // Радиальная координата (предполагаем y >= 0)
            double el_r = centroid.y;
            double el_z = centroid.x;
            bool isInner = el_r*el_r + (el_z - 0.5)*(el_z - 0.5) <= R_sqr;
            if(isInner){ //внутренний регион давления
                ini[4] = p_inner;
            }else {
                ini[4] = p_ambient;
            }
        
            initElemUs[elInd] = state_from_primitive_vars(ini);
        }

        initGhostElemUs.resize(geometryWorld.ghostElemCount, std::vector<double>(8, 0.0));
        for(const auto& [boundary, ghost] : geometryWorld.getNeighbourService().boundaryToGhostElements){
            int ghostInd = ghost - innerElemCount;
            initGhostElemUs[ghostInd] = initElemUs[boundary];
        }

        initBns.resize(edgp.edgeCount, 0.0);
        initGhostBNs.resize(geometryWorld.ghostElemCount*2, 0.0);
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        for (int i = 0; i < edgp.edgeCount; ++i) {
            Edge edge = edgp.edges[i];
            Vec2 midP = edge.midPoint;
            int edgeInd = edge.ind;
            Vec2 normal = edge.normalVector;
                                /*0rho      1v_z  2v_r    3v_phi   4p   5Bz   6Br   7Bphi 8gam_hcr*/
            std::vector<double> ini{rho,    0.0,  0.0,    0.0,    0.0,  Bz,  0.0,  0.0,  gam_hcr};
           
            // Радиальная координата (предполагаем y >= 0)
            double el_r = midP.y;
            double el_z = midP.x;
            bool isInner = el_r*el_r + (el_z - 0.5)*(el_z - 0.5) <= R_sqr;
            if(isInner){ //внутренний регион давления
                ini[4] = p_inner;
            }else {
                ini[4] = p_ambient;
            }

            initEdgeUs[edgeInd] = state_from_primitive_vars(ini);
            
            double Bn = initEdgeUs[edgeInd][5] * normal.x +initEdgeUs[edgeInd][6] * normal.y;
            if(i < innerEdgeCount) {
                initBns[edgeInd] = Bn;
            }else{
                int ghostInd = edgeInd - innerEdgeCount;
                initGhostBNs[ghostInd] = Bn;
            }
        }
        
    }
}

void MHDSolver2D::applyBoundaryConditions(NeighbourService& ns){
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
    else if(freeFLowBoundaries){
        // копируем значения в соответствующие фантомные ячейки (условие free flow)
        for(const auto& [boundary, ghost] : ns.boundaryToGhostElements){
            int ghostInd = ghost - innerElemCount;
            ghostElemUs[ghostInd] = elemUs[boundary];
        }
        // for(const auto& [top, bot]: ns.boundaryElemTopToBottom){
        //     int ghostIndTop = ns.boundaryToGhostElements[top] - innerElemCount;
        //     int ghostIndBot = ns.boundaryToGhostElements[bot] - innerElemCount;
        //     ghostElemUs[ghostIndTop][2] *=-1;
        //     ghostElemUs[ghostIndBot][6] *=-1;
        // }
        // for(const auto& [left, right]: ns.boundaryElemLeftToRight){
        //     int ghostIndLeft = ns.boundaryToGhostElements[left] - innerElemCount;
        //     int ghostIndRight = ns.boundaryToGhostElements[right] - innerElemCount;
        //     ghostElemUs[ghostIndLeft][1] *= -1;
        //     ghostElemUs[ghostIndRight][5] *= -1;
        // }
    }
//    else if(freeFLowBoundaries2) {
//        for (const auto &[top, bot]: ns.boundaryElemTopToBottom) {
//            Element elGhostTop = elPool.elements[ns.boundaryToGhostElements[top]];
//            for (const auto &edgeInd: elGhostTop.edgeIndexes) {
//                Edge edge = edgePool.edges[edgeInd];
//                if (edge.is_ghost) {
//                    int corspEdge = ns.edgeToGhostEdges[edge.ind];
//                    unrotated_fluxes[edge.ind] = (-1.0) * unrotated_fluxes[corspEdge];
//                    fluxes[edge.ind] = rotateStateFromAxisToNormal(unrotated_fluxes[edge.ind], edge.normalVector);
//                }
//            }
//            Element elGhostBot = elPool.elements[ns.boundaryToGhostElements[bot]];
//            for (const auto &edgeInd: elGhostBot.edgeIndexes) {
//                Edge edge = edgePool.edges[edgeInd];
//                if (edge.is_ghost) {
//                    int corspEdge = ns.edgeToGhostEdges[edge.ind];
//                    unrotated_fluxes[edge.ind] = (-1.0) * unrotated_fluxes[corspEdge];
//                    fluxes[edge.ind] = rotateStateFromAxisToNormal(unrotated_fluxes[edge.ind], edge.normalVector);
//                }
//            }
//        }
//    }
}

void MHDSolver2D::checkNan(bool& foundNan){
    //#pragma parallel for
    for(const auto& u: elemUs){
        for(const auto& val: u){
            if(std::isnan(val)){
                foundNan = true;
                break;
            }
        }
    }
}

void MHDSolver2D::applyZeroRConditions(const ElementPool& elPool, const EdgePool& edgePool, const NodePool& nodePool, const std::vector<std::vector<double>>& elemUs_prev){

    // v_r = v_phi = B_r = B_phi = 0 , r=0:
    for (int i = 0; i < innerElemCount; ++i) {
        Element elem = elPool.elements[i];
        bool isZeroR = false;
        if(elem.is_boundary){
            for(int& j: elem.edgeIndexes){
                Edge edge = edgePool.edges[j];
                if(nodePool.getNode(edge.nodeInd1).pos.y < 1e-15 && nodePool.getNode(edge.nodeInd2).pos.y < 1e-15){
                    isZeroR = true;
                }
            }
            if(isZeroR){
                elemUs[i][2] = 0.0;  // ρv_r = 0
                elemUs[i][6] = 0.0;  // B_r = 0
                elemUs[i][3] = elemUs_prev[i][3];  // Симметрия v_φ
                elemUs[i][7] = elemUs_prev[i][7];  // Симметрия B_φ
            }
        }
    }
}

void MHDSolver2D::runCylindricSolver() {
    /**
     * (z, r, \phi) <-> (x, y, z) (переобозначения в коде, не связь декартовых и цск)
     * (u, v, w) <-> (v_z, v_r, v_{\phi})
     * (Bx, By, Bz) <-> (B_z, B_r, B_{\phi})
     */
    // service
    EdgePool edgePool = geometryWorld.getEdgePool();
    ElementPool elPool = geometryWorld.getElementPool();
    NodePool nodePool = geometryWorld.getNodePool();
    NeighbourService ns = geometryWorld.getNeighbourService();

    // инициализируем состояния
    //task_type = 0; //TODO
    setInitCylindricElemUs();
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
        double norm = std::sqrt(edge.normalVector*edge.normalVector);
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

    std::vector<std::vector<double>> fluxes(edgePool.edges.size(), std::vector<double>(8, 0.0));
    std::vector<std::vector<double>> unrotated_fluxes(edgePool.edges.size(), std::vector<double>(8, 0.0));
    //MAX_ITERATIONS = 10;
    // основной цикл по времени
    while(currentTime < finalTime) {

        if (iterations >= MAX_ITERATIONS) {
            std::cout << "iterations limit!" << std::endl;
            break;
        }

        elemUs_prev.swap(elemUs);

        tau = std::max(min_tau, tau_from_cfl2D(cflNum, h,  elemUs_prev, gam_hcr));
        currentTime += tau;
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
            std::cout << std::setprecision(10) << "t = " << currentTime  << " (iteration " << iterations << ")"<< std::endl;
            writeVTU("OutputData/tmpres_" + std::to_string(iterations) + ".vtu", ghostOutput);
        }

        // вычисляем потоки, проходящие через каждое ребро
        // инициализируем вектор потоков через рёбра // MHD (HLLD) fluxes (from one element to another "<| -> |>")
#pragma parallel for
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
                continue; //TODO
                std::vector<double> U1 = rotateStateFromAxisToNormal(ghostElemUs[neighbour1.ind - innerElemCount], edge.normalVector);
                fluxes[edge.ind] = HLLD_flux(U1, U1, gam_hcr);
                unrotated_fluxes[edge.ind] =  fluxes[edge.ind];
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

        // по явной схеме обновляем газовые величины
        //#pragma omp parallel for
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = elPool.elements[i];
            if(elem.area < 1e-12){
                std::cout << "ZERO ELEM AREA!!!" << std::endl;
                continue;
            }
            
            if(elem.edgeIndexes.size() != 3){
                std::cout << "Bad edge vector size != 3" << std::endl;
                continue;
            }
                double r1 = nodePool.getNode(elem.nodeIndexes[0]).pos.y;
                double r2 = nodePool.getNode(elem.nodeIndexes[1]).pos.y;
                double r3 = nodePool.getNode(elem.nodeIndexes[2]).pos.y;
                double baricent = (r1 + r2 + r3) / 3.0;
                baricent = std::max(baricent, 1e-16);
                //double scale = tau / (baricent * elem.area);
                double volume = baricent * elem.area;
                double scale = tau / volume;
                std::vector<double> fluxSum(8, 0.0);
                for (int edgeIndex: elem.edgeIndexes) {
                    Edge edge_j = edgePool.edges[edgeIndex];
                    if (std::abs(edge_j.length) < 1e-12) {
                        std::cout << "BAD EDGE LENGTH = 0!!" << edge_j.length << std::endl;
                    }
                    double _r1 = nodePool.getNode(edge_j.nodeInd1).pos.y;
                    double _r2 = nodePool.getNode(edge_j.nodeInd2).pos.y;
                    double r_mid = (_r1 + _r2)/2.0;
                    r_mid = std::max(r_mid, 1e-16);
                    
                    if(r_mid < 1e-14){
                        continue;
                    }

                    double S_j = r_mid * edge_j.length;

                    Vec2 cTe = edge_j.midPoint - elem.centroid2D;
                    double scMult = cTe * edge_j.normalVector;
                    if (edge_j.neighbourInd1 == i) {
                        if(scMult <= 0.0){
                            std::cout << "WRONG NORMAL ORIENTATION!" << std::endl;
                        }
                        fluxSum = fluxSum + S_j * fluxes[edgeIndex];
                    } else if (edge_j.neighbourInd2 == i) {
                        if(scMult >= 0.0){
                            std::cout << "WRONG NORMAL ORIENTATION!" << std::endl;
                        }
                        fluxSum = fluxSum + (-1.0) * S_j * fluxes[edgeIndex];
                    } else {
                        std::cout << "No matching edge..."  << std::endl;
                    }
                }
                
                elemUs[i] = elemUs_prev[i] - scale * fluxSum;

               // Источниковые члены
                double rho_prev = elemUs_prev[i][0];
                double vr_prev = elemUs_prev[i][2] / rho_prev;  // v_r
                double vphi_prev = elemUs_prev[i][3] / rho_prev;  // v_φ
                double Br_prev = elemUs_prev[i][6];  // B_r
                double Bphi_prev = elemUs_prev[i][7];  // B_φ
                double pres = pressure(elemUs_prev[i], gam_hcr);
                double pSt = ptotal(pres, elemUs_prev[i][5], Br_prev, Bphi_prev);

                elemUs[i][2] += tau * (rho_prev * vphi_prev * vphi_prev - Bphi_prev * Bphi_prev + pSt)/baricent;
                elemUs[i][3] += tau * (-rho_prev * vr_prev * vphi_prev + Br_prev * Bphi_prev)/baricent;
                elemUs[i][7] += tau * (Bphi_prev * vr_prev - Br_prev * vphi_prev)/baricent;
            /*else{ //TODO: проверить потоки на рёбрах второй половины (возможно)
                int ghostInd = ns.boundaryToGhostElements[i];
                Element ghostEl = elPool.elements[ghostInd];
                int ghostNodeInd = *std::max_element(ghostEl.nodeIndexes.begin(),
                                                     ghostEl.nodeIndexes.end()); //тк при генерации гост элементы уже после основных создаются => их индексы всегда больше
                int ghostIndForState = ns.boundaryToGhostElements[i] - innerElemCount;
                double r_ghost = nodePool.getNode(ghostNodeInd).pos.y;
                if(r_ghost < 0.0){
                    continue;
                }else {
                    std::vector<int> commonEdges = findCommonElements(ghostEl.edgeIndexes, elem.edgeIndexes);
                    if (commonEdges.size() != 1) {
                        std::cout << "NO COMMON EDGES!" << std::endl;
                        std::cout << elem.edgeIndexes[0] << ",  " << elem.edgeIndexes[1] << ",  " << elem.edgeIndexes[2]
                                  << "; vs/s ghost: " << ghostEl.edgeIndexes[0] << ",  " << ghostEl.edgeIndexes[1]
                                  << ",  " << ghostEl.edgeIndexes[2] << std::endl;
                    }
                    int commonEdgeInd = commonEdges[0];
                    double area = elem.area * 2.0;
                    double baricent = (nodePool.getNode(elem.nodeIndexes[0]).pos.y + nodePool.getNode(elem.nodeIndexes[1]).pos.y\
                                        + nodePool.getNode(elem.nodeIndexes[2]).pos.y + nodePool.getNode(ghostNodeInd).pos.y) / 4.0;
                    if (baricent < 1e-14) {
                        baricent = 1e-14;
                    }
                    double scale = tau / (baricent * area);
                    for (const auto &edgeInd: elem.edgeIndexes) {
                        if (edgeInd == commonEdgeInd) {
                            continue;
                        }
                        Edge edge_j = edgePool.edges[edgeInd];
                        if (std::abs(edge_j.length) < 1e-16) {
                            std::cout << "BAD EDGE LENGTH = 0!! " << edge_j.length << std::endl;
                        }
                        double r1 = nodePool.getNode(edge_j.nodeInd1).pos.y;
                        double r2 = nodePool.getNode(edge_j.nodeInd2).pos.y;
                        double r_mid = (r1 + r2) / 2.0;
                        if (r_mid < 0) {
                            std::cout << "r_mid < 0! in boundary elem edge" << std::endl;
                        }
                        if (r_mid < 1e-14) {
                            r_mid = 1e-14;
                            continue;
                        }
                        if (edge_j.neighbourInd1 == i) {
                            fluxSum = fluxSum + r_mid * edge_j.length * fluxes[edgeInd];
                        } else if (edge_j.neighbourInd2 == i) {
                            fluxSum = fluxSum - r_mid * edge_j.length * fluxes[edgeInd];
                        } else {
                            std::cerr << "No matching edge..." << std::endl;
                        }
                    }
                    for (const auto &ghostEdgeInd: ghostEl.edgeIndexes) {
                        if (ghostEdgeInd == commonEdgeInd) {
                            continue;
                        }
                        Edge edge_j = edgePool.edges[ghostEdgeInd];
                        if (std::abs(edge_j.length) < 1e-16) {
                            std::cout << "BAD EDGE LENGTH = 0!! " << edge_j.length << std::endl;
                        }
                        double r1 = nodePool.getNode(edge_j.nodeInd1).pos.y;
                        double r2 = nodePool.getNode(edge_j.nodeInd2).pos.y;
                        double r_mid = (r1 + r2) / 2.0;
                        if (r_mid < 0) {
                            std::cout << r_mid << " r_mid < 0! in ghost elem edge" << std::endl;
                            continue;
                        }
                        if (r_mid < 1e-14) {
                            continue;
                            r_mid = 1e-14;
                        }
                        if (edge_j.neighbourInd1 == ghostEl.ind) {
                            fluxSum = fluxSum + r_mid * edge_j.length * fluxes[ghostEdgeInd];
                        } else if (edge_j.neighbourInd2 == ghostEl.ind) {
                            fluxSum = fluxSum - r_mid * edge_j.length * fluxes[ghostEdgeInd];
                        } else {
                            std::cerr << "No matching edge..." << std::endl;
                        }
                    }
                    elemUs[i] = elemUs_prev[i] - scale * fluxSum;
                }
            }*/
        }
        //applyZeroRConditions(elPool, edgePool, nodePool, elemUs_prev);
        //applyBoundaryConditions(ns);
        

//         //корректируем магнитные величины
//         //находим узловые значения нужных магнитных разностей //(v x B)z в узлах
        std::vector<double> nodeMagDiffs(nodePool.nodeCount, 0.0); //(v x B)z в узлах
//#pragma omp parallel for
        for (const auto &node: nodePool.nodes) {
            int tmp_count = 0;
            for (const auto &neighbourEdgeInd: ns.getEdgeNeighborsOfNode(node.ind)) {
                // добавить проверку на совпадения направлений в-ра скорости и нарпавляющего ребра
                Edge neighbourEdge = edgePool.edges[neighbourEdgeInd];
                // if(elPool.elements[neighbourEdge.neighbourInd1].is_ghost && neighbourEdge.neighbourInd2 == -1){
                //     //std::cout << "SEE A GHOST!" << std::endl;
                //     //continue;
                // }
                // Node node1 = nodePool.getNode(neighbourEdge.nodeInd1);
                // Node node2 = nodePool.getNode(neighbourEdge.nodeInd2);
                // int inflowOrientation = 1;
                // double u = edgeUs[neighbourEdgeInd][1]/edgeUs[neighbourEdgeInd][0];
                // double v = edgeUs[neighbourEdgeInd][2]/edgeUs[neighbourEdgeInd][0];
                // if(neighbourEdge.nodeInd1 == node.ind){
                //     //inflowOrientation = sgn((u*(node1.pos.x - node2.pos.x) + v*(node1.pos.y - node2.pos.y)));
                //     //std::cout << "scalar mult = " << (u*(node1.pos.x - node2.pos.x) + v*(node1.pos.y - node2.pos.y)) <<" , orient = " << inflowOrientation <<std::endl;
                //     //std::cin.get();
                // }
                // else if(neighbourEdge.nodeInd2 == node.ind){
                //     //inflowOrientation =  sgn((u*(node2.pos.x - node1.pos.x) + v*(node2.pos.y - node1.pos.y)));
                //     //std::cout << "scalar mult = " <<  (u*(node2.pos.x - node1.pos.x) + v*(node2.pos.y - node1.pos.y))  <<" , orient = " << inflowOrientation <<std::endl;
                //     //std::cin.get();
                // }
                // if(inflowOrientation > 0) {
                //     nodeMagDiffs[node.ind] += unrotated_fluxes[neighbourEdgeInd][6];
                //     ++tmp_count;
                // }
                nodeMagDiffs[node.ind] += unrotated_fluxes[neighbourEdgeInd][6];
                ++tmp_count;
            }
            if (tmp_count) {
                nodeMagDiffs[node.ind] /= tmp_count;
            }
        }

        //находим новое значение Bn в ребре (constrained transport)
        std::vector<double> bNs_prev(bNs);
        for (int i = 0; i < edgePool.edgeCount; ++i) {
            Edge edge = edgePool.edges[i];
            double r1 = nodePool.getNode(edge.nodeInd1).pos.y;
            double r2 = nodePool.getNode(edge.nodeInd2).pos.y;
            double r_mid = (r1 + r2)/2.0;
            if(r_mid < 1e-14){
                r_mid = 1e-14;
            }
            bNs[edge.ind] = bNs_prev[edge.ind] + (tau / (edge.length * r_mid)) * (r2 * nodeMagDiffs[edge.nodeInd2] -
                                                                        r1 * nodeMagDiffs[edge.nodeInd1]);
        }

        //сносим Bn в центр элемента
//#pragma omp parallel for
        for (int i = 0; i < innerElemCount; ++i) {
            Element elem = elPool.elements[i];
            Vec2 centroid = getElementCentroid2D(elem, nodePool);
            if(elem.edgeIndexes.empty()){
                std::cerr << "Empty element! (no edges)" << std::endl;
            }
            double temp_sum_Bx = 0.0;
            double temp_sum_By = 0.0;
            for (const auto &edgeInd: elem.edgeIndexes) {
                Edge edge = edgePool.edges[edgeInd];
                Vec2 cTe = edge.midPoint - centroid;
                double scmult = cTe*edge.normalVector;
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
                        temp_sum_Bx += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        temp_sum_By += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
                    }
                    else{
                        continue;
                        temp_sum_Bx += ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        temp_sum_By += ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
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
                        temp_sum_Bx -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        temp_sum_By -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
                    }
                    else{
                        continue;
                        temp_sum_Bx -= ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        temp_sum_By -= ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
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
        //applyZeroRConditions(elPool, edgePool, nodePool, elemUs_prev);
        //applyBoundaryConditions(ns);
        //

        checkNan(foundNan);
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



/**
 * GPU REALISATION OF CORTESIAN (X,Y,Z) 2D MHD SOLVER
 */
void MHDSolver2D::runGPUSolver() {

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

    // проверка на нормировку нормалей + предобработка для gpu
    for(const auto edge: edgePool.edges){
        double norm = std::sqrt(edge.normalVector*edge.normalVector);
        if(std::abs(1.0-norm) > 1e-15){
            std::cout << "bad normal! |1 - norm| = " << std::abs(1.0-norm) << " edge.is_ghost = " << edge.is_ghost <<std::endl;
        }
//        Element neighbour1 = elPool.elements[edge.neighbourInd1];
//        if(neighbour1.is_boundary && edge.neighbourInd2 == -1){
//            int ghostInd = ns.boundaryToGhostElements[neighbour1.ind] - innerElemCount;
//        }
//        else if(neighbour1.is_ghost && edge.neighbourInd2 == -1){}
//        else{}
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

        elemUs_prev = elemUs; //TODO: implement swap and test it

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
        std::vector<double> fluxes_flat(edgePool.edges.size() * 8, 0.0);
        std::vector<double> unrotated_fluxes_flat(edgePool.edges.size() * 8, 0.0);

        //computeHLLDFluxesGPU(elPool.elements, edgePool.edges, elemUs_prev, ghostElemUs, innerElemCount, ns.boundaryToGhostElements, gam_hcr, fluxes_flat, unrotated_fluxes_flat);

        std::vector<std::vector<double>> fluxes(edgePool.edges.size(), std::vector<double>(8, 0.0));
        std::vector<std::vector<double>> unrotated_fluxes(edgePool.edges.size(), std::vector<double>(8, 0.0));
        #pragma  parallel for
        for(int ii = 0; ii < fluxes.size(); ++ii){
            for(int kk =0 ; kk < 8; ++kk){
                fluxes[ii][kk] = fluxes_flat[ii*8 + kk];
                unrotated_fluxes[ii][kk] = unrotated_fluxes_flat[ii*8 + kk];
            }
        }

        //boundary conditions
        if(freeFLowBoundaries) {
            #pragma  parallel for
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
                        std::cout << "node1 # " << edge_j.nodeInd1 << " { " << node1.pos.x << ", " << node1.pos.y << " } \n";
                        std::cout << "node2 # " << edge_j.nodeInd2 << " { " << node2.pos.x << ", " << node2.pos.y << " } \n";
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
                        std::cout << "node1 # " << edge_j.nodeInd1 << " { " << node1.pos.x << ", " << node1.pos.y << " } \n";
                        std::cout << "node2 # " << edge_j.nodeInd2 << " { " << node2.pos.x << ", " << node2.pos.y << " } \n";
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
                        std::cout << "node1 # " << edge_j.nodeInd1 << " { " << node1.pos.x << ", " << node1.pos.y << " } \n";
                        std::cout << "node2 # " << edge_j.nodeInd2 << " { " << node2.pos.x << ", " << node2.pos.y << " } \n";
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
            #pragma  parallel for
            for(const auto& [top, bot]: ns.boundaryElemTopToBottom){
                int ghostIndTop = ns.boundaryToGhostElements[top] - innerElemCount;
                int ghostIndBot = ns.boundaryToGhostElements[bot] - innerElemCount;
                ghostElemUs[ghostIndTop] = elemUs[bot];
                ghostElemUs[ghostIndBot] = elemUs[top];
            }
            #pragma  parallel for
            for(const auto& [left, right]: ns.boundaryElemLeftToRight){
                int ghostIndLeft = ns.boundaryToGhostElements[left] - innerElemCount;
                int ghostIndRight = ns.boundaryToGhostElements[right] - innerElemCount;
                ghostElemUs[ghostIndLeft] =  elemUs[right];
                ghostElemUs[ghostIndRight] = elemUs[left];
            }
        }
        else{
            #pragma  parallel for
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
                    //inflowOrientation = sgn((u*(node1.pos.x - node2.pos.x) + v*(node1.pos.y - node2.pos.y)));
                    //std::cout << "scalar mult = " << (u*(node1.pos.x - node2.pos.x) + v*(node1.pos.y - node2.pos.y)) <<" , orient = " << inflowOrientation <<std::endl;
                    //std::cin.get();
                }
                else if(neighbourEdge.nodeInd2 == node.ind){
                    //inflowOrientation =  sgn((u*(node2.pos.x - node1.pos.x) + v*(node2.pos.y - node1.pos.y)));
                    //std::cout << "scalar mult = " <<  (u*(node2.pos.x - node1.pos.x) + v*(node2.pos.y - node1.pos.y))  <<" , orient = " << inflowOrientation <<std::endl;
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
            Vec2 centroid = getElementCentroid2D(elem, nodePool);
            if(elem.edgeIndexes.empty()){
                std::cerr << "Empty element! (no edges)" << std::endl;
            }
            double temp_sum_Bx = 0.0;
            double temp_sum_By = 0.0;
            for (const auto &edgeInd: elem.edgeIndexes) {
                Edge edge = edgePool.edges[edgeInd];
                Vec2 cTe = edge.midPoint - centroid;
                double scmult = cTe*edge.normalVector;
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
                        temp_sum_Bx += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        temp_sum_By += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
                    }
                    else{
                        continue;
                        temp_sum_Bx += ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        temp_sum_By += ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
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
                        temp_sum_Bx -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        temp_sum_By -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
                    }
                    else{
                        continue;
                        temp_sum_Bx -= ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.x - node_before.pos.x);
                        temp_sum_By -= ghostBNs[edgeInd - innerEdgeCount] * edge.length / (2 * elem.area) * (centroid.y - node_before.pos.y);
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
            #pragma  parallel for
            for(const auto& [top, bot]: ns.boundaryElemTopToBottom){
                int ghostIndTop = ns.boundaryToGhostElements[top] - innerElemCount;
                int ghostIndBot = ns.boundaryToGhostElements[bot] - innerElemCount;
                ghostElemUs[ghostIndTop] =  elemUs[bot];
                ghostElemUs[ghostIndBot] = elemUs[top];
            }
            #pragma  parallel for
            for(const auto& [left, right]: ns.boundaryElemLeftToRight){
                int ghostIndLeft = ns.boundaryToGhostElements[left] - innerElemCount;
                int ghostIndRight = ns.boundaryToGhostElements[right] - innerElemCount;
                ghostElemUs[ghostIndLeft] =  elemUs[right];
                ghostElemUs[ghostIndRight] = elemUs[left];
            }
        }
        else{
            // копируем значения в соответствующие фантомные ячейки (условие free flow)
            #pragma  parallel for
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

