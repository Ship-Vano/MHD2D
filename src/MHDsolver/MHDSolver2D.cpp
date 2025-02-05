//
// Created by Иван on 10/21/2024.
//

#include "MHDSolver2D.h"

// конструктор решателя
MHDSolver2D::MHDSolver2D(const World &world): geometryWorld(world), nodeUs(world.getNodePool().nodeCount),
    elemUs(world.getElementPool().elCount), edgeUs(world.getEdgePool().edgeCount), initElemUs(world.getElementPool().elCount),
    initBns(world.getEdgePool().edgeCount), bNs(world.getEdgePool().edgeCount){}

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
std::vector<double> MHDSolver2D::rotateStateFromAxisToNormal(vector<double> &U, const vector<double>& n) {
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
    double p = pressure(res, gam_hcr);
    double energy_rotated = energy(gam_hcr, p, res[0], res[1]/res[0], res[2]/res[0], res[3]/res[0], res[5], res[6], res[7]);
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
    double p = pressure(res, gam_hcr);
    double energy_rotated = energy(gam_hcr, p, res[0], res[1]/res[0], res[2]/res[0], res[3]/res[0], res[5], res[6], res[7]);
    //res[4] = energy_rotated;
    return res;
}

// условие Куранта (имплементация №1)
double MHDSolver2D::tau_from_cfl2D(const double& sigma, const double& hx, const double& hy, const std::vector<std::vector<double>>& states, const double& gam_hcr) {
    double max_speed = 0.0;

    for (const auto& state : states) {
        // для cfast брать модули в-в v и B
        double u = std::fabs(state[1] / state[0]);
        double v = std::fabs(state[2] / state[0]);
        double cf = cfast(state, gam_hcr);
        double local_speed = (u + cf) / hx + (v + cf) / hy;
        max_speed = std::max(max_speed, local_speed);
    }
    if(max_speed < 1e-16){
        max_speed = 1e-14;
    }
    double tau = sigma / max_speed;
    const double max_tau = 1e-2; // Define maximum allowable time step
    return std::min(tau, max_tau);
}

// условие Куранта (имплементация №2)
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

    // Brio-Wu problem
    if(task_type == 1) {
        cflNum = 0.1;
        freeFLowBoundaries = true;
        periodicBoundaries = false;
        std::cout << "SOLVING TASKTYPE 1 (BRIO-WU TEST)" << std::endl;
        /*rho  u   v   w   p   Bx   By   Bz gam_hcr*/
        std::vector<double> BrioWu_L1{1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0, gam_hcr};
        std::vector<double> BrioWu_R1{0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0, gam_hcr};
        ElementPool ep = geometryWorld.getElementPool();
        initElemUs.resize(ep.elCount, std::vector<double>(8, 0.0));
        for (const auto &elem: ep.elements) {
            std::vector<double> centroid = elem.centroid2D;
            int elInd = elem.ind;
            if (centroid[0] < 0.5) {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_L1);
            } else {
                initElemUs[elInd] = state_from_primitive_vars(BrioWu_R1);
            }
        }
        EdgePool edgp = geometryWorld.getEdgePool();
        initBns.resize(edgp.edgeCount, 0.0);
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        for (const auto &edge: edgp.edges) {
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
            initEdgeUs[edgeInd] = state;
        }
    }

    // alfven wave test
    else if(task_type == 2){
        std::cout << "SOLVING TASKTYPE 2 (ALFVEN WAVE TEST)" << std::endl;
        double r = 6.0;
        double nx = 1.0 / std::sqrt(r * r + 1);
        double ny = r / std::sqrt(r * r + 1);
        double rho0 = 1.0;
        double p0 = 1.0;
        double v0 = 0.0;
        double B0 = 1.0;
        double xi = 0.2;
        gam_hcr = 5.0/3.0;
        cflNum = 0.2;
        periodicBoundaries = true;
        ElementPool ep = geometryWorld.getElementPool();
        EdgePool edgp = geometryWorld.getEdgePool();
        initElemUs.resize(ep.elCount, std::vector<double>(8, 0.0));
        initEdgeUs.resize(edgp.edgeCount, std::vector<double>(8, 0.0));
        initBns.resize(edgp.edgeCount, 0.0);

        for (const auto &elem: ep.elements) {
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
            initElemUs[elem.ind] = state_from_primitive_vars(rho0, u_0, v_0, w_0, p0, Bx_0, By_0, Bz_0, gam_hcr);
        }

        for(const auto &edge: edgp.edges){
            if(edge.ind >= initEdgeUs.size() ){
                std::cout << edge.ind << std::endl;
            }
            double x = edge.midPoint[0];
            double y = edge.midPoint[1];
            double phase = 2.0 * M_PI / ny * (nx * x + ny * y);
            double u_0 = v0 * nx - xi * ny * std::cos(phase);
            double v_0 = v0 * ny + xi * nx * std::cos(phase);
            double w_0 = xi * std::sin(phase);
            double Bx_0 = B0 * nx + xi * ny * std::sqrt(4 * M_PI * rho0) * std::cos(phase);
            double By_0 = B0 * ny - xi * nx * std::sqrt(4 * M_PI * rho0) * std::cos(phase);
            double Bz_0 = - xi * std::sqrt(4 * M_PI * rho0) * std::sin(phase);
            initEdgeUs[edge.ind] = state_from_primitive_vars(rho0, u_0, v_0, w_0, p0, Bx_0, By_0, Bz_0, gam_hcr);
            double Bn = Bx_0 * edge.normalVector[0] + By_0 * edge.normalVector[1];
            initBns[edge.ind] = Bn;
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

    // сделать старые дубликаты состояний (предыдущие состояния) чтобы в новые записывать расчёты
    std::vector<std::vector<double>> elemUs_prev(elemUs);
    std::vector<std::vector<double>> nodeUs_prev(nodeUs);
    std::vector<std::vector<double>> edgeUs_prev(edgeUs);

    // проверка на нормировку нормалей
    for(const auto edge: edgePool.edges){
        double norm = std::sqrt(edge.normalVector[0]*edge.normalVector[0] + edge.normalVector[1]*edge.normalVector[1]);
        if(std::abs(1.0-norm) > 1e-15){
            std::cout << "bad normal! " << std::abs(1-norm) << std::endl;
        }
    }

    // дивергенция магнитного поля
    double divergence = computeDivergence(elemUs, edgePool);
    if (divergence > 1e-10) {
        std::cout << "Max divergence: " << divergence << std::endl;
    }

    double h = edgePool.minEdgeLen;
    std::cout << "Min h = " << h << std::endl;

    double currentTime = startTime;
    bool foundNan = false; // флаг для поиска NaN-значений
    int iterations = 0; // число текущих итераций

    // соновной цикл по времени
    while(currentTime < finalTime) {
        //std::cout << "iteration #" << iterations << std::endl;
        elemUs_prev.swap(elemUs);
        nodeUs_prev.swap(nodeUs);
        edgeUs_prev.swap(edgeUs);

        tau = std::max(min_tau, tau_from_cfl2D(cflNum, h, edgeUs_prev, gam_hcr, edgePool));

        // ограничение: в tauCfl идут только граничные состояния (требуется убрать после отладки гран условий)
        /*std::vector<std::vector<double>> usForCFL;
        for(const auto & elem: elPool.elements){
            bool is_boundary = false;
            for(const auto& edgind: elem.edgeIndexes){
                if(edgePool.edges[edgind].neighbourInd2 == -1){
                    is_boundary = true;
                    break;
                }
            }
            if(is_boundary){
                continue;
            }
            usForCFL.push_back(elemUs_prev[elem.ind]);
        }

        // выбор оптимального шага по времени
        tau = std::max(min_tau, tau_from_cfl2D(cflNum, h, h, usForCFL, gam_hcr));*/

        currentTime += tau;

        // подбор времени (если "перепрыгнули" за финальное время)
        if (currentTime > finalTime) {
            tau -= (currentTime - finalTime);
            currentTime = finalTime;
        }

        // запись в файл временного результата
        if (iterations % iterationsPerFrame == 0) {
            std::cout << std::setprecision(10) << "t = " << currentTime << std::endl;
            writeVTU("OutputData/tmpres_" + std::to_string(iterations) + ".vtu", geometryWorld, elemUs);
        }

        // вывд максимальной дивергенции
        if (debugDivergence) {
            divergence = computeDivergence(elemUs, edgePool);
            if (divergence > 1e-10) {
                std::cout << "Max divergence: " << divergence << std::endl;
            }
        }

        // г.у (фикт ячейки) поставить перед выч-м потоков

        // вычисляем потоки, проходящие через каждое ребро
        // инициализируем вектор потоков через рёбра // MHD (HLLD) fluxes (from one element to another "<| -> |>")
        std::vector<std::vector<double>> fluxes(edgePool.edgeCount, std::vector<double>(8, 0.0));
        std::vector<std::vector<double>> unrotated_fluxes(edgePool.edgeCount, std::vector<double>(8, 0.0));
//#pragma parallel for
        for (const auto &edge: edgePool.edges) {
            int neighbour1 = edge.neighbourInd1;
            int neighbour2 = edge.neighbourInd2;
            std::vector<double> U1 = rotateStateFromAxisToNormal(elemUs_prev[neighbour1], edge.normalVector);
            if (neighbour2 != -1) {
                std::vector<double> U2 = rotateStateFromAxisToNormal(elemUs_prev[neighbour2], edge.normalVector);
                fluxes[edge.ind] = HLLD_flux(U1, U2, gam_hcr);
                unrotated_fluxes[edge.ind] = HLLD_flux(U1, U2, gam_hcr);
                fluxes[edge.ind] = rotateStateFromNormalToAxisX(fluxes[edge.ind], edge.normalVector);
            } else { //здесь задавать гран условие ещё мб на поток
                fluxes[edge.ind] = HLLD_flux(U1, U1, gam_hcr);
                unrotated_fluxes[edge.ind] = HLLD_flux(U1, U1, gam_hcr);
                fluxes[edge.ind] = rotateStateFromNormalToAxisX(fluxes[edge.ind], edge.normalVector);
            }
        }


        // по явной схеме обновляем газовые величины
        //#pragma omp parallel for
        for (const auto &elem: elPool.elements) {
            bool is_boundary = false;
            int i = elem.ind;
            std::vector<double> fluxSum(8, 0.0);
            for (int edgeIndex: elem.edgeIndexes) {
                Edge edge_j = edgePool.edges[edgeIndex];
                if (edge_j.neighbourInd1 == i) {
                    fluxSum = fluxSum + edge_j.length * fluxes[edgeIndex];
                } else if (edge_j.neighbourInd2 == i) {
                    fluxSum = fluxSum - edge_j.length * fluxes[edgeIndex];
                } else {
                        std::cerr << "No matching edge..." << std::endl;
                }
            }
            if(elem.is_boundary){ //граничная ячейка : либо периодические, либо ничего не делаем - ИСТОРИЧЕСКИЕ ГУ
                if(periodicBoundaries){
                    if(ns.boundaryElemTopToBottom.count(i) != 0){
                        elemUs[i].swap(elemUs_prev[ns.boundaryElemTopToBottom[i]]);
                        elemUs[ns.boundaryElemTopToBottom[i]].swap(elemUs_prev[i]);
                    }
                    else if(ns.boundaryElemLeftToRight.count(i) != 0){
                        elemUs[i].swap(elemUs_prev[ns.boundaryElemLeftToRight[i]]);
                        elemUs[ns.boundaryElemLeftToRight[i]].swap(elemUs_prev[i]);
                    }
                }
            }
            else{
                elemUs[i] = elemUs_prev[i] - tau / elem.area * fluxSum;
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
                nodeMagDiffs[node.ind] += unrotated_fluxes[neighbourEdgeInd][6];
                ++tmp_count;
            }
            if (tmp_count) {
                nodeMagDiffs[node.ind] /= tmp_count;
            }
        }

        //находим новое значение Bn в ребре
        std::vector<double> bNs_prev(bNs);
        for (const auto& edge: edgePool.edges) {
            bNs[edge.ind] = bNs_prev[edge.ind] - (tau / edge.length) * (nodeMagDiffs[edge.nodeInd1] -
                                                        nodeMagDiffs[edge.nodeInd2]);
        }

        //сносим Bn в центр элемента
//#pragma omp parallel for
        for (const auto &elem: elPool.elements) {
            std::vector<double> centroid = getElementCentroid2D(elem, nodePool);
            if(elem.edgeIndexes.empty()){
                std::cout << "Empty element! (no edges)" << std::endl;
            }
            double temp_sum_Bx = 0.0;
            double temp_sum_By = 0.0;
            for (const auto &edgeInd: elem.edgeIndexes) {
                Edge edge = edgePool.edges[edgeInd];
                if (edge.neighbourInd1 == elem.ind) {
                    // у первого соседа в эдже заданы ноды в порядке полодительного обхода и нормаль тоже
                    const auto nodeInElemInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(),
                                                         edge.nodeInd1);
                    int node_before_ind =
                            nodeInElemInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim - 1] : *(
                                    nodeInElemInd - 1);
                    //std::cout << "current ind = " << edge.nodeInd1 <<" Elem's node indexes: " << elem.nodeIndexes[0] << " "<< elem.nodeIndexes[1] << " "<< elem.nodeIndexes[2] << ", node_beforeInd = " << node_before_ind << std::endl;
                    //std::cin.get();
                    Node node_before = nodePool.getNode(node_before_ind);
                    /*elemUs[elem.ind][5]*/temp_sum_Bx += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[0] - node_before.x);
                    /*elemUs[elem.ind][6]*/temp_sum_By += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[1] - node_before.y);
                } else {
                    // а вот для второго нужно умножать на -1 и в обратном порядке
                    const auto nodeInElemInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(),
                                                         edge.nodeInd2);
                    int node_before_ind =
                            nodeInElemInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim - 1] : *(
                                    nodeInElemInd - 1);
                    //        auto itNode1 = std::find(elementNodes.begin(), elementNodes.end(), node2);
                    Node node_before = nodePool.getNode(node_before_ind);
                    /*elemUs[elem.ind][5]*/temp_sum_Bx -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[0] - node_before.x);
                    /*elemUs[elem.ind][6]*/temp_sum_By -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[1] - node_before.y);
                }
            }if(elem.is_boundary){
               //do smth, например ничего - ИСТОРИЧЕСКИЕ ГУ
               // периодические г.у.
               if(periodicBoundaries){
                   if(ns.boundaryElemTopToBottom.count(elem.ind) != 0){
                       std::swap(elemUs[elem.ind][5], elemUs_prev[ns.boundaryElemTopToBottom[elem.ind]][5]);
                       std::swap(elemUs[ns.boundaryElemTopToBottom[elem.ind]][5], elemUs_prev[elem.ind][5]);

                       std::swap(elemUs[elem.ind][6], elemUs_prev[ns.boundaryElemTopToBottom[elem.ind]][6]);
                       std::swap(elemUs[ns.boundaryElemTopToBottom[elem.ind]][6], elemUs_prev[elem.ind][6]);
                   }
                   else if(ns.boundaryElemLeftToRight.count(elem.ind) != 0){
                       std::swap(elemUs[elem.ind][5], elemUs_prev[ns.boundaryElemLeftToRight[elem.ind]][5]);
                       std::swap(elemUs[ns.boundaryElemLeftToRight[elem.ind]][5], elemUs_prev[elem.ind][5]);
                       std::swap(elemUs[elem.ind][6], elemUs_prev[ns.boundaryElemLeftToRight[elem.ind]][6]);
                       std::swap(elemUs[ns.boundaryElemLeftToRight[elem.ind]][6], elemUs_prev[elem.ind][6]);
                   }
               }
            }else{
                elemUs[elem.ind][5] = temp_sum_Bx;
                elemUs[elem.ind][6] = temp_sum_By;
            }
        }

        for(const auto& edge: edgePool.edges){
            vector<double> state1 = elemUs[edge.neighbourInd1];
            if(edge.neighbourInd2 > -1) {
                vector<double> state2 = elemUs[edge.neighbourInd2];
                edgeUs[edge.ind] = 0.5 * (state1 + state2);
            }
            else {
                edgeUs[edge.ind] = state1;
            }
        }

        ++iterations;

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
        if(iterations > MAX_ITERATIONS){
            std::cout << "iterations limit!" << std::endl;
            break;
        }
    }

    std::cout << "Final time = " << currentTime << "; iterations = "<< iterations<<  std::endl;
}


void writeVTU(const std::string& filename, const World& geometryWorld, const std::vector<std::vector<double>>& elemUs) {
    const NodePool& np = geometryWorld.getNodePool();
    const ElementPool& ep = geometryWorld.getElementPool();

    std::ofstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << np.nodeCount << "\" NumberOfCells=\"" << ep.elCount << "\">\n";

    // Write points (nodes)
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& node : np.nodes) {
        file << "          " << node.x << " " << node.y << " " << node.z << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // Write cells (elements)
    file << "      <Cells>\n";

    // Connectivity
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& element : ep.elements) {
        for (const auto& nodeIndex : element.nodeIndexes) {
            file << "          " << nodeIndex << " ";
        }
        file << "\n";
    }
    file << "        </DataArray>\n";

    // Offsets
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (const auto& element : ep.elements) {
        offset += element.dim;
        file << "          " << offset << "\n";
    }
    file << "        </DataArray>\n";

    // Types
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (const auto& element : ep.elements) {
        if (element.dim == 3) {
            file << "          5\n"; // Triangle
        } else if (element.dim == 4) {
            file << "          9\n"; // Quadrilateral
        }
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";

    // Write solution data (elemUs)
    file << "      <CellData Scalars=\"elemUs\">\n";
    file << R"(        <DataArray type="Float64" NumberOfComponents=")" << elemUs[0].size() << "\" Name=\"elemUs\" format=\"ascii\">\n";
    for (const auto& U : elemUs) {
        for (const auto& value : U) {
            file << "          " << value << " ";
        }
        file << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </CellData>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
}


double computeDivergence(const std::vector<std::vector<double>>& elemUs, const EdgePool& edgePool) {
    double max_divergence = 0.0;

    for (const auto& edge : edgePool.edges) {
        double Bn1 = elemUs[edge.neighbourInd1][5] * edge.normalVector[0] +
                     elemUs[edge.neighbourInd1][6] * edge.normalVector[1];
        double Bn2 = (edge.neighbourInd2 > -1) ?
                     elemUs[edge.neighbourInd2][5] * edge.normalVector[0] +
                     elemUs[edge.neighbourInd2][6] * edge.normalVector[1]
                                               : Bn1;
        max_divergence = std::max(max_divergence, std::fabs(Bn1 - Bn2));
    }

    return max_divergence;
}




