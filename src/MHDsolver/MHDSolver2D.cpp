//
// Created by Иван on 10/21/2024.
//

#include "MHDSolver2D.h"


MHDSolver2D::MHDSolver2D(const World &world): geometryWorld(world), nodeUs(),
    elemUs(), edgeUs(), initElemUs(), initBns(), bNs(){}

double tau_from_cfl2D(const double& sigma, const double& hx, const double& hy, const std::vector<std::vector<double>>& states, const double& gam_hcr) {
    double max_speed = 0.0;

    for (const auto& state : states) {
        double u = std::fabs(state[1] / state[0]);
        double v = std::fabs(state[2] / state[0]);
        double cf = cfast(state, gam_hcr);
        double local_speed = (u + cf) / hx + (v + cf) / hy;
        max_speed = std::max(max_speed, local_speed);
    }
    if(max_speed < 1e-16){
        max_speed = 1e-8;
    }
    double tau = sigma / max_speed;
    const double max_tau = 1e-3; // Define maximum allowable time step
    return std::min(tau, max_tau);
}

/*вращения*/
// n = {cos(j), sin(j)} = {n.x, n.y}
// ({cos(j), -sin(j), 0}, {sin(j), cos(j), 0}, {0,0,1}) --- around OZ counterclockwise
// rotate 1: from normal to OX:      {v.x * n.x + vy * n.y, - v.x * n.y + v.y * n.x, v.z}
// вращаем {u,v,w} и {Bx, By, Bz}, остальные остаются на месте
std::vector<double> MHDSolver2D::rotateStateFromAxisToNormal(vector<double> &U, const vector<double>& n) {
    std::vector<double> res(U);
    res[1] =  U[1]*n[0] + U[2]*n[1];
    res[2] = -U[1]*n[1] + U[2]*n[0];
    res[5] =  U[5]*n[0] + U[6]*n[1];
    res[6] = -U[5]*n[1] + U[6]*n[0];
    return res;
}

// ({cos(j), sin(j), 0}, {-sin(j), cos(j), 0}, {0,0,1}) --- around OZ clockwise
// rotate 2: from OX to normal:      {v.x * n.x - vy * n.y,  v.x * n.y + v.y * n.x, v.z}
std::vector<double> MHDSolver2D::rotateStateFromNormalToAxisX(vector<double> &U, const vector<double>& n) {
    std::vector<double> res(U);
    res[1] =  U[1]*n[0] - U[2]*n[1];
    res[2] =  U[1]*n[1] + U[2]*n[0];
    res[5] =  U[5]*n[0] - U[6]*n[1];
    res[6] =  U[5]*n[1] + U[6]*n[0];
    return res;
}

void MHDSolver2D::setInitElemUs() {
    /*rho  u   v   w   p   Bx   By   Bz gam_hcr*/
    std::vector<double> BrioWu_L1{1.,   0., 0., 0., 1., 0.75, 1., 0., gam_hcr};
    std::vector<double> BrioWu_R1{0.125, 0., 0., 0., 0.1, 0.75, -1., 0., gam_hcr};
    ElementPool ep = geometryWorld.getElementPool();
    initElemUs.resize(ep.elCount, std::vector<double>(8, 0.0));
    for(const auto& elem: ep.elements){
        std::vector<double> centroid = elem.centroid2D;
        int elInd = elem.ind;
        if(centroid[0] < 0.5){
            initElemUs[elInd] = state_from_primitive_vars(BrioWu_L1);
        }
        else {
            initElemUs[elInd] = state_from_primitive_vars(BrioWu_R1);
        }
    }
    EdgePool edgp = geometryWorld.getEdgePool();
    initBns.resize(edgp.edgeCount, 0.0);
    for(const auto& edge: edgp.edges){
        std::vector<double> midP = edge.midPoint;
        int edgeInd = edge.ind;
        std::vector<double> normal = edge.normalVector;
        std::vector<double> state;
        if(midP[0] < 0.5){
            state = state_from_primitive_vars(BrioWu_L1);
        }
        else{
            state = state_from_primitive_vars(BrioWu_R1);
        }
        double Bn = state[5] * normal[0] + state[6] * normal[1];
        initBns[edgeInd] = Bn;
    }
}

void MHDSolver2D::runSolver() {

    omp_set_num_threads(omp_get_max_threads());

    setInitElemUs();

    // service
    EdgePool edgePool = geometryWorld.getEdgePool();
    ElementPool elPool = geometryWorld.getElementPool();
    NodePool nodePool = geometryWorld.getNodePool();
    NeighbourService ns = geometryWorld.getNeighbourService();

    // инициализируем состояния
    elemUs = initElemUs;
    bNs = initBns;

    // сделать старые дубликаты состояний (предыдущие состояния) чтобы в новые записывать расчёты
    std::vector<std::vector<double>> elemUs_prev(elemUs);
    std::vector<std::vector<double>> nodeUs_prev(nodeUs);
    std::vector<std::vector<double>> edgeUs_prev(edgeUs);

    // инициализируем вектор потоков через рёбра // MHD (HLLD) fluxes (from one element to another "<| -> |>")
    std::vector<std::vector<double>> fluxes(edgePool.edgeCount, std::vector<double>(8, 0.0));
    std::vector<std::vector<double>> unrotated_fluxes(edgePool.edgeCount, std::vector<double>(8, 0.0));

    double h = edgePool.minEdgeLen;
    std::cout << "Min h = " << h << std::endl;

    double currentTime = startTime;

    int iterations = 0;

    while(currentTime < finalTime) {
        elemUs_prev.swap(elemUs);
        nodeUs_prev.swap(nodeUs);
        edgeUs_prev.swap(edgeUs);

        tau = tau_from_cfl2D(cflNum, h, h, elemUs, gam_hcr);
        currentTime += tau;
        if(currentTime > finalTime){
            tau -= (currentTime - finalTime);
            currentTime = finalTime;
        }

        if(iterations % 100 == 0)
            std::cout << std::setprecision(10) << "t = "<< currentTime << std::endl;

        //(2) вычисляем потоки, проходящие через каждое ребро
#pragma parallel for
        for (const auto &edge: edgePool.edges) {
            Node node1 = nodePool.getNode(edge.nodeInd1);
            Node node2 = nodePool.getNode(edge.nodeInd2);
            int neighbour1 = edge.neighbourInd1;
            int neighbour2 = edge.neighbourInd2;
            std::vector<double> U1 = rotateStateFromNormalToAxisX(elemUs[neighbour1], edge.normalVector);
            if (neighbour2 > 0) {
                std::vector<double> U2 = rotateStateFromNormalToAxisX(elemUs[neighbour2], edge.normalVector);
                fluxes[edge.ind] = HLLD_flux(U1, U2, gam_hcr);
                unrotated_fluxes[edge.ind] = HLLD_flux(U1, U2, gam_hcr);
                fluxes[edge.ind] = rotateStateFromAxisToNormal(fluxes[edge.ind], edge.normalVector);
            } else { //здесь задавать гран условие ещё мб на поток
                fluxes[edge.ind] = HLLD_flux(U1, U1, gam_hcr);
                unrotated_fluxes[edge.ind] = HLLD_flux(U1, U1, gam_hcr);
                fluxes[edge.ind] = rotateStateFromAxisToNormal(fluxes[edge.ind], edge.normalVector);
            }
        }


        //по явной схеме обновляем газовые величины
        #pragma omp parallel for
        for (const auto &elem: elPool.elements) {
            int i = elem.ind;
            std::vector<double> fluxSum(8, 0.0);
            for (int edgeIndex: elem.edgeIndexes) {
                Edge edge_j = edgePool.edges[edgeIndex];
                if (edge_j.neighbourInd1 == i) {
                    fluxSum = fluxSum + edge_j.length * fluxes[edgeIndex];
                } else if (edge_j.neighbourInd2 == i){
                    fluxSum = fluxSum - edge_j.length * fluxes[edgeIndex];
                }
            }
            elemUs[i] = elemUs_prev[i] - tau / elem.area * fluxSum;
        }

        //корректируем магитные величины
        //находим узловые значения нужных магнитных разностей //(v x B)z в узлах
        std::vector<double> nodeMagDiffs(nodePool.nodeCount, 0.0); //(v x B)z в узлах
#pragma omp parallel for
        for (const auto &node: nodePool.nodes) {
            int tmp_count = 0;
            for (const auto &neighbourEdgeInd: ns.getEdgeNeighborsOfNode(node.ind)) {
                nodeMagDiffs[node.ind] += unrotated_fluxes[neighbourEdgeInd][6];
                ++tmp_count;
            }
            if (tmp_count) {
                nodeMagDiffs[node.ind] /= tmp_count;
            }
        }

        //находим новое значение Bn в ребре
        std::vector<double> bNs_prev(bNs);
        for (int i = 0; i < edgePool.edgeCount; ++i) {
            Edge edge = edgePool.edges[i];
            bNs[i] = bNs_prev[i] + tau / edge.length * (nodeMagDiffs[edge.nodeInd1] -
                                                        nodeMagDiffs[edge.nodeInd2]);
        }

        //сносим Bn в центр элемента
#pragma omp parallel for
        for (const auto &elem: elPool.elements) {
            elemUs[elem.ind][5] = 0.0;
            elemUs[elem.ind][6] = 0.0;
            std::vector<double> centroid = getElementCentroid2D(elem, nodePool);
            for (const auto &edgeInd: elem.edgeIndexes) {
                Edge edge = edgePool.edges[edgeInd];
                if (edge.neighbourInd1 == elem.ind) {
                    // у первого соседа в эдже заданы ноды в порядке полодительного обхода и нормаль тоже
                    const auto nodeInElemInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(),
                                                         edge.nodeInd1);
                    int node_before_ind =
                            nodeInElemInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim - 1] : *(
                                    nodeInElemInd - 1);
                    //        auto itNode1 = std::find(elementNodes.begin(), elementNodes.end(), node1);
                    Node node_before = nodePool.getNode(node_before_ind);
                    elemUs[elem.ind][5] += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[0] - node_before.x);
                    elemUs[elem.ind][6] += bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[1] - node_before.y);
                } else {
                    // а вот для второго нужно умножать на -1 и в обратном порядке
                    const auto nodeInElemInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(),
                                                         edge.nodeInd2);
                    int node_before_ind =
                            nodeInElemInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim - 1] : *(
                                    nodeInElemInd - 1);
                    //        auto itNode1 = std::find(elementNodes.begin(), elementNodes.end(), node2);
                    Node node_before = nodePool.getNode(node_before_ind);
                    elemUs[elem.ind][5] -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[0] - node_before.x);
                    elemUs[elem.ind][6] -= bNs[edgeInd] * edge.length / (2 * elem.area) * (centroid[1] - node_before.y);
                }
            }
        }
        ++iterations;
        //int elemNum = 0;
#pragma parallel for
        for(const auto& u: elemUs){
            for(const auto& val: u){
                if(std::isnan(val)){
                    std::cout << "found a nan value!" << std::endl;
                    //std::cout << "ElemNum = " << elemNum << std::endl;
                    std::cin.get();
                }
            }
           // ++elemNum;
        }

        if(iterations > MAX_ITERATIONS){
            writeVTU("OutputData/2D/output001.vtu", geometryWorld, elemUs);
            std::cout << "iterations limit!" << std::endl;
            break;
        }

        double divergence = computeDivergence(elemUs, edgePool);
        std::cout << "Max divergence: " << divergence << std::endl;
    }

    //std::cin.get();
    std::cout << "Final time = " << finalTime << std::endl;
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
        double Bn2 = (edge.neighbourInd2 >= 0) ?
                     elemUs[edge.neighbourInd2][5] * edge.normalVector[0] +
                     elemUs[edge.neighbourInd2][6] * edge.normalVector[1]
                                               : 0.0;
        max_divergence = std::max(max_divergence, std::fabs(Bn1 - Bn2));
    }

    return max_divergence;
}




