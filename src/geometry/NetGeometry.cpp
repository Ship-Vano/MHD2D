//
// Created by Иван on 10/21/2024.
//


#include "NetGeometry.h"

// узел: конструктор по умолчанию
Node::Node(int index, double xCoord, double yCoord, double zCoord)
        : ind(index), x(xCoord), y(yCoord), z(zCoord) {}

// элемент: конструктор по умолчанию
Element::Element(const int index, const std::vector<int> &nIndexes, int size)
        : ind(index), nodeIndexes(nIndexes), dim(size), edgeIndexes(), centroid2D() {
}

// ребро: конструктор по умолчанию
Edge::Edge(int index, int node1, int node2, int neighbor1, int neighbor2, double len, const std::vector<double>& normalVec, const std::vector<double>& midP)
        : ind(index), nodeInd1(node1), nodeInd2(node2),
          neighbourInd1(neighbor1), neighbourInd2(neighbor2), length(len), normalVector(normalVec), midPoint(midP) {}

// подсчёт площади элемента
double areaCalc(const Element& poly, const NodePool& nPool) {
    int dim = poly.dim;

    // получаем узлы элемента
    std::vector<Node> polyNodes;
    for(int i = 0; i < dim; ++i){
        polyNodes.push_back(nPool.nodes[poly.nodeIndexes[i]]);
    }

    // вычисляем площадь элемента
    double res = 0.0;
    double iSum = 0.0;
    double jSum = 0.0;
    double kSum = 0.0;
    for(int k = 1; k < dim-1; ++k){
        double yk_y1 = polyNodes[k].y - polyNodes[0].y;  //y_k - y_1
        double zk1_z1 = polyNodes[k+1].z - polyNodes[0].z; // z_{k+1} - z_1
        double zk_z1 = polyNodes[k].z - polyNodes[0].z;
        double yk1_y1 = polyNodes[k+1].y - polyNodes[0].y;
        double xk1_x1 = polyNodes[k+1].x - polyNodes[0].x;
        double xk_x1 = polyNodes[k].x - polyNodes[0].x;
        iSum += (yk_y1 * zk1_z1 - zk_z1 * yk1_y1);
        jSum += (zk_z1 * xk1_x1 - xk_x1 * zk1_z1);
        kSum += (xk_x1 * yk1_y1 - yk_y1 * xk1_x1);
    }
    res = 0.5 * std::sqrt(iSum*iSum + jSum*jSum + kSum*kSum);
    return res;
}

// точка середины элемента
std::vector<double> getElementCentroid2D(const Element &poly, const NodePool &nPool) {
    int dim = poly.dim;
    std::vector<double> centroid(2, 0.0);
    for(int i = 0; i < dim; ++i){
        Node node = nPool.getNode(poly.nodeIndexes[i]);
        centroid[0] += node.x;
        centroid[1] += node.y;
    }
    centroid[0] /= dim;
    centroid[1] /= dim;
    return centroid;
}

std::vector<double> getElementCentroid2D(const std::vector<Node> &nodes) {
    int dim = nodes.size();
    std::vector<double> centroid(2, 0.0);
    for(int i = 0; i < dim; ++i){
        Node node = nodes[i];
        centroid[0] += node.x;
        centroid[1] += node.y;
    }
    centroid[0] /= dim;
    centroid[1] /= dim;
    return centroid;
}

// точка середины отрезка между двумя узлами
std::vector<double> getMidPoint2D(const int nodeInd1, const int nodeInd2, const NodePool &nPool) {
    Node node1 = nPool.getNode(nodeInd1);
    Node node2 = nPool.getNode(nodeInd2);
    std::vector<double> mid{(node1.x + node2.x)/2.0, (node1.y + node2.y)/2.0};
    return mid;
}

// расстояние между узлами
double getDistance(const int nodeInd1, const int nodeInd2, const NodePool& nPool){
    Node node1 = nPool.getNode(nodeInd1);
    Node node2 = nPool.getNode(nodeInd2);
    return std::sqrt( (node1.x - node2.x)*(node1.x - node2.x) + (node1.y - node2.y)*(node1.y - node2.y) + (node1.z - node2.z)*(node1.z - node2.z) );
}

// набор узлов: конструктор по умолчанию
NodePool::NodePool(int size, const std::vector<Node>& nodeVec)
        : nodeCount(size), nodes(nodeVec) {}

// получить узел из набора по индексу
Node NodePool::getNode(int ind) const{
    return nodes[ind];
}

// набор рёбер: конструктор из вектора рёбер
EdgePool::EdgePool(int size, const std::vector<Edge>& edgeVec)
        : edgeCount(size), edges(edgeVec) {
    minEdgeLen = edgeVec[0].length;
    for (const auto& edge : edgeVec) {
        if(edge.length < minEdgeLen && edge.length > 1e-16){
            minEdgeLen = edge.length;
        }
    }
}


// нормаль к вектору между двумя узлами
std::vector<double> calculateNormalVector2D(const Node& node1, const Node& node2) {
    // координаты вектора
    double dx = node2.x - node1.x;
    double dy = node2.y - node1.y;

    // нормаль к вектору (поворот на 90 градусов против часовой стрелки)
    std::vector<double> normal = {-dy, dx};

    // нормировка
    double length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
    if (length > 0) {
        normal[0] /= length;
        normal[1] /= length;
    }
    return normal;
}

                                                                                                                                                                                        ////ЖЁТСКИЙ КОСТЫЛЬ, УБЕРУ ПОТОМ , КОГДА БУДЕТ МЕСТО НА BOOST
    // хэширование пары (для создания map из пар)
    namespace std {
        template <typename T1, typename T2>
        struct hash<std::pair<T1, T2>> {
            size_t operator()(const std::pair<T1, T2>& p) const {
                auto h1 = std::hash<T1>{}(p.first);
                auto h2 = std::hash<T2>{}(p.second);
                return h1 ^ (h2 << 1); // Combine the two hash values
            }
        };
    }

// набор рёбер: конструктор из набора узлов и набора элементов
EdgePool::EdgePool(const NodePool& np, ElementPool& ep) {
    std::unordered_map<std::pair<int, int>, std::unordered_set<int>, std::hash<std::pair<int, int>>> edgeMap;
    int edgeIndex = 0;
    minEdgeLen = 10000000.0;

    // проходимся по каждому элементу
    for (const auto& element : ep.elements) {
        int dim = element.dim;
        /*if(dim != 3){
            std::cout << "dim != 3" << std::endl;
        }*/
        for (int i = 0; i < dim; ++i) {
            int node1 = element.nodeIndexes[i];
            int node2 = element.nodeIndexes[(i + 1) % dim]; // для цикличности связей (0->1; 1->2; 2->0)
            if (node1 > node2) std::swap(node1, node2); // упорядочивание по возрастанию
            auto edgeKey = std::make_pair(node1, node2); // ключ для ребра (по двум узлам находим элементы-соседей)
            edgeMap[edgeKey].insert(element.ind); // вставляем элемент в map
        }
    }

    // проходимся по map из рёбер
    for (const auto& edgeEntry : edgeMap) {
        int node1 = edgeEntry.first.first; // первый узел ребра
        int node2 = edgeEntry.first.second; // второй узел ребра
        const auto& neighbors = edgeEntry.second; // соседние элементы
        int neighbor1 = -1, neighbor2 = -1; // индексы соседних элементов

        //проходимся по set из элементов, получаем индексы (если они есть)
        auto it = neighbors.begin();
        if (it != neighbors.end()) {
            neighbor1 = *it;
            ++it;
            if (it != neighbors.end()) {
                neighbor2 = *it;
            }
        }

        // делаем индекс первого соседа граничного ребра существующим (2-й сосед = -1, тк его нет)
        if(neighbor1 == -1){
            std::swap(neighbor1,neighbor2);
        }

        // нормаль к ребру
        std::vector<double> normalVector = calculateNormalVector2D(np.getNode(node1), np.getNode(node2));

        // середина ребра
        std::vector<double> edgeMid = getMidPoint2D(node1, node2, np);
        //std::cout << "CENTROID OF THE EL numb " << neighbor1 << " , centroid = (" << getElementCentroid2D(ep.elements[neighbor1], np)[0] << ", " <<getElementCentroid2D(ep.elements[neighbor1], np)[1]  <<") " << std::endl;

        // вектор, соединяющий середину элемента и середину ребра
        std::vector<double> neighbour1ToEdgeMidVector =
                    edgeMid - ep.elements[neighbor1].centroid2D;

        // ориентвция ребра
        int orientation = 1;
        orientation = normalVector * neighbour1ToEdgeMidVector > 0 ? 1 : -1; // скалярное произведение нормали и вектора центр-середина
                                                                                                            //     nei1->nei2
        // корректируем нормаль  (нормаль будет идти от первого соседа ко второму)                                  <| -> |>
        if(orientation < 0){
            normalVector = normalVector * (-1);
            orientation = 1;
        }

        // меняем порядок узлов, составляющих ребро для обеспечения правильной ориентации
        const auto &elementNodes = ep.elements[neighbor1].nodeIndexes; // узлы первого соседа
        auto itNode1 = std::find(elementNodes.begin(), elementNodes.end(), node1); // позиция первого узла в векторе индексов элемента (первый сосед)
        auto itNode2 = std::find(elementNodes.begin(), elementNodes.end(), node2); // позиция второго узла в векторе индексов элемента (первый сосед)

        // индексы узлов в элементе пронумерованы против часовой стрелки.
        // делаем проверку: первый узел ребра идёт после второго узла ребра в элементе (первом соседе), то меняем их местами в ребре
        if ((itNode1 > itNode2) && !(itNode1 == elementNodes.end()-1 && itNode2 == elementNodes.begin()) ) {
            std::swap(node1, node2); // Swap nodes to ensure counterclockwise order
            //std::cout << "swap1 EdgeIndex is " << edgeIndex << std::endl;
        }
        else if(itNode2 == elementNodes.end()-1 && itNode1 == elementNodes.begin()){
            std::swap(node1, node2);
            //std::cout << "swap2 EdgeIndex is " << edgeIndex << std::endl;
        }

        // вычисляем длину ребра
        double len = getDistance(node1, node2, np);
        if(len < minEdgeLen && len > 1e-16){
            minEdgeLen = len;
        }

        // создаём ребро и помещаем в общий набор
        edges.emplace_back(edgeIndex, node1, node2, neighbor1, neighbor2, len, normalVector, edgeMid);

        ++edgeIndex;
    }

    edgeCount = edges.size();  // задаём общее количество рёбер
}

// набор элементов: конструктор из вектора элементов
ElementPool::ElementPool(int nodesPerElement, int elCnt, const std::vector<Element>& elems)
        : elCount(elCnt), isSquare(nodesPerElement == SQUARE_ELEMENT_NODE_COUNT),
          isTriangular(nodesPerElement == TRIANGULAR_ELEMENT_NODE_COUNT), elements(elems) {
}

void World::setNodePool(const NodePool& np) {
    this->np = np;
}

NodePool World::getNodePool() const {
    return np;
}

void World::setElementPool(const ElementPool& ep) {
    this->ep = ep;
}

ElementPool World::getElementPool() const {
    return ep;
}


/*NEIGHBOUR SERVICE*/
NeighbourService::NeighbourService(const NodePool& np, const ElementPool& ep, const EdgePool& edgePool) {
    // выделяем место под maps
    nodeToElements.reserve(np.nodeCount);
    edgeToElements.reserve(edgePool.edgeCount);
    elementToElements.reserve(ep.elCount);
    nodeToEdgesMap.reserve(np.nodeCount);
    boundaryNodeLeftToRight.reserve(np.nodeCount / 4);
    boundaryNodeTopToBottom.reserve(np.nodeCount / 4);

    // map узел -> соседи-элементы
#pragma omp parallel for
    for (size_t i = 0; i < ep.elements.size(); ++i) {
        const auto& element = ep.elements[i];

        // thread-safe update
        for (int nodeIndex : element.nodeIndexes) {
#pragma omp critical(nodeToElements)
            {
                nodeToElements[nodeIndex].insert(element.ind);
            }
        }

    }

    // map узел -> соседи-рёбра
#pragma omp parallel for
    for (size_t edgeIndex = 0; edgeIndex < edgePool.edges.size(); ++edgeIndex) {
        const auto& edge = edgePool.edges[edgeIndex];
#pragma omp critical(nodeToEdgesMap)
        {
            nodeToEdgesMap[edge.nodeInd1].push_back(edgeIndex);
            nodeToEdgesMap[edge.nodeInd2].push_back(edgeIndex);
        }
    }
}

// найти ребро по индексам узлов
int NeighbourService::findEdgeByNodes(int node1Index, int node2Index, const EdgePool& edgePool) const {
    for (int i = 0; i < edgePool.edges.size(); ++i) {
        const auto& edge = edgePool.edges[i];
        if ((edge.nodeInd1 == node1Index && edge.nodeInd2 == node2Index) ||
            (edge.nodeInd1 == node2Index && edge.nodeInd2 == node1Index)) {
            return i; // индекс подходящего ребра
        }
    }
    return -1; // индекс -1 при отсутствии результата поиска
}

// получить рёбра, выходящие из данного узла (по индексу узла)
std::vector<int> NeighbourService::getEdgeNeighborsOfNode(int nodeIndex) const {
    auto it = nodeToEdgesMap.find(nodeIndex);
    if (it != nodeToEdgesMap.end()) {
        return it->second; // вовзращаем рёбра
    }
    return {};
}

std::unordered_set<int> NeighbourService::getNodeNeighbours(int nodeIndex) const {
    return nodeToElements.at(nodeIndex);
}

std::unordered_set<int> NeighbourService::getEdgeNeighbours(int edgeIndex) const {
    return edgeToElements.at(edgeIndex);
}

std::unordered_set<int> NeighbourService::getElementNeighbours(int elementIndex) const {
    return elementToElements.at(elementIndex);
}

//получить рёбра элемента
std::unordered_set<int> NeighbourService::getEdgesOfElement(int elementIndex) const {
    return elementToEdges.at(elementIndex);
}

// получить соседей ребра (элементы)
std::unordered_set<int> NeighbourService::getElementsOfEdge(int edgeIndex) const {
    return edgeToElementsMap.at(edgeIndex);
}


// печать связей-соседей
void NeighbourService::displayNeighbours() const {
    std::cout << "\n--- Neighbours ---" << std::endl;

    // Node to Element Neighbours
    std::cout << "\nNode to Element Neighbours:" << std::endl;
    for (const auto& [node, elements] : nodeToElements) {
        std::cout << "Node " << node << ": ";
        for (int elem : elements) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    // Element to Element Neighbours
    std::cout << "\nElement to Element Neighbours:" << std::endl;
    for (const auto& [element, neighbours] : elementToElements) {
        std::cout << "Element " << element << ": ";
        for (int neigh : neighbours) {
            std::cout << neigh << " ";
        }
        std::cout << std::endl;
    }

    // Edge to Element Neighbours
    std::cout << "\nEdge to Element Neighbours:" << std::endl;
    for (const auto& [edge, elements] : edgeToElements) {
        std::cout << "Edge " << edge << ": ";
        for (int elem : elements) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    // Node to Edges
    for(const auto& [node, edges]: nodeToEdgesMap) {
        std::cout << "Node " << node << " is connected to edges: ";
        for(const auto& edge: edges) {
            std::cout << edge << " ";
        }
        std::cout << std::endl;
    }
}


/*WORLD*/
NeighbourService& World::getNeighbourService() {
    return ns;
}

World::World(const std::string &fileName, const bool isRenderedBin) : np(), ep(), edgp(), ns(np, ep, edgp) {

    // Готовый пререндеренный мир (всё посчитано и записано в bin-файле)
    if(isRenderedBin) {
        importFromFile(fileName);
    }
    else {
        std::ifstream file(fileName);
        assert(file.is_open());

        std::vector<Node> nodes;
        std::vector<Element> elements;
        std::string tmp_line;

        while (std::getline(file, tmp_line)) {
            if (tmp_line == "$Nodes") {
                std::getline(file, tmp_line);  // размер
                std::getline(file, tmp_line);  // первое вхождение
                while (tmp_line != "$EndNodes") {
                    int ind;
                    double x, y, z;
                    std::istringstream ss(tmp_line);
                    ss >> ind >> x >> y >> z;  // считываем значения узла
                    nodes.emplace_back(ind - 1, x, y, z);
                    std::getline(file, tmp_line);
                }
            } else if (tmp_line == "$Elements") {
                std::getline(file, tmp_line); // размер
                std::getline(file, tmp_line); // первое вхождение
                while (tmp_line != "$EndElements") {
                    int ind, count;
                    std::istringstream ss(tmp_line);
                    ss >> ind >> count; // индекс и количество узлов

                    std::vector<int> indexes(count);
                    for (int i = 0; i < count; ++i) {
                        ss >> indexes[i]; // считываем индексы узлов
                        indexes[i] -= 1;
                    }
                    elements.emplace_back(ind-1, indexes, count);
                    std::getline(file, tmp_line);
                }
            }
        }

        std::cout << "Creating node and element pools..." << std::endl;
        np = NodePool(nodes.size(), nodes);
        ep = ElementPool(nodes[0].ind, elements.size(), elements); // предполагаем, что вс элементы имеют одинаковую размерность

        // вычисляем площади и центры для элементов
        for (int i = 0; i < ep.elCount; ++i) {
            ep.elements[i].area = areaCalc(ep.elements[i], np);
            ep.elements[i].centroid2D = getElementCentroid2D(ep.elements[i], np);
        }

        std::cout << "Creating edgepool..." << std::endl;
        edgp = EdgePool(np, ep);  // создаём набор рёбер
        std::cout << "Creating neighbour service..." << std::endl;
        ns = NeighbourService(np, ep, edgp); // создаём сервис соседей

        // индексы рёбер, образующих элемент
        for (auto &element: ep.elements) {
            for (int i = 0; i < element.dim; ++i) {
                int nodeInd = element.nodeIndexes[i];
                int nodeInd_after = element.nodeIndexes[(i + 1) % element.dim];
                int edgeInd = ns.findEdgeByNodes(nodeInd, nodeInd_after, edgp);
                if (edgeInd > -1) {
                    element.edgeIndexes.push_back(edgeInd);
                }
            }
        }

       // связь ребро - соседи
#pragma omp parallel for
        for (int i = 0; i < ep.elements.size(); ++i) {
            const auto &element = ep.elements[i];
            for (int edgeIndex: element.edgeIndexes) {
#pragma omp critical(edgeToElements)
                {
                    ns.edgeToElements[edgeIndex].insert(element.ind);
                }
            }
        }

        // связь элемент - элемент
#pragma omp parallel for
        for (int i = 0; i < ep.elCount; ++i) {
            const auto& element = ep.elements[i];
            std::unordered_set<int> connectedElements;

            for (int edgeIndex : element.edgeIndexes) {
                connectedElements.insert(ns.edgeToElements[edgeIndex].begin(),ns.edgeToElements[edgeIndex].end());
            }
            connectedElements.erase(connectedElements.find(element.ind)); // убираем сам элемент
            // thread-safe
#pragma omp critical(elementToElements)
            {
                ns.elementToElements[element.ind] = std::move(connectedElements);
            }
        }

        // анализ границ
        std::cout << "Analyzing boundaries...." << std::endl;
        minX = np.getNode(0).x;
        maxX = np.getNode(0).x;
        maxY = np.getNode(0).y;
        minY = np.getNode(0).y;
        // определение [minX, maxX] x [minY, maxY]
        for (const auto& node: np.nodes) {
            if(node.x > maxX){
                maxX = node.x;
            }
            else if(node.x < minX){
                maxX = node.x;
            }
            if(node.y > maxY){
                maxY = node.y;
            }
            else if(node.y < minY) {
                minY = node.y;
            }
        }

        // собираем граничные узлы
        for(const auto& node: np.nodes){
            if(std::abs(node.x - maxX) < 1e-15){
                boundaryRightNodes.push_back(node.ind);
            }
            else if(std::abs(node.x - minX) < 1e-15){
                boundaryLeftNodes.push_back(node.ind);
            }
            if(std::abs(node.y - maxY) < 1e-15){
                boundaryTopNodes.push_back(node.ind);
            }
            else if(std::abs(node.y - minY) < 1e-15){
                boundaryBottomNodes.push_back(node.ind);
            }
        }

        // граничные соответствия между узлами
        if(boundaryTopNodes.size() == boundaryBottomNodes.size()){
            std::cout << "Top corresponds to Bot"<< std::endl;
            for(const auto& indTop: boundaryTopNodes){
                const auto nodeTop = np.getNode(indTop);
                ns.boundaryNodeTopToBottom[nodeTop.ind] = -1;
                for(const auto& indBot: boundaryBottomNodes){
                    const auto nodeBot = np.getNode(indBot);
                    if(std::abs(nodeTop.x - nodeBot.x) < 1e-10){
                        ns.boundaryNodeTopToBottom[nodeTop.ind] = nodeBot.ind;
                    }
                }
                if( ns.boundaryNodeTopToBottom[nodeTop.ind] == -1){
                    std::cerr << "Couldn't find the correspondig node!" << std::endl;
                }
            }
        }
        else{
            std::cerr << "Top doesn't correspond to Bot"<< std::endl;
        }
        if(boundaryLeftNodes.size() == boundaryRightNodes.size()){
            std::cout << "Left corresponds to Right"<< std::endl;
            for(const auto& indLeft: boundaryLeftNodes){
                const auto nodeLeft = np.getNode(indLeft);
                ns.boundaryNodeLeftToRight[nodeLeft.ind] = -1;
                for(const auto& indRight: boundaryRightNodes){
                    const auto nodeRight = np.getNode(indRight);
                    if(std::abs(nodeLeft.y - nodeRight.y) < 1e-8){
                        ns.boundaryNodeLeftToRight[nodeLeft.ind] = nodeRight.ind;
                    }
                }
                if( ns.boundaryNodeLeftToRight[nodeLeft.ind] == -1){
                    std::cerr << "Couldn't find the correspondig node!" << std::endl;
                }
            }
        }
        else{
            std::cerr << "Left doesn't correspond to Right!"<< std::endl;
        }

        //gather boundary edges and elems
        for(const auto& edge: edgp.edges){
            if(edge.neighbourInd2 == -1) {
                ep.elements[edge.neighbourInd1].is_boundary = true;
                if (ns.boundaryNodeTopToBottom.count(edge.nodeInd1) == 1  && \
                    ns.boundaryNodeTopToBottom.count(edge.nodeInd2) == 1){
                    //соответствующее ребро
                    ns.boundaryEdgeTopToBottom[edge.ind] = -1;
                    for(const auto& ind1:ns.nodeToEdgesMap[ns.boundaryNodeTopToBottom[edge.nodeInd1]]){
                        for(const auto& ind2:ns.nodeToEdgesMap[ns.boundaryNodeTopToBottom[edge.nodeInd2]]){
                            if(ind1 == ind2){
                                ns.boundaryEdgeTopToBottom[edge.ind] = ind1;
                                ns.boundaryElemTopToBottom[edge.neighbourInd1] = ep.elements[edgp.edges[ind1].neighbourInd1].ind;
                            }
                        }
                    }
                    if(ns.boundaryEdgeTopToBottom[edge.ind] == -1){
                        std::cerr << "Couldn't find a corresponding edge top-bot!" << std::endl;
                    }
                }
                else if(ns.boundaryNodeLeftToRight.count(edge.nodeInd1) == 1 &&  \
                    ns.boundaryNodeLeftToRight.count(edge.nodeInd2) == 1){
                    ns.boundaryEdgeLeftToRight[edge.ind] = -1;
                    for(const auto& ind1:ns.nodeToEdgesMap[ns.boundaryNodeLeftToRight[edge.nodeInd1]]){
                        for(const auto& ind2:ns.nodeToEdgesMap[ns.boundaryNodeLeftToRight[edge.nodeInd2]]){
                            if(ind1 == ind2){
                                ns.boundaryEdgeLeftToRight[edge.ind] = ind1;
                                ns.boundaryElemLeftToRight[edge.neighbourInd1] = ep.elements[edgp.edges[ind1].neighbourInd1].ind;
                            }
                        }
                    }
                    if(ns.boundaryEdgeLeftToRight[edge.ind] == -1){
                        std::cerr << "Couldn't find a corresponding edge left-right!" << std::endl;
                    }
                }
            }
        }

        // ghost cells generation
        int ghostElemInd = ep.elCount;
        int ghostNodeInd = np.nodeCount;
        int ghostEdgeInd = edgp.edgeCount;
        std::vector<Node> ghostNodes;           // фантомные узлы
        std::vector<Element> ghostElements;        // фантомные ячейки
        std::vector<Edge> ghostEdges;         // фантомные ребра
        for(const auto& elem: ep.elements){
            // определяем граничный элемент
            if(elem.is_boundary){
                // определяем граничное ребро
                for(const auto& edgeInd: elem.edgeIndexes){
                    Edge edge = edgp.edges[edgeInd];
                    if(edge.neighbourInd2 == -1){
                        // отражаем 3й узел
                        const auto nodeInElemInd = std::find(elem.nodeIndexes.begin(), elem.nodeIndexes.end(),
                                                             edge.nodeInd1);
                        int node_before_ind =
                                nodeInElemInd == elem.nodeIndexes.begin() ? elem.nodeIndexes[elem.dim - 1] : *(
                                        nodeInElemInd - 1);
                        Node node1st = np.getNode(edge.nodeInd1);
                        Node node2nd = np.getNode(edge.nodeInd2);
                        Node node3rd = np.getNode(node_before_ind);
                        std::vector<double> reflPoint = reflectNodeOverVector(node3rd, node1st, node2nd);
                        Node reflNode(ghostNodeInd, reflPoint[0], reflPoint[1], 0.0);
                        reflNode.is_ghost = true;

                        // строим элемент
                        std::vector<int> ghostElNodes{node2nd.ind, node1st.ind, ghostNodeInd};
                        if(!isCounterClockwise(ghostElNodes)){
                            std::swap(ghostElNodes[0], ghostElNodes[1]);
                        }
                        Element ghostEl(ghostElemInd, ghostElNodes, 3);
                        ghostEl.is_ghost = true;

                        // строим рёбра
                        Edge ghostEdge1(ghostEdgeInd, edge.nodeInd1, ghostNodeInd, ghostElemInd, -1,
                                        getDistance(edge.nodeInd1, ghostNodeInd, np),
                                        calculateNormalVector2D(node1st, reflNode),
                                        getMidPoint2D(edge.nodeInd1, ghostNodeInd, np));
                        if(std::abs(1.0-std::sqrt(ghostEdge1.normalVector[0]*ghostEdge1.normalVector[0] + ghostEdge1.normalVector[1]*ghostEdge1.normalVector[1])) > 1e-15){
                            std::cout << "bad normal while creating a ghost edge! normal = { " << ghostEdge1.normalVector[0] << " , " << ghostEdge1.normalVector[1] << " }"<<std::endl;
                        }

                        ghostEdge1.is_ghost = true;
                        Edge ghostEdge2(ghostEdgeInd + 1, ghostNodeInd, edge.nodeInd2, ghostElemInd, -1,
                                        getDistance(ghostNodeInd, edge.nodeInd2, np),
                                        calculateNormalVector2D(reflNode, node2nd),
                                        getMidPoint2D(ghostNodeInd, edge.nodeInd2, np));
                        if(std::abs(1.0-std::sqrt(ghostEdge2.normalVector[0]*ghostEdge2.normalVector[0] + ghostEdge2.normalVector[1]*ghostEdge2.normalVector[1])) > 1e-15){
                            std::cout << "bad normal while creating a ghost edge!" <<std::endl;
                        }

                        ghostEdge2.is_ghost = true;
                        ghostEl.edgeIndexes = std::vector<int>{edge.ind, ghostEdgeInd, ghostEdgeInd + 1};
                        ghostEl.area = elem.area;
                        ghostEl.centroid2D = getElementCentroid2D({node2nd, node1st, reflNode});

                        // neighbour service connection
                        ns.boundaryToGhostElements[elem.ind] = ghostEl.ind;

                        // Update indices
                        ghostNodes.push_back(reflNode);
                        ghostElements.push_back(ghostEl);
                        ghostEdges.push_back(ghostEdge1);
                        ghostEdges.push_back(ghostEdge2);
                        ++ghostElemInd;
                        ++ghostNodeInd;
                        ghostEdgeInd += 2;
                    }
                }

            }
        }

        // Append ghost nodes, elements, and edges to pools
        np.nodes.insert(np.nodes.end(), ghostNodes.begin(), ghostNodes.end());
        np.nodeCount = np.nodes.size();
        ep.elements.insert(ep.elements.end(), ghostElements.begin(), ghostElements.end());
        ep.elCount = ep.elements.size();
        edgp.edges.insert(edgp.edges.end(), ghostEdges.begin(), ghostEdges.end());
        edgp.edgeCount = edgp.edges.size();

    } /*end: 2d variant constructor*/
}

void World::display() const {
    // Display Node Pool
    std::cout << "Node Pool:" << std::endl;
    std::cout << "Total Nodes: " << np.nodeCount << std::endl;
    for (const auto& node : np.nodes) {
        std::cout << "Node Index: " << node.ind << ", Coordinates: ("
                  << node.x << ", " << node.y << ", " << node.z << ")" << std::endl;
    }

    // Display Element Pool
    std::cout << "\nElement Pool:" << std::endl;
    std::cout << "Total Elements: " << ep.elCount << std::endl;
    for (const auto& element : ep.elements) {
        std::cout << "Element Index: " << element.ind << ", Node Count: " << element.dim
                  << ", Node Indices: ";
        for (const auto& index : element.nodeIndexes) {
            std::cout << index << " ";
        }
        std::cout << ", Area: " << element.area << ", Edges: (";
        for(const auto& index: element.edgeIndexes){
            std::cout << index << " ";
        }
        std::cout << ")"<<std::endl;
    }

    // Display Edge Pool and Edge Neighbors
    std::cout << "\nEdge Pool:" << std::endl;
    std::cout << "Total Edges: " << edgp.edgeCount << std::endl;
    for (const auto& edge : edgp.edges) {
        std::cout << "Edge Index: " << edge.ind << ", Nodes: ("
                  << edge.nodeInd1 << ", " << edge.nodeInd2 << "), "
                  << "Neighbors: (" << edge.neighbourInd1 << ", " << edge.neighbourInd2 << ")"
                  << ", Normal: ("
                  << edge.normalVector[0] << ", " << edge.normalVector[1] << "), len = "
                  << edge.length << ", MidPoint = ("<< edge.midPoint[0] << ", "<< edge.midPoint[1] << ")"<<std::endl;
    }

    // Display Neighbor Information for Nodes, Edges, and Elements
    ns.displayNeighbours();
}

EdgePool World::getEdgePool() const {
    return edgp;
}

void World::setEdgePool(const EdgePool &edgp) {
    this->edgp = edgp;
}

void World::exportToFile(const string &filename) const{
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    // Write NodePool
    int nodeCount = np.nodeCount;
    file.write(reinterpret_cast<const char*>(&nodeCount), sizeof(nodeCount));
    for (const auto& node : np.nodes) {
        file.write(reinterpret_cast<const char*>(&node.ind), sizeof(node.ind));
        file.write(reinterpret_cast<const char*>(&node.x), sizeof(node.x));
        file.write(reinterpret_cast<const char*>(&node.y), sizeof(node.y));
        file.write(reinterpret_cast<const char*>(&node.z), sizeof(node.z));
        file.write(reinterpret_cast<const char*>(&node.is_ghost), sizeof(node.is_ghost));
    }

    // Write ElementPool
    int elementCount = ep.elCount;
    file.write(reinterpret_cast<const char*>(&elementCount), sizeof(elementCount));
    for (const auto& element : ep.elements) {
        file.write(reinterpret_cast<const char*>(&element.ind), sizeof(element.ind));
        int dim = element.dim;
        file.write(reinterpret_cast<const char*>(&dim), sizeof(dim));
        for (const auto& nodeIndex : element.nodeIndexes) {
            file.write(reinterpret_cast<const char*>(&nodeIndex), sizeof(nodeIndex));
        }
        if(element.dim != element.edgeIndexes.size()){
            std::cout << element.dim << " vs " << element.edgeIndexes.size() << std::endl;
            std::cout << element.ind << " element has area = " << element.area << std::endl;
            std::cout << "Elem's nodes are: ";
            for(const auto& ind: element.nodeIndexes){
                Node curnode = np.getNode(ind);
                std::cout << ind << "( "<< curnode.x << ", " << curnode.y << ", " << curnode.z << " ); ";
            }

            std::cout << std::endl;
        }
        for(const auto& edgeIndex: element.edgeIndexes){
            file.write(reinterpret_cast<const char*>(&edgeIndex), sizeof(edgeIndex));
        }
        double area = element.area;
        file.write(reinterpret_cast<const char*>(&area), sizeof(area));
        for (const auto& coord : element.centroid2D) {
            file.write(reinterpret_cast<const char*>(&coord), sizeof(coord));
        }
        bool is_boundary = element.is_boundary;
        file.write(reinterpret_cast<const char*>(&is_boundary), sizeof(is_boundary));
        bool is_ghost = element.is_ghost;
        file.write(reinterpret_cast<const char*>(&is_ghost), sizeof(is_ghost));
    }

    // Write EdgePool
    int edgeCount = edgp.edgeCount;
    file.write(reinterpret_cast<const char*>(&edgeCount), sizeof(edgeCount));
    for (const auto& edge : edgp.edges) {
        file.write(reinterpret_cast<const char*>(&edge.ind), sizeof(edge.ind));
        file.write(reinterpret_cast<const char*>(&edge.nodeInd1), sizeof(edge.nodeInd1));
        file.write(reinterpret_cast<const char*>(&edge.nodeInd2), sizeof(edge.nodeInd2));
        file.write(reinterpret_cast<const char*>(&edge.neighbourInd1), sizeof(edge.neighbourInd1));
        file.write(reinterpret_cast<const char*>(&edge.neighbourInd2), sizeof(edge.neighbourInd2));
        file.write(reinterpret_cast<const char*>(&edge.length), sizeof(edge.length));
        for (const auto& val : edge.normalVector) {
            file.write(reinterpret_cast<const char*>(&val), sizeof(val));
        }
        for (const auto& val : edge.midPoint) {
            file.write(reinterpret_cast<const char*>(&val), sizeof(val));
        }
        file.write(reinterpret_cast<const char*>(&edge.is_ghost), sizeof(edge.is_ghost));
    }

    // Serialize NeighbourService
    // Export node-to-elements map
    int nodeToElementsSize = ns.nodeToElements.size();
    file.write(reinterpret_cast<const char*>(&nodeToElementsSize), sizeof(nodeToElementsSize));
    for (const auto& [node, elements] : ns.nodeToElements) {
        file.write(reinterpret_cast<const char*>(&node), sizeof(node));
        int elementCount = elements.size();
        file.write(reinterpret_cast<const char*>(&elementCount), sizeof(elementCount));
        for (int elem : elements) {
            file.write(reinterpret_cast<const char*>(&elem), sizeof(elem));
        }
    }

    // Export element-to-elements map
    int elementToElementsSize = ns.elementToElements.size();
    file.write(reinterpret_cast<const char*>(&elementToElementsSize), sizeof(elementToElementsSize));
    for (const auto& [element, neighbors] : ns.elementToElements) {
        file.write(reinterpret_cast<const char*>(&element), sizeof(element));
        int neighborCount = neighbors.size();
        file.write(reinterpret_cast<const char*>(&neighborCount), sizeof(neighborCount));
        for (int neighbor : neighbors) {
            file.write(reinterpret_cast<const char*>(&neighbor), sizeof(neighbor));
        }
    }

    // Export edge-to-elements map
    int edgeToElementsSize = ns.edgeToElements.size();
    file.write(reinterpret_cast<const char*>(&edgeToElementsSize), sizeof(edgeToElementsSize));
    for (const auto& [edge, elements] : ns.edgeToElements) {
        file.write(reinterpret_cast<const char*>(&edge), sizeof(edge));
        int elementCount = elements.size();
        file.write(reinterpret_cast<const char*>(&elementCount), sizeof(elementCount));
        for (int elem : elements) {
            file.write(reinterpret_cast<const char*>(&elem), sizeof(elem));
        }
    }

    // Export node-to-edges map
    int nodeToEdgesMapSize = ns.nodeToEdgesMap.size();
    file.write(reinterpret_cast<const char*>(&nodeToEdgesMapSize), sizeof(nodeToEdgesMapSize));
    for (const auto& [node, edges] : ns.nodeToEdgesMap) {
        file.write(reinterpret_cast<const char*>(&node), sizeof(node));
        int edgeCount = edges.size();
        file.write(reinterpret_cast<const char*>(&edgeCount), sizeof(edgeCount));
        for (int edge : edges) {
            file.write(reinterpret_cast<const char*>(&edge), sizeof(edge));
        }
    }

    // Export boundaries
    int boundaryNodeLeftToRight_size =  ns.boundaryNodeLeftToRight.size();
    file.write(reinterpret_cast<const char*>(&boundaryNodeLeftToRight_size), sizeof(boundaryNodeLeftToRight_size));
    for (const auto& [nodeLeft, nodeRight] : ns.boundaryNodeLeftToRight) {
        file.write(reinterpret_cast<const char*>(&nodeLeft), sizeof(nodeLeft));
        file.write(reinterpret_cast<const char*>(&nodeRight), sizeof(nodeRight));
    }
    int boundaryNodeTopToBottom_size =  ns.boundaryNodeTopToBottom.size();
    file.write(reinterpret_cast<const char*>(&boundaryNodeTopToBottom_size), sizeof(boundaryNodeTopToBottom_size));
    for (const auto& [nodeTop, nodeBot] : ns.boundaryNodeTopToBottom) {
        file.write(reinterpret_cast<const char*>(&nodeTop), sizeof(nodeTop));
        file.write(reinterpret_cast<const char*>(&nodeBot), sizeof(nodeBot));
    }
    // edges
    int boundaryEdgeLeftToRight_size =  ns.boundaryEdgeLeftToRight.size();
    file.write(reinterpret_cast<const char*>(&boundaryEdgeLeftToRight_size), sizeof(boundaryEdgeLeftToRight_size));
    for (const auto& [edgeLeft, edgeRight] : ns.boundaryEdgeLeftToRight) {
        file.write(reinterpret_cast<const char*>(&edgeLeft), sizeof(edgeLeft));
        file.write(reinterpret_cast<const char*>(&edgeRight), sizeof(edgeRight));
    }
    int boundaryEdgeTopToBottom_size =  ns.boundaryEdgeTopToBottom.size();
    file.write(reinterpret_cast<const char*>(&boundaryEdgeTopToBottom_size), sizeof(boundaryEdgeTopToBottom_size));
    for (const auto& [edgeTop, edgeBot] : ns.boundaryEdgeTopToBottom) {
        file.write(reinterpret_cast<const char*>(&edgeTop), sizeof(edgeTop));
        file.write(reinterpret_cast<const char*>(&edgeBot), sizeof(edgeBot));
    }
    // elems
    int boundaryElemLeftToRight_size =  ns.boundaryElemLeftToRight.size();
    file.write(reinterpret_cast<const char*>(&boundaryElemLeftToRight_size), sizeof(boundaryElemLeftToRight_size));
    for (const auto& [elemLeft, elemRight] : ns.boundaryElemLeftToRight) {
        file.write(reinterpret_cast<const char*>(&elemLeft), sizeof(elemLeft));
        file.write(reinterpret_cast<const char*>(&elemRight), sizeof(elemRight));
    }
    int boundaryElemTopToBottom_size =  ns.boundaryElemTopToBottom.size();
    file.write(reinterpret_cast<const char*>(&boundaryElemTopToBottom_size), sizeof(boundaryElemTopToBottom_size));
    for (const auto& [elemTop, elemBot] : ns.boundaryElemTopToBottom) {
        file.write(reinterpret_cast<const char*>(&elemTop), sizeof(elemTop));
        file.write(reinterpret_cast<const char*>(&elemBot), sizeof(elemBot));
    }

    // Export boundaryELem to ghostELem
    int boundaryToGhost_size = ns.boundaryToGhostElements.size();
    file.write(reinterpret_cast<const char*>(&boundaryToGhost_size), sizeof(boundaryToGhost_size));
    for (const auto& [boundary, ghost] : ns.boundaryToGhostElements) {
        file.write(reinterpret_cast<const char*>(&boundary), sizeof(boundary));
        file.write(reinterpret_cast<const char*>(&ghost), sizeof(ghost));
    }

    file.close();
    std::cout << "World exported to " << filename << std::endl;
}

void World::importFromFile(const string &filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    // Read NodePool
    int nodeCount;
    file.read(reinterpret_cast<char*>(&nodeCount), sizeof(nodeCount));
    std::vector<Node> nodes(nodeCount);
    for (auto& node : nodes) {
        file.read(reinterpret_cast<char*>(&node.ind), sizeof(node.ind));
        file.read(reinterpret_cast<char*>(&node.x), sizeof(node.x));
        file.read(reinterpret_cast<char*>(&node.y), sizeof(node.y));
        file.read(reinterpret_cast<char*>(&node.z), sizeof(node.z));
        file.read(reinterpret_cast<char*>(&node.is_ghost), sizeof(node.is_ghost));
    }
    np = NodePool(nodeCount, nodes);

    // Read ElementPool
    int elementCount;
    file.read(reinterpret_cast<char*>(&elementCount), sizeof(elementCount));
    std::vector<Element> elements(elementCount);
    for (auto& element : elements) {
        file.read(reinterpret_cast<char*>(&element.ind), sizeof(element.ind));
        int dim;
        file.read(reinterpret_cast<char*>(&dim), sizeof(dim));
        element.dim = dim;
        element.nodeIndexes.resize(dim);
        for (auto& nodeIndex : element.nodeIndexes) {
            file.read(reinterpret_cast<char*>(&nodeIndex), sizeof(nodeIndex));
        }
        element.edgeIndexes.resize(dim);
        for (auto& edgeIndex : element.edgeIndexes) {
            file.read(reinterpret_cast<char*>(&edgeIndex), sizeof(edgeIndex));
        }
        file.read(reinterpret_cast<char*>(&element.area), sizeof(element.area));
        element.centroid2D.resize(2);
        for (auto& coord : element.centroid2D) {
            file.read(reinterpret_cast<char*>(&coord), sizeof(coord));
        }
        file.read(reinterpret_cast<char*>(&element.is_boundary), sizeof(element.is_boundary));
        file.read(reinterpret_cast<char*>(&element.is_ghost), sizeof(element.is_ghost));
    }
    ep = ElementPool(elements[0].dim, elementCount, elements);

    // Read EdgePool
    int edgeCount;
    file.read(reinterpret_cast<char*>(&edgeCount), sizeof(edgeCount));
    std::vector<Edge> edges(edgeCount);
    for (auto& edge : edges) {
        file.read(reinterpret_cast<char*>(&edge.ind), sizeof(edge.ind));
        file.read(reinterpret_cast<char*>(&edge.nodeInd1), sizeof(edge.nodeInd1));
        file.read(reinterpret_cast<char*>(&edge.nodeInd2), sizeof(edge.nodeInd2));
        file.read(reinterpret_cast<char*>(&edge.neighbourInd1), sizeof(edge.neighbourInd1));
        file.read(reinterpret_cast<char*>(&edge.neighbourInd2), sizeof(edge.neighbourInd2));
        file.read(reinterpret_cast<char*>(&edge.length), sizeof(edge.length));
        edge.normalVector.resize(2);
        for (auto& val : edge.normalVector) {
            file.read(reinterpret_cast<char*>(&val), sizeof(val));
        }
        edge.midPoint.resize(2);
        for (auto& val : edge.midPoint) {
            file.read(reinterpret_cast<char*>(&val), sizeof(val));
        }
        file.read(reinterpret_cast<char*>(&edge.is_ghost), sizeof(edge.is_ghost));
    }
    edgp = EdgePool(edgeCount, edges);

    // Deserialize NeighbourService
    // Import node-to-elements map
    int nodeToElementsSize;
    file.read(reinterpret_cast<char*>(&nodeToElementsSize), sizeof(nodeToElementsSize));
    ns.nodeToElements.clear();
    for (int i = 0; i < nodeToElementsSize; ++i) {
        int node;
        file.read(reinterpret_cast<char*>(&node), sizeof(node));
        int elementCount;
        file.read(reinterpret_cast<char*>(&elementCount), sizeof(elementCount));
        std::unordered_set<int> elements;
        for (int j = 0; j < elementCount; ++j) {
            int elem;
            file.read(reinterpret_cast<char*>(&elem), sizeof(elem));
            elements.insert(elem);
        }
        ns.nodeToElements[node] = elements;
    }

    // Import element-to-elements map
    int elementToElementsSize;
    file.read(reinterpret_cast<char*>(&elementToElementsSize), sizeof(elementToElementsSize));
    ns.elementToElements.clear();
    for (int i = 0; i < elementToElementsSize; ++i) {
        int element;
        file.read(reinterpret_cast<char*>(&element), sizeof(element));
        int neighborCount;
        file.read(reinterpret_cast<char*>(&neighborCount), sizeof(neighborCount));
        std::unordered_set<int> neighbors;
        for (int j = 0; j < neighborCount; ++j) {
            int neighbor;
            file.read(reinterpret_cast<char*>(&neighbor), sizeof(neighbor));
            neighbors.insert(neighbor);
        }
        ns.elementToElements[element] = neighbors;
    }

    // Import edge-to-elements map
    int edgeToElementsSize;
    file.read(reinterpret_cast<char*>(&edgeToElementsSize), sizeof(edgeToElementsSize));
    ns.edgeToElements.clear();
    for (int i = 0; i < edgeToElementsSize; ++i) {
        int edge;
        file.read(reinterpret_cast<char*>(&edge), sizeof(edge));
        int elementCount;
        file.read(reinterpret_cast<char*>(&elementCount), sizeof(elementCount));
        std::unordered_set<int> elements;
        for (int j = 0; j < elementCount; ++j) {
            int elem;
            file.read(reinterpret_cast<char*>(&elem), sizeof(elem));
            elements.insert(elem);
        }
        ns.edgeToElements[edge] = elements;
    }

    // Import node-to-edges map
    int nodeToEdgesMapSize;
    file.read(reinterpret_cast<char*>(&nodeToEdgesMapSize), sizeof(nodeToEdgesMapSize));
    ns.nodeToEdgesMap.clear();
    for (int i = 0; i < nodeToEdgesMapSize; ++i) {
        int node;
        file.read(reinterpret_cast<char*>(&node), sizeof(node));
        int edgeCount;
        file.read(reinterpret_cast<char*>(&edgeCount), sizeof(edgeCount));
        std::vector<int> edges(edgeCount);
        for (int j = 0; j < edgeCount; ++j) {
            file.read(reinterpret_cast<char*>(&edges[j]), sizeof(edges[j]));
        }
        ns.nodeToEdgesMap[node] = edges;
    }

    // Import boundaries
    int boundaryNodeLeftToRight_size;
    file.read(reinterpret_cast<char*>(&boundaryNodeLeftToRight_size), sizeof(boundaryNodeLeftToRight_size));
    ns.boundaryNodeLeftToRight.clear();
    for (int i = 0; i < boundaryNodeLeftToRight_size; ++i) {
        int nodeLeft;
        int nodeRight;
        file.read(reinterpret_cast<char*>(&nodeLeft), sizeof(nodeLeft));
        file.read(reinterpret_cast<char*>(&nodeRight), sizeof(nodeRight));
        ns.boundaryNodeLeftToRight[nodeLeft] = nodeRight;
    }
    int boundaryNodeTopToBottom_size;
    file.read(reinterpret_cast<char*>(&boundaryNodeTopToBottom_size), sizeof(boundaryNodeTopToBottom_size));
    ns.boundaryNodeTopToBottom.clear();
    for (int i = 0; i < boundaryNodeTopToBottom_size; ++i) {
        int nodeTop;
        int nodeBot;
        file.read(reinterpret_cast<char*>(&nodeTop), sizeof(nodeTop));
        file.read(reinterpret_cast<char*>(&nodeBot), sizeof(nodeBot));
        ns.boundaryNodeTopToBottom[nodeTop] = nodeBot;
    }
    // edges
    int boundaryEdgeLeftToRight_size;
    file.read(reinterpret_cast<char*>(&boundaryEdgeLeftToRight_size), sizeof(boundaryEdgeLeftToRight_size));
    ns.boundaryEdgeLeftToRight.clear();
    for (int i = 0; i < boundaryEdgeLeftToRight_size; ++i) {
        int EdgeLeft;
        int EdgeRight;
        file.read(reinterpret_cast<char*>(&EdgeLeft), sizeof(EdgeLeft));
        file.read(reinterpret_cast<char*>(&EdgeRight), sizeof(EdgeRight));
        ns.boundaryEdgeLeftToRight[EdgeLeft] = EdgeRight;
    }
    int boundaryEdgeTopToBottom_size;
    file.read(reinterpret_cast<char*>(&boundaryEdgeTopToBottom_size), sizeof(boundaryEdgeTopToBottom_size));
    ns.boundaryEdgeTopToBottom.clear();
    for (int i = 0; i < boundaryEdgeTopToBottom_size; ++i) {
        int EdgeTop;
        int EdgeBot;
        file.read(reinterpret_cast<char*>(&EdgeTop), sizeof(EdgeTop));
        file.read(reinterpret_cast<char*>(&EdgeBot), sizeof(EdgeBot));
        ns.boundaryEdgeTopToBottom[EdgeTop] = EdgeBot;
    }
    // elems
    int boundaryElemLeftToRight_size;
    file.read(reinterpret_cast<char*>(&boundaryElemLeftToRight_size), sizeof(boundaryElemLeftToRight_size));
    ns.boundaryElemLeftToRight.clear();
    for (int i = 0; i < boundaryElemLeftToRight_size; ++i) {
        int ElemLeft;
        int ElemRight;
        file.read(reinterpret_cast<char*>(&ElemLeft), sizeof(ElemLeft));
        file.read(reinterpret_cast<char*>(&ElemRight), sizeof(ElemRight));
        ns.boundaryElemLeftToRight[ElemLeft] = ElemRight;
    }
    int boundaryElemTopToBottom_size;
    file.read(reinterpret_cast<char*>(&boundaryElemTopToBottom_size), sizeof(boundaryElemTopToBottom_size));
    ns.boundaryElemTopToBottom.clear();
    for (int i = 0; i < boundaryElemTopToBottom_size; ++i) {
        int ElemTop;
        int ElemBot;
        file.read(reinterpret_cast<char*>(&ElemTop), sizeof(ElemTop));
        file.read(reinterpret_cast<char*>(&ElemBot), sizeof(ElemBot));
        ns.boundaryElemTopToBottom[ElemTop] = ElemBot;
    }

    int boundaryToGhost_size;
    file.read(reinterpret_cast<char*>(&boundaryToGhost_size), sizeof(boundaryToGhost_size));
    ns.boundaryToGhostElements.clear();
    for(int i = 0; i < boundaryToGhost_size; ++i){
        int boundaryElem;
        int ghostElem;
        file.read(reinterpret_cast<char*>(&boundaryElem), sizeof(boundaryElem));
        file.read(reinterpret_cast<char*>(&ghostElem), sizeof(ghostElem));
        ns.boundaryToGhostElements[boundaryElem] = ghostElem;
    }

    file.close();
    std::cout << "World imported from " << filename << std::endl;
}

void setNeighbourEdge(Element& el, const int edgeInd){
    el.nodeIndexes.push_back(edgeInd);
}

std::vector<double> reflectNodeOverVector(const Node& nodeToReflect, const Node& node1, const Node& node2){

    double x1 = node1.x;
    double y1 = node1.y;
    double x2 = node2.x;
    double y2 = node2.y;
    double x3 = nodeToReflect.x; // point to reflect
    double y3 = nodeToReflect.y; // point to reflect

    double ab_x = x2 - x1;
    double ab_y = y2 - y1;
    double ac_x = x3 - x1;
    double ac_y = y3 - y1;

    // AC proj on AB
    double ab_sqr = ab_x * ab_x + ab_y * ab_y;
    double ac_cdot_ab = ac_x * ab_x + ac_y * ab_y;
    double coef = ac_cdot_ab / ab_sqr;
    double proj_x = coef * ab_x;
    double proj_y = coef * ab_y;

    // perpendicular part
    double perp_x = ac_x - proj_x;
    double perp_y = ac_y - proj_y;

    // The reflected vector
    double reflv_x = proj_x - perp_x;
    double reglv_y = proj_y - perp_y;

    return std::vector<double>{ x1 + reflv_x, y1 + reglv_y };

}

bool World::isCounterClockwise(const std::vector<int> &nodeIndices) {
    double sum = 0;
    for (size_t i = 0; i < nodeIndices.size(); ++i) {
        Node n1 = np.getNode(nodeIndices[i]);
        Node n2 = np.getNode(nodeIndices[(i + 1) % nodeIndices.size()]);
        sum += (n2.x - n1.x) * (n2.y + n1.y);
    }
    return sum > 0; // CCW if sum > 0
}