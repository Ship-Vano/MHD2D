//
// Created by Иван on 10/21/2024.
//

#ifndef MAGNETTOPRJCT_NETGEOMETRY_H
#define MAGNETTOPRJCT_NETGEOMETRY_H

#include "NamesNConstants.h"
#include "service/LinOp.h"
#include <omp.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <utility>
#include <iostream>
#include <cmath>
#include <cassert>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

// Узел
class Node {
public:
    int ind; // порядковый номер
    double x;
    double y;
    double z;
    bool is_ghost = false;
    Node() : ind(0), x(0.0), y(0.0), z(0.0) {}
    Node(int index, double xCoord, double yCoord, double zCoord);
};

// Ребро
class Edge{
public:
    int ind; // порядковый номер
    int nodeInd1; // номер первого узла
    int nodeInd2; // номер второго узла
    int neighbourInd1; //номер первого соседнего элемента
    int neighbourInd2; //номер второго соседнего элемента
    double length;  // длина
    bool is_ghost = false;
    std::vector<double> normalVector; // компоненты вектора нормали
    std::vector<double> midPoint; // центр ребра
    Edge(): ind(0), nodeInd1(0), nodeInd2(0),
              neighbourInd1(-1), neighbourInd2(-1), length(0.0),
              normalVector(2, 0.0), midPoint(2, 0.0) {}
    Edge(int index, int node1, int node2, int neighbor1, int neighbor2,
         double len, const std::vector<double>& normalVec, const std::vector<double>& midPoint);
};

// Элемент
class Element {
public:
    int ind; // порядковый номер
    int dim; // размерность (кол-во узлов)
    std::vector<int> nodeIndexes; //номера узлов
    std::vector<int> edgeIndexes; //номера рёбер
    std::vector<double> centroid2D; // центр элемента
    bool is_boundary = false;
    bool is_ghost = false;
    double area = 0.0; //площадь
    Element() : ind(0), dim(0), area(0.0), centroid2D(2, 0.0) {}
    Element(const int index, const std::vector<int> &nIndexes, int size);
};

// Набор узлов
class NodePool {
public:
    int nodeCount; // количество узлов в наборе
    std::vector<Node> nodes; // узлы

    NodePool(int size, const std::vector<Node>& nodeVec);
    NodePool() : nodeCount(0), nodes() {} // Default constructor

    Node getNode(int ind) const; // получить узел по индексу
};

// Набор элементов
class ElementPool {
public:
    int elCount;  // количество элементов в наборе
    bool isSquare; // флаг: черырёхугольные элементы
    bool isTriangular; // флаг: треугольные элементы
    std::vector<Element> elements; // элементы

    ElementPool(int nodesPerElement, int elCnt, const std::vector<Element>& elements);
    ElementPool() : elCount(0), isSquare(false), isTriangular(false), elements() {} // Default constructor
};

// Набор рёбер
class EdgePool{
public:
    int edgeCount; // количество рёбер
    std::vector<Edge> edges; // рёбра
    double minEdgeLen; // минимальная длина ребра
    EdgePool(int size, const std::vector<Edge>& edgeVec);
    EdgePool(const NodePool& np, ElementPool& ep);
    EdgePool() : edgeCount(0), edges(), minEdgeLen(1.0) {} // Default constructor
};

class NeighbourService {
public:
    std::unordered_map<int, std::unordered_set<int>> nodeToElements; // Nodes -> Elements (exp+)
    std::unordered_map<int, std::vector<int>> nodeToEdgesMap;       //Nodes -> Edges  (exp+)
    std::unordered_map<int, std::unordered_set<int>> edgeToElements; // Edges -> Elements (exp+)
    std::unordered_map<int, std::unordered_set<int>> elementToElements; // Elements -> Elements (exp+)
    std::unordered_map<int, std::unordered_set<int>> elementToEdges; // Elements -> Edges
    std::unordered_map<int, std::unordered_set<int>> edgeToElementsMap; // Edges -> Elements
    std::unordered_map<int, int> boundaryNodeLeftToRight; // left boundary nodes -> right boundary nodes
    std::unordered_map<int, int> boundaryNodeTopToBottom; // top boundary nodes -> bot boundary nodes
    std::unordered_map<int, int> boundaryEdgeLeftToRight; // left boundary edges -> right boundary edges
    std::unordered_map<int, int> boundaryEdgeTopToBottom; // top boundary edges -> bot boundary edges
    std::unordered_map<int, int> boundaryElemLeftToRight; // left boundary elements -> right boundary elements
    std::unordered_map<int, int> boundaryElemTopToBottom; // top boundary elements -> bot boundary elements
    std::unordered_map<int, int> boundaryToGhostElements; // boundary element -> ghost element
    NeighbourService(const NodePool& np, const ElementPool& ep, const EdgePool& edgePool);
    std::vector<int> getEdgeNeighborsOfNode(int nodeIndex) const;
    std::unordered_set<int> getNodeNeighbours(int nodeIndex) const;
    std::unordered_set<int> getEdgeNeighbours(int edgeIndex) const;
    std::unordered_set<int> getElementNeighbours(int elementIndex) const;
    std::unordered_set<int> getEdgesOfElement(int elementIndex) const;
    std::unordered_set<int> getElementsOfEdge(int edgeIndex) const;
    int findEdgeByNodes(int node1Index, int node2Index, const EdgePool& edgePool) const;

    void displayNeighbours() const;
};


// Problem's World
class World {
private:
    NodePool np; // набор узлов
    ElementPool ep; // набор элементов
    EdgePool edgp; // набор рёбер
    NeighbourService ns; // сервис соседей
public:
    double maxX; // максимальное значение абсциссы
    double minX; // минимальное значение абсциссы
    double maxY; // максимальное значение ординаты
    double minY; // минимальное значение ординаты
    std::vector<int> boundaryLeftNodes; // индексы левых граничных узлов
    std::vector<int> boundaryLeftElems; // индексы левых граничных элементов
    std::vector<int> boundaryLeftEdges; // индексы левых граничных рёбер
    std::vector<int> boundaryRightNodes; // индексы правых граничных узлов
    std::vector<int> boundaryRightElems; // индексы правых граничных элементов
    std::vector<int> boundaryRightEdges; // индексы правых граничных рёбер
    std::vector<int> boundaryTopNodes; // индексы верхних граничных узлов
    std::vector<int> boundaryTopElems; // индексы верхних граничных элементов
    std::vector<int> boundaryTopEdges; // индексы верхних граничных рёбер
    std::vector<int> boundaryBottomNodes; // индексы нижних граничных узлов
    std::vector<int> boundaryBottomElems; // индексы нижних граничных элементов
    std::vector<int> boundaryBottomEdges; // индексы нижних граничных рёбер
    void setNodePool(const NodePool& np);
    NodePool getNodePool() const;
    void setElementPool(const ElementPool& ep);
    ElementPool getElementPool() const;
    EdgePool getEdgePool() const;
    void setEdgePool(const EdgePool& edgp);
    World(const std::string &fileName, const bool isRenderedBin);
    void display() const;
    NeighbourService& getNeighbourService();
    void exportToFile(const std::string& filename) const; // экспорт в бинарный файл
    void importFromFile(const std::string& filename); // импорт из файла узлов и элементов
    bool isCounterClockwise(const std::vector<int> &nodeIndices);
};

double areaCalc(const Element& poly, const NodePool& nPool); // подсчёт площади элемента

std::vector<double> getElementCentroid2D(const Element& poly, const NodePool& nPool); // подсчёт середины элемента
std::vector<double> getElementCentroid2D(const std::vector<Node> &nodes);

std::vector<double> getMidPoint2D(const int nodeInd1, const int nodeInd2, const NodePool& nPool); // подсчёт середины отрезка, соединяющего два узла

double getDistance(const int nodeInd1, const int nodeInd2, const NodePool& nPool); // подсчёт расстояния между двумя узлами

void setNeighbourEdge(Element& el, const int edgeInd);

std::vector<double> reflectNodeOverVector(const Node& nodeToReflect, const Node& node1, const Node& node2);

#endif // MAGNETTOPRJCT_NETGEOMETRY_H
