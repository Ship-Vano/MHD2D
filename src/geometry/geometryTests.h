//
// Created by Иван on 1/18/2025.
//

#ifndef MHD2D_GEOMETRYTESTS_H
#define MHD2D_GEOMETRYTESTS_H

#include "NetGeometry.h"
#include <climits>

/* ТЕСТЫ ГЕОМЕТРИИ*/
void testNode() {
    Node n(0, 1.0, 2.0, 3.0);
    assert(n.ind == 0);
    assert(n.pos.x == 1.0);
    assert(n.pos.y == 2.0);
    assert(n.pos.z == 3.0);

    Node n1(-100, -1.0e-9, 2.0e9, -3.0e9);
    assert(n1.ind == -100);
    assert(n1.pos.x == -1.0e-9);
    assert(n1.pos.y == 2.0e9);
    assert(n1.pos.z == -3.0e9);

    Node n2(INT_MAX, 1.0, 2.0, 3.0);
    assert(n2.ind == INT_MAX);

    Node n3(INT_MIN, -1.0, -2.0, -3.0);
    assert(n3.ind == INT_MIN);

    std::cout << "Node test passed!" << std::endl;
}
void testElement() {
    std::vector<int> nodeIndexes = {0, 1, 2};
    Element e(0, nodeIndexes, 3);
    assert(e.ind == 0);
    assert(e.nodeIndexes == nodeIndexes);
    assert(e.dim == 3);

    // Test with repeated indices
    std::vector<int> repeatedIndexes = {0, 0, 0};
    Element e1(0, repeatedIndexes, 3);
    assert(e1.ind == 0);
    assert(e1.nodeIndexes == repeatedIndexes);
    assert(e1.dim == 3);

    // Degenerate case with collinear points
    std::vector<int> degenerateIndexes = {0, 1, 2};
    Element e2(1, degenerateIndexes, 3);
    assert(e2.ind == 1);
    assert(e2.nodeIndexes == degenerateIndexes);
    assert(e2.dim == 3);

    std::cout << "Element test passed!" << std::endl;
}
void testEdge() {
    Vec2 normalVec = {0.0, 1.0};
    Vec2 midP{0.5, 1.5};
    Edge e(0, 0, 1, -1, -1, 1, normalVec, midP);
    assert(e.ind == 0);
    assert(e.nodeInd1 == 0);
    assert(e.nodeInd2 == 1);
    assert(e.neighbourInd1 == -1);
    assert(e.neighbourInd2 == -1);
    assert(e.length == 1.0);
    assert(e.normalVector == normalVec);
    assert(e.midPoint == midP);

    // Zero-length edge
    Vec2 zeroNormal = {0.0, 0.0};
    Vec2 zeroMid = {0.0, 0.0};
    Edge e1(0, 0, 0, -1, -1, 0, zeroNormal, zeroMid);
    assert(e1.length == 0.0);

    // Test normal vector precision
    Vec2 normalVec1 = {1.0 / std::sqrt(2), 1.0 / std::sqrt(2)};
    Vec2 midP1 = {0.5, 0.5};
    Edge e2(1, 0, 1, -1, -1, 1, normalVec1, midP1);
    assert(std::abs(e2.normalVector.x - 0.7071) < 1e-4);
    assert(std::abs(e2.normalVector.y - 0.7071) < 1e-4);

    std::cout << "Edge test passed!" << std::endl;
}
void testAreaCalc() {
    std::vector<Node> nodes1 = {
            Node(0, 0.0, 0.0, 0.0),
            Node(1, 1.0, 0.0, 0.0),
            Node(2, 0.0, 1.0, 0.0)
    };
    std::vector<int> nodeIndexes1 = {0, 1, 2};
    Element e1(0, nodeIndexes1, 3);
    NodePool np1(3, nodes1);
    double area1 = areaCalc(e1, np1);
    assert(area1 == 0.5); // Area of a triangle with base and height of 1

    // Degenerate triangle (all points lie on a line)
    std::vector<Node> degenerateNodes = {
            Node(0, 0.0, 0.0, 0.0),
            Node(1, 1.0, 1.0, 0.0),
            Node(2, 2.0, 2.0, 0.0)
    };
    std::vector<int> nodeIndexes = {0, 1, 2};
    Element e(0, nodeIndexes, 3);
    NodePool np(3, degenerateNodes);
    double area = areaCalc(e, np);
    assert(area == 0.0);

    // Large triangle
    std::vector<Node> largeNodes = {
            Node(0, 1e9, 0.0, 0.0),
            Node(1, 1e9, 1e9, 0.0),
            Node(2, 0.0, 0.0, 0.0)
    };
    nodeIndexes = {0, 1, 2};
    Element e2(0, nodeIndexes, 3);
    NodePool np2(3, largeNodes);
    double area2 = areaCalc(e2, np2);
    assert(area2 == 5.0e17);

    std::cout << "Area calculation test passed!" << std::endl;
}
void testCentroidCalc() {
    std::vector<Node> nodes1 = {
            Node(0, 0.0, 0.0, 0.0),
            Node(1, 1.0, 0.0, 0.0),
            Node(2, 0.0, 1.0, 0.0)
    };
    std::vector<int> nodeIndexes1 = {0, 1, 2};
    Element e1(0, nodeIndexes1, 3);
    NodePool np1(3, nodes1);
    Vec2 centroid1 = getElementCentroid2D(e1, np1);
    assert(centroid1.x == 1.0 / 3.0);
    assert(centroid1.y == 1.0 / 3.0);

    // Degenerate triangle (all points on a line)
    std::vector<Node> degenerateNodes = {
            Node(0, 0.0, 0.0, 0.0),
            Node(1, 1.0, 1.0, 0.0),
            Node(2, 2.0, 2.0, 0.0)
    };
    std::vector<int> nodeIndexes = {0, 1, 2};
    Element e(0, nodeIndexes, 3);
    NodePool np(3, degenerateNodes);
    Vec2 centroid = getElementCentroid2D(e, np);
    assert(centroid.x == 1.0);
    assert(centroid.y == 1.0);

    // Large triangle
    std::vector<Node> largeNodes = {
            Node(0, 1e9, 0.0, 0.0),
            Node(1, 1e9, 1e9, 0.0),
            Node(2, 0.0, 0.0, 0.0)
    };
    nodeIndexes = {0, 1, 2};
    Element e2(0, nodeIndexes, 3);
    NodePool np2(3, largeNodes);
    centroid = getElementCentroid2D(e2, np2);
    assert(std::abs(centroid.x - 6.6667e8) < 1e4);
    assert(std::abs(centroid.y - 3.3333e8) < 1e4);


    std::cout << "Centroid calculation test passed!" << std::endl;
}
void testEdgePool() {
    std::vector<Node> nodes1 = {
            Node(0, 0.0, 0.0, 0.0),
            Node(1, 1.0, 0.0, 0.0),
            Node(2, 0.0, 1.0, 0.0),
            Node(3, 1.0, 1.0, 0.0)
    };
    std::vector<Element> elements1 = {
            Element(0, {0, 1, 2}, 3),
            Element(1, {0, 2, 3}, 3)
    };
    NodePool np1(4, nodes1);
    ElementPool ep1(3, 2, elements1);
    EdgePool edgePool1(np1, ep1);
    assert(edgePool1.edges.size() == 5);  // Check if the expected number of edges is created
    std::cout << "EdgePool test passed!" << std::endl;
}
void testNeighbourService() {
    std::vector<Node> nodes1 = {
            Node(0, 0.0, 0.0, 0.0),
            Node(1, 1.0, 0.0, 0.0),
            Node(2, 0.0, 1.0, 0.0),
            Node(3, 1.0, 1.0, 0.0)
    };
    std::vector<Element> elements1 = {
            Element(0, {0, 1, 2}, 3),
            Element(1, {1, 2, 3}, 3)
    };
    NodePool np1(4, nodes1);
    ElementPool ep1(3, 2, elements1);
    EdgePool edgePool1(np1, ep1);
    NeighbourService ns1(np1, ep1, edgePool1);

    std::unordered_set<int> nodeNeighbours1 = ns1.getNodeNeighbours(1);
    assert(nodeNeighbours1.size() == 2);  // Node 1 should be shared by two elements

    std::vector<Node> nodes = {
            Node(0, 0.0, 0.0, 0.0),
            Node(1, 1.0, 0.0, 0.0),
            Node(2, 0.0, 1.0, 0.0)
    };
    std::vector<Element> elements = {Element(0, {0, 1, 2}, 3)};
    NodePool np(3, nodes);
    ElementPool ep(3, 1, elements);
    EdgePool edgePool(np, ep);
    NeighbourService ns(np, ep, edgePool);


    std::cout << "NeighbourService test passed!" << std::endl;
}

void testReflectPointOverVector(){
    Node A(0,0.0,0.0,0.0);
    Node B(1, 2.0,2.0,0.0);
    Node C(2, 3.0,1.0,0.0);
    std::vector<double> resPoint = reflectNodeOverVector(C,A,B);
    assert(resPoint[0] == 1.0);
    assert(resPoint[1] == 3.0);

    std::cout << "ReflectPointOverVector test passed!" << std::endl;
}
#endif //MHD2D_GEOMETRYTESTS_H
