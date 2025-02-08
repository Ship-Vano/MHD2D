//
// Created by Иван on 1/18/2025.
//

#include "NetGeometry.h"
#include "geometryTests.h"
#include "json/json.h"

#include <string>

int main() {
    omp_set_num_threads(1);
    //TODO: move it to the global config file
    std::string configPath = "InputData/netConfig.json";
    std::ifstream in(configPath, std::ios::in);

    Json::Reader json_reader;
    Json::Value json_root;
    bool read_succeeded = json_reader.parse(in, json_root);
    assert(read_succeeded);

    bool executeTests = json_root.get("executeTests", false).asBool();
    if(executeTests){
        std::cout << "Executing program tests..." << std::endl;
        testNode();
        testElement();
        testEdge();
        testAreaCalc();
        testCentroidCalc();
        testEdgePool();
        testNeighbourService();
        testReflectPointOverVector();
    }

    std::string dataFileName = json_root.get("fileName", "").asString();
    std::cout << "fileName is " << dataFileName << std::endl;
    omp_set_num_threads(omp_get_max_threads());
    World world(dataFileName, false);
    std::cout << "World has been generated successfully!" << std::endl;
    std::cout << "[INFO]: NNodes = " << world.getNodePool().nodeCount << ", NElems = " << world.getElementPool().elCount << std::endl;

    for (const auto &elem: world.getElementPool().elements) {
        if (elem.edgeIndexes.empty()) {
            std::cout << "Empty element! (no edges)" << std::endl;
        }
        if(elem.area <= 0.0){
            std::cout << "Zero-area element! (elem.area = 0)" << std::endl;
        }
    }
    bool displayWorld = json_root.get("displayWorld", false).asBool();
    if(displayWorld){
        world.display();
    }
    bool exportWorld = json_root.get("exportWorld", false).asBool();
    if(exportWorld){
        std::cout << "Exporting world..." << std::endl;
        std::string exportFileName = json_root.get("exportFileName", "OutputData/unnamedWorld").asString();
        world.exportToFile(exportFileName);
    }
    return 0;
}
