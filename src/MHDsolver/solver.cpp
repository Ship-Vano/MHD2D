//
// Created by Иван on 1/18/2025.
//
#include "NetGeometry.h"
#include "MHDSolver2D.h"
#include "json/json.h"

int main(){

    //TODO: move it to the global config file
    std::string configPath = "InputData/solverConfig.json";
    std::ifstream in(configPath, std::ios::in);

    Json::Reader json_reader;
    Json::Value json_root;
    bool read_succeeded = json_reader.parse(in, json_root);
    assert(read_succeeded);

    std::string importFileName = json_root.get("importFileName", "").asString();
    std::cout << "Generating the world from " << importFileName << " file by reading it..." << std::endl;
    World world(importFileName, true);
    std::cout << "MINlen = "<< world.getEdgePool().minEdgeLen << std::endl;
    int taskType = json_root.get("taskType", 1).asInt();
    // world.display();
    //std::cin.get();
    omp_set_num_threads(omp_get_max_threads());
    MHDSolver2D solver(world);
    solver.task_type = taskType;
    solver.runSolver();
    std::cout << "solver complete" << std::endl;

    std::string exportFileName = json_root.get("exportFileName", "OutputData/unnamed_res.vtu").asString();
    std::cout << "Writing vtu output to " << exportFileName << std::endl;
    writeVTU(exportFileName, world, solver.elemUs);

    return 0;
}