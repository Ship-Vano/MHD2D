//
// Created by Иван on 1/18/2025.
//
#include "NetGeometry.h"
#include "MHDSolver2D.h"
#include "json/json.h"

int main(){
    std::cout << "running solver..." << std::endl;
    omp_set_num_threads(1);
    //TODO: move it to the global config file
    std::string configPath = "InputData/solverConfig.json";
    std::ifstream in(configPath, std::ios::in);

    Json::Reader json_reader;
    Json::Value json_root;
    bool read_succeeded = json_reader.parse(in, json_root);
    assert(read_succeeded);

    // std::string importFileName = json_root.get("importFileName", "").asString();
    // std::cout << "Generating the world from " << importFileName << " file by reading it..." << std::endl;
    // bool generatedMesh = json_root.get("generatedMesh", false).asBool();

    // World world = World(importFileName, true);
    // if(!generatedMesh) {
    //     world = World("InputData/mesh.txt", false);
    // }

    World world = World("InputData/mesh.txt", false);
    
    std::cout << "MINlen = "<< world.getEdgePool().minEdgeLen << std::endl;
    int taskType = json_root.get("taskType", 1).asInt();
    double finalTime = json_root.get("finalTime", 0.1).asDouble();
    bool debugDivergence = json_root.get("debugDivergence", false).asBool();
    bool ghostOutput = json_root.get("ghostOutput", false).asBool();
    int iterationsPerFrame = json_root.get("iterationsPerFrame", 10).asInt();
    // world.display();
    //std::cin.get();
    omp_set_num_threads(omp_get_max_threads());
    MHDSolver2D solver(world);
    solver.task_type = taskType;
    solver.finalTime = finalTime;
    solver.debugDivergence = debugDivergence;
    solver.iterationsPerFrame = iterationsPerFrame;
    solver.ghostOutput = ghostOutput;

    bool cylindrical= json_root.get("cylindrical", false).asBool();
    bool gpu = json_root.get("gpu", false).asBool();

    if(cylindrical){
//        float time;
//        cudaEvent_t start, stop;
//        cudaEventCreate(&start);
//        cudaEventCreate(&stop);
//        cudaEventRecord(start, 0);
       solver.runCylindricSolver();
//        cudaEventRecord(stop, 0);
//        cudaEventSynchronize(stop);
//        cudaEventElapsedTime(&time, start, stop);
//        std::cout << "Time to generate: " << time <<  " ms \n";
    }
    else if(gpu){
//        float time;
//        cudaEvent_t start, stop;
//        cudaEventCreate(&start);
//        cudaEventCreate(&stop);
//        cudaEventRecord(start, 0);
//        solver.runGPUSolver();
//        cudaEventRecord(stop, 0);
//        cudaEventSynchronize(stop);
//        cudaEventElapsedTime(&time, start, stop);
//        std::cout << "Time to generate: " << time <<  " ms \n";
    }
    else{
        solver.runSolver();
    }
    std::cout << "solver complete" << std::endl;

    std::string exportFileName = json_root.get("exportFileName", "OutputData/unnamed_res.vtu").asString();
    std::cout << "Writing vtu output to " << exportFileName << std::endl;
    solver.writeVTU(exportFileName, ghostOutput);

    return 0;
}