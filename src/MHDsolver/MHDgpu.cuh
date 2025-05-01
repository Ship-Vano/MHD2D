#ifndef GPU_HLLDFLUX_CUH
#define GPU_HLLDFLUX_CUH

#include "MHDSolver2D.h"      // для определения HLLD_flux и типов
#include "GPUFunctions.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include<vector>
#include <algorithm>
#include <math.h>

void computeHLLDFluxesGPU(const std::vector<Element> &elemHost,
                          const std::vector<Edge>    &edgeHost,
                          const std::vector<std::vector<double>> &elemUs,
                          const std::vector<std::vector<double>> &ghostElemUs,
                          const int innerElemCount,
                          unordered_map<int, int> &bTg_map,
                          const double gam_hcr,
                          std::vector<std::vector<double>> &fluxHost,
                          std::vector<std::vector<double>> &unrotated_fluxHost);
// Структуры данных (должны совпадать с хостовыми)
struct GPU_Element {
    double U[8];      // вектор состояния: rho, rho*u, rho*v, Bx, By, Bz, E, psi
};

struct GPU_Edge {
    int left;         // индекс левого элемента
    int right;        // индекс правого элемента (или -1 для границы)
    double nx, ny;    // компоненты нормали
    bool neighbour1_is_boundary;
    bool neighbour1_is_ghost;
    int ghostInd;
};

// Результирующий поток
struct GPU_Flux {
    double F[8];
};


// CUDA-ядро: вычисление HLLD-потоков для каждого ребра
__global__ void kernelComputeHLLD(const GPU_Element *elems,
                                  const GPU_Element *ghost_elems,
                                  const GPU_Edge    *edges,
                                  GPU_Flux          *fluxes,
                                  GPU_Flux          *unrotated_fluxes,
                                  int               numEdges,
                                  const double gam_hcr);

#endif // GPU_HLLDFLUX_CUH