#include "MHDgpu.cuh"
#include <cuda_runtime.h>

// Хостовая обёртка для запуска ядра
void computeHLLDFluxesGPU(const std::vector<Element> &elemHost,
                          const std::vector<Edge>    &edgeHost,
                          const std::vector<std::vector<double>> &elemUs,
                          const std::vector<std::vector<double>> &ghostElemUs,
                          const int innerElemCount,
                          unordered_map<int, int> bTg_map,
                          const double gam_hcr,
                          std::vector<double> &fluxHost,
                          std::vector<double> &unrotated_fluxHost) {

    int N = (int)elemHost.size();
    int ghostN = N - innerElemCount;
    int M = (int)edgeHost.size();

    //  1. Выделяем device-память
    GPU_Element *d_elems;
    GPU_Element *d_ghost_elems;
    GPU_Edge    *d_edges;
    double    *d_fluxes;
    double    *d_unrotated_fluxes;
    cudaMalloc(&d_elems, sizeof(GPU_Element)*innerElemCount);                    //  2
    cudaMalloc(&d_ghost_elems, sizeof(GPU_Element)*ghostN);
    cudaMalloc(&d_edges, sizeof(GPU_Edge)*M);                       //  3
    cudaMalloc(&d_fluxes, sizeof(double)*M*8);                      //  4
    cudaMalloc(&d_unrotated_fluxes, sizeof(double)*M*8);

    //  2. Подготавливаем и копируем данные
    std::vector<GPU_Element> elems(innerElemCount);
    for (int i = 0; i < innerElemCount; ++i) {
        for (int k = 0; k < 8; ++k) elems[i].U[k] = elemUs[i][k];
    }
    std::vector<GPU_Element> ghost_elems(ghostN);
    for(int i = 0; i < ghostN; ++i){
        for(int k = 0; k < 8; ++k) ghost_elems[i].U[k] = ghostElemUs[i][k];
    }

    std::vector<GPU_Edge> edges(M);
    for (int i = 0; i < M; ++i) {
        edges[i].left  = edgeHost[i].neighbourInd1;
        edges[i].right = edgeHost[i].neighbourInd2;
        edges[i].nx    = edgeHost[i].normalVector[0];
        edges[i].ny    = edgeHost[i].normalVector[1];
        edges[i].neighbour1_is_boundary = elemHost[edgeHost[i].neighbourInd1].is_boundary;
        edges[i].neighbour1_is_ghost = elemHost[edgeHost[i].neighbourInd1].is_ghost;
        if(edges[i].neighbour1_is_boundary && edges[i].right == -1){
            int ghostInd = bTg_map[edges[i].left] - innerElemCount;
            edges[i].ghostInd = ghostInd;
        }
        else if(edges[i].neighbour1_is_ghost && edges[i].right == -1){
            int ghostInd = edges[i].left - innerElemCount;
            edges[i].ghostInd = ghostInd;
        }
        else{
            edges[i].ghostInd = -1;
        }
    }
    cudaMemcpy(d_elems, elems.data(), sizeof(GPU_Element)*innerElemCount, cudaMemcpyHostToDevice);   //  7
    cudaMemcpy(d_ghost_elems, ghost_elems.data(), sizeof(GPU_Element)*ghostN, cudaMemcpyHostToDevice);
    cudaMemcpy(d_edges, edges.data(), sizeof(GPU_Edge)*M,    cudaMemcpyHostToDevice);   //  8

    //  3. Запуск kernel'а
    int threadsPerBlock = 256;                                     //  9
    int blocks = (M + threadsPerBlock - 1) / threadsPerBlock;      // 10
    kernelComputeHLLD<<<blocks, threadsPerBlock>>>(d_elems, d_ghost_elems, d_edges, d_fluxes, d_unrotated_fluxes, M, gam_hcr); // 11

    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess){
        std::cerr << "Kernel launch failed!\n";
    }
    cudaDeviceSynchronize();
    //  4. Копируем результат
    cudaMemcpy(fluxHost.data(), d_fluxes, sizeof(double)*8*M, cudaMemcpyDeviceToHost); // 13
    cudaMemcpy(unrotated_fluxHost.data(), d_unrotated_fluxes, sizeof(double)*8*M, cudaMemcpyDeviceToHost);

    //  5. Освобождаем
    cudaFree(d_elems);
    cudaFree(d_ghost_elems),
    cudaFree(d_edges);
    cudaFree(d_fluxes);
    cudaFree(d_unrotated_fluxes);
}

// CUDA-ядро: вычисление HLLD-потоков для каждого ребра
__global__ void kernelComputeHLLD(const GPU_Element *elems,
                                  const GPU_Element *ghost_elems,
                                  const GPU_Edge    *edges,
                                  double          *fluxes,
                                  double          *unrotated_fluxes,
                                  int               numEdges,
                                  const double gam_hcr) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;                    //  1
    if (i >= numEdges) return;                                        //  2

    double *F = fluxes + i*8;
    double *Fu = unrotated_fluxes + i*8;

    const GPU_Edge &e = edges[i];                                     //  3

    if(e.neighbour1_is_boundary && e.right == -1){
        const GPU_Element &EL = elems[e.left];
        const GPU_Element &ER = ghost_elems[e.ghostInd];
        device_HLLD_flux(EL.U, ER.U, e.nx, e.ny, F, Fu, gam_hcr);
    }
    else if(e.neighbour1_is_ghost && e.right == -1){
        const GPU_Element &EL = ghost_elems[e.ghostInd];
        const GPU_Element &ER = ghost_elems[e.ghostInd];
        device_HLLD_flux(EL.U, ER.U, e.nx, e.ny, F, Fu, gam_hcr);
    }
    else{
        const GPU_Element &EL = elems[e.left];
        const GPU_Element &ER = elems[e.right];
        device_HLLD_flux(EL.U, ER.U, e.nx, e.ny, F, Fu, gam_hcr);
    }
}