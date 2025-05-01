#ifndef GPU_FUNCTIONS_H
#define GPU_FUNCTIONS_H

#include <cuda_runtime.h>

__device__ double signum_gpu(double x);
__device__ double pressure_gpu(const double &gam_hcr, const double &e, const double &rho, const double &u, const double &v, const double &w, const double &Bx, const double &By, const double &Bz);

__device__ double ptotal_gpu(const double &p, const double &Bx, const double &By, const double &Bz);

__device__ double cfast_gpu(const double *U, const double& gam_hcr);

// Определяем МГД-поток
__device__ void compute_MHD_flux(const double &rho, const double &u, const double &v, const double &w, const double &e, const double &Bx, const double &By, const double &Bz, const double &pT, const double &gam_hcr, double *Fout);

__device__ void compute_MHD_flux(const double *U, const double &gam_hcr, double *Fout);

__device__ void compute_HLLD_flux(const double *U_L, const double *U_R, const double gam_hcr, double *Fout);

__device__ void device_HLLD_flux(const double *UL, const double *UR,
                                 const double nx, const double ny,
                                 double *Fout, double *unrotatedFout, const double gam_hcr);
#endif