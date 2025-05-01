#include "GPUFunctions.cuh"

__device__ double signum_gpu(double x) {
    if (x == 0.0) return 0.0;
    return signbit(x) ? -1.0 : 1.0;
}

__device__ double pressure_gpu(const double &gam_hcr, const double &e, const double &rho, const double &u, const double &v, const double &w, const double &Bx, const double &By, const double &Bz){
    return (gam_hcr - 1.0) * (e - 0.5 * rho * (u * u + v * v + w * w) - 0.5*(Bx * Bx + By * By + Bz * Bz));
}

// Суммарное давление
__device__ double ptotal_gpu(const double &p, const double &Bx, const double &By, const double &Bz) {
    return p + 0.5 * (Bx*Bx + By*By + Bz*Bz);
}

__device__ double cfast_gpu(const double *U, const double& gam_hcr) {
    double rho = U[0];
    double u = U[1]/rho;
    double v = U[2]/rho;
    double w = U[3]/rho;
    double e = U[4];
    double Bx = U[5];
    double By = U[6];
    double Bz = U[7];

    //|B|^2
    double BB = Bx * Bx + By * By + Bz * Bz;
    //p
    //double p = std::max(pressure(gam_hcr, e, rho, u, v, w, Bx, By, Bz), 1e-15);
    double p = pressure_gpu(gam_hcr, e, rho, u, v, w, Bx, By, Bz);

    double cfast_gpu = sqrt((gam_hcr * p + BB + sqrt((gam_hcr * p + BB) * (gam_hcr * p + BB) - 4.0 * gam_hcr * p * Bx * Bx)) / (2 * rho));

    return cfast_gpu;
}

// Определяем МГД-поток
__device__ void compute_MHD_flux(const double &rho, const double &u, const double &v, const double &w, const double &e, const double &Bx, const double &By, const double &Bz, const double &pT, const double &gam_hcr, double *Fout) {

    Fout[0] = rho * u;
    Fout[1] = rho * u * u + pT - Bx * Bx /*/(4*PI)*/;
    Fout[2] = rho * v * u - Bx * By /*/ (4 * PI)*/;
    Fout[3] = rho * w * u - Bx * Bz /*/ (4 * PI)*/;
    Fout[4] = (e + pT) * u - Bx * (u * Bx + v * By + w * Bz) /*/ (4 * PI)*/;
    Fout[5] = 0.0;
    Fout[6] = By * u - Bx * v;
    Fout[7] = Bz * u - Bx * w;

}

__device__ void compute_MHD_flux(const double *U, const double &gam_hcr, double *Fout) {
    double rho = U[0];
    double u = U[1] / rho;
    double v = U[2] / rho;
    double w = U[3] / rho;
    double e = U[4];
    double Bx = U[5];
    double By = U[6];
    double Bz = U[7];

    double p = pressure_gpu(gam_hcr, e, rho, u, v, w, Bx, By, Bz);
    double pT = ptotal_gpu(p, Bx, By, Bz);

    compute_MHD_flux(rho, u, v, w, e, Bx, By, Bz, pT, gam_hcr, Fout);
}

__device__ void compute_HLLD_flux(const double *U_L, const double *U_R, const double gam_hcr, double *Fout){
    double rho_L = U_L[0];
    double u_L = U_L[1]/rho_L;
    double v_L = U_L[2]/rho_L;
    double w_L = U_L[3]/rho_L;
    double e_L = U_L[4];
    double Bx_L = U_L[5];
    double By_L = U_L[6];
    double Bz_L = U_L[7];
    double p_L = pressure_gpu(gam_hcr, e_L, rho_L, u_L, v_L, w_L, Bx_L, By_L, Bz_L);
    double pT_L = ptotal_gpu(p_L, Bx_L, By_L, Bz_L);

    double rho_R = U_R[0];
    double u_R = U_R[1]/rho_R;
    double v_R = U_R[2]/rho_R;
    double w_R = U_R[3]/rho_R;
    double e_R = U_R[4];
    double Bx_R = U_R[5];
    double By_R = U_R[6];
    double Bz_R = U_R[7];
    double p_R = pressure_gpu(gam_hcr, e_R, rho_R, u_R, v_R, w_R, Bx_R, By_R, Bz_R);
    double pT_R = ptotal_gpu(p_R, Bx_R, By_R, Bz_R);

    //быстрые магнитозвуковые скорости на левом и правом концах
    double cf_L = cfast_gpu(U_L, gam_hcr);
    double cf_R = cfast_gpu(U_R, gam_hcr);

    double Bx = (sqrt(rho_L)*Bx_L + sqrt(rho_R)*Bx_R) / (sqrt(rho_L) + sqrt(rho_R));

    //скорость левого сигнала, рассчитываемая как минимальное значение скорости левого состояния (uL) и быстрой магнитозвуковой скорости (cfL).
    //double SL = std::min(u_L - cf_L, u_R - cf_R);
    double SL = fmin(u_L, u_R) - fmax(cf_L, cf_R);
    // Скорость правого сигнала, рассчитываемая как максимальное значение скорости правого состояния (uR) и быстрой магнитозвуковой скорости (cfR).ы
    //double SR = std::max(u_L + cf_L, u_R + cf_R);
    double SR = fmin(u_L, u_R) + fmax(cf_L, cf_R);

    double SR_m_uR = SR - u_R;
    double SL_m_uL = SL - u_L;
    double den_SM = SR_m_uR * rho_R - SL_m_uL * rho_L;

    // Скорость середины разрыва
    double SM = (SR_m_uR * rho_R * u_R - SL_m_uL * rho_L * u_L - pT_R + pT_L) / den_SM;

    double SM_m_uR = SM - u_R;
    double SM_m_uL = SM - u_L;

    double rho_L_star = rho_L;
    if(fabs(SL-SM) > 1e-15) {
        rho_L_star = rho_L * SL_m_uL / (SL - SM);
    }

    double rho_R_star = rho_R;
    if(fabs(SR-SM) > 1e-15){
        rho_R_star = rho_R * SR_m_uR / (SR - SM);
    }

    double SL_star = SM - std::fabs(Bx)/sqrt(rho_L_star);
    double SR_star = SM + std::fabs(Bx)/sqrt(rho_R_star);

    double pT_star = (SR_m_uR * rho_R * pT_L - SL_m_uL * rho_L * pT_R + rho_L * rho_R * SR_m_uR * SL_m_uL * (u_R - u_L)) / den_SM;

    // звёздочка слева
    double u_L_star = u_L;
    double pT_L_star = pT_L;
    double v_L_star = v_L;
    double w_L_star = w_L;
    double By_L_star = 0.0;
    double Bz_L_star = 0.0;
    double denom_L = rho_L * SL_m_uL * (SL - SM) - Bx * Bx;
    if(fabs(denom_L) > 1e-15){
        pT_L_star = pT_star;
        u_L_star = SM;
        v_L_star = v_L - Bx * By_L * SM_m_uL/denom_L;
        w_L_star = w_L - Bx * Bz_L * SM_m_uL/denom_L;
        double factor_L = rho_L * SL_m_uL * SL_m_uL - Bx * Bx;
        By_L_star = By_L * factor_L/denom_L;
        Bz_L_star = Bz_L * factor_L/denom_L;
    }
    double e_L_star = (SL_m_uL *e_L - pT_L*u_L + pT_star*SM + Bx*(u_L * Bx_L + v_L * By_L + w_L * Bz_L - u_L_star*Bx - v_L_star*By_L_star - w_L_star*Bz_L_star))/(SL-SM);

    // звёздочка справа
    double u_R_star = u_R;
    double pT_R_star = pT_R;
    double v_R_star = v_R;
    double w_R_star = w_R;
    double By_R_star = 0.0;
    double Bz_R_star = 0.0;
    double denom_R = rho_R * SR_m_uR * (SR - SM) - Bx * Bx;
    if(fabs(denom_R) > 1e-15){
        pT_R_star = pT_star;
        u_R_star = SM;
        v_R_star = v_R - Bx * By_R * SM_m_uR/denom_R;
        w_R_star = w_R - Bx * Bz_R * SM_m_uR/denom_R;
        double factor_R = rho_R * SR_m_uR * SR_m_uR - Bx * Bx;
        By_R_star = By_R * factor_R/denom_R;
        Bz_R_star = Bz_R * factor_R/denom_R;
    }
    double e_R_star = (SR_m_uR *e_R - pT_R*u_R + pT_star*SM + Bx*(u_R * Bx_R + v_R * By_R + w_R * Bz_R - u_R_star*Bx - v_R_star*By_R_star - w_R_star*Bz_R_star))/(SR-SM);

    if (SL > 0.0) {
//!!!   //FL
        compute_MHD_flux(U_L, gam_hcr, Fout);
    }
    else if (/*SL <= 0 &&*/ SL_star >= 0.0) {
//!!!   //F*L
        compute_MHD_flux(rho_L_star, u_L_star, v_L_star, w_L_star, e_L_star, Bx, By_L_star, Bz_L_star, pT_L_star, gam_hcr, Fout);
    }
    else if (/*SL_star <= 0 &&*/ SM >= 0.0) {
//!!!   //F**L
        double u_L_2star = SM;
        double pT_L_2star = pT_L_star;
        double rho_star_sum_of_sqrts = sqrt(rho_L_star)+ sqrt(rho_R_star);
        double sign_Bx = signum_gpu(Bx_L);
        double v_L_2star = (sqrt(rho_L_star)*v_L_star + sqrt(rho_R_star)*v_R_star + (By_R_star-By_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double w_L_2star = (sqrt(rho_L_star)*w_L_star + sqrt(rho_R_star)*w_R_star + (Bz_R_star-Bz_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double By_L_2star = (sqrt(rho_L_star)*By_R_star + sqrt(rho_R_star)*By_L_star + sqrt(rho_L_star*rho_R_star)*(v_R_star-v_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double Bz_L_2star = (sqrt(rho_L_star)*Bz_R_star + sqrt(rho_R_star)*Bz_L_star + sqrt(rho_L_star*rho_R_star)*(w_R_star-w_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double e_L_2star = e_L_star - sqrt(rho_L_star)*(u_L_star * Bx + v_L_star * By_L_star + w_L_star * Bz_L_star - u_L_2star * Bx - v_L_2star * By_L_2star - w_L_2star * Bz_L_2star)*sign_Bx;
        compute_MHD_flux(rho_L_star, u_L_2star, v_L_2star, w_L_2star, e_L_2star, Bx, By_L_2star, Bz_L_2star, pT_L_2star, gam_hcr, Fout);
    }
    else if (/*SM <= 0 &&*/ SR_star >= 0.0){
//!!!   //F**R
        double u_R_2star = SM;
        double pT_R_2star = pT_R_star;
        double rho_star_sum_of_sqrts = sqrt(rho_L_star)+ sqrt(rho_R_star);
        double sign_Bx = signum_gpu(Bx_L);
        double v_R_2star = (sqrt(rho_L_star)*v_L_star + sqrt(rho_R_star)*v_R_star + (By_R_star-By_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double w_R_2star = (sqrt(rho_L_star)*w_L_star + sqrt(rho_R_star)*w_R_star + (Bz_R_star-Bz_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double By_R_2star = (sqrt(rho_L_star)*By_R_star + sqrt(rho_R_star)*By_L_star + sqrt(rho_L_star*rho_R_star)*(v_R_star-v_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double Bz_R_2star = (sqrt(rho_L_star)*Bz_R_star + sqrt(rho_R_star)*Bz_L_star + sqrt(rho_L_star*rho_R_star)*(w_R_star-w_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double e_R_2star = e_R_star + sqrt(rho_R_star)*(u_R_star * Bx + v_R_star * By_R_star + w_R_star * Bz_R_star - u_R_2star * Bx - v_R_2star * By_R_2star - w_R_2star * Bz_R_2star)*sign_Bx;
        compute_MHD_flux(rho_R_star, u_R_2star, v_R_2star, w_R_2star, e_R_2star, Bx, By_R_2star, Bz_R_2star, pT_R_2star, gam_hcr, Fout);
    }
    else if (/*SR_star <= 0 &&*/ SR >= 0.0){
//!!!   //F*R
        compute_MHD_flux(rho_R_star, u_R_star, v_R_star, w_R_star, e_R_star, Bx, By_R_star, Bz_R_star, pT_R_star, gam_hcr, Fout);
    }
        //if (SR < 0)
    else  {
//!!!   //FR
        compute_MHD_flux(U_R, gam_hcr, Fout);
    }

}

__device__ void device_HLLD_flux(const double *UL, const double *UR,
                                 const double nx, const double ny,
                                 double *Fout, double *unrotatedFout, const double gam_hcr) {

    //  1. поворот в нормаль
    double uL[8];
    uL[0] = UL[0];
    uL[1] = UL[1]*nx + UL[2]*ny;
    uL[2] = -UL[1]*ny + UL[2]*nx;
    uL[3] = UL[3];
    uL[4] = UL[4];
    uL[5] = UL[5]*nx + UL[6]*ny;
    uL[6] = -UL[5]*ny + UL[6]*nx;
    uL[7] = UL[7];

    // аналогично для UR
    double uR[8];
    uR[0] = UR[0];
    uR[1] = UR[1]*nx + UR[2]*ny;
    uR[2] = -UR[1]*ny + UR[2]*nx;
    uR[3] = UR[3];
    uR[4] = UR[4];
    uR[5] = UR[5]*nx + UR[6]*ny;
    uR[6] = -UR[5]*ny + UR[6]*nx;
    uR[7] = UR[7];

    compute_HLLD_flux(uL, uR, gam_hcr, unrotatedFout);

    // поворот для unrotated
    Fout[0] = unrotatedFout[0];
    Fout[1] = unrotatedFout[1]*nx - unrotatedFout[2]*ny;
    Fout[2] = unrotatedFout[1]*ny + unrotatedFout[2]*nx;
    Fout[3] = unrotatedFout[3];
    Fout[4] = unrotatedFout[4];
    Fout[5] = unrotatedFout[5]*nx - unrotatedFout[6]*ny;
    Fout[6] = unrotatedFout[5]*ny + unrotatedFout[6]*nx;
    Fout[7] = unrotatedFout[7];
}