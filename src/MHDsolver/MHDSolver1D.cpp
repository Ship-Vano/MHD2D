//
// Created by Иван on 4/29/2024.
//
#include "MHDSolver1D.h"

/*          0      1      2      3    4   5   6   7
 * state:  rho,  rho*u, rho*v, rho*w, e, Bx, Bz, By
 * */

// Давление
double pressure(const double &gam_hcr, const double &e, const double &rho, const double &u, const double &v, const double &w, const double &Bx, const double &By, const double &Bz){
    return (gam_hcr - 1) * (e - 0.5 * rho * (u * u + v * v + w * w) - 0.5*(Bx * Bx + By * By + Bz * Bz));
}

double pressure(const std::vector<double> &U, const double &gam_hcr){
    double rho = U[0];
    double u = U[1]/rho;
    double v = U[2]/rho;
    double w = U[3]/rho;
    double e = U[4];
    double Bx = U[5];
    double By = U[6];
    double Bz = U[7];
    return (gam_hcr - 1.0) * (e - 0.5 * rho * (u * u + v * v + w * w) - 0.5*(Bx * Bx + By * By + Bz * Bz));
}

// Энергия
double energy(const double &gam_hcr, const double &p, const double &rho, const double &u, const double &v, const double &w, const double &Bx, const double &By, const double &Bz){
    return p/(gam_hcr-1.0) + 0.5*rho*(u * u + v * v + w * w) + 0.5 * (Bx * Bx + By * By + Bz * Bz);
}

// Суммарное давление
double ptotal(const double &p, const double &Bx, const double &By, const double &Bz) {
    return p + 0.5 * (Bx*Bx + By*By + Bz*Bz);
}

// Вектор состояния из параметров
std::vector<double>
state_from_primitive_vars(const double &rho, const double &u, const double &v, const double &w, const double &p,
                          const double &Bx, const double &By, const double &Bz, const double &gam_hcr) {
    std::vector<double> U(8,0.0);

    double mx = rho * u;
    double my = rho * v;
    double mz = rho * w;
    double e = energy(gam_hcr, p, rho, u, v, w, Bx, By, Bz);

    U[0] = rho;
    U[1] = mx;
    U[2] = my;
    U[3] = mz;
    U[4] = e;
    U[5] = Bx;
    U[6] = By;
    U[7] = Bz;

    return U;
}

std::vector<double> state_from_primitive_vars(const std::vector<double>& primitiveVars){
    /*rho  u   v   w   p   Bx   By   Bz  gam_hcr*/
    return state_from_primitive_vars(primitiveVars[0], primitiveVars[1], primitiveVars[2],
                                     primitiveVars[3], primitiveVars[4], primitiveVars[5],
                                     primitiveVars[6], primitiveVars[7], primitiveVars[8]);
}

// Определяем МГД-поток (из состояния)
std::vector<double> MHD_flux(const std::vector<double>& U, const double &gam_hcr) {
    double rho = U[0];
    double u = U[1]/rho;
    double v = U[2]/rho;
    double w = U[3]/rho;
    double e = U[4];
    double Bx = U[5];
    double By = U[6];
    double Bz = U[7];

    double p = pressure(gam_hcr, e, rho, u, v, w, Bx, By, Bz);
    double pT = ptotal(p, Bx, By, Bz);

    std::vector<double> F(8, 0.0);
    F[0] = rho * u;
    F[1] = rho * u * u + pT - Bx * Bx /*/(4*PI)*/;
    F[2] = rho * v * u - Bx * By /*/ (4 * PI)*/;
    F[3] = rho * w * u - Bx * Bz /*/ (4 * PI)*/;
    F[4] = (e + pT) * u - Bx * (u * Bx + v * By + w * Bz) /*/ (4 * PI)*/;
    F[5] = 0.0;
    F[6] = By * u - Bx * v;
    F[7] = Bz * u - Bx * w;

    return F;
}

// Определяем МГД-поток
std::vector<double> MHD_flux(const double &rho, const double &u, const double &v, const double &w, const double &e, const double &Bx, const double &By, const double &Bz, const double &pT, const double &gam_hcr) {

    std::vector<double> F(8, 0);
    F[0] = rho * u;
    F[1] = rho * u * u + pT - Bx * Bx /*/(4*PI)*/;
    F[2] = rho * v * u - Bx * By /*/ (4 * PI)*/;
    F[3] = rho * w * u - Bx * Bz /*/ (4 * PI)*/;
    F[4] = (e + pT) * u - Bx * (u * Bx + v * By + w * Bz) /*/ (4 * PI)*/;
    F[5] = 0.0;
    F[6] = By * u - Bx * v;
    F[7] = Bz * u - Bx * w;

    return F;
}

// Скорость быстрой волны в точке
double cfast(const std::vector<double>& U, const double& gam_hcr) {
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
    double p = pressure(gam_hcr, e, rho, u, v, w, Bx, By, Bz);
    /*if (rho <= 0.0 || p <= 0.0) {
        std::cerr << "Warning: Unphysical state detected in cfast. "
                  << "rho = " << rho << ", p = " << p << std::endl;
        std::cin.get();
    }*/

    double cfast = std::sqrt((gam_hcr * p + BB + std::sqrt((gam_hcr * p + BB) * (gam_hcr * p + BB) - 4.0 * gam_hcr * p * Bx * Bx)) / (2 * rho));

    if (cfast > 1e10 || std::isnan(cfast)) {
        std::cerr << "Unrealistic cfast detected: " << cfast
                  << ", rho = " << rho << ", p = " << p
                  << ", BB = " << BB << std::endl;
        std::cin.get();
    }

    return cfast;
}

// Определяем HLL поток F
std::vector<double> HLL_flux(const std::vector<double>& U_L, const std::vector<double>& U_R, const double &gam_hcr) {
    /*          0      1      2      3    4   5   6   7
     * state:  rho,  rho*u, rho*v, rho*w, e, Bx, Bz, By
     */
    double rho_L = U_L[0];
    double u_L = U_L[1]/rho_L;
    double v_L = U_L[2]/rho_L;
    double w_L = U_L[3]/rho_L;
    double e_L = U_L[4];
    double Bx_L = U_L[5];
    double By_L = U_L[6];
    double Bz_L = U_L[7];

    double rho_R = U_R[0];
    double u_R = U_R[1]/rho_R;
    double v_R = U_R[2]/rho_R;
    double w_R = U_R[3]/rho_R;
    double e_R = U_R[4];
    double Bx_R = U_R[5];
    double By_R = U_R[6];
    double Bz_R = U_R[7];

    //быстрые магнитозвуковые скорости на левом и правом концах
    double cf_L = cfast(U_L, gam_hcr);
    double cf_R = cfast(U_R, gam_hcr);

    //скорость левого сигнала, рассчитываемая как минимальное значение скорости левого состояния (uL) и быстрой магнитозвуковой скорости (cfL).
    double SL = std::min(u_L - cf_L, u_R - cf_R);//std::min(uL, uR) - std::max(cfL, cfR);//
    // Скорость правого сигнала, рассчитываемая как максимальное значение скорости правого состояния (uR) и быстрой магнитозвуковой скорости (cfR).ы
    double SR = std::max(u_L + cf_L, u_R + cf_R);

    if (SL > 0.0) {
        return MHD_flux(U_L, gam_hcr);
    }
    else if (SL <= 0.0 && SR >= 0) {
        return 1 / (SR - SL) * (SR * MHD_flux(U_L, gam_hcr) - SL * MHD_flux(U_R, gam_hcr) + SL * SR * (U_R - U_L));
    }
        //if (SR < 0)
    else  {
        return MHD_flux(U_R, gam_hcr);
    }
}

// Определяем HLLC поток F
std::vector<double> HLLC_flux(const std::vector<double>& U_L, const std::vector<double>& U_R, const double &gam_hcr) {
    /*          0      1      2      3    4   5   6   7
     * state:  rho,  rho*u, rho*v, rho*w, e, Bx, Bz, By
     */
    double rho_L = U_L[0];
    double u_L = U_L[1]/rho_L;
    double v_L = U_L[2]/rho_L;
    double w_L = U_L[3]/rho_L;
    double e_L = U_L[4];
    double Bx_L = U_L[5];
    double By_L = U_L[6];
    double Bz_L = U_L[7];
    double p_L = pressure(gam_hcr, e_L, rho_L, u_L, v_L, w_L, Bx_L, By_L, Bz_L);

    double rho_R = U_R[0];
    double u_R = U_R[1]/rho_R;
    double v_R = U_R[2]/rho_R;
    double w_R = U_R[3]/rho_R;
    double e_R = U_R[4];
    double Bx_R = U_R[5];
    double By_R = U_R[6];
    double Bz_R = U_R[7];
    double p_R = pressure(gam_hcr, e_R, rho_R, u_R, v_R, w_R, Bx_R, By_R, Bz_R);

    //быстрые магнитозвуковые скорости на левом и правом концах
    double cf_L = cfast(U_L, gam_hcr);
    double cf_R = cfast(U_R, gam_hcr);

    //скорость левого сигнала, рассчитываемая как минимальное значение скорости левого состояния (uL) и быстрой магнитозвуковой скорости (cfL).
    double SL = std::min(u_L - cf_L, u_R - cf_R);//std::min(uL, uR) - std::max(cfL, cfR);//
    // Скорость правого сигнала, рассчитываемая как максимальное значение скорости правого состояния (uR) и быстрой магнитозвуковой скорости (cfR).ы
    double SR = std::max(u_L + cf_L, u_R + cf_R);

    double SR_m_uR = SR - u_R;
    double SL_m_uL = SL - u_L;

    /*          0      1      2      3    4   5   6   7
     * state:  rho,  rho*u, rho*v, rho*w, e, Bx, By, Bz
     * */
    //std::vector<double> U_star = 1/(SR-SL) * (SR * U_R - SL * U_L - F_R + F_L);

    // Скорость середины разрыва
    double SM = ((SR_m_uR)*rho_R*u_R - (SL_m_uL)*rho_L*u_L - p_R + p_L)/((SR_m_uR)*rho_R - (SL_m_uL)*rho_L);
    //double SM = U_star[1]/U_star[0];

    if (SL > 0) {
        return MHD_flux(U_L, gam_hcr);
    }
    else if (SL <= 0 && SM >= 0) {
        std::vector<double> F_L = MHD_flux(U_L, gam_hcr);
        std::vector<double> F_R = MHD_flux(U_R, gam_hcr);
        std::vector<double> U_star = 1/(SR-SL) * (SR * U_R - SL * U_L - F_R + F_L);
        double rho_L_star = rho_L * (SL - u_L)/(SL-SM);
        //double rho_L_star = U_star[0];
        double u_L_star = SM;
        double Bx_L_star = U_star[5];
        double By_L_star = U_star[6];
        double Bz_L_star = U_star[7];
        //double v_L_star = U_star[1]/rho_L_star;
        //double w_L_star = U_star[2]/rho_L_star;
        double v_L_star = v_L + Bx_L*(By_L-By_L_star)/(rho_L*SL_m_uL);
        double w_L_star = w_L + Bx_L*(Bz_L-Bz_L_star)/(rho_L*SL_m_uL);
        if(Bx_L == 0){
            By_L_star = By_L * SL_m_uL / (SL - SM);
            Bz_L_star = Bz_L * SL_m_uL / (SL - SM);
            v_L_star = v_L;
            w_L_star = w_L;
        }
        double pT_L = ptotal(p_L, Bx_L, By_L, Bz_L);
        double pT_star = pT_L + rho_L*SL_m_uL*(SM - u_L);
        double e_L_star = (SL_m_uL*e_L - pT_L*u_L + pT_star*SM + Bx_L *(u_L*Bx_L + v_L*By_L + w_L*Bz_L - u_L_star*Bx_L_star - v_L_star*By_L_star - w_L_star*Bz_L_star))/(SL-SM);
        //MHD_flux(rho, u, v, w, e, Bx, By, Bz, pT, gam_hcr)
        return MHD_flux(rho_L_star, u_L_star, v_L_star, w_L_star, e_L_star, Bx_L, By_L_star, Bz_L_star, pT_star, gam_hcr);
    }
    else if (SM <= 0 && SR >= 0){
/*          0      1      2      3    4   5   6   7
        * state:  rho,  rho*u, rho*v, rho*w, e, Bx, By, Bz
        * */
        std::vector<double> F_L = MHD_flux(U_L, gam_hcr);
        std::vector<double> F_R = MHD_flux(U_R, gam_hcr);
        std::vector<double> U_star = 1/(SR-SL) * (SR * U_R - SL * U_L - F_R + F_L);
        double rho_R_star = rho_R * (SR - u_R)/(SR-SM);
        //double rho_R_star = U_star[0];
        double u_R_star = SM;
        double Bx_R_star = U_star[5];
        double By_R_star = U_star[6];
        double Bz_R_star = U_star[7];
        //double v_R_star = U_star[1]/rho_R_star;
        //double w_R_star = U_star[2]/rho_R_star;
        double v_R_star = v_R + Bx_R*(By_R-By_R_star)/(rho_R*(SR-u_R));
        double w_R_star = w_R + Bx_R*(Bz_R-Bz_R_star)/(rho_R*(SR-u_R));
        if(Bx_R == 0){
            By_R_star = By_R * SR_m_uR / (SR - SM);
            Bz_R_star = Bz_R * SR_m_uR / (SR - SM);
            v_R_star = v_R;
            w_R_star = w_R;
        }
        double pT_R = ptotal(p_R, Bx_R, By_R, Bz_R);
        double pT_star = pT_R + rho_R*SR_m_uR*(SM - u_R);
        double e_R_star = (SR_m_uR*e_R - pT_R*u_R + pT_star*SM + Bx_R *(u_R*Bx_R + v_R*By_R + w_R*Bz_R - u_R_star*Bx_R_star - v_R_star*By_R_star - w_R_star*Bz_R_star))/(SR-SM);
        //MHD_flux(rho, u, v, w, e, Bx, By, Bz, pT, gam_hcr)
        return MHD_flux(rho_R_star, u_R_star, v_R_star, w_R_star, e_R_star, Bx_R, By_R_star, Bz_R_star, pT_star, gam_hcr);
    }
        //if (SR < 0)
    else  {
        //return F_R;
        return MHD_flux(U_R, gam_hcr);
    }
}

// Определяем HLLD поток F
std::vector<double> HLLD_flux(const std::vector<double>& U_L, const std::vector<double>& U_R, const double &gam_hcr) {
    /*          0      1      2      3    4   5   6   7
     * state:  rho,  rho*u, rho*v, rho*w, e, Bx, Bz, By
     */
    double rho_L = U_L[0];
    double u_L = U_L[1]/rho_L;
    double v_L = U_L[2]/rho_L;
    double w_L = U_L[3]/rho_L;
    double e_L = U_L[4];
    double Bx_L = U_L[5];
    double By_L = U_L[6];
    double Bz_L = U_L[7];
    double p_L = pressure(gam_hcr, e_L, rho_L, u_L, v_L, w_L, Bx_L, By_L, Bz_L);
    double pT_L = ptotal(p_L, Bx_L, By_L, Bz_L);

    double rho_R = U_R[0];
    double u_R = U_R[1]/rho_R;
    double v_R = U_R[2]/rho_R;
    double w_R = U_R[3]/rho_R;
    double e_R = U_R[4];
    double Bx_R = U_R[5];
    double By_R = U_R[6];
    double Bz_R = U_R[7];
    double p_R = pressure(gam_hcr, e_R, rho_R, u_R, v_R, w_R, Bx_R, By_R, Bz_R);
    double pT_R = ptotal(p_R, Bx_R, By_R, Bz_R);

    //быстрые магнитозвуковые скорости на левом и правом концах
    double cf_L = cfast(U_L, gam_hcr);
    double cf_R = cfast(U_R, gam_hcr);

    //double Bx = (Bx_L + Bx_R)/2;
    // или взять среднее по Роу
    double Bx = (std::sqrt(rho_L)*Bx_L + std::sqrt(rho_R)*Bx_R) / (std::sqrt(rho_L) + std::sqrt(rho_R));

    //скорость левого сигнала, рассчитываемая как минимальное значение скорости левого состояния (uL) и быстрой магнитозвуковой скорости (cfL).
    //double SL = std::min(u_L - cf_L, u_R - cf_R);
    double SL = std::min(u_L, u_R) - std::max(cf_L, cf_R);
    // Скорость правого сигнала, рассчитываемая как максимальное значение скорости правого состояния (uR) и быстрой магнитозвуковой скорости (cfR).ы
    //double SR = std::max(u_L + cf_L, u_R + cf_R);
    double SR = std::max(u_L, u_R) + std::max(cf_L, cf_R);

    double SR_m_uR = SR - u_R;
    double SL_m_uL = SL - u_L;
    double den_SM = SR_m_uR * rho_R - SL_m_uL * rho_L;

    // Скорость середины разрыва
    double SM = (SR_m_uR * rho_R * u_R - SL_m_uL * rho_L * u_L - pT_R + pT_L) / den_SM;

    double SM_m_uR = SM - u_R;
    double SM_m_uL = SM - u_L;

    double rho_L_star = rho_L;
    if(std::fabs(SL-SM) > 1e-16) {
        rho_L_star = rho_L * SL_m_uL / (SL - SM);
    }

    double rho_R_star = rho_R;
    if(std::fabs(SR-SM) > 1e-16){
        rho_R_star = rho_R * SR_m_uR / (SR - SM);
    }

    double SL_star = SM - std::fabs(Bx)/std::sqrt(rho_L_star);
    double SR_star = SM + std::fabs(Bx)/std::sqrt(rho_R_star);

    double pT_star = (SR_m_uR * rho_R * pT_L - SL_m_uL * rho_L * pT_R + rho_L * rho_R * SR_m_uR * SL_m_uL * (u_R - u_L)) / den_SM;

    // звёздочка слева
    double u_L_star = u_L;
    double pT_L_star = pT_L;
    double v_L_star = v_L;
    double w_L_star = w_L;
    double By_L_star = 0.0;
    double Bz_L_star = 0.0;
    double denom_L = rho_L * SL_m_uL * (SL - SM) - Bx * Bx;
    if(std::fabs(denom_L) > 1e-16){
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
    if(std::fabs(denom_R) > 1e-16){
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
        return MHD_flux(U_L, gam_hcr);
    }
    else if (/*SL <= 0 &&*/ SL_star >= 0.0) {
//!!!   //F*L
        return MHD_flux(rho_L_star, u_L_star, v_L_star, w_L_star, e_L_star, Bx, By_L_star, Bz_L_star, pT_L_star, gam_hcr);
    }
    else if (/*SL_star <= 0 &&*/ SM >= 0.0) {
//!!!   //F**L
        double u_L_2star = SM;
        double pT_L_2star = pT_L_star;
        double rho_star_sum_of_sqrts = std::sqrt(rho_L_star)+ std::sqrt(rho_R_star);
        double sign_Bx = signum(Bx_L);
        double v_L_2star = (std::sqrt(rho_L_star)*v_L_star + std::sqrt(rho_R_star)*v_R_star + (By_R_star-By_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double w_L_2star = (std::sqrt(rho_L_star)*w_L_star + std::sqrt(rho_R_star)*w_R_star + (Bz_R_star-Bz_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double By_L_2star = (std::sqrt(rho_L_star)*By_R_star + std::sqrt(rho_R_star)*By_L_star + std::sqrt(rho_L_star*rho_R_star)*(v_R_star-v_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double Bz_L_2star = (std::sqrt(rho_L_star)*Bz_R_star + std::sqrt(rho_R_star)*Bz_L_star + std::sqrt(rho_L_star*rho_R_star)*(w_R_star-w_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double e_L_2star = e_L_star - std::sqrt(rho_L_star)*(u_L_star * Bx + v_L_star * By_L_star + w_L_star * Bz_L_star - u_L_2star * Bx - v_L_2star * By_L_2star - w_L_2star * Bz_L_2star)*sign_Bx;
        return MHD_flux(rho_L_star, u_L_2star, v_L_2star, w_L_2star, e_L_2star, Bx, By_L_2star, Bz_L_2star, pT_L_2star, gam_hcr);
    }
    else if (/*SM <= 0 &&*/ SR_star >= 0.0){
//!!!   //F**R
        double u_R_2star = SM;
        double pT_R_2star = pT_R_star;
        double rho_star_sum_of_sqrts = std::sqrt(rho_L_star)+ std::sqrt(rho_R_star);
        double sign_Bx = signum(Bx_L);
        double v_R_2star = (std::sqrt(rho_L_star)*v_L_star + std::sqrt(rho_R_star)*v_R_star + (By_R_star-By_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double w_R_2star = (std::sqrt(rho_L_star)*w_L_star + std::sqrt(rho_R_star)*w_R_star + (Bz_R_star-Bz_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double By_R_2star = (std::sqrt(rho_L_star)*By_R_star + std::sqrt(rho_R_star)*By_L_star + std::sqrt(rho_L_star*rho_R_star)*(v_R_star-v_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double Bz_R_2star = (std::sqrt(rho_L_star)*Bz_R_star + std::sqrt(rho_R_star)*Bz_L_star + std::sqrt(rho_L_star*rho_R_star)*(w_R_star-w_L_star)*sign_Bx)/rho_star_sum_of_sqrts;
        double e_R_2star = e_R_star + std::sqrt(rho_R_star)*(u_R_star * Bx + v_R_star * By_R_star + w_R_star * Bz_R_star - u_R_2star * Bx - v_R_2star * By_R_2star - w_R_2star * Bz_R_2star)*sign_Bx;
        return MHD_flux(rho_R_star, u_R_2star, v_R_2star, w_R_2star, e_R_2star, Bx, By_R_2star, Bz_R_2star, pT_R_2star, gam_hcr);
    }
    else if (/*SR_star <= 0 &&*/ SR >= 0.0){
//!!!   //F*R
        return MHD_flux(rho_R_star, u_R_star, v_R_star, w_R_star, e_R_star, Bx, By_R_star, Bz_R_star, pT_R_star, gam_hcr);
    }
        //if (SR < 0)
    else  {
//!!!   //FR
        return MHD_flux(U_R, gam_hcr);
    }
}


/* Инициализация начального состояния системы (каждой точке пространства ставим в соотвествие величину начального отклонения)
 * args: problem
 * return: initial state
 * */
std::vector<std::vector<double>> initializeState(const MHDProblem1D &problem){
    std::vector<std::vector<double>> new_state(problem.num_space_steps+1, std::vector<double>(8,0));
    double x_i = problem.x0;

    for(int i = 0; i <= problem.num_space_steps; ++i){
        new_state[i] = problem.initStateFunc(x_i);
        x_i += problem.h;
    }

    return new_state;
}

double tau_from_cfl(const double& sigma, const double& h, const std::vector<std::vector<double>>& states, const int& num_space_steps, const double& gam_hcr){
    double u = std::fabs(states[0][1]/states[0][0]);
    double current_speed = u + cfast(states[0], gam_hcr);
    double max_speed = current_speed;
    for(int i = 1; i < num_space_steps; ++i){
        u = std::fabs(states[i][1]/states[i][0]);
        current_speed = u + cfast(states[i], gam_hcr);
        if (current_speed > max_speed){
            max_speed = current_speed;
        }
    }
    double optimalTau = sigma * h;
    if(max_speed > 0.0){
        optimalTau /= max_speed;
    }
    return optimalTau;
}

bool HLLScheme(const MHDProblem1D &problem, const std::string &filename) {

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting HLL Scheme..." << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open())
    {
        // Физические параметры
        double gam_hcr = problem.gam_hcr;
        double x0 = problem.x0;    // начало отсчёта по пространству
        double L = problem.L;      // характерный пространственный размер (длина струны)
        double X = x0 + L;         // координата правого края
        double t0 = problem.t0;    // время отсчёта
        double T = problem.T;      // время окончания
        double tau = problem.tau;  // шаг по времени
        double h = problem.h;      // шаг по пространству
        double gam_courant = problem.gam_courant;  // число Куранта

        // Шаги по времени и пространству
        int num_time_steps = problem.num_time_steps;
        int num_space_steps = problem.num_space_steps;

        double t_i = t0;
        double x_i = x0;

        std::vector<std::vector<double>> state_j = initializeState(problem);
        std::vector<std::vector<double>> state_jp(state_j);
        std::vector<std::vector<double>> fluxes(num_space_steps, std::vector<double>(8,0));

        // Запись первого слоя в файл
        fpoints << t_i << std::endl;
        for (int i = 0; i <= num_space_steps; ++i)
        {
            writeVectorToFile(fpoints, state_j[i], pressure(state_j[i], gam_hcr));
        }
        fpoints << "--------" << std::endl;

        //tau = tau_from_cfl(gam_courant, h, state_j, num_space_steps, gam_hcr);
        // Эволюция системы во времени
        while(t_i < T){
            //for(int j = 1; j <= num_time_steps; ++j) {
            tau = tau_from_cfl(gam_courant, h, state_j, num_space_steps, gam_hcr);
            t_i += tau;

            // HLL-потоки
#pragma parallel for
            for(int i = 0; i < num_space_steps; ++i){
                fluxes[i] = HLL_flux(state_j[i], state_j[i+1], gam_hcr);
            }

            // Обход пространства
            for (int i = 1; i < num_space_steps; ++i) {
                x_i += h;
                state_jp[i] = state_j[i] - tau/h * (fluxes[i]-fluxes[i-1]);
            }
            // Граничные условия слева
            if(!problem.periodic_boundaries) {
                state_jp[0] = problem.leftBoundaryFunction(t_i);
            }
            else{
                state_jp[0] = state_jp[num_space_steps-1];
            }

            // Граничные условия справа
            if(!problem.periodic_boundaries) {
                state_jp[num_space_steps] = problem.rightBoundaryFunction(t_i);
            }
            else{
                state_jp[num_space_steps] = state_jp[1];
            }
            // Запись в файл
            fpoints << t_i << std::endl;
            for (int i = 0; i <= num_space_steps; ++i)
            {
                writeVectorToFile(fpoints, state_jp[i], pressure(state_jp[i], gam_hcr));
            }
            fpoints << "--------" << std::endl;

            // j  jp   -> jp j
            state_j.swap(state_jp);
        }
        fpoints.close();

        std::cout << "LOG[INFO]: HLL-calculation is completed." << std::endl;

        return true;

    } else {
        std::cout << "log[ERROR]: Couldn't open or create a file." << std::endl;
        return false;
    }
};

bool HLLCScheme(const MHDProblem1D &problem, const std::string &filename) {

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting HLLC Scheme..." << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open())
    {
        // Физические параметры
        double gam_hcr = problem.gam_hcr;
        double x0 = problem.x0;    // начало отсчёта по пространству
        double L = problem.L;      // характерный пространственный размер (длина струны)
        double X = x0 + L;         // координата правого края
        double t0 = problem.t0;    // время отсчёта
        double T = problem.T;      // время окончания
        double tau = problem.tau;  // шаг по времени
        double h = problem.h;      // шаг по пространству
        double gam_courant = problem.gam_courant;  // число Куранта

        // Шаги по времени и пространству
        int num_time_steps = problem.num_time_steps;
        int num_space_steps = problem.num_space_steps;

        double t_i = t0;
        double x_i = x0;

        std::vector<std::vector<double>> state_j = initializeState(problem);
        std::vector<std::vector<double>> state_jp(state_j);
        std::vector<std::vector<double>> fluxes(num_space_steps, std::vector<double>(8,0));

        // Запись первого слоя в файл
        fpoints << t_i << std::endl;
        for (int i = 0; i <= num_space_steps; ++i)
        {
            writeVectorToFile(fpoints, state_j[i], pressure(state_j[i], gam_hcr));
        }
        fpoints << "--------" << std::endl;

        //tau = tau_from_cfl(gam_courant, h, state_j, num_space_steps, gam_hcr);
        // Эволюция системы во времени
        while(t_i < T){
            //for(int j = 1; j <= num_time_steps; ++j) {
            tau = tau_from_cfl(gam_courant, h, state_j, num_space_steps, gam_hcr);
            t_i += tau;

            // HLL-потоки
#pragma parallel for
            for(int i = 0; i < num_space_steps; ++i){
                fluxes[i] = HLLC_flux(state_j[i], state_j[i+1], gam_hcr);
            }

            // Обход пространства
            for (int i = 1; i < num_space_steps; ++i) {
                x_i += h;
                state_jp[i] = state_j[i] - tau/h * (fluxes[i]-fluxes[i-1]);
            }

            // Граничные условия слева
            if(!problem.periodic_boundaries) {
                state_jp[0] = problem.leftBoundaryFunction(t_i);
            }
            else{
                state_jp[0] = state_jp[num_space_steps-1];
            }

            // Граничные условия справа
            if(!problem.periodic_boundaries) {
                state_jp[num_space_steps] = problem.rightBoundaryFunction(t_i);
            }
            else{
                state_jp[num_space_steps] = state_jp[1];
            }

            // Запись в файл
            fpoints << t_i << std::endl;
            for (int i = 0; i <= num_space_steps; ++i)
            {
                writeVectorToFile(fpoints, state_jp[i], pressure(state_jp[i], gam_hcr));
            }
            fpoints << "--------" << std::endl;

            // j  jp   -> jp j
            state_j.swap(state_jp);
        }
        fpoints.close();

        std::cout << "LOG[INFO]: HLLC-calculation is completed." << std::endl;

        return true;

    } else {
        std::cout << "log[ERROR]: Couldn't open or create a file." << std::endl;
        return false;
    }
};

bool HLLDScheme(const MHDProblem1D &problem, const std::string &filename) {

    // Создание файла
    std::string path = "./OutputData/" + filename;
    std::ofstream fpoints(path);
    std::cout << "log[INFO]: Starting HLLD Scheme..." << std::endl;
    std::cout << "log[INFO]: Opening a file \"" << filename << "\" to write..." << std::endl;
    if (fpoints.is_open())
    {
        // Физические параметры
        double gam_hcr = problem.gam_hcr;
        double x0 = problem.x0;    // начало отсчёта по пространству
        double L = problem.L;      // характерный пространственный размер (длина струны)
        double X = x0 + L;         // координата правого края
        double t0 = problem.t0;    // время отсчёта
        double T = problem.T;      // время окончания
        double tau = problem.tau;  // шаг по времени
        double h = problem.h;      // шаг по пространству
        double gam_courant = problem.gam_courant;  // число Куранта

        // Шаги по времени и пространству
        int num_time_steps = problem.num_time_steps;
        int num_space_steps = problem.num_space_steps;

        double t_i = t0;
        double x_i = x0;

        std::vector<std::vector<double>> state_j = initializeState(problem);
        std::vector<std::vector<double>> state_jp(state_j);
        std::vector<std::vector<double>> fluxes(num_space_steps, std::vector<double>(8,0));

        // Запись первого слоя в файл
        fpoints << t_i << std::endl;
        for (int i = 0; i <= num_space_steps; ++i)
        {
            writeVectorToFile(fpoints, state_j[i], pressure(state_j[i], gam_hcr));
        }
        fpoints << "--------" << std::endl;

        //tau = tau_from_cfl(gam_courant, h, state_j, num_space_steps, gam_hcr);
        // Эволюция системы во времени
        while(t_i < T){
            //for(int j = 1; j <= num_time_steps; ++j) {
            tau = tau_from_cfl(gam_courant, h, state_j, num_space_steps, gam_hcr);
            t_i += tau;
            std::cout << t_i << std::endl;
            // HLLD-потоки
#pragma parallel for
            for(int i = 0; i < num_space_steps; ++i){
                fluxes[i] = HLLD_flux(state_j[i], state_j[i+1], gam_hcr);
            }

            // Обход пространства
            for (int i = 1; i < num_space_steps; ++i) {
                x_i += h;
                state_jp[i] = state_j[i] - tau/h * (fluxes[i]-fluxes[i-1]);
            }

            // Граничные условия слева
            if(!problem.periodic_boundaries) {
                state_jp[0] = problem.leftBoundaryFunction(t_i);
            }
            else{
                state_jp[0] = state_jp[num_space_steps-1];
            }

            // Граничные условия справа
            if(!problem.periodic_boundaries) {
                state_jp[num_space_steps] = problem.rightBoundaryFunction(t_i);
            }
            else{
                state_jp[num_space_steps] = state_jp[1];
            }

            // Запись в файл
            fpoints << t_i << std::endl;
            for (int i = 0; i <= num_space_steps; ++i)
            {
                writeVectorToFile(fpoints, state_jp[i], pressure(state_jp[i], gam_hcr));
            }
            fpoints << "--------" << std::endl;

            // j  jp   -> jp j
            state_j.swap(state_jp);
        }
        fpoints.close();

        std::cout << "LOG[INFO]: HLLD-calculation is completed." << std::endl;

        return true;

    } else {
        std::cout << "log[ERROR]: Couldn't open or create a file." << std::endl;
        return false;
    }
};