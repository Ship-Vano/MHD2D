//
// Created by Иван on 4/30/2024.
//

#include "ShockTube1D.h"

ShockTube1D::ShockTube1D(const double &gam_hcr_init, const std::vector<double> &L_state, const std::vector<double> &R_state, const double &mid_point) {
    if(L_state.size() == 8 && R_state.size() == 8){
        gam_hcr = gam_hcr_init;

        rhoL = L_state[0];
        uL = L_state[1];
        vL = L_state[2];
        wL = L_state[3];
        pL = L_state[4];
        BxL = L_state[5];
        ByL = L_state[6];
        BzL = L_state[7];

        x_mid = mid_point;

        rhoR = R_state[0];
        uR = R_state[1];
        vR = R_state[2];
        wR = R_state[3];
        pR = R_state[4];
        BxR = R_state[5];
        ByR = R_state[6];
        BzR = R_state[7];

        initDistrib = ([&](double x) {
            if(x < x_mid){
                return state_from_primitive_vars(rhoL, uL, vL, wL, pL, BxL, ByL, BzL, gam_hcr);
            }
            else{
                return state_from_primitive_vars(rhoR, uR, vR, wR, pR, BxR, ByR, BzR, gam_hcr);
            }
        });
        initDistrib_is_set = true;

        leftBound =([&](double t){
            return state_from_primitive_vars(rhoL, uL, vL, wL, pL, BxL, ByL, BzL, gam_hcr);
        });
        leftBound_is_set = true;

        rightBound = ([&](double t){
            return state_from_primitive_vars(rhoR, uR, vR, wR, pR, BxR, ByR, BzR, gam_hcr);
        });
        rightBound_is_set = true;
    }
    else{
        std::cout << "LOG[ERROR]: inappropriate size of vectors! The size of both must be 8." << std::endl;
    }
}
