#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include<vector>
#include <algorithm>

void computeHLLDFluxes(std::vector<std::vector<double>>& fluxes, std::vector<std::vector<double>>& unrotated_fluxes);