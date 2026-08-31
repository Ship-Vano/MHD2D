#include "MHDSolver1D.h"
#include "MHDSolver2D.h"
#include "NetGeometry.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(const bool condition, const std::string& message) {
    if (!condition) fail(message);
}

void require_near(const double actual, const double expected, const double tolerance,
                  const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        fail(message + ": actual=" + std::to_string(actual) +
             ", expected=" + std::to_string(expected));
    }
}

double max_abs_difference(const std::vector<double>& first, const std::vector<double>& second) {
    require(first.size() == second.size(), "state-vector sizes differ");
    double result = 0.0;
    for (std::size_t component = 0; component < first.size(); ++component) {
        result = std::max(result, std::abs(first[component] - second[component]));
    }
    return result;
}

void require_finite(const std::vector<double>& values, const std::string& label) {
    for (const double value : values) {
        require(std::isfinite(value), label + " contains a non-finite value");
    }
}

double fast_speed_for_face(const std::vector<double>& state, const double normalB,
                           const double gamma) {
    const double rho = state[0];
    const double pressureValue = pressure(state, gamma);
    const double totalB2 = state[5] * state[5] + state[6] * state[6] + state[7] * state[7];
    const double sound2 = gamma * pressureValue / rho;
    const double alfven2 = totalB2 / rho;
    const double normalAlfven2 = normalB * normalB / rho;
    const double sum = sound2 + alfven2;
    return std::sqrt(0.5 * (sum + std::sqrt(std::max(0.0, sum * sum -
                                                               4.0 * sound2 * normalAlfven2))));
}

std::vector<double> discrete_curl_from_nodes(const World& world,
                                             const double azX, const double azY) {
    const NodePool nodes = world.getNodePool();
    const EdgePool edges = world.getEdgePool();
    const int physicalEdges = edges.edgeCount - world.ghostElemCount * 2;
    std::vector<double> normalField(physicalEdges, 0.0);
    for (int edgeIndex = 0; edgeIndex < physicalEdges; ++edgeIndex) {
        const Edge& edge = edges.edges[edgeIndex];
        const auto az = [&](const Node& node) {
            return azX * node.pos.x + azY * node.pos.y;
        };
        normalField[edgeIndex] = edgeCtOrientation(edge, nodes) *
            (az(nodes.nodes[edge.nodeInd2]) - az(nodes.nodes[edge.nodeInd1])) / edge.length;
    }
    return normalField;
}

struct OneDimensionalCpResult {
    double relativeL1 = 0.0;
    std::size_t fallbacks = 0;
};

std::vector<double> cp_alfven_state_1d(const double x, const double gamma) {
    const double phase = 2.0 * std::acos(-1.0) * x;
    const double transverse = 0.1 * std::sin(phase);
    const double vertical = 0.1 * std::cos(phase);
    return state_from_primitive_vars(1.0, 0.0, transverse, vertical, 0.1,
                                     1.0, transverse, vertical, gamma);
}

OneDimensionalCpResult evolve_cp_alfven_1d(const int cellCount) {
    constexpr double gamma = 5.0 / 3.0;
    require(cellCount > 0, "1D CP test needs cells");
    std::vector<std::vector<double>> state(cellCount);
    std::vector<std::vector<double>> next(cellCount, std::vector<double>(8, 0.0));
    std::vector<std::vector<double>> flux(cellCount, std::vector<double>(8, 0.0));
    for (int cell = 0; cell < cellCount; ++cell) {
        state[cell] = cp_alfven_state_1d((static_cast<double>(cell) + 0.5) / cellCount, gamma);
    }

    // Ten substeps per cell over a unit period: the test deliberately keeps
    // forward Euler/first-order space, but isolates the corrected HLLD flux
    // from the unstructured CT coupling.
    const double dt = 0.1 / static_cast<double>(cellCount);
    OneDimensionalCpResult result;
    for (int step = 0; step < 10 * cellCount; ++step) {
        for (int face = 0; face < cellCount; ++face) {
            flux[face] = HLLD_flux_corrected(
                state[face], state[(face + 1) % cellCount], gamma, &result.fallbacks);
        }
        for (int cell = 0; cell < cellCount; ++cell) {
            const int leftFace = (cell + cellCount - 1) % cellCount;
            for (int component = 0; component < 8; ++component) {
                next[cell][component] = state[cell][component] -
                    dt * cellCount * (flux[cell][component] - flux[leftFace][component]);
            }
            require(next[cell][0] > 0.0 && pressure(next[cell], gamma) > 0.0,
                    "1D CP evolution produced an inadmissible state");
        }
        state.swap(next);
    }

    double accumulatedRelative = 0.0;
    for (const int component : {2, 3, 6, 7}) {
        double numerator = 0.0;
        double denominator = 0.0;
        for (int cell = 0; cell < cellCount; ++cell) {
            const std::vector<double> exact =
                cp_alfven_state_1d((static_cast<double>(cell) + 0.5) / cellCount, gamma);
            const double observed = component < 4 ? state[cell][component] / state[cell][0]
                                                  : state[cell][component];
            const double reference = component < 4 ? exact[component] / exact[0]
                                                    : exact[component];
            numerator += std::abs(observed - reference);
            denominator += std::abs(reference);
        }
        accumulatedRelative += numerator / denominator;
    }
    result.relativeL1 = accumulatedRelative / 4.0;
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: legacy_corrected_kernel_test <irregular_mesh> <sliver_mesh> <cw_mesh>\n";
        return 2;
    }
    try {
        constexpr double gamma = 5.0 / 3.0;
        // Conversion is a physics invariant shared by initialisation, HLLD
        // and VTU diagnostics: pressure recovered from the conservative
        // state must include both kinetic and magnetic energy exactly once.
        constexpr double conversionRho = 1.7;
        constexpr double conversionU = -0.31;
        constexpr double conversionV = 0.27;
        constexpr double conversionW = 0.13;
        constexpr double conversionP = 0.83;
        constexpr double conversionBx = -0.62;
        constexpr double conversionBy = 0.44;
        constexpr double conversionBz = 0.19;
        const std::vector<double> converted = state_from_primitive_vars(
            conversionRho, conversionU, conversionV, conversionW, conversionP,
            conversionBx, conversionBy, conversionBz, gamma);
        require_near(converted[0], conversionRho, 2.0e-15,
                     "primitive-to-conservative conversion changed density");
        require_near(converted[1] / converted[0], conversionU, 2.0e-15,
                     "primitive-to-conservative conversion changed x velocity");
        require_near(converted[2] / converted[0], conversionV, 2.0e-15,
                     "primitive-to-conservative conversion changed y velocity");
        require_near(converted[3] / converted[0], conversionW, 2.0e-15,
                     "primitive-to-conservative conversion changed z velocity");
        require_near(pressure(converted, gamma), conversionP, 2.0e-15,
                     "conservative-to-primitive pressure round trip failed");

        const std::vector<double> constant =
            state_from_primitive_vars(1.2, 0.3, -0.2, 0.1, 0.9, 0.7, -0.4, 0.2, gamma);
        std::size_t fallbackCount = 0;
        const std::vector<double> equalFlux =
            HLLD_flux_corrected(constant, constant, gamma, &fallbackCount);
        const std::vector<double> physicalFlux = MHD_flux(constant, gamma);
        require(max_abs_difference(equalFlux, physicalFlux) <= 3.0e-12,
                "HLLD consistency F(U,U)=F_phys failed");
        require(fallbackCount == 0, "HLLD fell back on an admissible equal state");
        require_near(equalFlux[5], 0.0, 3.0e-14, "HLLD must not transport normal B");

        const std::vector<double> leftSupersonic =
            state_from_primitive_vars(1.0, 10.0, 0.2, 0.0, 1.0, 0.5, 0.1, 0.0, gamma);
        const std::vector<double> rightSupersonic =
            state_from_primitive_vars(0.8, 9.0, -0.1, 0.0, 0.7, 0.5, -0.2, 0.0, gamma);
        require(max_abs_difference(HLLD_flux_corrected(leftSupersonic, rightSupersonic, gamma),
                                   MHD_flux(leftSupersonic, gamma)) <= 3.0e-12,
                "HLLD positive supersonic branch is not the left physical flux");
        const std::vector<double> leftNegative =
            state_from_primitive_vars(1.0, -10.0, 0.2, 0.0, 1.0, 0.5, 0.1, 0.0, gamma);
        const std::vector<double> rightNegative =
            state_from_primitive_vars(0.8, -9.0, -0.1, 0.0, 0.7, 0.5, -0.2, 0.0, gamma);
        require(max_abs_difference(HLLD_flux_corrected(leftNegative, rightNegative, gamma),
                                   MHD_flux(rightNegative, gamma)) <= 3.0e-12,
                "HLLD negative supersonic branch is not the right physical flux");

        const std::vector<double> brioLeft =
            state_from_primitive_vars(1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0, 2.0);
        const std::vector<double> brioRight =
            state_from_primitive_vars(0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0, 2.0);
        const std::vector<double> brioFlux = HLLD_flux_corrected(brioLeft, brioRight, 2.0);
        require_finite(brioFlux, "Brio--Wu HLLD flux");
        require_near(brioFlux[5], 0.0, 3.0e-14, "Brio--Wu HLLD changed normal B flux");

        const std::vector<double> zeroNormalLeft =
            state_from_primitive_vars(1.0, 0.0, 0.4, 0.0, 1.0, 0.0, 0.2, -0.1, gamma);
        const std::vector<double> zeroNormalRight =
            state_from_primitive_vars(0.9, 0.1, -0.3, 0.2, 0.8, 0.0, -0.2, 0.3, gamma);
        require_finite(HLLD_flux_corrected(zeroNormalLeft, zeroNormalRight, gamma),
                       "HLLD Bn=0 degeneracy flux");

        // A strictly positive near-vacuum state is deliberately not floored.
        // The HLLD intermediate denominator is ill-conditioned at this
        // scale, therefore the observable HLLE fallback is the safe result.
        const std::vector<double> dilute =
            state_from_primitive_vars(1.0e-30, 0.0, 0.0, 0.0, 1.0e-30,
                                      0.0, 0.0, 0.0, gamma);
        fallbackCount = 0;
        require_finite(HLLD_flux_corrected(dilute, dilute, gamma, &fallbackCount),
                       "HLLE fallback flux for positive near-vacuum state");
        require(fallbackCount == 1,
                "HLLD did not expose the HLLE fallback for a degenerate positive state");

        bool rejectedInvalidState = false;
        std::vector<double> invalid = constant;
        invalid[4] = 0.0;
        try {
            (void)HLLD_flux_corrected(invalid, constant, gamma);
        } catch (const std::domain_error&) {
            rejectedInvalidState = true;
        }
        require(rejectedInvalidState, "HLLD accepted an invalid state instead of reporting it");

        const OneDimensionalCpResult cp32 = evolve_cp_alfven_1d(32);
        const OneDimensionalCpResult cp64 = evolve_cp_alfven_1d(64);
        const OneDimensionalCpResult cp128 = evolve_cp_alfven_1d(128);
        require(cp32.fallbacks == 0 && cp64.fallbacks == 0 && cp128.fallbacks == 0,
                "smooth 1D CP wave unexpectedly used HLLE fallback");
        require(cp64.relativeL1 < cp32.relativeL1 && cp128.relativeL1 < cp64.relativeL1,
                "corrected HLLD does not converge for the isolated 1D CP wave");
        std::cout << "1D CP relative L1: N=32 " << cp32.relativeL1
                  << ", N=64 " << cp64.relativeL1
                  << ", N=128 " << cp128.relativeL1 << '\n';

        World irregular(argv[1], false);
        NeighbourService& irregularNeighbours = irregular.getNeighbourService();
        require_near(irregular.minX, 0.0, 1.0e-15,
                     "World minX depends on the first (interior) mesh node");
        require_near(irregular.maxX, 1.0, 1.0e-15, "World maxX is wrong");
        require_near(irregular.minY, 0.0, 1.0e-15, "World minY is wrong");
        require_near(irregular.maxY, 1.0, 1.0e-15, "World maxY is wrong");
        require(irregular.ghostElemCount > 4, "irregular mesh did not create boundary ghosts");
        require(static_cast<int>(irregularNeighbours.boundaryEdgeToGhostElement.size()) ==
                    irregular.ghostElemCount,
                "boundary-edge to ghost mapping is incomplete");
        require(static_cast<int>(irregularNeighbours.ghostElementToBoundaryEdge.size()) ==
                    irregular.ghostElemCount,
                "ghost to boundary-edge mapping is incomplete");
        for (const auto& [edge, ghost] : irregularNeighbours.boundaryEdgeToGhostElement) {
            const auto inverse = irregularNeighbours.ghostElementToBoundaryEdge.find(ghost);
            require(inverse != irregularNeighbours.ghostElementToBoundaryEdge.end() &&
                        inverse->second == edge,
                    "boundary/ghost maps are not inverse");
        }

        MHDSolver2D divergenceSolver(irregular);
        const EdgePool irregularEdges = irregular.getEdgePool();
        const int irregularPhysicalEdges =
            irregularEdges.edgeCount - irregular.ghostElemCount * 2;
        divergenceSolver.innerElemCount =
            irregular.getElementPool().elCount - irregular.ghostElemCount;
        divergenceSolver.innerEdgeCount = irregularPhysicalEdges;
        divergenceSolver.bNs.assign(irregularPhysicalEdges, 0.0);
        divergenceSolver.bNs = discrete_curl_from_nodes(irregular, -0.3, 0.75);
        const NodePool irregularNodes = irregular.getNodePool();
        const DivergenceMetrics initialCurl =
            divergenceSolver.computeDivergenceMetrics(1.0);
        require(initialCurl.maxFlux <= 5.0e-14,
                "nodal vector potential is not discretely divergence-free");
        require(initialCurl.maxAbs <= 5.0e-12,
                "discrete curl produced non-roundoff divergence");

        // RT0 must reproduce a constant field exactly from its normal face
        // fluxes on an actually irregular triangle mesh.  This is stronger
        // than checking div B alone: a telescoping but incorrectly oriented
        // reconstruction could still have roundoff divergence.
        constexpr double uniformBx = 0.73;
        constexpr double uniformBy = -0.41;
        std::vector<double> constantFaceField(irregularPhysicalEdges, 0.0);
        for (int edgeIndex = 0; edgeIndex < irregularPhysicalEdges; ++edgeIndex) {
            const Edge& edge = irregularEdges.edges[edgeIndex];
            constantFaceField[edgeIndex] =
                uniformBx * edge.normalVector.x + uniformBy * edge.normalVector.y;
        }
        const ElementPool irregularElements = irregular.getElementPool();
        for (int elementIndex = 0; elementIndex < divergenceSolver.innerElemCount; ++elementIndex) {
            const Vec2 reconstructed = reconstructCellMagneticFieldRT0(
                irregularElements.elements[elementIndex], irregularNodes, irregularEdges,
                constantFaceField, irregularPhysicalEdges);
            require_near(reconstructed.x, uniformBx, 2.0e-12,
                         "RT0 reconstruction does not recover constant Bx");
            require_near(reconstructed.y, uniformBy, 2.0e-12,
                         "RT0 reconstruction does not recover constant By");
        }

        // Input element winding is not a physical degree of freedom.  A
        // clockwise copy of the same mesh must retain the physical B.n and
        // the discrete curl invariant, including after one shared CT update.
        World clockwise(argv[3], false);
        const EdgePool clockwiseEdges = clockwise.getEdgePool();
        const ElementPool clockwiseElements = clockwise.getElementPool();
        const NodePool clockwiseNodes = clockwise.getNodePool();
        const int clockwisePhysicalEdges = clockwiseEdges.edgeCount - clockwise.ghostElemCount * 2;
        const int clockwisePhysicalElements = clockwiseElements.elCount - clockwise.ghostElemCount;
        require(clockwisePhysicalEdges == irregularPhysicalEdges,
                "CW mesh changed the physical edge count");
        for (int edgeIndex = 0; edgeIndex < clockwisePhysicalEdges; ++edgeIndex) {
            const Edge& edge = clockwiseEdges.edges[edgeIndex];
            const double orientation = edgeCtOrientation(edge, clockwiseNodes);
            require(orientation > 0.0,
                    "CW text mesh was not canonicalized before CT topology construction");
        }

        // The production loader canonicalizes text cells, but the CT helper
        // must still retain its explicit orientation factor: old binary meshes
        // and future topology sources may carry the opposite edge tangent.
        Edge oppositeTangent = clockwiseEdges.edges[0];
        std::swap(oppositeTangent.nodeInd1, oppositeTangent.nodeInd2);
        require(edgeCtOrientation(oppositeTangent, clockwiseNodes) < 0.0,
                "negative CT orientation branch was not detected");
        constexpr double previousBn = 0.37;
        constexpr double firstNodeG = 0.2;
        constexpr double secondNodeG = -0.1;
        require_near(
            advanceFaceNormalBFromCt(oppositeTangent, clockwiseNodes, previousBn, 0.125,
                                     firstNodeG, secondNodeG),
            previousBn - 0.125 / oppositeTangent.length * (secondNodeG - firstNodeG),
            2.0e-14, "negative CT orientation update has the wrong sign");
        std::vector<double> clockwiseCurl =
            discrete_curl_from_nodes(clockwise, -uniformBy, uniformBx);
        for (int edgeIndex = 0; edgeIndex < clockwisePhysicalEdges; ++edgeIndex) {
            const Edge& edge = clockwiseEdges.edges[edgeIndex];
            const double expectedBn =
                uniformBx * edge.normalVector.x + uniformBy * edge.normalVector.y;
            require_near(clockwiseCurl[edgeIndex], expectedBn, 3.0e-12,
                         "potential-derived Bn changes under CW input winding");
        }
        MHDSolver2D clockwiseSolver(clockwise);
        clockwiseSolver.innerElemCount = clockwisePhysicalElements;
        clockwiseSolver.innerEdgeCount = clockwisePhysicalEdges;
        clockwiseSolver.bNs = clockwiseCurl;
        require(clockwiseSolver.computeDivergenceMetrics(1.0).maxFlux <= 5.0e-14,
                "CW potential field is not discretely divergence-free");
        for (int edgeIndex = 0; edgeIndex < clockwisePhysicalEdges; ++edgeIndex) {
            const Edge& edge = clockwiseEdges.edges[edgeIndex];
            const double psiFirst = std::sin(1.7 * clockwiseNodes.nodes[edge.nodeInd1].pos.x) +
                                    0.4 * clockwiseNodes.nodes[edge.nodeInd1].pos.y;
            const double psiSecond = std::sin(1.7 * clockwiseNodes.nodes[edge.nodeInd2].pos.x) +
                                     0.4 * clockwiseNodes.nodes[edge.nodeInd2].pos.y;
            clockwiseSolver.bNs[edgeIndex] = advanceFaceNormalBFromCt(
                edge, clockwiseNodes, clockwiseSolver.bNs[edgeIndex], 0.125,
                psiFirst, psiSecond);
        }
        require(clockwiseSolver.computeDivergenceMetrics(1.0).maxFlux <= 8.0e-14,
                "CW CT update violates the discrete divergence invariant");

        // The CT Faraday increment itself is another edge-wise nodal
        // difference, so its divergence must telescope on every triangle.
        for (int edgeIndex = 0; edgeIndex < irregularPhysicalEdges; ++edgeIndex) {
            const Edge& edge = irregularEdges.edges[edgeIndex];
            const double psiFirst = std::sin(1.7 * irregularNodes.nodes[edge.nodeInd1].pos.x) +
                                    0.4 * irregularNodes.nodes[edge.nodeInd1].pos.y;
            const double psiSecond = std::sin(1.7 * irregularNodes.nodes[edge.nodeInd2].pos.x) +
                                     0.4 * irregularNodes.nodes[edge.nodeInd2].pos.y;
            divergenceSolver.bNs[edgeIndex] = advanceFaceNormalBFromCt(
                edge, irregularNodes, divergenceSolver.bNs[edgeIndex], 0.125,
                psiFirst, psiSecond);
        }
        const DivergenceMetrics updatedCurl =
            divergenceSolver.computeDivergenceMetrics(1.0);
        require(updatedCurl.maxFlux <= 8.0e-14,
                "edge update from nodal EMF violated discrete divergence");

        int interiorEdge = -1;
        for (int edgeIndex = 0; edgeIndex < irregularPhysicalEdges; ++edgeIndex) {
            if (irregularEdges.edges[edgeIndex].neighbourInd2 >= 0) {
                interiorEdge = edgeIndex;
                break;
            }
        }
        require(interiorEdge >= 0, "irregular mesh has no interior edge");
        divergenceSolver.bNs[interiorEdge] += 1.0e-3;
        const DivergenceMetrics perturbed =
            divergenceSolver.computeDivergenceMetrics(1.0);
        require(perturbed.maxAbs > 1.0e-4,
                "max(abs(div B)) did not detect a signed interior face perturbation");

        World sliver(argv[2], false);
        MHDSolver2D sliverSolver(sliver);
        sliverSolver.cflNum = 0.4;
        const EdgePool sliverEdges = sliver.getEdgePool();
        const ElementPool sliverElements = sliver.getElementPool();
        const int sliverPhysicalElements = sliverElements.elCount - sliver.ghostElemCount;
        const int sliverPhysicalEdges = sliverEdges.edgeCount - sliver.ghostElemCount * 2;
        const std::vector<double> uniform =
            state_from_primitive_vars(1.0, 0.0, 0.0, 0.0, 1.0, 0.5, -0.2, 0.1, gamma);
        std::vector<std::vector<double>> uniformStates(sliverPhysicalElements, uniform);
        std::vector<double> sliverBn(sliverPhysicalEdges, 0.0);
        double maxSignal = 0.0;
        for (int edgeIndex = 0; edgeIndex < sliverPhysicalEdges; ++edgeIndex) {
            const Edge& edge = sliverEdges.edges[edgeIndex];
            sliverBn[edgeIndex] =
                uniform[5] * edge.normalVector.x + uniform[6] * edge.normalVector.y;
            maxSignal = std::max(maxSignal, fast_speed_for_face(uniform, sliverBn[edgeIndex], gamma));
        }
        double expectedCflDt = std::numeric_limits<double>::infinity();
        for (int elementIndex = 0; elementIndex < sliverPhysicalElements; ++elementIndex) {
            double spectralPerimeter = 0.0;
            for (const int edgeIndex : sliverElements.elements[elementIndex].edgeIndexes) {
                spectralPerimeter += sliverEdges.edges[edgeIndex].length *
                                     fast_speed_for_face(uniform, sliverBn[edgeIndex], gamma);
            }
            expectedCflDt = std::min(expectedCflDt,
                                     0.4 * sliverElements.elements[elementIndex].area / spectralPerimeter);
        }
        const double computedCflDt =
            sliverSolver.stableTimeStepUnstructured(uniformStates, sliverBn, gamma);
        require_near(computedCflDt, expectedCflDt,
                     2.0e-13 * std::max(1.0, expectedCflDt),
                     "unstructured CFL bound differs from A/sum(lambda*edge)");
        const double historicalMinEdgeDt = 0.4 * sliverEdges.minEdgeLen / maxSignal;
        require(computedCflDt < historicalMinEdgeDt,
                "sliver-triangle CFL did not become stricter than min-edge CFL");
        sliverSolver.cflNum = 1.0e-8;
        require(sliverSolver.stableTimeStepUnstructured(uniformStates, sliverBn, gamma) < 1.0e-7,
                "unstructured CFL still has a hidden lower time-step floor");

        std::cout << "legacy_corrected kernel and geometry regressions passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "legacy_corrected regression failure: " << error.what() << '\n';
        return 1;
    }
}
