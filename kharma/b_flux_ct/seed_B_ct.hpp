// Seed a torus of some type with a magnetic field according to its density
#pragma once

#include "decs.hpp"
#include "types.hpp"

namespace B_FluxCT
{

/**
 * Seed a divergence-free magnetic field of user's choice, optionally
 * proportional to existing fluid density.
 * Updates primitive and conserved variables.
 */
TaskStatus SeedBField(MeshBlockData<Real> *rc, ParameterInput *pin);

KOKKOS_INLINE_FUNCTION void averaged_curl_3D(const GRCoordinates& G, const GridVector& A, const GridVector& B_U,
                                             const int& k, const int& j, const int& i)
{
    // Take a flux-ct step from the corner potentials.
    // This needs to be 3D because post-tilt A may not point in the phi direction only

    // A3,2 derivative
    const Real A3c2f = (A(V3, k, j + 1, i)     + A(V3, k, j + 1, i + 1) + 
                        A(V3, k + 1, j + 1, i) + A(V3, k + 1, j + 1, i + 1)) / 4;
    const Real A3c2b = (A(V3, k, j, i)     + A(V3, k, j, i + 1) +
                        A(V3, k + 1, j, i) + A(V3, k + 1, j, i + 1)) / 4;
    // A2,3 derivative
    const Real A2c3f = (A(V2, k + 1, j, i)     + A(V2, k + 1, j, i + 1) +
                        A(V2, k + 1, j + 1, i) + A(V2, k + 1, j + 1, i + 1)) / 4;
    const Real A2c3b = (A(V2, k, j, i)     + A(V2, k, j, i + 1) +
                        A(V2, k, j + 1, i) + A(V2, k, j + 1, i + 1)) / 4;
    B_U(V1, k, j, i) = (A3c2f - A3c2b) / G.Dxc<2>(j) - (A2c3f - A2c3b) / G.Dxc<3>(k);

    // A1,3 derivative
    const Real A1c3f = (A(V1, k + 1, j, i)     + A(V1, k + 1, j, i + 1) + 
                        A(V1, k + 1, j + 1, i) + A(V1, k + 1, j + 1, i + 1)) / 4;
    const Real A1c3b = (A(V1, k, j, i)     + A(V1, k, j, i + 1) +
                        A(V1, k, j + 1, i) + A(V1, k, j + 1, i + 1)) / 4;
    // A3,1 derivative
    const Real A3c1f = (A(V3, k, j, i + 1)     + A(V3, k + 1, j, i + 1) +
                        A(V3, k, j + 1, i + 1) + A(V3, k + 1, j + 1, i + 1)) / 4;
    const Real A3c1b = (A(V3, k, j, i)     + A(V3, k + 1, j, i) +
                        A(V3, k, j + 1, i) + A(V3, k + 1, j + 1, i)) / 4;
    B_U(V2, k, j, i) = (A1c3f - A1c3b) / G.Dxc<3>(k) - (A3c1f - A3c1b) / G.Dxc<1>(i);

    // A2,1 derivative
    const Real A2c1f = (A(V2, k, j, i + 1)     + A(V2, k, j + 1, i + 1) + 
                        A(V2, k + 1, j, i + 1) + A(V2, k + 1, j + 1, i + 1)) / 4;
    const Real A2c1b = (A(V2, k, j, i)     + A(V2, k, j + 1, i) +
                        A(V2, k + 1, j, i) + A(V2, k + 1, j + 1, i)) / 4;
    // A1,2 derivative
    const Real A1c2f = (A(V1, k, j + 1, i)     + A(V1, k, j + 1, i + 1) +
                        A(V1, k + 1, j + 1, i) + A(V1, k + 1, j + 1, i + 1)) / 4;
    const Real A1c2b = (A(V1, k, j, i)     + A(V1, k, j, i + 1) +
                        A(V1, k + 1, j, i) + A(V1, k + 1, j, i + 1)) / 4;
    B_U(V3, k, j, i) = (A2c1f - A2c1b) / G.Dxc<1>(i) - (A1c2f - A1c2b) / G.Dxc<2>(j);
}

KOKKOS_INLINE_FUNCTION void averaged_curl_2D(const GRCoordinates& G, const GridVector& A, const GridVector& B_U,
                                             const int& k, const int& j, const int& i)
{
    // A3,2 derivative
    const Real A3c2f = (A(V3, k, j + 1, i) + A(V3, k, j + 1, i + 1)) / 2;
    const Real A3c2b = (A(V3, k, j, i)     + A(V3, k, j, i + 1)) / 2;
    B_U(V1, k, j, i) = (A3c2f - A3c2b) / G.Dxc<2>(j);

    // A3,1 derivative
    const Real A3c1f = (A(V3, k, j, i + 1) + A(V3, k, j + 1, i + 1)) / 2;
    const Real A3c1b = (A(V3, k, j, i)     + A(V3, k, j + 1, i)) / 2;
    B_U(V2, k, j, i) = - (A3c1f - A3c1b) / G.Dxc<1>(i);

    B_U(V3, k, j, i) = 0;
}

} // namespace B_FluxCT
