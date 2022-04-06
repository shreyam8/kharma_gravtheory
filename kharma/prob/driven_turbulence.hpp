#pragma once

#include "decs.hpp"
#include "types.hpp"

#include <parthenon/parthenon.hpp>

using namespace std;
using namespace parthenon;

TaskStatus InitializeDrivenTurbulence(MeshBlockData<Real> *rc, ParameterInput *pin)
{
    Flag(rc, "Initializing Driven Turbulence problem");
    auto pmb = rc->GetBlockPointer();
    GridScalar rho = rc->Get("prims.rho").data;
    GridScalar u = rc->Get("prims.u").data;

    const Real gam = pmb->packages.Get("GRMHD")->Param<Real>("gamma");
    const Real rho0 = pin->GetOrAddReal("driven_turbulence", "rho", 1.0);
    const Real cs0 = pin->GetOrAddReal("driven_turbulence", "cs0", 8.6e-4);
    const Real edot_frac = pin->GetOrAddReal("driven_turbulence", "edot_frac", 0.5);
    const Real x1min = pin->GetOrAddReal("parthenon/mesh", "x1min", 0);
    const Real x1max = pin->GetOrAddReal("parthenon/mesh", "x1max",  1);
    const Real x2min = pin->GetOrAddReal("parthenon/mesh", "x2min", 0);
    const Real x2max = pin->GetOrAddReal("parthenon/mesh", "x2max",  1);
    const Real x3min = pin->GetOrAddReal("parthenon/mesh", "x3min", -1);
    const Real x3max = pin->GetOrAddReal("parthenon/mesh", "x3max",  1);

    const Real edot = edot_frac * rho0 * pow(cs0, 3);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("drive_edot")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("drive_edot", edot);
    const Real lx1 = x1max-x1min;   const Real lx2 = x2max-x2min;
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("lx1")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("lx1", lx1);
    if(! (pmb->packages.Get("GRMHD")->AllParams().hasKey("lx2")))
        pmb->packages.Get("GRMHD")->AddParam<Real>("lx2", lx2);
    //adding for later use in create_grf

    const Real u0 = cs0 * cs0 * rho0 / (gam - 1) / gam; //from flux_functions.hpp
    IndexRange ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
    IndexRange jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
    IndexRange kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
    pmb->par_for("driven_turb_rho_u_init", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA_3D {
            rho(k, j, i) = rho0;
            u(k, j, i) = u0;
        }
    );
    return TaskStatus::complete;
}
