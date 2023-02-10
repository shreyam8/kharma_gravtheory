/* 
 *  File: bondi.cpp
 *  
 *  BSD 3-Clause License
 *  
 *  Copyright (c) 2020, AFD Group at UIUC
 *  All rights reserved.
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1. Redistributions of source code must retain the above copyright notice, this
 *     list of conditions and the following disclaimer.
 *  
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  
 *  3. Neither the name of the copyright holder nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "bondi.hpp"

#include "flux_functions.hpp"

/**
 * Initialization of a Bondi problem with specified sonic point, BH mdot, and horizon radius
 * TODO mdot and rs are redundant and should be merged into one parameter. Uh, no.
 */
TaskStatus InitializeBondi(std::shared_ptr<MeshBlockData<Real>>& rc, ParameterInput *pin)
{
    Flag(rc, "Initializing Bondi problem");
    auto pmb = rc->GetBlockPointer();

    const Real mdot = pin->GetOrAddReal("bondi", "mdot", 1.0);
    const Real rs = pin->GetOrAddReal("bondi", "rs", 8.0);

    // By default, stay away from the EH if BH is spinning
    const Real rin_bondi_default = (pin->GetReal("coordinates", "a") > 0.1) ? 2.0 : 0.0;
    const Real rin_bondi = pin->GetOrAddReal("bondi", "rin", rin_bondi_default);

    // Whether to set magnetic field as a part of boundary conditions
    const Real set_b = pin->GetOrAddBoolean("bondi", "set_b", false);

    // Add these to package properties, since they continue to be needed on boundaries
    // TODO Problems need params
    if(! pmb->packages.Get("GRMHD")->AllParams().hasKey("mdot"))
        pmb->packages.Get("GRMHD")->AddParam<Real>("mdot", mdot);
    if(! pmb->packages.Get("GRMHD")->AllParams().hasKey("rs"))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rs", rs);
    if(! pmb->packages.Get("GRMHD")->AllParams().hasKey("rin_bondi"))
        pmb->packages.Get("GRMHD")->AddParam<Real>("rin_bondi", rin_bondi);

    if(! pmb->packages.Get("GRMHD")->AllParams().hasKey("set_b"))
        pmb->packages.Get("GRMHD")->AddParam<bool>("set_b", set_b);

    // Set this problem to control the outer X1 boundary
    // remember to disable inflow_check in parameter file!
    auto bound_pkg = static_cast<KHARMAPackage*>(pmb->packages.Get("Boundaries").get());
    bound_pkg->KHARMAOuterX1Boundary = SetBondi;

    // Set the interior domain to the analytic solution to begin
    // This tests that PostInitialize will correctly fill ghost zones with the boundary we set
    SetBondi(rc, IndexDomain::interior);

    Flag(rc, "Initialized");
    return TaskStatus::complete;
}

TaskStatus SetBondi(std::shared_ptr<MeshBlockData<Real>>& rc, IndexDomain domain, bool coarse)
{
    Flag(rc, "Setting Bondi zones");
    auto pmb = rc->GetBlockPointer();

    PackIndexMap prims_map, cons_map;
    auto P = GRMHD::PackMHDPrims(rc.get(), prims_map);
    auto U = GRMHD::PackMHDCons(rc.get(), cons_map);
    const VarMap m_u(cons_map, true), m_p(prims_map, false);

    const Real mdot = pmb->packages.Get("GRMHD")->Param<Real>("mdot");
    const Real rs = pmb->packages.Get("GRMHD")->Param<Real>("rs");
    const Real rin_bondi = pmb->packages.Get("GRMHD")->Param<Real>("rin_bondi");
    const Real gam = pmb->packages.Get("GRMHD")->Param<Real>("gamma");

    const EMHD::EMHD_parameters& emhd_params = EMHD::GetEMHDParameters(pmb->packages);

    // Just the X1 right boundary
    GRCoordinates G = pmb->coords;
    SphKSCoords ks = mpark::get<SphKSCoords>(G.coords.base);
    SphBLCoords bl = SphBLCoords(ks.a);
    CoordinateEmbedding cs = G.coords;

    // Set the Bondi conditions wherever we're asked
    auto bounds = coarse ? pmb->c_cellbounds : pmb->cellbounds;
    IndexRange ib = bounds.GetBoundsJ(domain);
    IndexRange jb = bounds.GetBoundsJ(domain);
    IndexRange kb = bounds.GetBoundsK(domain);
    pmb->par_for("bondi_boundary", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
            get_prim_bondi(G, cs, P, m_p, gam, bl, ks, mdot, rs, rin_bondi, k, j, i);
        }
    );

    if (pmb->packages.Get("GRMHD")->Param<bool>("set_b")) {
        pmb->par_for("bondi_boundary", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
            KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
                // This ignores rin_bondi to keep divB consistent
                // B \prop r^-3
                GReal Xembed[GR_DIM];
                G.coord_embed(k, j, i, Loci::center, Xembed);
                P(m_p.B1, k, j, i) = m::pow(Xembed[1], -3);
            }
        );
    }

    Flag(rc, "Set");
    return TaskStatus::complete;
}
