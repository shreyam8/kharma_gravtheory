/* 
 *  File: inverter.hpp
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
#pragma once

#include "decs.hpp"
#include "types.hpp"

// Implementation of u_to_p must be defined before instantiation below
// Additionally, invert_template contains the Type and Status enums
#include "invert_template.hpp"
#include "onedw.hpp"

#include "pack.hpp"

using namespace parthenon;

/**
 * Recover primitive variables from conserved forms.
 * Currently can only use the 1D_W scheme of Noble et al. (2006),
 * but this is the spot for alternate implementations
 */
namespace Inverter {

std::shared_ptr<KHARMAPackage> Initialize(ParameterInput *pin, std::shared_ptr<Packages_t>& packages);

/**
 * Get the primitive variables
 * This just computes P, and only for the GRHD fluid varaibles rho, u, uvec.
 *
 * Defaults to entire domain, as the KHARMA algorithm relies on applying UtoP over ghost zones.
 * 
 * input: U, whatever form
 * output: U and P match down to inversion errors
 */
void BlockUtoP(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse);

/**
 * Smooth over inversion failures, usually by averaging values of the primitive variables from each neighboring zone
 * a.k.a. Diffusion?  What diffusion?  There is no diffusion here.
 * 
 * LOCKSTEP: this function expects and should preserve P<->U
 */
TaskStatus FixUtoP(MeshBlockData<Real> *rc);

/**
 * Print details of any inversion failures or fixed zones
 */
TaskStatus PostStepDiagnostics(const SimTime& tm, MeshData<Real> *md);

/**
 * Internal inversion fn, templated on inverter type.  Calls through to templated u_to_p
 * This is called with the correct template argument from BlockUtoP
 */
template<Type inverter>
inline void BlockInvert(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse)
{
    Flag(rc, "Filling Primitives");
    auto pmb = rc->GetBlockPointer();
    const auto& G = pmb->coords;

    PackIndexMap prims_map, cons_map;
    auto U = GRMHD::PackMHDCons(rc, cons_map);
    auto P = GRMHD::PackHDPrims(rc, prims_map);
    const VarMap m_u(cons_map, true), m_p(prims_map, false);

    GridScalar pflag = rc->Get("pflag").data;

    const Real gam = pmb->packages.Get("GRMHD")->Param<Real>("gamma");

    const Real err_tol = pmb->packages.Get("Inverter")->Param<Real>("err_tol");
    const int iter_max = pmb->packages.Get("Inverter")->Param<int>("iter_max");
    const Real stepsize = pmb->packages.Get("Inverter")->Param<Real>("stepsize");

    // Get the primitives from our conserved versions
    // Currently this runs over *all* zones, including all ghosts, even
    // uninitialized zones which are still zero.  We select for initialized
    // zones only in the loop below, to avoid failures to converge while
    // calculating primtive vars over as much of the domain as possible
    // We could (did formerly) save some time here by running over
    // only zones with initialized conserved variables, but the domain
    // of such values is not rectangular in the current handling
    auto bounds = coarse ? pmb->c_cellbounds : pmb->cellbounds;
    const IndexRange ib = bounds.GetBoundsI(domain);
    const IndexRange jb = bounds.GetBoundsJ(domain);
    const IndexRange kb = bounds.GetBoundsK(domain);
    const IndexRange ib_b = bounds.GetBoundsI(IndexDomain::interior);
    const IndexRange jb_b = bounds.GetBoundsJ(IndexDomain::interior);
    const IndexRange kb_b = bounds.GetBoundsK(IndexDomain::interior);

    pmb->par_for("U_to_P", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
            if (inside(k, j, i, kb_b, jb_b, ib_b) ||
                m::abs(P(m_p.RHO, k, j, i)) > SMALL || m::abs(P(m_p.UU, k, j, i)) > SMALL) {
                // Run over all interior zones and any initialized ghosts
                pflag(k, j, i) = (double) Inverter::u_to_p<inverter>(G, U, m_u, gam, k, j, i, P, m_p, Loci::center);
            } else {
                // Don't *use* un-initialized zones for fixes, but also don't *fix* them
                pflag(k, j, i) = -1;
            }
        }
    );
    Flag(rc, "Filled");
}

}