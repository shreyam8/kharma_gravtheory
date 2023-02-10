/* 
 *  File: inverter.cpp
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
#include "inverter.hpp"

// This will include headers in the correct order
#include "invert_template.hpp"

#include "reductions.hpp"

std::shared_ptr<KHARMAPackage> Inverter::Initialize(ParameterInput *pin, std::shared_ptr<Packages_t>& packages)
{
    auto pkg = std::make_shared<KHARMAPackage>("Inverter");
    Params &params = pkg->AllParams();

    // TODO TODO THESE ARE NO-OPS
    Real err_tol = pin->GetOrAddReal("inverter", "err_tol", 1e-8);
    params.Add("err_tol", err_tol);
    int iter_max = pin->GetOrAddInteger("inverter", "iter_max", 8);
    params.Add("iter_max", iter_max);
    Real stepsize = pin->GetOrAddReal("inverter", "stepsize", 1e-5);
    params.Add("stepsize", stepsize);

    std::string inverter_name = pin->GetOrAddString("inverter", "type", "onedw");
    if (inverter_name == "onedw") {
        params.Add("inverter_type", Type::onedw);
    } else if (inverter_name == "none") {
        params.Add("inverter_type", Type::none);
    }

    // Flag denoting UtoP inversion failures
    // Only needed if we're actually calling UtoP, but always allocated as it's retrieved often
    // Needs boundary sync if treating primitive variables as fundamental
    bool sync_prims = packages->Get("Driver")->Param<bool>("sync_prims");
    bool implicit_grmhd = packages->Get("GRMHD")->Param<bool>("implicit");
    Metadata m;
    if (sync_prims && !implicit_grmhd) {
        m = Metadata({Metadata::Real, Metadata::Cell, Metadata::Derived, Metadata::OneCopy, Metadata::FillGhost});
    } else {
        m = Metadata({Metadata::Real, Metadata::Cell, Metadata::Derived, Metadata::OneCopy});
    }
    pkg->AddField("pflag", m);

    // Don't operate if GRMHD variables are being evolved implicitly
    // This package is still loaded because fixes
    if (!implicit_grmhd) {
        pkg->BlockUtoP = Inverter::BlockUtoP;
    }

    pkg->PostStepDiagnosticsMesh = Inverter::PostStepDiagnostics;

    return pkg;
}

void Inverter::BlockUtoP(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse)
{
    auto& type = rc->GetBlockPointer()->packages.Get("Inverter")->Param<Type>("inverter_type");
    switch(type) {
    case Type::onedw:
        BlockInvert<Type::onedw>(rc, domain, coarse);
        break;
    case Type::none:
        break;
    }
}

TaskStatus Inverter::PostStepDiagnostics(const SimTime& tm, MeshData<Real> *md)
{
    Flag("Printing Floor diagnostics");
    auto pmesh = md->GetMeshPointer();
    auto pmb0 = md->GetBlockData(0)->GetBlockPointer();
    // Options
    const auto& pars = pmesh->packages.Get("Globals")->AllParams();
    const int flag_verbose = pars.Get<int>("flag_verbose");

    // Debugging/diagnostic info about floor and inversion flags
    if (flag_verbose >= 1) {
        Flag("Printing flags");
        Reductions::CountFlags(md, "pflag", Inverter::status_names, IndexDomain::interior, flag_verbose, false);
    }
    return TaskStatus::complete;
}
