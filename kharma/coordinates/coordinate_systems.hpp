/* 
 *  File: coordinate_systems.hpp
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

// See note in coordintate_embedding
#include <mpark/variant.hpp>
//#include <variant>
//namespace mpark = std;

#include "decs.hpp"

#include "matrix.hpp"
#include "kharma_utils.hpp"
#include "root_find.hpp"

#define LEGACY_TH 1

/**
 * Embedding/Base systems implemented:
 * Minkowski space: Cartesian and Spherical coordinates
 * Kerr Space: Spherical KS and BL coordinates
 * 
 * Transformations:
 * Nulls in Cartesian and Spherical coordinates
 * "Modified": r=exp(x1), th=pi*x2 + (1-hslope)*etc
 * "Funky" modified: additional non-invertible cylindrization of th
 * 
 * TODO Cartesian KS base
 * TODO snake coordinate transform for Cartesian Minkowski
 * TODO CMKS, MKS3 transforms, proper Cartesian<->Spherical functions stolen from e.g. coordinate_utils.hpp
 * TODO overhaul the LEGACY_TH stuff
 * TODO currently avoids returning gcov which might be singular,
 *      is this the correct play vs handling in inversions?
 */

/**
 * EMBEDDING SYSTEMS:
 * These are the usual systems of coordinates for different spacetimes.
 * Each system/class must define at least gcov_embed, returning the metric in terms of their own coordinates Xembed
 * Some extra convenience classes have been defined for some systems.
 */
// ____________________________________________________________________________________________________________________________

/**
 * Cartesian Coordinates over flat space.
 */
class CartMinkowskiCoords {
    public:
        static constexpr char name[] = "CartMinkowskiCoords";
        static constexpr bool spherical = false;
        static constexpr GReal a = 0.0;
        KOKKOS_INLINE_FUNCTION void gcov_embed(const GReal Xembed[GR_DIM], Real gcov[GR_DIM][GR_DIM]) const
        {
            DLOOP2 gcov[mu][nu] = (mu == nu) - 2*(mu == 0 && nu == 0);
        }
};
// ____________________________________________________________________________________________________________________________

/**
 * Spherical coordinates for flat space
 */
class SphMinkowskiCoords {
    public:
        static constexpr char name[] = "SphMinkowskiCoords";
        static constexpr bool spherical = true;
        static constexpr GReal a = 0.0;
        KOKKOS_INLINE_FUNCTION void gcov_embed(const GReal Xembed[GR_DIM], Real gcov[GR_DIM][GR_DIM]) const
        {
            const GReal r = m::max(Xembed[1], SMALL);
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);
            const GReal sth = m::sin(th);

            gzero2(gcov);
            gcov[0][0] = 1.;
            gcov[1][1] = 1.;
            gcov[2][2] = r*r;
            gcov[3][3] = sth*sth*r*r;
        }
};
// ____________________________________________________________________________________________________________________________

/**
 * Spherical Kerr-Schild coordinates
 */
class SphKSCoords {
    public:
        static constexpr char name[] = "SphKSCoords";
        // BH Spin is a property of KS
        const GReal a;
        static constexpr bool spherical = true;

        KOKKOS_FUNCTION SphKSCoords(GReal spin): a(spin) {};

        KOKKOS_INLINE_FUNCTION void gcov_embed(const GReal Xembed[GR_DIM], Real gcov[GR_DIM][GR_DIM]) const
        {
            const GReal r = Xembed[1];
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);

            const GReal cth = m::cos(th);
            const GReal sth = m::sin(th);
            const GReal sin2 = sth*sth;
            const GReal rho2 = r*r + a*a*cth*cth;

            gcov[0][0] = -1. + 2.*r/rho2;
            gcov[0][1] = 2.*r/rho2;
            gcov[0][2] = 0.;
            gcov[0][3] = -2.*a*r*sin2/rho2;

            gcov[1][0] = 2.*r/rho2;
            gcov[1][1] = 1. + 2.*r/rho2;
            gcov[1][2] = 0.;
            gcov[1][3] = -a*sin2*(1. + 2.*r/rho2);

            gcov[2][0] = 0.;
            gcov[2][1] = 0.;
            gcov[2][2] = rho2;
            gcov[2][3] = 0.;

            gcov[3][0] = -2.*a*r*sin2/rho2;
            gcov[3][1] = -a*sin2*(1. + 2.*r/rho2);
            gcov[3][2] = 0.;
            gcov[3][3] = sin2*(rho2 + a*a*sin2*(1. + 2.*r/rho2));
        }

        // For converting from BL
        KOKKOS_INLINE_FUNCTION void vec_from_bl(const GReal Xembed[GR_DIM], const Real vcon_bl[GR_DIM], Real vcon[GR_DIM]) const
        {
            GReal r = Xembed[1];
            Real trans[GR_DIM][GR_DIM];
            DLOOP2 trans[mu][nu] = (mu == nu);
            trans[0][1] = 2.*r/(r*r - 2.*r + a*a);
            trans[3][1] = a/(r*r - 2.*r + a*a);

            gzero(vcon);
            DLOOP2 vcon[mu] += trans[mu][nu]*vcon_bl[nu];
        }

        KOKKOS_INLINE_FUNCTION void vec_to_bl(const GReal Xembed[GR_DIM], const Real vcon_bl[GR_DIM], Real vcon[GR_DIM]) const
        {
            GReal r = Xembed[1];
            GReal rtrans[GR_DIM][GR_DIM], trans[GR_DIM][GR_DIM];
            DLOOP2 rtrans[mu][nu] = (mu == nu);
            rtrans[0][1] = 2.*r/(r*r - 2.*r + a*a);
            rtrans[3][1] = a/(r*r - 2.*r + a*a);

            invert(&rtrans[0][0], &trans[0][0]);

            gzero(vcon);
            DLOOP2 vcon[mu] += trans[mu][nu]*vcon_bl[nu];
        }
};
// ____________________________________________________________________________________________________________________________

// Question : Could I potentially replace all the places that SphKSExtG is used with DCSKSCoords? 
// ____________________________________________________________________________________________________________________________

/**
 * Spherical Kerr-Schild coordinates w/ external gravity term
 */
class SphKSExtG {
    public:
        static constexpr char name[] = "SphKSExtG";
        // BH Spin is a property of KS
        const GReal a;
        static constexpr bool spherical = true;

        static constexpr GReal A = 1.46797639e-8;
        static constexpr GReal B = 1.29411117;

        KOKKOS_FUNCTION SphKSExtG(GReal spin): a(spin) {};

        KOKKOS_INLINE_FUNCTION void gcov_embed(const GReal Xembed[GR_DIM], Real gcov[GR_DIM][GR_DIM]) const
        {
            const GReal r = Xembed[1];
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);

            const GReal cth = m::cos(th);
            const GReal sth = m::sin(th);
            const GReal sin2 = sth*sth;
            const GReal rho2 = r*r + a*a*cth*cth;

            const GReal Phi_g = (A / (B - 1.)) * (m::pow(r, B-1.) - m::pow(2, B-1.));

            gcov[0][0] = -1. + 2.*r/rho2 - 2. * Phi_g;
            gcov[0][1] = 2.*r/rho2 - 2. * Phi_g;
            gcov[0][2] = 0.;
            gcov[0][3] = -2.*a*r*sin2/rho2;

            gcov[1][0] = 2.*r/rho2 - 2. * Phi_g;
            gcov[1][1] = 1. + 2.*r/rho2 - 2. * Phi_g;
            gcov[1][2] = 0.;
            gcov[1][3] = -a*sin2*(1. + 2.*r/rho2);

            gcov[2][0] = 0.;
            gcov[2][1] = 0.;
            gcov[2][2] = rho2;
            gcov[2][3] = 0.;

            gcov[3][0] = -2.*a*r*sin2/rho2;
            gcov[3][1] = -a*sin2*(1. + 2.*r/rho2);
            gcov[3][2] = 0.;
            gcov[3][3] = sin2*(rho2 + a*a*sin2*(1. + 2.*r/rho2));
        }

        // For converting from BL
        // TODO will we ever need a from_ks?
        KOKKOS_INLINE_FUNCTION void vec_from_bl(const GReal Xembed[GR_DIM], const Real vcon_bl[GR_DIM], Real vcon[GR_DIM]) const
        {
            GReal r = Xembed[1];
            Real trans[GR_DIM][GR_DIM];
            DLOOP2 trans[mu][nu] = (mu == nu);

            // external gravity from GIZMO
            const GReal Phi_g = (A/(B-1.)) * (m::pow(r,B-1.)-m::pow(2,B-1.));

            trans[0][1] = (2./r - 2.*Phi_g)/(1. - 2./r + 2.*Phi_g);
            trans[3][1] = a/(r*r - 2.*r + a*a);

            gzero(vcon);
            DLOOP2 vcon[mu] += trans[mu][nu]*vcon_bl[nu];
        }

        KOKKOS_INLINE_FUNCTION void vec_to_bl(const GReal Xembed[GR_DIM], const Real vcon_bl[GR_DIM], Real vcon[GR_DIM]) const
        {
            GReal r = Xembed[1];
            GReal rtrans[GR_DIM][GR_DIM], trans[GR_DIM][GR_DIM];
            DLOOP2 rtrans[mu][nu] = (mu == nu);

            const GReal Phi_g = (A / (B-1.)) * (m::pow(r, B-1.) - m::pow(2, B-1.));

            rtrans[0][1] = (2./r - 2.*Phi_g)/(1. - 2./r + 2.*Phi_g);
            rtrans[3][1] = a/(r*r - 2.*r + a*a);

            invert(&rtrans[0][0], &trans[0][0]);

            gzero(vcon);
            DLOOP2 vcon[mu] += trans[mu][nu]*vcon_bl[nu];
        }
};

// ____________________________________________________________________________________________________________________________

/**
 * Boyer-Lindquist coordinates as an embedding system
 */
class SphBLCoords {
    public:
        static constexpr char name[] = "SphBLCoords";
        // BH Spin is a property of BL
        const GReal a;
        static constexpr bool spherical = true;

        KOKKOS_FUNCTION SphBLCoords(GReal spin): a(spin) {}

        KOKKOS_INLINE_FUNCTION void gcov_embed(const GReal Xembed[GR_DIM], Real gcov[GR_DIM][GR_DIM]) const
        {
            const GReal r = Xembed[1];
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);
            const GReal cth = m::cos(th), sth = m::sin(th);

            const GReal sin2 = sth*sth;
            const GReal a2 = a*a;
            const GReal r2 = r*r;
            // TODO(BSP) this and gcov_embed for KS should look more similar...
            const GReal mmu = 1. + a2*cth*cth/r2; // mu is taken as an index

            gzero2(gcov);
            gcov[0][0]  = -(1. - 2./(r*mmu));
            gcov[0][3]  = -2.*a*sin2/(r*mmu);
            gcov[1][1]   = mmu/(1. - 2./r + a2/r2);
            gcov[2][2]   = r2*mmu;
            gcov[3][0]  = -2.*a*sin2/(r*mmu);
            gcov[3][3]   = sin2*(r2 + a2 + 2.*a2*sin2/(r*mmu));
        }

        // TODO(BSP) vec to/from ks, put guaranteed ks/bl fns into embedding

};
// ____________________________________________________________________________________________________________________________

/**
 * Boyer-Lindquist coordinates as an embedding system
 */
class SphBLExtG {
    public:
        static constexpr char name[] = "SphBLExtG";
        // BH Spin is a property of BL
        const GReal a;
        static constexpr bool spherical = true;

        static constexpr GReal A = 1.46797639e-8;
        static constexpr GReal B = 1.29411117;

        KOKKOS_FUNCTION SphBLExtG(GReal spin): a(spin) {}

        KOKKOS_INLINE_FUNCTION void gcov_embed(const GReal Xembed[GR_DIM], Real gcov[GR_DIM][GR_DIM]) const
        {
            const GReal r = Xembed[1];
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);
            const GReal cth = m::cos(th), sth = m::sin(th);

            const GReal sin2 = sth*sth;
            const GReal a2 = a*a;
            const GReal r2 = r*r;
            const GReal mmu = 1. + a2*cth*cth/r2; // mu is taken as an index

            const GReal Phi_g = (A / (B-1.)) * (m::pow(r, B-1.) - m::pow(2, B-1.));

            gzero2(gcov);
            gcov[0][0]  = -(1. - 2./(r*mmu)) - 2. * Phi_g;;
            gcov[0][3]  = -2.*a*sin2/(r*mmu);
            gcov[1][1]   = mmu / (1. - 2./r + 2.*Phi_g);
            gcov[2][2]   = r2*mmu;
            gcov[3][0]  = -2.*a*sin2/(r*mmu);
            gcov[3][3]   = sin2*(r2 + a2 + 2.*a2*sin2/(r*mmu));
        }
};
// Write a class DCSKSCoords with a new gcov definition. 
// ____________________________________________________________________________________________________________________________
/**
 * Changes Made. DCS theory KS coordinates. 
 */
class DCSKSCoords {
    public:
        static constexpr char name[] = "DCSKSCoords";
        // BH Spin is a property of KS
        const GReal a;
        const GReal zeta;
        static constexpr bool spherical = true;

        KOKKOS_FUNCTION DCSKSCoords(GReal spin, GReal z): a(spin), zeta(z) {} //semicolon here ?

        KOKKOS_INLINE_FUNCTION void gcov_embed(const GReal Xembed[GR_DIM], Real gcov[GR_DIM][GR_DIM]) const
        {
            const GReal r = Xembed[1];
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);
            
            // Assign gcov matrix to zero. 

            const GReal cth = m::cos(th);
            const GReal sth = m::sin(th);
            const GReal s2t = sth*sth;
            const GReal c2t = cth*cth ; 
            const GReal c4t = c2t*c2t ;
            const GReal rho2 = r*r + a*a*cth*cth;
            const GReal ep2 = a*a ;
            const GReal ep4 = ep2*ep2 ;

            gzero2(gcov);

            gcov[0][0] = -1. + 2.*r/rho2 + zeta*((ep2*(-338688.*(-1. + 3.*c2t) + 4221.*pow(r,6.)*(-1. + 3.*c2t) + 4221.*pow(r,7.)*(-1. + 3.*c2t) - 420.*pow(r,2.)*(263. + 321.*c2t) + 20.*pow(r,3.)*(-5428. + 2025.*c2t) - 20.*pow(r,4.)*(1523. + 2781.*c2t) + 
            168.*r*(-275. + 3471.*c2t) + 2.*pow(r,5.)*(-2482. + 6711.*c2t)))/(37632.*pow(r,10.)) + (ep4*(-436560.*pow(r,10.)*(-1. + 3.*c2t) - 
            436560.*pow(r,11.)*(-1. + 3.*c2t) + pow(r,5.)*(-27367138. + 338395440.*c2t - 79842960.*c4t) + pow(r,7.)*(5994779. + 9180120.*c2t - 23096700.*c4t) + pow(r,8.)*(1835988. + 14588424.*c2t - 20269020.*c4t) + 
            pow(r,9.)*(609036. + 6695880.*c2t - 14028580.*c4t) - 3048192.*(299. - 1440.*c2t + 720.*c4t) + 1512.*r*(200377. - 262200.*c2t + 546900.*c4t) + 540.*pow(r,2.)*(209695. - 1641704.*c2t + 3850652.*c4t) + 36.*pow(r,4.)*(-1309617. + 14023212.*c2t + 
            5328620.*c4t) + pow(r,6.)*(1493875. + 66961128.*c2t + 66934880.*c4t) - 12.*pow(r,3.)*(-7413845. - 9229392.*c2t + 156818740.*c4t)))/(13547520.*pow(r,14.))) ; 
            
            gcov[0][1] = 2.*r/rho2 + zeta*((ep2*(8442.*pow(r,7.) + pow(r,4.)*(504690. - 543240.*c2t) + pow(r,6.)*(55419. - 50652.*c2t) - 338688.*(-1. + 3.*c2t) - 2016.*pow(r,2.)*(-747. + 1291.*c2t) - 168.*r*(-7789. + 20721.*c2t) - 20.*pow(r,3)*(-50957. + 73764.*c2t) 
            - 2*pow(r,5.)*(-94018. + 97629.*c2t)))/(37632*pow(r,10)) + (ep4*(-1746240*pow(r,11) + 240*pow(r,10)*(-48727 + 43656*c2t) - 6096384*(299 - 1440*c2t + 720*c4t) + 8*pow(r,9)*(-4902855 + 3171874*c2t + 332305*c4t) - 3456*pow(r,2)*(5955596 - 12196010*c2t + 3150155*c4t) - 
            3024*r*(3416327 - 17156040*c2t + 8162220*c4t) + 96*pow(r,3)*(-80010145 + 69893298*c2t + 34916840*c4t) + 8*pow(r,7)*(-55039417 - 5982330*c2t + 69327885*c4t) + pow(r,8)*(-164813529 + 
            45059288*c2t + 116931460*c4t) + 36*pow(r,4)*(-27676299 - 98275560*c2t + 239844820*c4t) + 8*pow(r,6)*(-76467017 - 93284976*c2t + 253740770*c4t) + 4*pow(r,5)*(-83689817 - 868048872*c2t + 1345549080*c4t)))/(27095040*pow(r,14))); 
            
            gcov[0][3] = -2.*a*r*s2t/rho2 + zeta*(-1./112.*((189. + 120.*r + 70.*pow(r,2.))*a*(-1. + c2t))/pow(r,6.) + (pow(a,3.)*(-1. + c2t)*(-10160640.*(-1. + 3.*c2t) + 105828.*pow(r,7.)*(-1. + 5.*c2t) + 30.*pow(r,5.)*(2247. + 95.*c2t) + 
            55440.*r*(-25. + 453.*c2t) + 15.*pow(r,6.)*(-271. + 24875.*c2t) + 360.*pow(r,2.)*(-9205. + 64257.*c2t) + 16.*pow(r,4.)*(-40883. + 112735.*c2t) + 8.*pow(r,3.)*(-382792. + 1712045.*c2t)))/(1128960.*pow(r,10.)));


            gcov[1][0] = gcov[0][1] ;
            gcov[1][1] = 1. + 2.*r/rho2 + zeta*(-1/37632.*(ep2*(4221.*pow(r,7.)*(-5. + 3.*c2t) + 338688.*(-1. + 3.*c2t) + 3820.*pow(r,3.)*(-562. + 783.*c2t) + 21.*pow(r,6.)*(-5479. + 5427.*c2t) + 168.*r*(-15853. + 44913.*c2t) + 
            20.*pow(r,4.)*(-51992. + 51543.*c2t) + 84.*pow(r,2.)*(-37171. + 60363.*c2t) + 6.*pow(r,5.)*(-63506. + 67323.*c2t)))/pow(r,10.) + (ep4*(436560.*pow(r,11.)*(-5. + 3.*c2t) + 240.*pow(r,10.)*(-50546. + 49113.*c2t) - 3048192.*(299. - 1440.*c2t + 720.*c4t) + 
            4.*pow(r,9.)*(-9957969. + 4669778.*c2t + 4171755.*c4t) - 1512.*r*(7033031. - 34574280.*c2t + 16871340.*c4t) + 5.*pow(r,7.)*(-89262023. - 11407752.*c2t + 115543956.*c4t) + 72.*pow(r,4.)*(-13183341. - 56149386.*c2t + 117258100.*c4t) - 
            108.*pow(r,2.)*(191627547. - 398480840.*c2t + 120058220.*c4t) + pow(r,8.)*(-166649517. + 30470864.*c2t + 137200480.*c4t) + 12.*pow(r,3.)*(-647495005. + 549916992.*c2t + 436153460.*c4t) + pow(r,6.)*(-613230011. - 813240936.*c2t + 1962991280.*c4t) + 
            pow(r,5.)*(-307392130. - 3810590928.*c2t + 5462039280.*c4t)))/(13547520.*pow(r,14.)));  
    
            gcov[1][3] = -a*s2t*(1. + 2.*r/rho2) + zeta*(((571536 + 393540*r + 219380*pow(r,2) + 64660*pow(r,3) + 19880*pow(r,4) + 4221*pow(r,5))*a*(-1 + c2t))/(12544*pow(r,6)) - (pow(a,3)*(-1 + c2t)*(1570950*pow(r,8) + 327420*pow(r,9) - 
            60963840*(-49 + 48*c2t) - 7560*r*(-133087 + 6297*c2t) + 1620*pow(r,2)*(47721 + 449461*c2t) + 25*pow(r,6)*(711795 + 481847*c2t) + 30*pow(r,5)*(1256960 + 2241169*c2t) + pow(r,7)*(4312974 + 2699780*c2t) + 120*pow(r,3)*(-870409 + 5949739*c2t) + 
            12*pow(r,4)*(379503 + 21372890*c2t)))/(3386880*pow(r,10))); 

   
            gcov[2][2] = rho2 - (zeta*(1080.*pow(r,4.)*(338688. + 80808.*r + 67380.*pow(r,2.) + 10360.*pow(r,3.) + 18908.*pow(r,4.) + 9940.*pow(r,5.) + 4221.*pow(r,6.))*ep2*(1. - 3.*c2t) + ep4*(3141900*pow(r,9)*(-1 + 3*c2t) + 1309680*pow(r,10)*
            (-1 + 3*c2t) - 9144576*(299 - 1440*c2t + 720*c4t) - 4536*r*(41543 - 947400*c2t + 420780*c4t) + 1620*pow(r,2)*(-93041 - 85688*c2t + 1454300*c4t) + 12*pow(r,8)*(-535899 - 116700*c2t + 2367475*c4t) + 18*pow(r,5)*(-7227915 + 2913384*c2t + 17254600*c4t) + 
            36*pow(r,3)*(4260605 - 30961368*c2t + 27684860*c4t) + 3*pow(r,6)*(-20890109 + 18841800*c2t + 46796480*c4t) + pow(r,7)*(-20549463 + 6081360*c2t + 59683820*c4t) + 4*pow(r,4)*(-33335811 - 120639012*c2t + 227696380*c4t))))/(40642560*pow(r,11)) ; 
            gcov[2][3] = 0.;

            gcov[3][0] = gcov[0][3] ; 
            gcov[3][1] = gcov[1][3] ;
  
            gcov[3][3] = s2t*(rho2 + a*a*s2t*(1. + 2.*r/rho2)) + zeta*(-1./37632.*((338688. + 80808.*r + 67380.*pow(r,2.) + 10360.*pow(r,3.) + 18908.*pow(r,4.) + 9940.*pow(r,5.) + 4221.*pow(r,6.))*ep2*(-1. + c2t)*(-1. + 3.*c2t))/pow(r,7.) + (ep4*(-1. + c2t)*(3141900.*pow(r,9.)*(-1. + 3.*c2t) + 
            1309680.*pow(r,10.)*(-1. + 3.*c2t) - 9144576.*(299. - 1440.*c2t + 720.*c4t) + 4536.*r*(119737. + 302280.*c2t + 63060.*c4t) + 1620.*pow(r,2.)*(240495. - 1419832.*c2t + 2454908.*c4t) + 12.*pow(r,8.)*(-156009. - 1636260.*c2t + 
            3507145.*c4t) + 18.*pow(r,5.)*(-4337355. - 8648856.*c2t + 25926280.*c4t) + 36.*pow(r,3.)*(10727645. - 56829528.*c2t + 47085980.*c4t) + 3.*pow(r,6.)*(-6926429. - 37012920.*c2t + 88687520.*c4t) + pow(r,7.)*(-696903. - 73328880.*c2t + 119241500.*
            c4t) + 4.*pow(r,4.)*(-9548811. - 215787012.*c2t + 299057380.*c4t)))/(40642560.*pow(r,11.))) ;

        }

    // New Transformation matrix here !!! 

        KOKKOS_INLINE_FUNCTION void vec_from_bl(const GReal Xembed[GR_DIM], const Real vcon_bl[GR_DIM], Real vcon[GR_DIM]) const
        {
            GReal r = Xembed[1];
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);
            Real trans[GR_DIM][GR_DIM];
            DLOOP2 trans[mu][nu] = (mu == nu);

            const GReal cth = m::cos(th);
            const GReal sth = m::sin(th);
            const GReal s2t = sth*sth;
            const GReal c2t = cth*cth ; 
            const GReal c4t = c2t*c2t ;
            const GReal ep2 = a*a ;
            const GReal ep4 = ep2*ep2 ;

            trans[0][0] = 1;
            trans[0][1] = 2.*r/(r*r - 2.*r + a*a) + zeta*(((4741632 - 469728*r - 320544*pow(r,2) - 551200*pow(r,3) - 7190*pow(r,4) + 128*pow(r,5) + 4767*pow(r,6) + 8442*pow(r,7))*ep2)/(37632*pow((-2 + r),2)*pow(r,8)) + (ep4*(2275440*pow(r,11) - 1746240*pow(r,12) + 
            134120448*(1141 - 1440*c2t + 720*c4t) + 1632*pow(r,3)*(4447141 - 2110000*c2t + 1055000*c4t) + 8*pow(r,10)*(693009 - 3839450*c2t + 1919725*c4t) + 10368*pow(r,2)*(3412537 - 7500540*c2t + 3750270*c4t) - 16*pow(r,7)*(-15723816 - 16243430*c2t + 8121715*c4t) - 5*pow(r,9)*(10616445 - 18370504*c2t + 
            9185252*c4t) - 12096*r*(13535587 - 22161480*c2t + 11080740*c4t) - 16*pow(r,6)*(38276799 - 63133580*c2t + 31566790*c4t) + 4*pow(r,8)*(52778599 - 79370940*c2t + 39685470*c4t) + 192*pow(r,4)*(12270982 - 88661110*c2t + 44330555*c4t) - 4*pow(r,5)*(615913343 - 1841915640*c2t + 920957820*c4t)))/
            (27095040*pow((-2 + r),3)*pow(r,12))) ;
            trans[0][2] = 0;
            trans[0][3] = 0;

            trans[1][0] = 0 ;
            trans[1][1] = 1 ;
            trans[1][2] = 0 ;
            trans[1][3] = 0 ;

            trans[2][0] = 0 ;
            trans[2][1] = 0 ;
            trans[2][2] = 1 ;
            trans[2][3] = 0 ;

            trans[3][0] = 0;
            trans[3][1] = a/(r*r - 2.*r + a*a) + zeta*(((-1185408 - 242424*r - 60900*pow(r,2) + 90060*pow(r,3) + 24900*pow(r,4) + 11438*pow(r,5) + 4221*pow(r,6))*a)/(12544*(-2 + r)*pow(r,8)) + (pow(a,3)*(-261270*pow(r,10) - 327420*pow(r,11) + pow(r,6)*(61143744 - 31673420*c2t) + 12070840320*
            (-1 + c2t) - 997920*r*(-7591 + 10741*c2t) - 6480*pow(r,2)*(-54353 + 77003*c2t) - 34*pow(r,9)*(2901 + 112925*c2t) - 30*pow(r,7)*(-535687 + 1069338*c2t) - 24*pow(r,3)*(1595277 + 1427905*c2t) + 24*pow(r,5)*(732089 + 2621680*c2t) + pow(r,8)*(-8895627 + 3802665*c2t) + 
            36*pow(r,4)*(-12720549 + 29077735*c2t)))/(3386880*pow((-2 + r),2)*pow(r,12)));
            trans[3][2] = 0;
            trans[3][3] = 1;

            gzero(vcon);
            DLOOP2 vcon[mu] += trans[mu][nu]*vcon_bl[nu];
        }

        KOKKOS_INLINE_FUNCTION void vec_to_bl(const GReal Xembed[GR_DIM], const Real vcon_bl[GR_DIM], Real vcon[GR_DIM]) const
        {
            GReal r = Xembed[1];
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);
            GReal rtrans[GR_DIM][GR_DIM], trans[GR_DIM][GR_DIM];
            DLOOP2 rtrans[mu][nu] = (mu == nu);

            const GReal cth = m::cos(th);
            const GReal sth = m::sin(th);
            const GReal s2t = sth*sth;
            const GReal c2t = cth*cth ; 
            const GReal c4t = c2t*c2t ;
            const GReal ep2 = a*a ;
            const GReal ep4 = ep2*ep2 ;

            trans[0][0] = 1;
            trans[0][1] = 2.*r/(r*r - 2.*r + a*a) + zeta*(((4741632 - 469728*r - 320544*pow(r,2) - 551200*pow(r,3) - 7190*pow(r,4) + 128*pow(r,5) + 4767*pow(r,6) + 8442*pow(r,7))*ep2)/(37632*pow((-2 + r),2)*pow(r,8)) + (ep4*(2275440*pow(r,11) - 1746240*pow(r,12) + 
            134120448*(1141 - 1440*c2t + 720*c4t) + 1632*pow(r,3)*(4447141 - 2110000*c2t + 1055000*c4t) + 8*pow(r,10)*(693009 - 3839450*c2t + 1919725*c4t) + 10368*pow(r,2)*(3412537 - 7500540*c2t + 3750270*c4t) - 16*pow(r,7)*(-15723816 - 16243430*c2t + 8121715*c4t) - 5*pow(r,9)*(10616445 - 18370504*c2t + 
            9185252*c4t) - 12096*r*(13535587 - 22161480*c2t + 11080740*c4t) - 16*pow(r,6)*(38276799 - 63133580*c2t + 31566790*c4t) + 4*pow(r,8)*(52778599 - 79370940*c2t + 39685470*c4t) + 192*pow(r,4)*(12270982 - 88661110*c2t + 44330555*c4t) - 4*pow(r,5)*(615913343 - 1841915640*c2t + 920957820*c4t)))/
            (27095040*pow((-2 + r),3)*pow(r,12))) ;
            trans[0][2] = 0;
            trans[0][3] = 0;

            trans[1][0] = 0 ;
            trans[1][1] = 1 ;
            trans[1][2] = 0 ;
            trans[1][3] = 0 ;

            trans[2][0] = 0 ;
            trans[2][1] = 0 ;
            trans[2][2] = 1 ;
            trans[2][3] = 0 ;

            trans[3][0] = 0;
            trans[3][1] = a/(r*r - 2.*r + a*a) + zeta*(((-1185408 - 242424*r - 60900*pow(r,2) + 90060*pow(r,3) + 24900*pow(r,4) + 11438*pow(r,5) + 4221*pow(r,6))*a)/(12544*(-2 + r)*pow(r,8)) + (pow(a,3)*(-261270*pow(r,10) - 327420*pow(r,11) + pow(r,6)*(61143744 - 31673420*c2t) + 12070840320*
            (-1 + c2t) - 997920*r*(-7591 + 10741*c2t) - 6480*pow(r,2)*(-54353 + 77003*c2t) - 34*pow(r,9)*(2901 + 112925*c2t) - 30*pow(r,7)*(-535687 + 1069338*c2t) - 24*pow(r,3)*(1595277 + 1427905*c2t) + 24*pow(r,5)*(732089 + 2621680*c2t) + pow(r,8)*(-8895627 + 3802665*c2t) + 
            36*pow(r,4)*(-12720549 + 29077735*c2t)))/(3386880*pow((-2 + r),2)*pow(r,12)));
            trans[3][2] = 0;
            trans[3][3] = 1;

            invert(&rtrans[0][0], &trans[0][0]); // INVERTING BECAUSE IT IS NOW FROM KS TO BL!!! 

            gzero(vcon);
            DLOOP2 vcon[mu] += trans[mu][nu]*vcon_bl[nu];
        }
};

// ____________________________________________________________________________________________________________________________

/**
 * Changes Made. DCS theory BL coordinates. 
 */
class DCSBLCoords {
    public:
        static constexpr char name[] = "DCSBLCoords";

        // BH Spin is a property of KS
        const GReal a;
        const GReal zeta;
        static constexpr bool spherical = true;

        KOKKOS_FUNCTION DCSBLCoords(GReal spin, GReal z): a(spin), zeta(z) {}

        KOKKOS_INLINE_FUNCTION void gcov_embed(const GReal Xembed[GR_DIM], Real gcov[GR_DIM][GR_DIM]) const
        {
            const GReal r = Xembed[1];
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);
            
            // Assign gcov matrix to zero.
            gzero2(gcov); 

            const GReal cth = m::cos(th);
            const GReal sth = m::sin(th);
            const GReal s2t = sth*sth;
            const GReal c2t = cth*cth; 
            const GReal c4t = c2t*c2t;
            const GReal rho2 = r*r + a*a*cth*cth;
            const GReal ep2 = a*a;
            const GReal ep4 = ep2*ep2;
            const GReal DD = 1. - 2./r + ep2/(r*r);
            const GReal mu = 1. + ep2*cth*cth/(r*r);


            gcov[0][0] = -(1. - 2./(r*mu))+ zeta*(-1./13547520.*(ep2*(-121927680.*pow(r,4.) + 16632000.*pow(r,5.) + 39765600.*pow(r,6.) + 39081600.*pow(r,7.) + 10965600.*pow(r,8.) + 1787040.*pow(r,9.) + 1519560.*pow(r,10.) + 1519560.*pow(r,11.) + 365783040.*pow(r,4.)*c2t - 
            209926080.*pow(r,5.)*c2t + 48535200.*pow(r,6.)*c2t - 14580000.*pow(r,7.)*c2t + 20023200.*pow(r,8.)*c2t - 4831920.*pow(r,9.)*c2t - 4558680.*pow(r,10.)*c2t - 4558680.*pow(r,11.)*c2t))/pow(r,14.) - (ep4*(911409408. - 302970024.*r - 113235300.*pow(r,2.) - 88966140.*pow(r,3.) + 
            47146212.*pow(r,4) + 27367138.*pow(r,5.) - 1493875.*pow(r,6.) - 5994779.*pow(r,7.) - 1835988.*pow(r,8.) - 609036.*pow(r,9.) - 436560.*pow(r,10.) - 436560.*pow(r,11.) - 4389396480.*c2t + 396446400.*r*c2t + 886520160.*pow(r,2.)*c2t - 110752704.*pow(r,3.)*c2t - 504835632.*pow(r,4.)*c2t - 
            338395440.*pow(r,5.)*c2t - 66961128.*pow(r,6.)*c2t - 9180120.*pow(r,7.)*c2t - 14588424.*pow(r,8.)*c2t - 6695880.*pow(r,9.)*c2t + 1309680.*pow(r,10.)*c2t + 1309680.*pow(r,11.)*c2t + 2194698240.*c4t - 826912800.*r*c4t - 2079352080.*pow(r,2.)*c4t + 1881824880.*pow(r,3.)*c4t - 191830320.*pow(r,4.)*c4t + 
            79842960.*pow(r,5.)*c4t - 66934880.*pow(r,6.)*c4t + 23096700.*pow(r,7.)*c4t + 20269020.*pow(r,8.)*c4t + 14028580.*pow(r,9.)*c4t))/(13547520.*pow(r,14.))) ;

            gcov[0][3] = -2.*a*s2t/(r*mu) + zeta*(((10561985280*pow(r,8) + 6706022400*pow(r,9) + 3911846400*pow(r,10))*a*s2t)/(6258954240*pow(r,14)) + (pow(a,3)*(-56330588160*pow(r,4) + 7683984000*pow(r,5) + 18371707200*pow(r,6) + 16977590784*pow(r,7) + 3626485632*pow(r,8) - 373721040*pow(r,9) + 22536360*pow(r,10) + 
            586710432*pow(r,11) + 168991764480*pow(r,4)*c2t - 139233790080*pow(r,5)*c2t - 128246690880*pow(r,6)*c2t - 75932619840*pow(r,7)*c2t - 10000045440*pow(r,8)*c2t - 15800400*pow(r,9)*c2t - 2068605000*pow(r,10)*c2t - 2933552160*pow(r,11)*c2t)*s2t)/(6258954240*pow(r,14)) + 
            (pow(a,5)*(421071146496 - 139972151088*r - 52314708600*pow(r,2) - 41102356680*pow(r,3) + 21781549944*pow(r,4) + 12643617756*pow(r,5) - 1006758858*pow(r,6) - 3329787258*pow(r,7) - 1451349732*pow(r,8) - 764978904*pow(r,9) - 477908475*pow(r,10) - 253500126*pow(r,11) - 2027901173760*c2t + 
            521141765760*r*c2t + 560625468480*pow(r,2)*c2t - 16204624128*pow(r,3)*c2t - 235727442720*pow(r,4)*c2t - 134472345504*pow(r,5)*c2t - 3369309552*pow(r,6)*c2t + 11440144800*pow(r,7)*c2t - 1729175952*pow(r,8)*c2t - 1325370726*pow(r,9)*c2t + 1656071175*pow(r,10)*c2t + 1267500630*pow(r,11)*c2t + 
            1013950586880*c4t - 720017242560*r*c4t - 1006093962720*pow(r,2)*c4t + 1068136866720*pow(r,3)*c4t + 695430723680*pow(r,4)*c4t + 291985421280*pow(r,5)*c4t + 8459512320*pow(r,6)*c4t - 2700008360*pow(r,7)*c4t + 10394560932*pow(r,8)*c4t + 8193773574*pow(r,9)*c4t)*s2t)/(6258954240*pow(r,14))) ;

            gcov[1][1] = mu/DD + zeta*((ep2*(-3657830400*pow(r,4) + 4212190080*pow(r,5) - 1144039680*pow(r,6) + 243781920*pow(r,7) - 330004800*pow(r,8) + 93121920*pow(r,9) - 546480*pow(r,10) - 1519560*pow(r,11) + 1519560*pow(r,12) + 10973491200*pow(r,4)*c2t - 8933016960*pow(r,5)*c2t + 1390556160*pow(r,6)*c2t - 
            846015840*pow(r,7)*c2t + 753991200*pow(r,8)*c2t - 124839360*pow(r,9)*c2t + 4285440*pow(r,10)*c2t + 7204680*pow(r,11)*c2t - 4558680*pow(r,12)*c2t))/(13547520*pow((2 - r),3)*pow(r,12)) + (ep4*(41924832768 - 55560983184*r + 19661682528*pow(r,2) - 3892185204*pow(r,3) + 4726045380*pow(r,4) - 
            1468894340*pow(r,5) - 249401400*pow(r,6) - 11023875*pow(r,7) + 66718519*pow(r,8) - 15329955*pow(r,9) - 1396596*pow(r,10) + 436560*pow(r,11) - 436560*pow(r,12) - 201912238080*c2t + 222952383360*r*c2t - 54077483520*pow(r,2)*c2t - 4105960416*pow(r,3)*c2t - 6745675680*pow(r,4)*c2t + 3910163568*pow(r,5)*c2t + 968164320*pow(r,6)*c2t - 
            287437416*pow(r,7)*c2t - 290871480*pow(r,8)*c2t + 76131840*pow(r,9)*c2t - 3237000*pow(r,10)*c2t - 1838880*pow(r,11)*c2t + 1309680*pow(r,12)*c2t + 100956119040*c4t - 137881154880*r*c4t + 37014140160*pow(r,2)*c4t + 6359274480*pow(r,3)*c4t + 5547126000*pow(r,4)*c4t - 3376291680*pow(r,5)*c4t - 
            527042000*pow(r,6)*c4t + 55986760*pow(r,7)*c4t + 159391380*pow(r,8)*c4t - 56199680*pow(r,9)*c4t + 9469900*pow(r,10)*c4t))/(13547520*pow((2 - r),3)*pow(r,12)));

            gcov[2][2] = r*r*mu + zeta*((ep2*(-365783040*pow(r,4) - 87272640*pow(r,5) - 72770400*pow(r,6) - 11188800*pow(r,7) - 20420640*pow(r,8) - 10735200*pow(r,9) - 4558680*pow(r,10) + 1097349120*pow(r,4)*c2t + 261817920*pow(r,5)*c2t + 218311200*pow(r,6)*c2t + 33566400*pow(r,7)*c2t + 61261920*pow(r,8)*c2t + 32205600*pow(r,9)*c2t + 
            13676040*pow(r,10)*c2t))/(40642560*pow(r,11)) + (ep4*(2734228224 + 188439048*r + 150726420*pow(r,2) - 153381780*pow(r,3) + 133343244*pow(r,4) + 130102470*pow(r,5) + 62670327*pow(r,6) + 20549463*pow(r,7) + 6430788*pow(r,8) + 3141900*pow(r,9) + 1309680*pow(r,10) - 13168189440*c2t - 4297406400*r*c2t + 138814560*pow(r,2)*c2t + 1114609248*pow(r,3)*c2t + 
            482556048*pow(r,4)*c2t - 52440912*pow(r,5)*c2t - 56525400*pow(r,6)*c2t - 6081360*pow(r,7)*c2t + 1400400*pow(r,8)*c2t - 9425700*pow(r,9)*c2t - 3929040*pow(r,10)*c2t + 6584094720*c4t + 1908658080*r*c4t - 2355966000*pow(r,2)*c4t - 996654960*pow(r,3)*c4t - 910785520*pow(r,4)*c4t - 
            310582800*pow(r,5)*c4t - 140389440*pow(r,6)*c4t - 59683820*pow(r,7)*c4t - 28409700*pow(r,8)*c4t))/(40642560*pow(r,11)));

            gcov[3][0] = gcov[0][3]; 

            gcov[3][3] = r*r*sth*sth*(1. + ep2/(r*r) + 2.*ep2*s2t/(r*r*r*mu)) + zeta*((ep2*(-365783040*pow(r,4) - 87272640*pow(r,5) - 72770400*pow(r,6) - 11188800*pow(r,7) - 20420640*pow(r,8) - 10735200*pow(r,9) - 4558680*pow(r,10) + 1097349120*pow(r,4)*c2t + 261817920*pow(r,5)*c2t + 218311200*pow(r,6)*c2t + 33566400*pow(r,7)*c2t + 
            61261920*pow(r,8)*c2t + 32205600*pow(r,9)*c2t + 13676040*pow(r,10)*c2t)*s2t)/(40642560*pow(r,11)) + (ep4*(2734228224 - 543127032*r - 389601900*pow(r,2) - 386195220*pow(r,3) + 38195244*pow(r,4) + 78072390*pow(r,5) + 20779287*pow(r,6) + 696903*pow(r,7) + 1872108*pow(r,8) + 3141900*pow(r,9) + 1309680*pow(r,10) - 
            13168189440*c2t - 1371142080*r*c2t + 2300127840*pow(r,2)*c2t + 2045863008*pow(r,3)*c2t + 863148048*pow(r,4)*c2t + 155679408*pow(r,5)*c2t + 111038760*pow(r,6)*c2t + 73328880*pow(r,7)*c2t + 19635120*pow(r,8)*c2t - 9425700*pow(r,9)*c2t - 3929040*pow(r,10)*c2t + 6584094720*c4t - 286040160*r*c4t - 3976950960*pow(r,2)*c4t - 
            1695095280*pow(r,3)*c4t - 1196229520*pow(r,4)*c4t - 466673040*pow(r,5)*c4t - 266062560*pow(r,6)*c4t - 119241500*pow(r,7)*c4t - 42085740*pow(r,8)*c4t)*s2t)/(40642560*pow(r,11)));

        }

};
// ____________________________________________________________________________________________________________________________
/**
 * Changes Made. EDGB theory KS coordinates. 
 */
class EDGBKSCoords {
    public:
        static constexpr char name[] = "EDGBKSCoords";

        // BH Spin is a property of KS
        const GReal a;
        const GReal zeta;
        static constexpr bool spherical = true;

        KOKKOS_FUNCTION EDGBKSCoords(GReal spin, GReal z): a(spin), zeta(z) {}

        KOKKOS_INLINE_FUNCTION void gcov_embed(const GReal Xembed[GR_DIM], Real gcov[GR_DIM][GR_DIM]) const
        {
            const GReal r = Xembed[1];
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);
            
            // Assign gcov matrix to zero.
            gzero2(gcov); 

            const GReal cth = m::cos(th);
            const GReal sth = m::sin(th);
            const GReal s2 = sth*sth;
            const GReal c2 = cth*cth; 
            const GReal c4 = c2*c2;
            const GReal rho2 = r*r + a*a*cth*cth;
            const GReal a2 = a*a;
            const GReal ep2 = a2;
            const GReal ep4 = ep2*ep2;
            const GReal a4 = a2*a2;
            const GReal DD = 1. - 2./r + a2/(r*r);
            const GReal mu = 1. + a2*cth*cth/(r*r);

            gcov[0][0] = -1. + 2.*r/rho2+ zeta*(-1./15.*(-400.+96.*r+66.*pow(r,2)+130.*pow(r,3)+5.*pow(r,4))/pow(r,7)+ (a2*(pow(r,7)*(444696.-562338.*c2)  + 8820000.*(-1. + 3.*c2) - 19600.*r*(-467. + 1251.*c2) - 63.*pow(r,8)*(-3267. + 8926.*c2) 
            + 1050.*pow(r,3)*(-1465. + 11997.*c2) - 2100.*pow(r,2)*(-955. + 22577.*c2) - 6.*pow(r,6)*(-59329. + 82437.*c2) + 15.*pow(r,5)*(-52533. + 455029.*c2) + 10.*pow(r,4)*(-281221. + 1218513.*c2)))/(110250.*pow(r,11))+ (a4*(675.*pow(r,12)*(-19717. + 67726.*c2) + 450.*pow(r,11)*(43312. + 101589.*c2) 
            - 164640000.*(-70. + 585.*c2 + 156.*c4) - 8232000.*r*(3949. - 11058.*c2 + 6948.*c4) - 39200.*pow(r,2)*(189191. - 824825.*c2 + 972045.*c4) + 60.*pow(r,9)*(717867. - 13885852.*c2 + 12733507.*c4) + 30.*pow(r,10)*(209773. - 9090216.*c2 + 16888370.*c4) + 1400.*pow(r,3)*(648009. - 14691730.*c2 
            + 26074500.*c4) + 420.*pow(r,4)*(553219. - 32471380.*c2 + 222891320.*c4) - 14.*pow(r,7)*(-11393603. + 38599350.*c2 + 359928985.*c4) + 2.*pow(r,8)*(59964679. - 491173980.*c2 + 452824800.*c4) - 28.*pow(r,5)*(21852807. + 12094180.*c2 + 762315255.*c4) 
            - 14.*pow(r,6)*(42004007. - 226439060.*c2 + 1041312310.*c4)))/(30870000.*pow(r,15)));

            gcov[0][1] = 2.*r/rho2+ zeta*(-1./30.*(-800. + 912.*r + 516.*pow(r,2) + 470.*pow(r,3) + 50.*pow(r,4) + 15.*pow(r,5))/ pow(r,7) + (a2*(-566034.*pow(r,8) + 55125.*pow(r,9) +  17640000.*(-1. + 3.*c2) + 12600.*pow(r,3)*(4006. + 1877.*c2) + 78400.*r*(-779. + 2412.*c2) + 42.*pow(r,7)*(-94339. 
            + 107112.*c2) - 1400.*pow(r,2)*(20431. + 132243.*c2) + 36.*pow(r,6)*(-436917. + 491281.*c2) + 20.*pow(r,4)*(-875941. + 2263053.*c2) + 20.*pow(r,5)*(-1122937. + 2632446.*c2)))/(220500.*pow(r,11))- (a4*(-80247150.*pow(r,12) - 5788125.*pow(r,13) + 450.*pow(r,11)*(-1196839. + 812712.*c2) 
            + 329280000.*(-70. + 585.*c2 + 156.*c4) + 49392000.*r*(-1717. + 21664.*c2 + 9076.*c4) + 60.*pow(r,10)*(-23601289. + 13406112.*c2 + 2187250.*c4) + 78400.*pow(r,2)*(-2206069. - 11318105.*c2 + 13657725.*c4) + 14000.*pow(r,3)*(22540153. - 88517480.*c2 + 62290230.*c4) 
            + 30.*pow(r,9)*(-145291221. + 30934768.*c2 + 146683252.*c4) + 280.*pow(r,4)*(793805393. - 2014699860.*c2 + 507428040.*c4) + 28.*pow(r,6)*(-711534337. - 1340707620.*c2 + 4191169150.*c4) + 28.*pow(r,5)*(2484408549. - 2498856340.*c2 + 4680981810.*c4) 
            + 14.*pow(r,7)*(-1576790517. - 1991632680.*c2 + 4890821060.*c4) + 4.*pow(r,8)*(-2741883352. - 2122491690.*c2 + 5145464715.*c4)))/(61740000.*pow(r,15)));

            gcov[0][3] = -2.*a*r*s2/rho2+ zeta*(-1./15.*((-400. + 144.*r + 90.*pow(r,2) + 140.*pow(r,3) + 9.*pow(r,4))*a*(-1. + c2))/pow(r,7) - (pow(a,3)*(-1. + c2)* (pow(r,4)*(2736210. - 4410530.*c2) + pow(r,5)*(766015. - 3620183.*c2) - 
            8820000.*(-1. + 3.*c2) + 19600.*r*(-467. + 1551.*c2) - 12.*pow(r,6)*(26511. +  6310.*c2) + 750.*pow(r,3)*(2051. + 8733.*c2) +  2100.*pow(r,2)*(-955. + 21233.*c2) + 3.*pow(r,8)*(-63529. + 262520.*c2) + pow(r,7)*(-406611. + 563055.*c2)))/(110250.*pow(r,11)));

            gcov[1][0] = gcov[0][1] ;

            gcov[1][1] = 1. + 2.*r/rho2 + zeta*(-1./15.*(-400. + 816.*r + 450.*pow(r,2) + 340.*pow(r,3) + 45.*pow(r,4) + 15.*pow(r,5))/pow(r,7) + (a2*(55125.*pow(r,9) + 8820000.*(-1. + 3.*c2) + 1050.*pow(r,3)*(49537. + 10527.*c2) + 19600.*r*(-3583. + 10899.*c2) + 21.*pow(r,8)*(-36755. + 26778.*c2) 
            + 42.*pow(r,7)*(-104927. + 120501.*c2) - 700.*pow(r,2)*(43727. + 196755.*c2) + 6.*pow(r,6)*(-2680831. + 3030123.*c2) + 10.*pow(r,4)*(-1470661. + 3307593.*c2) + 5.*pow(r,5)*(-4334149. + 9164697.*c2)))/(110250.*pow(r,11))- (a4*(-5788125.*pow(r,13) + 225.*pow(r,12)*(-415805. + 203178.*c2) 
            + 1350.*pow(r,11)*(-384509. + 304767.*c2) + 164640000.*(-70. + 585.*c2 + 156.*c4) + 8232000.*r*(-14251. + 141042.*c2 + 47508.*c4) + 30.*pow(r,10)*(-46992805. + 17722008.*c2 + 21262870.*c4) + 39200.*pow(r,2)*(-4601329. - 21811385.*c2 + 26343405.*c4) 
            + 30.*pow(r,9)*(-143855487. + 3163064.*c2 + 172150266.*c4) + 1400.*pow(r,3)*(226049539. - 899866530.*c2 + 648976800.*c4) + 140.*pow(r,4)*(1589270443. - 4126813860.*c2 + 1683530040.*c4) + 10.*pow(r,8)*(-1084760405. - 947231472.*c2 + 2148750846.*c4) 
            + 28.*pow(r,5)*(2462555742. - 2510950520.*c2 + 3918666555.*c4) + 14.*pow(r,7)*(-1565396914. - 2030232030.*c2 + 4530892075.*c4) + 14.*pow(r,6)*(-1465072681. - 2454976180.*c2 + 7341025990.*c4)))/(30870000.*pow(r,15)));  

            gcov[1][3] = -a *s2*(1. + 2.*r/rho2)+ zeta*(-1./36750.*((16660000. - 5350800.*r + 4797450.*pow(r,2) + 3526600.*pow(r,3) + 2965560.*pow(r,4) + 918855.*pow(r,5) + 187446.*pow(r,6))*a*(-1. + c2))/pow(r,7)+ (pow(a,3)*(-1. + c2)*(22085775.*pow(r,9) 
            + 4571505.*pow(r,10) + 49392000.*(580. + 327.*c2) + 548800.*r*(-23041. + 74715.*c2) + 6300.*pow(r,3)*(-446807. + 973501.*c2) + 126.*pow(r,7)*(1064483. + 1485790.*c2) + 9800.*pow(r,2)*(-1223527. + 1991748.*c2) + 12.*pow(r,8)*(4458631. + 3456783.*c2) 
            + 280.*pow(r,4)*(3773463. + 15733496.*c2) + 42.*pow(r,6)*(6767669. + 23527525.*c2) + 56.*pow(r,5)*(17855552. + 49207893.*c2)))/(3087000.*pow(r,11)));

            gcov[2][2] = rho2+ ((8820000. - 6213200.*r - 3416700.*pow(r,2) - 1855650.*pow(r,3) + 887110.*pow(r,4) + 800733.*pow(r,5) + 435540.*pow(r,6) + 187446.*pow(r,7))*a2*zeta*(1. - 3.*c2))/(110250.*pow(r,8))+ (a4*zeta*(45715050.*pow(r,11)*(-1. + 3.*c2) + 5625.*pow(r,10)*(-20749. + 58131.*c2) 
            + 493920000.*(-70. + 585.*c2 + 156.*c4) + 24696000.*r*(3049. - 10698.*c2 + 8868.*c4) + 117600.*pow(r,2)*(280331. - 1711445.*c2 + 1596165.*c4) + 180.*pow(r,9)*(-1286466. - 846865.*c2 + 5819941.*c4) + 4200.*pow(r,3)*(2362411. - 16650910.*c2 + 14489100.*c4) 
            - 1260.*pow(r,4)*(-3173281. - 5026080.*c2 + 26477920.*c4) + 42.*pow(r,8)*(-18071967. - 940590.*c2 + 54146980.*c4) + 42.*pow(r,6)*(-19116713. - 46592740.*c2 + 138130070.*c4) - 28.*pow(r,5)*(11804979. - 261030540.*c2 + 235282135.*c4) 
            + 6.*pow(r,7)*(-259078241. - 99440670.*c2 + 857000595.*c4)))/(92610000.*pow(r,12));

            gcov[3][0] = gcov[0][3];

            gcov[3][1] = gcov[1][3];

            gcov[3][3] = s2 * (rho2+a2*s2*(1. + 2. *r/rho2))+ zeta*(((8820000. - 6213200.*r - 3416700.*pow(r,2) - 1855650.*pow(r,3) + 887110.*pow(r,4) + 800733.*pow(r,5) + 435540.*pow(r,6) 
            + 187446.*pow(r,7))*a2* (-1. + c2)*(-1. + 3.*c2))/(110250.*pow(r,8)) -  (a4*(-1. + c2)*(45715050.*pow(r,11)*(-1. + 3.*c2) + 5625.*pow(r,10)*(-20749. + 58131.*c2) + 493920000.*(-70. + 585.*c2 + 156.*c4) + 24696000.*r*(3649. - 8958.*c2 + 6528.*c4) 
            + 352800.*pow(r,2)*(84857. - 350495.*c2 + 320655.*c4) + 12600.*pow(r,3)*(-82303. - 1443030.*c2 + 1592200.*c4) + 180.*pow(r,9)*(-411718. - 4345857.*c2 + 8444185.*c4) - 1260.*pow(r,4)*(1578719. - 11450880.*c2 + 28150720.*c4) 
            + 42.*pow(r,8)*(-1863327. - 67980150.*c2 + 104977900.*c4) + 28.*pow(r,5)*(-14247879. - 109560360.*c2 + 137751665.*c4) + 42.*pow(r,6)*(30654807. - 316973820.*c2 + 358739630.*c4) + 6.*pow(r,7)*(-25024421. - 1143700950.*c2 + 1667207055.*c4)))/(92610000.*pow(r,12)));

        }

        KOKKOS_INLINE_FUNCTION void vec_from_bl(const GReal Xembed[GR_DIM], const Real vcon_bl[GR_DIM], Real vcon[GR_DIM]) const
        {
            GReal r = Xembed[1];
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);
            Real trans[GR_DIM][GR_DIM];
            DLOOP2 trans[mu][nu] = (mu == nu);

            const GReal cth = m::cos(th);
            const GReal sth = m::sin(th);
            const GReal s2t = sth*sth;
            const GReal c2t = cth*cth ; 
            const GReal c4t = c2t*c2t ;
            const GReal ep2 = a*a ;
            const GReal ep4 = ep2*ep2 ;

            trans[0][0] = 1 ;  
            trans[0][1] = 2.*r/(r*r - 2.*r + a*a) + zeta*(-1./30.*(-2240. + 144.*r + 96.*pow(r,2.) + 390.*pow(r,3.) + 20.*pow(r,4.) + 15.*pow(r,5.))/
            (pow((-2. + r),2.)*pow(r,5.)) + ((564480000 - 1030019200*r + 964723200*pow(r,2) - 425087600*pow(r,3) + 76461320*pow(r,4) - 50376120*pow(r,5) + 25335172*pow(r,6) + 
            149172*pow(r,7) + 1595034*pow(r,8) - 786534*pow(r,9) + 55125*pow(r,10))*ep2)/(220500*pow((-2 + r),3)*pow(r,9)) + (ep4*(-109543950*pow(r,13) + 45518400*pow(r,14) + 5788125*pow(r,15) + pow(r,7)*(42017164588 + 
            92850567600*c2t - 46425283800*c4t) - 15805440000*(-515 - 312*c2t + 156*c4t) - 197568000*r*(89583 - 12472*c2t + 6236*c4t) + 28224000*pow(r,2)*(467681 - 364878*c2t + 182439*c4t) - 
            9576000*pow(r,3)*(308359 - 435916*c2t + 217958*c4t) - 360*pow(r,12)*(1391168 - 3179270*c2t + 1589635*c4t) + 330*pow(r,11)*(4960011 - 16855960*c2t + 8427980*c4t) - 672*pow(r,6)*(-422234373 - 
            329391250*c2t + 164695625*c4t) - 4480*pow(r,4)*(6618859 - 417919560*c2t + 208959780*c4t) - 16*pow(r,8)*(1064817521 - 1991920800*c2t + 995960400*c4t) + 
            448*pow(r,5)*(-1334378083 - 3683238390*c2t + 1841619195*c4t) - 4*pow(r,10)*(1003478993 - 4282320390*c2t + 2141160195*c4t) + pow(r,9)*(685673062 - 47384968800*c2t + 23692484400*
            c4t)))/(61740000*pow((-2 + r),4)*pow(r,13))) ;
            trans[0][2] = 0 ;
            trans[0][3] = 0 ;

            trans[1][0] = 0 ;
            trans[1][1] = 1 ; 
            trans[1][2] = 0 ;
            trans[1][3] = 0 ;


            trans[2][0] = 0 ;
            trans[2][1] = 0 ; 
            trans[2][2] = 1 ;
            trans[2][3] = 0 ;

   
            trans[3][0] = 0 ;
            trans[3][1] = a/(r*r - 2.*r + a*a) + zeta*(-1./36750.*((70560000. - 96902400.*r + 57429400.*pow(r,2.) - 11130000.*pow(r,3.) + 4106590.*pow(r,4.) - 4567120.*pow(r,5.) + 76674.*pow(r,6.) + 
            169071.*pow(r,7.) + 187446.*pow(r,8.))*a)/(pow((-2. + r),2.)*pow(r,9.)) + (a*ep2*(-5343255*pow(r,12) + 4571505*pow(r,13) + pow(r,8)*(28851876 - 1310998764*c2t) - 
            15410304000*(15 + 8*c2t) - 39513600*r*(-11162 + 3119*c2t) - 16800*pow(r,3)*(87866 + 383337*c2t) + 78400*pow(r,2)*(-2896207 + 2497914*c2t) + 18*pow(r,11)*(-672919 + 3179270*c2t) - 
            168*pow(r,7)*(8855021 + 6121791*c2t) + 2072*pow(r,5)*(7482122 + 7837239*c2t) - 6*pow(r,10)*(-8868567 + 27278270*c2t) + 6*pow(r,9)*(-13597657 + 88187473*c2t) - 
            560*pow(r,4)*(418996 + 89333967*c2t) + 28*pow(r,6)*(-15725399 + 92343093*c2t)))/(3087000*pow((-2 + r),3)*pow(r,13))) ; 
            trans[3][2] = 0 ; 
            trans[3][3] = 1 ;

            gzero(vcon);
            DLOOP2 vcon[mu] += trans[mu][nu]*vcon_bl[nu]; // CHANGES MADE, (DID NOT INCLUDE THIS BEFORE I THINK)

        }

        KOKKOS_INLINE_FUNCTION void vec_to_bl(const GReal Xembed[GR_DIM], const Real vcon_bl[GR_DIM], Real vcon[GR_DIM]) const
        {
            GReal r = Xembed[1];
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);
            GReal rtrans[GR_DIM][GR_DIM], trans[GR_DIM][GR_DIM];
            DLOOP2 rtrans[mu][nu] = (mu == nu);

            const GReal cth = m::cos(th);
            const GReal sth = m::sin(th);
            const GReal s2t = sth*sth;
            const GReal c2t = cth*cth ; 
            const GReal c4t = c2t*c2t ;
            const GReal ep2 = a*a ;
            const GReal ep4 = ep2*ep2 ;

            trans[0][0] = 1 ;  
            trans[0][1] = 2.*r/(r*r - 2.*r + a*a) + zeta*(-1./30.*(-2240. + 144.*r + 96.*pow(r,2.) + 390.*pow(r,3.) + 20.*pow(r,4.) + 15.*pow(r,5.))/
            (pow((-2. + r),2.)*pow(r,5.)) + ((564480000 - 1030019200*r + 964723200*pow(r,2) - 425087600*pow(r,3) + 76461320*pow(r,4) - 50376120*pow(r,5) + 25335172*pow(r,6) + 
            149172*pow(r,7) + 1595034*pow(r,8) - 786534*pow(r,9) + 55125*pow(r,10))*ep2)/(220500*pow((-2 + r),3)*pow(r,9)) + (ep4*(-109543950*pow(r,13) + 45518400*pow(r,14) + 5788125*pow(r,15) + pow(r,7)*(42017164588 + 
            92850567600*c2t - 46425283800*c4t) - 15805440000*(-515 - 312*c2t + 156*c4t) - 197568000*r*(89583 - 12472*c2t + 6236*c4t) + 28224000*pow(r,2)*(467681 - 364878*c2t + 182439*c4t) - 
            9576000*pow(r,3)*(308359 - 435916*c2t + 217958*c4t) - 360*pow(r,12)*(1391168 - 3179270*c2t + 1589635*c4t) + 330*pow(r,11)*(4960011 - 16855960*c2t + 8427980*c4t) - 672*pow(r,6)*(-422234373 - 
            329391250*c2t + 164695625*c4t) - 4480*pow(r,4)*(6618859 - 417919560*c2t + 208959780*c4t) - 16*pow(r,8)*(1064817521 - 1991920800*c2t + 995960400*c4t) + 
            448*pow(r,5)*(-1334378083 - 3683238390*c2t + 1841619195*c4t) - 4*pow(r,10)*(1003478993 - 4282320390*c2t + 2141160195*c4t) + pow(r,9)*(685673062 - 47384968800*c2t + 23692484400*
            c4t)))/(61740000*pow((-2 + r),4)*pow(r,13))) ;
            trans[0][2] = 0 ;
            trans[0][3] = 0 ;

            trans[1][0] = 0 ;
            trans[1][1] = 1 ; 
            trans[1][2] = 0 ;
            trans[1][3] = 0 ;


            trans[2][0] = 0 ;
            trans[2][1] = 0 ; 
            trans[2][2] = 1 ;
            trans[2][3] = 0 ;

   
            trans[3][0] = 0 ;
            trans[3][1] = a/(r*r - 2.*r + a*a) + zeta*(-1./36750.*((70560000. - 96902400.*r + 57429400.*pow(r,2.) - 11130000.*pow(r,3.) + 4106590.*pow(r,4.) - 4567120.*pow(r,5.) + 76674.*pow(r,6.) + 
            169071.*pow(r,7.) + 187446.*pow(r,8.))*a)/(pow((-2. + r),2.)*pow(r,9.)) + (a*ep2*(-5343255*pow(r,12) + 4571505*pow(r,13) + pow(r,8)*(28851876 - 1310998764*c2t) - 
            15410304000*(15 + 8*c2t) - 39513600*r*(-11162 + 3119*c2t) - 16800*pow(r,3)*(87866 + 383337*c2t) + 78400*pow(r,2)*(-2896207 + 2497914*c2t) + 18*pow(r,11)*(-672919 + 3179270*c2t) - 
            168*pow(r,7)*(8855021 + 6121791*c2t) + 2072*pow(r,5)*(7482122 + 7837239*c2t) - 6*pow(r,10)*(-8868567 + 27278270*c2t) + 6*pow(r,9)*(-13597657 + 88187473*c2t) - 
            560*pow(r,4)*(418996 + 89333967*c2t) + 28*pow(r,6)*(-15725399 + 92343093*c2t)))/(3087000*pow((-2 + r),3)*pow(r,13))) ; 
            trans[3][2] = 0 ; 
            trans[3][3] = 1 ;


            invert(&rtrans[0][0], &trans[0][0]); // INVERTING BECAUSE IT IS NOW FROM KS TO BL!!! 

            gzero(vcon);
            DLOOP2 vcon[mu] += trans[mu][nu]*vcon_bl[nu];
        }
};

// ____________________________________________________________________________________________________________________________
/**
 * Changes Made. EDGB theory KS coordinates. 
 */

class EDGBBLCoords {
    public:
        static constexpr char name[] = "EDGBBLCoords";

        // BH Spin is a property of KS
        const GReal a;
        const GReal zeta;
        static constexpr bool spherical = true;

        KOKKOS_FUNCTION EDGBBLCoords(GReal spin, GReal z): a(spin), zeta(z) {}

        KOKKOS_INLINE_FUNCTION void gcov_embed(const GReal Xembed[GR_DIM], Real gcov[GR_DIM][GR_DIM]) const
        {
            const GReal r = Xembed[1];
            const GReal th = excise(excise(Xembed[2], 0.0, SMALL), M_PI, SMALL);
            
            // Assign gcov matrix to zero.
            gzero2(gcov); 

            const GReal cth = m::cos(th);
            const GReal sth = m::sin(th);
            const GReal s2t = sth*sth;
            const GReal c2t = cth*cth; 
            const GReal c4t = c2t*c2t;
            const GReal rho2 = r*r + a*a*cth*cth;
            const GReal a2 = a*a;
            const GReal ep2 = a2;
            const GReal ep4 = ep2*ep2;
            const GReal a4 = a2*a2;
            const GReal r2 = r*r;
            const GReal r4 = r2*r2;
            const GReal DD = 1. - 2./r + ep2/(r*r);
            const GReal mu = 1. + ep2*cth*cth/(r*r);

            gcov[0][0] = -(1. - 2./(r*mu)) + zeta*(-1./30870000.*(-823200000.*pow(r,8.) + 197568000.*pow(r,9.) + 135828000.*pow(r,10.) + 267540000.*pow(r,11.) + 10290000.*pow(r,12.))/pow(r,15.) - (ep2*(2469600000*pow(r,4) - 2562896000*pow(r,5) - 561540000*pow(r,6) + 
            430710000*pow(r,7) + 787418800*pow(r,8) + 220638600*pow(r,9) - 99672720*pow(r,10) - 124514880*pow(r,11) - 57629880*pow(r,12) - 7408800000*pow(r,4)*c2t + 6865488000*pow(r,5)*c2t + 13275276000*pow(r,6)*c2t - 3527118000*pow(r,7)*c2t - 3411836400*pow(r,8)*c2t - 1911121800*pow(r,9)*c2t + 138494160*pow(r,10)*c2t + 
            157454640*pow(r,11)*c2t + 157454640*pow(r,12)*c2t))/(30870000*pow(r,15)) - (ep4*(-11524800000 + 32508168000*r + 7416287200*pow(r,2) - 907212600*pow(r,3) - 232351980*pow(r,4) + 611878596*pow(r,5) + 588056098*pow(r,6) - 159510442*pow(r,7) - 119929358*pow(r,8) - 43072020*pow(r,9) - 
            6293190*pow(r,10) - 19490400*pow(r,11) + 13308975*pow(r,12) + 96314400000*c2t - 91029456000*r*c2t - 32333140000*pow(r,2)*c2t + 20568422000*pow(r,3)*c2t + 13637979600*pow(r,4)*c2t + 338637040*pow(r,5)*c2t - 3170146840*pow(r,6)*c2t + 540390900*pow(r,7)*c2t + 
            982347960*pow(r,8)*c2t + 833151120*pow(r,9)*c2t + 272706480*pow(r,10)*c2t - 45715050*pow(r,11)*c2t - 45715050*pow(r,12)*c2t + 25683840000*c4t + 57195936000*r*c4t + 38104164000*pow(r,2)*c4t - 
            36504300000*pow(r,3)*c4t - 93614354400*pow(r,4)*c4t + 21344827140*pow(r,5)*c4t + 14578372340*pow(r,6)*c4t + 5039005790*pow(r,7)*c4t - 905649600*pow(r,8)*c4t - 764010420*pow(r,9)*c4t - 506651100*pow(r,10)*c4t))/(30870000*pow(r,15))) ;

            gcov[0][3] = -2.*a*s2t/(r*mu) + zeta*(((-27165600000*pow(r,8) + 9779616000*pow(r,9) + 6112260000*pow(r,10) + 9507960000*pow(r,11) + 611226000*pow(r,12))*a*s2t)/(1018710000*pow(r,15)) + (a*ep2*(81496800000*pow(r,4) - 
            84575568000*pow(r,5) - 18530820000*pow(r,6) + 14213430000*pow(r,7) + 25282580400*pow(r,8) + 7077978600*pow(r,9) - 2939539680*pow(r,10) - 3757085640*pow(r,11) - 1761023880*pow(r,12) - 244490400000*pow(r,4)*c2t + 280892304000*pow(r,5)*c2t + 412005132000*pow(r,6)*c2t + 60519690000*pow(r,7)*
            c2t - 40753297200*pow(r,8)*c2t - 33450490920*pow(r,9)*c2t - 699652800*pow(r,10)*c2t + 5202628200*pow(r,11)*c2t + 7277054400*pow(r,12)*c2t)*s2t)/(1018710000*pow(r,15)) + (a*ep2*ep2*(-380318400000 + 1072769544000*r + 244737477600*pow(r,2) - 
            29938015800*pow(r,3) - 7667615340*pow(r,4) + 20191993668*pow(r,5) + 19405851234*pow(r,6) - 5564719986*pow(r,7) - 4295984814*pow(r,8) - 1611965940*pow(r,9) - 289290630*pow(r,10) - 665034700*pow(r,11) + 430455575*pow(r,12) + 3178375200000*c2t - 3492952848000*r*c2t - 
            831196212000*pow(r,2)*c2t + 956920734000*pow(r,3)*c2t + 553351906800*pow(r,4)*c2t - 36480832080*pow(r,5)*c2t - 165779867400*pow(r,6)*c2t - 15014000340*pow(r,7)*c2t + 17039883000*pow(r,8)*c2t + 23469871860*pow(r,9)*c2t + 7463546250*pow(r,10)*c2t - 3041764000*pow(r,11)*c2t - 
            2725302250*pow(r,12)*c2t + 847566720000*c4t + 1626676128000*r*c4t + 487020996000*pow(r,2)*c4t - 1869733404000*pow(r,3)*c4t - 3210861931200*pow(r,4)*c4t - 623042369180*pow(r,5)*c4t + 77644523460*pow(r,6)*c4t + 80476531230*pow(r,7)*c4t - 2936254160*pow(r,8)*c4t - 
            30360437600*pow(r,9)*c4t - 21708712650*pow(r,10)*c4t)*s2t)/(1018710000*pow(r,15)));

            gcov[3][0] = gcov[0][3];

            gcov[1][1] = mu/DD + zeta*((15146880000*pow(r,8) - 15542016000*pow(r,9) + 3934896000*pow(r,10) - 1992144000*pow(r,11) + 1955100000*pow(r,12) - 535080000*pow(r,13) + 92610000*pow(r,14) - 30870000*pow(r,15))/(30870000*pow((2 - r),4)*pow(r,13)) + 
            (ep2*(167932800000*pow(r,4) - 410612160000*pow(r,5) + 337475152000*pow(r,6) - 118020616000*pow(r,7) + 40725512800*pow(r,8) - 29021241200*pow(r,9) + 10117034480*pow(r,10) - 389495960*pow(r,11) - 468335280*pow(r,12) + 17475360*pow(r,13) + 6179880*pow(r,14) + 15435000*pow(r,15) - 
            503798400000*pow(r,4)*c2t + 1182444480000*pow(r,5)*c2t - 1144005744000*pow(r,6)*c2t + 572127528000*pow(r,7)*c2t - 187670330400*pow(r,8)*c2t + 80877375600*pow(r,9)*c2t - 35105083440*pow(r,10)*c2t + 5644637880*pow(r,11)*c2t - 76754160*pow(r,12)*c2t + 
            503233920*pow(r,13)*c2t - 157454640*pow(r,14)*c2t))/(30870000*pow((2 - r),4)*pow(r,13)) + (ep4*(-1152480000000 + 4742521056000*r - 4815584995200*pow(r,2) + 1528575008800*pow(r,3) + 5539906960*pow(r,4) + 21883963080*pow(r,5) - 43044821876*pow(r,6) - 
            11107489456*pow(r,7) + 15782429458*pow(r,8) - 3057061468*pow(r,9) - 404570410*pow(r,10) + 111508170*pow(r,11) - 293151690*pow(r,12) + 107455050*pow(r,13) - 32602725*pow(r,14) + 5788125*pow(r,15) + 9631440000000*c2t - 23014498752000*r*c2t + 18254377680000*pow(r,2)*
            c2t - 4405182600000*pow(r,3)*c2t - 71403735200*pow(r,4)*c2t - 642983462800*pow(r,5)*c2t + 318392457040*pow(r,6)*c2t + 65421739600*pow(r,7)*c2t - 32610296680*pow(r,8)*c2t + 4098219700*pow(r,9)*c2t - 4796219280*pow(r,10)*c2t + 1774459200*pow(r,11)*c2t - 93636960*pow(r,12)*c2t - 
            152580150*pow(r,13)*c2t + 45715050*pow(r,14)*c2t + 2568384000000*c4t + 1358082432000*r*c4t - 5259111984000*pow(r,2)*c4t + 1943658192000*pow(r,3)*c4t + 813282640800*pow(r,4)*c4t - 486702829200*pow(r,5)*c4t - 9611772800*pow(r,6)*c4t + 28705236420*pow(r,7)*c4t + 14053790220*pow(r,8)*c4t - 
            18201492290*pow(r,9)*c4t + 9309896580*pow(r,10)*c4t - 2703498060*pow(r,11)*c4t + 349196460*pow(r,12)*c4t))/(30870000*pow((2 - r),4)*pow(r,13)));

            gcov[2][2] = r2*mu + zeta*((ep2*(7408800000*pow(r,4) - 5219088000*pow(r,5) - 2870028000*pow(r,6) - 1558746000*pow(r,7) + 745172400*pow(r,8) + 672615720*pow(r,9) + 365853600*pow(r,10) + 157454640*pow(r,11) - 22226400000*pow(r,4)*c2t + 15657264000*pow(r,5)*c2t + 8610084000*pow(r,6)*c2t + 
            4676238000*pow(r,7)*c2t - 2235517200*pow(r,8)*c2t - 2017847160*pow(r,9)*c2t - 1097560800*pow(r,10)*c2t - 472363920*pow(r,11)*c2t))/(92610000*pow(r,12)) + (ep4*(-34574400000 + 75298104000*r + 32966925600*pow(r,2) + 9922126200*pow(r,3) + 3998334060*pow(r,4) - 330539412*pow(r,5) - 802901946*pow(r,6) - 
            1554469446*pow(r,7) - 759022614*pow(r,8) - 231563880*pow(r,9) - 116713125*pow(r,10) - 45715050*pow(r,11) + 288943200000*c2t - 264197808000*r*c2t - 201265932000*pow(r,2)*c2t - 69933822000*pow(r,3)*c2t + 6332860800*pow(r,4)*c2t + 7308855120*pow(r,5)*c2t - 1956895080*pow(r,6)*c2t - 
            596644020*pow(r,7)*c2t - 39504780*pow(r,8)*c2t - 152435700*pow(r,9)*c2t + 326986875*pow(r,10)*c2t + 137145150*pow(r,11)*c2t + 77051520000*c4t + 219004128000*r*c4t + 187709004000*pow(r,2)*c4t + 
            60854220000*pow(r,3)*c4t - 33362179200*pow(r,4)*c4t - 6587899780*pow(r,5)*c4t + 5801462940*pow(r,6)*c4t + 5142003570*pow(r,7)*c4t + 2274173160*pow(r,8)*c4t + 1047589380*pow(r,9)*c4t))/(92610000*pow(r,12)));

            gcov[3][3] = r2*sth*sth*(1. + a2/r2 + 2.*a2*s2t/(r2*r*mu)) + zeta*((ep2*(7408800000*pow(r,4) - 5219088000*pow(r,5) - 2870028000*pow(r,6) - 1558746000*pow(r,7) + 745172400*pow(r,8) + 672615720*pow(r,9) + 
            365853600*pow(r,10) + 157454640*pow(r,11) - 22226400000*pow(r,4)*c2t + 15657264000*pow(r,5)*c2t + 8610084000*pow(r,6)*c2t + 4676238000*pow(r,7)*c2t - 2235517200*pow(r,8)*c2t - 2017847160*pow(r,9)*c2t - 1097560800*pow(r,10)*c2t - 472363920*pow(r,11)*c2t)*s2t)/(92610000*pow(r,12)) + 
            (ep4*(-34574400000 + 90115704000*r + 29937549600*pow(r,2) - 1037017800*pow(r,3) - 1989185940*pow(r,4) - 398940612*pow(r,5) + 1287501894*pow(r,6) - 150146526*pow(r,7) - 78259734*pow(r,8) - 74109240*pow(r,9) - 116713125*pow(r,10) - 45715050*pow(r,11) + 288943200000*c2t - 
            221226768000*r*c2t - 123654636000*pow(r,2)*c2t - 18182178000*pow(r,3)*c2t + 14428108800*pow(r,4)*c2t - 3067690080*pow(r,5)*c2t - 13312900440*pow(r,6)*c2t - 6862205700*pow(r,7)*c2t - 2855166300*pow(r,8)*c2t - 782254260*pow(r,9)*c2t + 326986875*pow(r,10)*c2t + 137145150*pow(r,11)*c2t + 77051520000*c4t + 
            161215488000*r*c4t + 113127084000*pow(r,2)*c4t + 20061720000*pow(r,3)*c4t - 35469907200*pow(r,4)*c4t + 3857046620*pow(r,5)*c4t + 15067064460*pow(r,6)*c4t + 10003242330*pow(r,7)*c4t + 4409071800*pow(r,8)*c4t + 1519953300*pow(r,9)*c4t)*s2t)/(92610000*pow(r,12)));

        }
};


// ____________________________________________________________________________________________________________________________
// ____________________________________________________________________________________________________________________________

/**
 * COORDINATE TRANSFORMS:
 * These are transformations which can be applied to base coordinates
 * Each class must define enough functions to apply the transform to coordinates and vectors,
 * both forward and in reverse.
 * That comes out to 4 functions: coord_to_embed, coord_to_native, dXdx, dxdX
 */

/**
 * This class represents a null transformation from the embedding cooridnates, i.e. just using them directly
 */
class NullTransform {
    public:
        static constexpr char name[] = "NullTransform";
        static constexpr GReal startx[3] = {-1, -1, -1};
        static constexpr GReal stopx[3] = {-1, -1, -1};
        // Coordinate transformations
        // Any coordinate value protections (th < 0, th > pi, phi > 2pi) should be in the base system
        KOKKOS_INLINE_FUNCTION void coord_to_embed(const GReal Xnative[GR_DIM], GReal Xembed[GR_DIM]) const
        {
            DLOOP1 Xembed[mu] = Xnative[mu];
        }
        KOKKOS_INLINE_FUNCTION void coord_to_native(const GReal Xembed[GR_DIM], GReal Xnative[GR_DIM]) const
        {
            DLOOP1 Xnative[mu] = Xembed[mu];
        }
        // Tangent space transformation matrices
        KOKKOS_INLINE_FUNCTION void dxdX(const GReal X[GR_DIM], Real dxdX[GR_DIM][GR_DIM]) const
        {
            DLOOP2 dxdX[mu][nu] = (mu == nu);
        }
        KOKKOS_INLINE_FUNCTION void dXdx(const GReal X[GR_DIM], Real dXdx[GR_DIM][GR_DIM]) const
        {
            DLOOP2 dXdx[mu][nu] = (mu == nu);
        }
};

/**
 * Just exponentiate the radial coordinate
 * Makes sense only for spherical base systems!
 */
class ExponentialTransform {
    public:
        static constexpr char name[] = "ExponentialTransform";
        static constexpr GReal startx[3] = {-1, 0., 0.};
        static constexpr GReal stopx[3] = {-1, M_PI, 2*M_PI};

        // Coordinate transformations
        KOKKOS_INLINE_FUNCTION void coord_to_embed(const GReal Xnative[GR_DIM], GReal Xembed[GR_DIM]) const
        {
            Xembed[0] = Xnative[0];
            Xembed[1] = m::exp(Xnative[1]);
#if LEGACY_TH
            Xembed[2] = excise(excise(Xnative[2], 0.0, SMALL), M_PI, SMALL);
#else
            Xembed[2] = Xnative[2];
#endif
            Xembed[3] = Xnative[3];
        }
        KOKKOS_INLINE_FUNCTION void coord_to_native(const GReal Xembed[GR_DIM], GReal Xnative[GR_DIM]) const
        {
            Xnative[0] = Xembed[0];
            Xnative[1] = m::log(Xembed[1]);
            Xnative[2] = Xembed[2];
            Xnative[3] = Xembed[3];
        }
        /**
         * Transformation matrix for contravariant vectors to embedding, or covariant vectors to native
         */
        KOKKOS_INLINE_FUNCTION void dxdX(const GReal Xnative[GR_DIM], Real dxdX[GR_DIM][GR_DIM]) const
        {
            gzero2(dxdX);
            dxdX[0][0] = 1.;
            dxdX[1][1] = m::exp(Xnative[1]);
            dxdX[2][2] = 1.;
            dxdX[3][3] = 1.;
        }
        /**
         * Transformation matrix for contravariant vectors to native, or covariant vectors to embedding
         */
        KOKKOS_INLINE_FUNCTION void dXdx(const GReal Xnative[GR_DIM], Real dXdx[GR_DIM][GR_DIM]) const
        {
            gzero2(dXdx);
            dXdx[0][0] = 1.;
            dXdx[1][1] = 1 / m::exp(Xnative[1]);
            dXdx[2][2] = 1.;
            dXdx[3][3] = 1.;
        }
};

/**
 * SuperExponential coordinates, for super simulations
 * Implementation follows HARMPI described in Tchekhovskoy+
 */
class SuperExponentialTransform {
    public:
        static constexpr char name[] = "SuperExponentialTransform";
        static constexpr GReal startx[3] = {-1, 0., 0.};
        static constexpr GReal stopx[3] = {-1, M_PI, 2*M_PI};

        const GReal xe1br, xn1br;
        const double npow2, cpow2;

        // Constructor
        KOKKOS_FUNCTION SuperExponentialTransform(GReal xe1br_in, double npow2_in, double cpow2_in):
            xe1br(xe1br_in), npow2(npow2_in), cpow2(cpow2_in), xn1br(m::log(xe1br_in)) {}

        // Coordinate transformations
        KOKKOS_INLINE_FUNCTION void coord_to_embed(const GReal Xnative[GR_DIM], GReal Xembed[GR_DIM]) const
        {
            Xembed[0] = Xnative[0];
            const GReal super_dist = Xnative[1] - xn1br;
            Xembed[1] = m::exp(Xnative[1] + (super_dist > 0) * cpow2 * m::pow(super_dist, npow2));
#if LEGACY_TH
            Xembed[2] = excise(excise(Xnative[2], 0.0, SMALL), M_PI, SMALL);
#else
            Xembed[2] = Xnative[2];
#endif
            Xembed[3] = Xnative[3];
        }
        KOKKOS_INLINE_FUNCTION void coord_to_native(const GReal Xembed[GR_DIM], GReal Xnative[GR_DIM]) const
        {
            Xnative[0] = Xembed[0];
            Xnative[2] = Xembed[2];
            Xnative[3] = Xembed[3];
            // TODO can just take log for x1 < xe1br
            ROOT_FIND_1
        }
        /**
         * Transformation matrix for contravariant vectors to embedding, or covariant vectors to native
         */
        KOKKOS_INLINE_FUNCTION void dxdX(const GReal Xnative[GR_DIM], Real dxdX[GR_DIM][GR_DIM]) const
        {
            gzero2(dxdX);
            dxdX[0][0] = 1.;
            const GReal super_dist = Xnative[1] - xn1br;
            dxdX[1][1] = m::exp(Xnative[1] + (super_dist > 0) * cpow2 * m::pow(super_dist, npow2))
                            * (1 + (super_dist > 0) * cpow2 * npow2 * m::pow(super_dist, npow2-1));
            dxdX[2][2] = 1.;
            dxdX[3][3] = 1.;
        }
        /**
         * Transformation matrix for contravariant vectors to native, or covariant vectors to embedding
         */
        KOKKOS_INLINE_FUNCTION void dXdx(const GReal Xnative[GR_DIM], Real dXdx[GR_DIM][GR_DIM]) const
        {
            gzero2(dXdx);
            dXdx[0][0] = 1.;
            const GReal super_dist = Xnative[1] - xn1br;
            dXdx[1][1] = 1 / (m::exp(Xnative[1] + (super_dist > 0) * cpow2 * m::pow(super_dist, npow2))
                              * (1 + (super_dist > 0) * cpow2 * npow2 * m::pow(super_dist, npow2-1)));
            dXdx[2][2] = 1.;
            dXdx[3][3] = 1.;
        }
};

/**
 * Modified Kerr-Schild coordinates "MKS"
 * Makes sense only for spherical base systems!
 */
class ModifyTransform {
    public:
        static constexpr char name[] = "ModifyTransform";
        static constexpr GReal startx[3] = {-1, 0., 0.};
        static constexpr GReal stopx[3] = {-1, 1., 2*M_PI};

        const GReal hslope;

        // Constructor
        KOKKOS_FUNCTION ModifyTransform(GReal hslope_in): hslope(hslope_in) {}

        // Coordinate transformations
        KOKKOS_INLINE_FUNCTION void coord_to_embed(const GReal Xnative[GR_DIM], GReal Xembed[GR_DIM]) const
        {
            Xembed[0] = Xnative[0];
            Xembed[1] = m::exp(Xnative[1]);
#if LEGACY_TH
            const GReal th = M_PI*Xnative[2] + ((1. - hslope)/2.)*m::sin(2.*M_PI*Xnative[2]);
            Xembed[2] = excise(excise(th, 0.0, SMALL), M_PI, SMALL);
#else
            Xembed[2] = M_PI*Xnative[2] + ((1. - hslope)/2.)*m::sin(2.*M_PI*Xnative[2]);
#endif
            Xembed[3] = Xnative[3];
        }
        KOKKOS_INLINE_FUNCTION void coord_to_native(const GReal Xembed[GR_DIM], GReal Xnative[GR_DIM]) const
        {
            Xnative[0] = Xembed[0];
            Xnative[1] = m::log(Xembed[1]);
            Xnative[3] = Xembed[3];
            // Treat the special case with a large macro
            ROOT_FIND
        }
        /**
         * Transformation matrix for contravariant vectors to embedding, or covariant vectors to native
         */
        KOKKOS_INLINE_FUNCTION void dxdX(const GReal Xnative[GR_DIM], Real dxdX[GR_DIM][GR_DIM]) const
        {
            gzero2(dxdX);
            dxdX[0][0] = 1.;
            dxdX[1][1] = m::exp(Xnative[1]);
            dxdX[2][2] = M_PI - (hslope - 1.)*M_PI*m::cos(2.*M_PI*Xnative[2]);
            dxdX[3][3] = 1.;
        }
        /**
         * Transformation matrix for contravariant vectors to native, or covariant vectors to embedding
         */
        KOKKOS_INLINE_FUNCTION void dXdx(const GReal Xnative[GR_DIM], Real dXdx[GR_DIM][GR_DIM]) const
        {
            gzero2(dXdx);
            dXdx[0][0] = 1.;
            dXdx[1][1] = 1 / m::exp(Xnative[1]);
            dXdx[2][2] = 1 / (M_PI - (hslope - 1.)*M_PI*m::cos(2.*M_PI*Xnative[2]));
            dXdx[3][3] = 1.;
        }
};

/**
 * "Funky" Modified Kerr-Schild coordinates
 * Make sense only for spherical base systems!
 */
class FunkyTransform {
    public:
        static constexpr char name[] = "FunkyTransform";
        static constexpr GReal startx[3] = {-1, 0., 0.};
        static constexpr GReal stopx[3] = {-1, 1., 2*M_PI};

        const GReal startx1;
        const GReal hslope, poly_xt, poly_alpha, mks_smooth;
        // Must be *defined* afterward to use constructor below
        const GReal poly_norm;

        // Constructor
        KOKKOS_FUNCTION FunkyTransform(GReal startx1_in, GReal hslope_in, GReal mks_smooth_in, GReal poly_xt_in, GReal poly_alpha_in):
            startx1(startx1_in), hslope(hslope_in), mks_smooth(mks_smooth_in), poly_xt(poly_xt_in), poly_alpha(poly_alpha_in),
            poly_norm(0.5 * M_PI * 1./(1. + 1./(poly_alpha + 1.) * 1./m::pow(poly_xt, poly_alpha))) {}

        // Coordinate transformations
        KOKKOS_INLINE_FUNCTION void coord_to_embed(const GReal Xnative[GR_DIM], GReal Xembed[GR_DIM]) const
        {
            Xembed[0] = Xnative[0];
            Xembed[1] = m::exp(Xnative[1]);

            const GReal thG = M_PI*Xnative[2] + ((1. - hslope)/2.)*m::sin(2.*M_PI*Xnative[2]);
            const GReal y = 2*Xnative[2] - 1.;
            const GReal thJ = poly_norm * y * (1. + m::pow(y/poly_xt,poly_alpha) / (poly_alpha + 1.)) + 0.5 * M_PI;
#if LEGACY_TH
            const GReal th = thG + m::exp(mks_smooth * (startx1 - Xnative[1])) * (thJ - thG);
            Xembed[2] = excise(excise(th, 0.0, SMALL), M_PI, SMALL);
#else
            Xembed[2] = thG + m::exp(mks_smooth * (startx1 - Xnative[1])) * (thJ - thG);
#endif
            Xembed[3] = Xnative[3];
        }
        KOKKOS_INLINE_FUNCTION void coord_to_native(const GReal Xembed[GR_DIM], GReal Xnative[GR_DIM]) const
        {
            Xnative[0] = Xembed[0];
            Xnative[1] = m::log(Xembed[1]);
            Xnative[3] = Xembed[3];
            // Treat the special case with a macro
            ROOT_FIND
        }
        /**
         * Transformation matrix for contravariant vectors to embedding, or covariant vectors to native
         */
        KOKKOS_INLINE_FUNCTION void dxdX(const GReal Xnative[GR_DIM], Real dxdX[GR_DIM][GR_DIM]) const
        {
            gzero2(dxdX);
            dxdX[0][0] = 1.;
            dxdX[1][1] = m::exp(Xnative[1]);
            dxdX[2][1] = -exp(mks_smooth * (startx1 - Xnative[1])) * mks_smooth
                * (
                M_PI / 2. -
                M_PI * Xnative[2]
                    + poly_norm * (2. * Xnative[2] - 1.)
                        * (1
                            + (m::pow((-1. + 2 * Xnative[2]) / poly_xt, poly_alpha))
                                / (1 + poly_alpha))
                    - 1. / 2. * (1. - hslope) * m::sin(2. * M_PI * Xnative[2]));
            dxdX[2][2] = M_PI + (1. - hslope) * M_PI * m::cos(2. * M_PI * Xnative[2])
                + m::exp(mks_smooth * (startx1 - Xnative[1]))
                    * (-M_PI
                        + 2. * poly_norm
                            * (1.
                                + m::pow((2. * Xnative[2] - 1.) / poly_xt, poly_alpha)
                                    / (poly_alpha + 1.))
                        + (2. * poly_alpha * poly_norm * (2. * Xnative[2] - 1.)
                            * m::pow((2. * Xnative[2] - 1.) / poly_xt, poly_alpha - 1.))
                            / ((1. + poly_alpha) * poly_xt)
                        - (1. - hslope) * M_PI * m::cos(2. * M_PI * Xnative[2]));
            dxdX[3][3] = 1.;
        }
        /**
         * Transformation matrix for contravariant vectors to native, or covariant vectors to embedding
         */
        KOKKOS_INLINE_FUNCTION void dXdx(const GReal Xnative[GR_DIM], Real dXdx[GR_DIM][GR_DIM]) const
        {
            // Okay this one should probably stay numerical
            Real dxdX_tmp[GR_DIM][GR_DIM];
            dxdX(Xnative, dxdX_tmp);
            invert(&dxdX_tmp[0][0],&dXdx[0][0]);
        }
};

// Bundle coordinates and transforms into umbrella variant types
// These act as a wannabe "interface" or "parent class" with the exception that access requires "mpark::visit"
// See coordinate_embedding.hpp
using SomeBaseCoords = mpark::variant<SphMinkowskiCoords, CartMinkowskiCoords, SphBLCoords, SphKSCoords, SphBLExtG, SphKSExtG, DCSKSCoords, DCSBLCoords, EDGBKSCoords, EDGBBLCoords>; // Changes Made. 
using SomeTransform = mpark::variant<NullTransform, ExponentialTransform, SuperExponentialTransform, ModifyTransform, FunkyTransform>;