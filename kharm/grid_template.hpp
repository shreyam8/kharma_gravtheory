/*
 * Class representing a logically Cartesian grid of points in a CoordinateSystem, including:
 * * Grid zone locations, start & end points
 * * Ghost or "halo" zones, iterators over grids with/without halos
 */
#pragma once

#include "decs.hpp"
#include "coordinates.hpp"

#include <vector>
#include <memory>

// TODO standardize namespaces
using namespace Kokkos;
using namespace std;

// TODO list
// Carrying geometry on grid loci
// Grid/geom/linalg operations? (considering spatially-symmetric geom?)
// MPI/grid less than full support
// Passive variables

/**
 * Struct holding all parameters related to the logically Cartesian grid.  Purposefully minimal 
 */
template<class Coords>
class Grid
{
protected:
    GeomTensor gcon_direct, gcov_direct;
    GeomScalar gdet_direct;
    GeomConn conn_direct;

public:
    // TODO if we ditch AMR we can probably scrape up time making grids static
    // That or *very* disciplined pass-by-ref
    int n1, n2, n3, nvar;
    int ng, gn1, gn2, gn3;
    int n1tot, n2tot, n3tot;
    int n1start, n2start, n3start;
    GReal startx1, startx2, startx3;
    GReal dx1, dx2, dx3;

    Coords *coords;

    // TODO see if these slow anything.  Leaves me the option to return them analytically for e.g. very fast Minkowski
    KOKKOS_INLINE_FUNCTION Real gcon(const Loci loc, const int i, const int j, const int mu, const int nu) const;
    KOKKOS_INLINE_FUNCTION Real gcov(const Loci loc, const int i, const int j, const int mu, const int nu) const;
    KOKKOS_INLINE_FUNCTION Real gdet(const Loci loc, const int i, const int j) const;
    KOKKOS_INLINE_FUNCTION Real conn(const int i, const int j, const int mu, const int nu, const int lam) const;

    // Constructors
    Grid(CoordinateSystem *coordinates, std::vector<int> shape, std::vector<GReal> startx, std::vector<GReal> endx, int ng_in=3, int nvar_in=8);
    Grid(CoordinateSystem *coordinates, std::vector<int> fullshape, std::vector<int> startn, std::vector<int> shape, std::vector<GReal> startx, std::vector<GReal> endx, int ng_in=3, int nvar_in=8);
    void init_grids();

    // Coordinates of the grid, i.e. "native"
    KOKKOS_INLINE_FUNCTION void coord(const int i, const int j, const int k, const Loci loc, Real X[NDIM]) const;
    // TODO think on this.  More the domain of coords, but on the other hand, real convenient
    KOKKOS_INLINE_FUNCTION void ks_coord(const int i, const int j, const int k, const Loci loc, Real &r, Real &th) const
    {
        Real X[NDIM];
        coord(i, j, k, loc, X);
        coords->ks_coord(X, r, th);
    }

    // Transformations using the cached geometry
    // TODO could dramatically expand with enough loop fission
    KOKKOS_INLINE_FUNCTION void lower(const Real vcon[NDIM], Real vcov[NDIM],
                                        const int i, const int j, const int k, const Loci loc) const;
    KOKKOS_INLINE_FUNCTION void raise(const Real vcov[NDIM], Real vcon[NDIM],
                                        const int i, const int j, const int k, const Loci loc) const;
    KOKKOS_INLINE_FUNCTION void lower_grid(const GridVector vcon, GridVector vcov,
                                        const int i, const int j, const int k, const Loci loc) const;
    KOKKOS_INLINE_FUNCTION void raise_grid(const GridVector vcov, GridVector vcon,
                                        const int i, const int j, const int k, const Loci loc) const;

    // Indexing
    KOKKOS_INLINE_FUNCTION int i4(const int i, const int j, const int k, const int p, const bool use_ghosts=true) const
    {
        if (use_ghosts) {
            return i*gn2*gn3*nvar + j*gn3*nvar + k*nvar + p;
        } else {
            return i*n2*n3*nvar + j*n3*nvar + k*nvar + p;
        }
    }
    KOKKOS_INLINE_FUNCTION int i3(const int i, const int j, const int k, const bool use_ghosts=true) const
    {
        if (use_ghosts) {
            return i*gn2*gn3 + j*gn3 + k;
        } else {
            return i*n2*n3 + j*n3 + k;
        }
    }

    // RangePolicies over grid.  Construct them on the fly to stay slim
    // TODO can I consolidate these?  Types are fun
    MDRangePolicy<Rank<3>> bulk_0() const {return MDRangePolicy<Rank<3>>({0, 0, 0}, {n1, n2, n3});}
    MDRangePolicy<Rank<3>> bulk_ng() const {return MDRangePolicy<Rank<3>>({ng, ng, ng}, {n1+ng, n2+ng, n3+ng});}
    MDRangePolicy<Rank<3>> all_0() const {return MDRangePolicy<Rank<3>>({0, 0, 0}, {n1+2*ng, n2+2*ng, n3+2*ng});}
    MDRangePolicy<OpenMP,Rank<3>> h_bulk_0() const {return MDRangePolicy<OpenMP,Rank<3>>({0, 0, 0}, {n1, n2, n3});}
    MDRangePolicy<OpenMP,Rank<3>> h_bulk_ng() const {return MDRangePolicy<OpenMP,Rank<3>>({ng, ng, ng}, {n1+ng, n2+ng, n3+ng});}
    MDRangePolicy<OpenMP,Rank<3>> h_all_0() const {return MDRangePolicy<OpenMP,Rank<3>>({0, 0, 0}, {n1+2*ng, n2+2*ng, n3+2*ng});}

    MDRangePolicy<Rank<3>> bulk_plus(const int i) const {return MDRangePolicy<Rank<3>>({ng-i, ng-i, ng-i}, {n1+ng+i, n2+ng+i, n3+ng+i});}

    // Versions with 4th index for fluid variables
    MDRangePolicy<Rank<4>> bulk_0_p() const {return MDRangePolicy<Rank<4>>({0, 0, 0, 0}, {n1, n2, n3, nvar});}
    MDRangePolicy<Rank<4>> bulk_ng_p() const {return MDRangePolicy<Rank<4>>({ng, ng, ng, 0}, {n1+ng, n2+ng, n3+ng, nvar});}
    MDRangePolicy<Rank<4>> all_0_p() const {return MDRangePolicy<Rank<4>>({0, 0, 0, 0}, {n1+2*ng, n2+2*ng, n3+2*ng, nvar});}
    MDRangePolicy<OpenMP,Rank<4>> h_bulk_0_p() const {return MDRangePolicy<OpenMP,Rank<4>>({0, 0, 0, 0}, {n1, n2, n3, nvar});}
    MDRangePolicy<OpenMP,Rank<4>> h_bulk_ng_p() const {return MDRangePolicy<OpenMP,Rank<4>>({ng, ng, ng, 0}, {n1+ng, n2+ng, n3+ng, nvar});}
    MDRangePolicy<OpenMP,Rank<4>> h_all_0_p() const {return MDRangePolicy<OpenMP,Rank<4>>({0, 0, 0, 0}, {n1+2*ng, n2+2*ng, n3+2*ng, nvar});}

    MDRangePolicy<Rank<4>> bulk_plus_p(const int i) const {return MDRangePolicy<Rank<4>>({ng-i, ng-i, ng-i, 0}, {n1+ng+i, n2+ng+i, n3+ng+i, nvar});}

    // TODO Versions for geometry stuff?

};

/**
 * Construct a grid which starts at 0 and covers the entire global space
 */
template<class Coords>
Grid<Coords>::Grid(CoordinateSystem *coordinates, std::vector<int> shape, std::vector<GReal> startx,
            std::vector<GReal> endx, const int ng_in, const int nvar_in)
{
    nvar = nvar_in;
    ng = ng_in;

    n1tot = n1 = shape[0];
    n2tot = n2 = shape[1];
    n3tot = n3 = shape[2];

    n1start = 0;
    n2start = 0;
    n3start = 0;

    gn1 = n1 + 2*ng;
    gn2 = n2 + 2*ng;
    gn3 = n3 + 2*ng;

    startx1 = startx[0];
    startx2 = startx[1];
    startx3 = startx[2];

    dx1 = (endx[0] - startx1) / n1;
    dx2 = (endx[1] - startx2) / n2;
    dx3 = (endx[2] - startx3) / n3;

    // TODO do I really want to manage this from here?
    coords = new Coords();

    // Allocate and initialize all the device-side caches of geometry data
    // Since it's axisymmetric and therefore reads are amortized, this is almost
    // certainly waaaay faster than computing all geometry on the fly
    gcon_direct = GeomTensor("gcon", NLOC, n1, n2);
    gcov_direct = GeomTensor("gcov", NLOC, n1, n2);
    gdet_direct = GeomScalar("gdet", NLOC, n1, n2);
    conn_direct = GeomConn("conn", n1, n2);
}

/**
 * Construct a sub-grid starting at some point in a global space
 */
template<class Coords>
Grid<Coords>::Grid(CoordinateSystem *coordinates, std::vector<int> fullshape, std::vector<int> startn,
            std::vector<int> shape, std::vector<GReal> startx, std::vector<GReal> endx,
            const int ng_in, const int nvar_in)
{
    nvar = nvar_in;
    ng = ng_in;

    n1tot = fullshape[0];
    n2tot = fullshape[1];
    n3tot = fullshape[2];

    n1start = startn[0];
    n2start = startn[1];
    n3start = startn[2];

    n1 = shape[0];
    n2 = shape[1];
    n3 = shape[2];

    gn1 = n1 + 2*ng;
    gn2 = n2 + 2*ng;
    gn3 = n3 + 2*ng;

    startx1 = startx[0];
    startx2 = startx[1];
    startx3 = startx[2];

    dx1 = (endx[0] - startx1) / n1;
    dx2 = (endx[1] - startx2) / n2;
    dx3 = (endx[2] - startx3) / n3;

    coords = new Coords();

    gcon_direct = GeomTensor("gcon", NLOC, n1, n2);
    gcov_direct = GeomTensor("gcov", NLOC, n1, n2);
    gdet_direct = GeomScalar("gdet", NLOC, n1, n2);
    conn_direct = GeomConn("conn", n1, n2);
}

// For most coordinate systems, cache the geometry and return the cached values when called
template<class Coords>
void Grid<Coords>::init_grids() {
    // Cache gcon, gcov, gdet, and the connection coeffs
    Kokkos::parallel_for("init_geom", MDRangePolicy<Rank<2>>({0, 0}, {n1+2*ng, n2+2*ng}),
        KOKKOS_LAMBDA (const int i, const int j) {
            GReal X[NDIM];
            Real gcov_loc[NDIM][NDIM], gcon_loc[NDIM][NDIM];
            for (int loc=0; loc < NLOC; ++loc) {
                coord(i, j, 0, (Loci)loc, X);
                coords->gcov_native(X, gcov_loc);
                gdet_direct((Loci)loc, i, j) = coords->gcon_native(gcov_loc, gcon_loc);
                DLOOP2 {
                    gcov_direct((Loci)loc, i, j, mu, nu) = gcov_loc[mu][nu];
                    gcon_direct((Loci)loc, i, j, mu, nu) = gcon_loc[mu][nu];
                }
            }
        }
    );
    Kokkos::parallel_for("init_conn", MDRangePolicy<Rank<2>>({0, 0}, {n1+2*ng, n2+2*ng}),
        KOKKOS_LAMBDA (const int i, const int j) {
            GReal X[NDIM];
            coord(i, j, 0, Loci::center, X);
            Real conn_loc[NDIM][NDIM][NDIM];
            coords->conn_func(X, conn_loc);
            DLOOP2 for(int kap=0; kap<NDIM; ++kap)
                conn_direct(i, j, mu, nu, kap) = conn_loc[mu][nu][kap];
        }
    );
}
template <class Coords>
KOKKOS_INLINE_FUNCTION Real Grid<Coords>::gcon(const Loci loc, const int i, const int j, const int mu, const int nu) const
{
    return gcon_direct(loc, i, j, mu, nu);
}
template <class Coords>
KOKKOS_INLINE_FUNCTION Real Grid<Coords>::gcov(const Loci loc, const int i, const int j, const int mu, const int nu) const
{
    return gcov_direct(loc, i, j, mu, nu);
}
template <class Coords>
KOKKOS_INLINE_FUNCTION Real Grid<Coords>::gdet(const Loci loc, const int i, const int j) const
{
    return gdet_direct(loc, i, j);
}
template <class Coords>
KOKKOS_INLINE_FUNCTION Real Grid<Coords>::conn(const int i, const int j, const int mu, const int nu, const int lam) const
{
    return conn_direct(i, j, mu, nu, lam);
}

// But in Minkowski, screw memory, we have the answers already
template <>
void Grid<Minkowski>::init_grids() {};
template <>
KOKKOS_INLINE_FUNCTION Real Grid<Minkowski>::gcon(const Loci loc, const int i, const int j, const int mu, const int nu) const
{
    if (mu == 0 && nu == 0 ) {
        return -1;
    } else if (mu == nu) {
        return 1;
    } else {
        return 0;
    }
}
template <>
KOKKOS_INLINE_FUNCTION Real Grid<Minkowski>::gcov(const Loci loc, const int i, const int j, const int mu, const int nu) const
{
    if (mu == 0 && nu == 0 ) {
        return -1;
    } else if (mu == nu) {
        return 1;
    } else {
        return 0;
    }
}
template <>
KOKKOS_INLINE_FUNCTION Real Grid<Minkowski>::gdet(const Loci loc, const int i, const int j) const
{
    return 1;
}
template <>
KOKKOS_INLINE_FUNCTION Real Grid<Minkowski>::conn(const int i, const int j, const int mu, const int nu, const int lam) const
{
    return 0;
}

/**
 * Function to return native coordinates of a grid
 * TODO is it more instruction-efficient to split this per location or have a separate one for centers?
 */
template <class Coords>
KOKKOS_INLINE_FUNCTION void Grid<Coords>::coord(const int i, const int j, const int k, Loci loc, Real *X) const
{
    X[0] = 0;
    switch (loc)
    {
    case face1:
        X[1] = startx1 + (n1start + i - ng) * dx1;
        X[2] = startx2 + (n2start + j + 0.5 - ng) * dx2;
        X[3] = startx3 + (n3start + k + 0.5 - ng) * dx3;
        break;
    case face2:
        X[1] = startx1 + (n1start + i + 0.5 - ng) * dx1;
        X[2] = startx2 + (n2start + j - ng) * dx2;
        X[3] = startx3 + (k + 0.5 - ng) * dx3;
        break;
    case face3:
        X[1] = startx1 + (n1start + i + 0.5 - ng) * dx1;
        X[2] = startx2 + (n2start + j + 0.5 - ng) * dx2;
        X[3] = startx3 + (n3start + k - ng) * dx3;
        break;
    case center:
        X[1] = startx1 + (n1start + i + 0.5 - ng) * dx1;
        X[2] = startx2 + (n2start + j + 0.5 - ng) * dx2;
        X[3] = startx3 + (n3start + k + 0.5 - ng) * dx3;
        break;
    case corner:
        X[1] = startx1 + (n1start + i - ng) * dx1;
        X[2] = startx2 + (n2start + j - ng) * dx2;
        X[3] = startx3 + (n3start + k - ng) * dx3;
        break;
    }
}

template <class Coords>
KOKKOS_INLINE_FUNCTION void Grid<Coords>::lower(const Real vcon[NDIM], Real vcov[NDIM],
                                        const int i, const int j, const int k, const Loci loc) const
{
    DLOOP2 vcov[mu] += gcov(loc, i, j, mu, nu) * vcon[nu];
}

template <class Coords>
KOKKOS_INLINE_FUNCTION void Grid<Coords>::raise(const Real vcov[NDIM], Real vcon[NDIM],
                                        const int i, const int j, const int k, const Loci loc) const
{
    DLOOP2 vcon[mu] += gcon(loc, i, j, mu, nu) * vcov[nu];
}

template <class Coords>
KOKKOS_INLINE_FUNCTION void Grid<Coords>::lower_grid(const GridVector vcon, GridVector vcov,
                                        const int i, const int j, const int k, const Loci loc) const
{
    DLOOP2 vcov(i, j, k, mu) += gcov(loc, i, j, mu, nu) * vcon(i, j, k, nu);
}

template <class Coords>
KOKKOS_INLINE_FUNCTION void Grid<Coords>::raise_grid(const GridVector vcov, GridVector vcon,
                                        const int i, const int j, const int k, const Loci loc) const
{
    DLOOP2 vcon(i, j, k, mu) += gcon(loc, i, j, mu, nu) * vcov(i, j, k, nu);
}