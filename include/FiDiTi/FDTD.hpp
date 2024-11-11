#ifndef FIDITI_FDTD_HPP
#define FIDITI_FDTD_HPP

#include "Math/Rect.hpp"

#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>
#include <filesystem>


namespace fidi
{

namespace constants
{
    inline constexpr double pi       = std::asin(1) * 2;
    inline constexpr double epsilon0 = 8.854187817e-12;
    inline constexpr double mu0      = 1.25663706212e-6;
    inline constexpr double c        = 2.99792458e8;  // Speed of light (vacuum) = 1 / sqrt(epsilon0 * mu0)

    //inline constexpr double imp0     = std::sqrt(mu0 / epsilon0);
}

namespace Std
{
    // C++20
    template <class T, class Alloc, class Pred>
    constexpr std::size_t
    erase_if(std::vector<T, Alloc>& c, Pred pred)
    {
        auto it = std::remove_if(c.begin(), c.end(), pred);
        auto n = std::distance(it, c.end());
        c.erase(it, c.end());
        return n;
    }
}

namespace fn
{
    template <class Scalar>
    Scalar gauss(Scalar x, Scalar d = 30., Scalar w = 10.)  // Gauss curve
    {
        Scalar a = (x - d) / w;
        return std::exp(-a * a);
    }

    // Ricker wavelet (2nd deriv of Gauss)
    template <class Scalar>
    Scalar ricker(Scalar t, Scalar d = 30., Scalar P = 10.)  // P: Peak frequency
    {
        Scalar X = constants::pi * P * (t - d);
        return (1 - 2*X*X) * std::exp(-X*X);
    }

    template <class Scalar>
    auto makePulse(Scalar delay, Scalar width)
    {
        return [delay, width](Scalar t) { return gauss(t, delay, width); };
    }
}


namespace fdtd
{

enum class Boundary
{
    XMin = 0, XMax,
    YMin, YMax,
    ZMin, ZMax
};

int toAxis(Boundary b)
{
    return int(b) / 2;
}

struct Material
{
    Scalar rA = 1;  // epsilon_r
    Scalar rB = 1;  // mu_r
    Scalar sA = 0;  // sigma_e (el. conductivity)
    Scalar sB = 0;  // sigma_m (magn. conductivity)
};

using MaterialMap = std::unordered_map<int, Material>;

struct ReferenceMaterial  // Vacuum (Air)
{
    Scalar mA = constants::epsilon0;
    Scalar mB = constants::mu0;
};

struct UpdateCoeffs  // Pre-compute at start of simulation, internal only
{
    Scalar cAA = 1;  // A[i] = A[i]*cAA + (B[i] - B[i-1])*cAB
    Scalar cBB = 1;  // B[i] = B[i]*cBB + (A[i+1] - A[i])*cBA
    Scalar cAB = 0;
    Scalar cBA = 0;

    UpdateCoeffs(Material m, ReferenceMaterial m0, Scalar h, Scalar k)
    {
        Scalar lossA = m.sA * k / (2 * m.rA * m0.mA);
        Scalar lossB = m.sB * k / (2 * m.rB * m0.mB);
        cAA = (1 - lossA) / (1 + lossA);
        cBB = (1 - lossB) / (1 + lossB);
        cAB = k / (h * m.rA * m0.mA) / (1 + lossA);
        cBA = k / (h * m.rB * m0.mB) / (1 + lossB);
    }

    UpdateCoeffs() = default;
};


// Point source (Hard source, or additive/current source)
struct Source
{
    int cell = 0;
    Scalar duration = 0;  // Replace source with a regular node after this time (0 = no expiry)
    // Field component ... (using Az for now)

    std::function<Scalar(Scalar)> f = fn::makePulse<Scalar>(30, 10);
};

// Material must not change in normal-direction at the boundary (order 1: 2 cells; order 2: 3 cells)
// Material must be lossless (no conductivity)
template <int D>
struct ABC
{
    VecNi<D> direction;  // Normal direction, into the domain
    int order = 1;

    std::vector<int> cells;
    std::vector<Scalar> oldValues;

    std::vector<Scalar>* field;
};

template <int D>
struct TFSF;  // Contains FDTD<1>, requires class definition -> See below


template <int D>
class FDTD
{
    VecNi<D> N;  // Grid: N.x * N.y * N.z cells

    std::vector<Scalar> Ax, Ay, Az, Bx, By, Bz;
    std::vector<Scalar> cAA, cAB, cBB, cBA;  // Interpolate: cAA[3], cAB[3], cBB[3], cBA[3]

    std::vector<int> matIds;
    MaterialMap mats;
    ReferenceMaterial refMat;

    Scalar S_c = 1 / std::sqrt(D);  // Stability constant (Courant number) = c0*dt/dx
    Scalar dx = 1;
    Scalar dt = 0;  // S_c * dx / c0

    std::vector<Source> hardSrcsA, hardSrcsB;
    std::vector<Source> additiveSrcsA, additiveSrcsB;
    std::vector<ABC<D>> abcsA, abcsB;
    std::vector<TFSF<D>> tfsfSrcs;

    int timeStep = -1;  // First step at 0; -1 means initial sources not applied yet

    int frame = 0;          // Output frame counter
    int frameInterval = 4;  // Output interval; 0 means no automatic printing

 public:
    FDTD(VecNi<D> gridSize)
      : FDTD(gridSize, 1 / std::sqrt(D))
    { }

    FDTD(VecNi<D> gridSize, Scalar stabilityConstant)
      : N(gridSize),
        S_c(stabilityConstant),
        dt(S_c * dx / constants::c)
    {
        if (S_c > 1 / std::sqrt(D))
            std::cout << "Warning: unstable Courant number" << std::endl;

        int numCells = RectNi<D>{{}, gridSize}.volume();
        matIds.resize(numCells);
        Az.resize(numCells);
        By.resize(numCells);
        if (D > 1)
            Bx.resize(numCells);
        if (D > 2)
        {
            Ax.resize(numCells);
            Ay.resize(numCells);
            Bz.resize(numCells);
        }
        cAA.resize(numCells);
        cAB.resize(numCells);
        cBB.resize(numCells);
        cBA.resize(numCells);

        setMaterials({{0, Material()}});
    }

    int step()
    {
        ensureInitialState();

        updateFieldB();
        applySourcesB(timeStep + 0.5);

        updateFieldA();
        applySourcesA(timeStep + 1);

        return ++timeStep;
    }

    void run(int numSteps)
    {
        print();  // Initial state

        for (int q = 0; q < numSteps; ++q)
        {
            step();
            print();
        }
    }

    // Call this only after adding all sources, to set the initial fields
    // (only needed if reading the field values at time = 0)
    void ensureInitialState()
    {
        if (timeStep == -1)
        {
            applySourcesA(0);
            timeStep = 0;
        }
    }

    // Call ensureInitialState() to get correct values at time = 0
    const std::vector<Scalar>& fieldA(unsigned i) const
    {
        assert(i < 3);
        switch(i)
        {
            case 0:  return Ax;
            case 1:  return Ay;
            default: return Az;
        }
    }

    // Call ensureInitialState() to get correct values at time = 0
    const std::vector<Scalar>& fieldB(unsigned i) const
    {
        assert(i < 3);
        switch(i)
        {
            case 0:  return Bx;
            case 1:  return By;
            default: return Bz;
        }
    }

    void setMaterials(MaterialMap materials)  // defineMaterials(), setMaterialParams()
    {
        mats = std::move(materials);
        applyMaterials();
    }

    void setCellMaterials(std::vector<int> materialIds)
    {
        matIds = std::move(materialIds);
        applyMaterials();
    }

    void addAbsorbingBoundary(int order = 1)
    {
        for(int i=0; i < 2*D; ++i)
            addAbsorbingBoundary(Boundary{i}, order);
    }

    void addAbsorbingBoundary(Boundary b, int order = 1)
    {
        ABC<D> abc;
        int I = toAxis(b);
        bool atEnd = int(b) % 2;

        abc.order = order;
        abc.cells = boundaryCells(b, 1);  // TODO: Exclude unneeded corner/edge
        abc.oldValues.resize(abc.cells.size() * (order == 2 ? 6 : 1));

        abc.direction = basisVec<D>(I, atEnd ? -1 : 1);

        for (auto* p : boundaryFields(b))
        {
            abc.field = p;
            if (atEnd)
                abcsB.push_back(abc);
            else
                abcsA.push_back(abc);
        }
    }

    void addHardSource(VecNi<D> pos)
    {
        hardSrcsA.push_back(Source{to_idx(pos)});
    }

    void addPlaneSource(Boundary b)
    {
        forEachIndex(boundaryRect(b), [&](int index) { addHardSource(to_vec(index)); });
    }

    std::vector<int> boundaryCells(Boundary b, int width = 1, int offsetMin = 0, int offsetMax = 0) const
    {
        auto rect = boundaryRect(b, width, offsetMin, offsetMax);
        return subgridCells(rect);
    }

    Rect<D, int> boundaryRect(Boundary b, int width = 1, int offsetMin = 0, int offsetMax = 0) const
    {
        int d = toAxis(b);
        bool atMax = int(b) % 2;
        VecNi<D> n = basisVec<D, int>(d);
        VecNi<D> t = 1 - n;
        VecNi<D> corner1 = t * offsetMin;
        VecNi<D> corner2 = N - t * offsetMax;

        if(atMax) corner1 += n * (N[d] - width);
        else      corner2 -= n * (N[d] - width);

        return clamp(RectNi<D>{corner1, corner2}, RectNi<D>{{}, N});
    }

    std::vector<int> subgridCells(const RectNi<D>& bounds) const
    {
        std::vector<int> cells;
        cells.reserve(bounds.volume());

        forEachIndex(bounds, [&](int index) { cells.push_back(index); });

        return cells;
    }

    // bounds: {0,N} -> N cells, N+1 points
    template <class F>
    void forEachIndex(const RectNi<D>& bounds, F f) const
    {
        forEachCell(bounds, [&](VecNi<D> p) { f(to_idx(p)); });
    }

    template <class F>
    void forEachIndex(const RectNi<D>& bounds, F f)
    {
        forEachCell(bounds, [&](VecNi<D> p) { f(to_idx(p)); });
    }

    // Grid position to cell index
    int to_idx(VecNi<D> v) const
    {
        // inner_product(v, Vec{1, N.x, N.x*N.y})
        int idx = 0, product = 1;
        for (int i = 0; i < D; ++i)
        {
            idx += v[i] * product;
            product *= N[i];
        }
        return idx;
    }

    // Cell index to grid position
    VecNi<D> to_vec(int c) const  // to_pos
    {
        VecNi<D> v;
        int product = 1;
        for (int i = 0; i < D; ++i)
        {
            v[i] = c % (product * N[i]) / product;
            product *= N[i];
        }
        return v;
    }


 private:
    void updateFieldA()
    {
        if constexpr(D == 1)
        {
            update1D(Az, By, cAA, cAB, 1, true, +1);
        }
        else if constexpr(D == 2)
        {
            update2D(Az, Bx, By, cAA, cAB, {1, 0}, {0, 1}, true);
        }
        else if constexpr(D == 3)
        {
            std::vector<Scalar>* A[3] = {&Ax, &Ay, &Az};
            std::vector<Scalar>* B[3] = {&Bx, &By, &Bz};

            for (int i = 0; i < 3; ++i)
            {
                int j = (i+1) % 3;
                int k = (i+2) % 3;
                update2D(*A[i], *B[j], *B[k], cAA, cAB,
                         basisVec<3>(j, 1), basisVec<3>(k, 1), true);
            }
        }
    }

    void updateFieldB()
    {
        if constexpr(D == 1)
        {
            update1D(By, Az, cBB, cBA, 1, false, +1);
        }
        else if constexpr(D == 2)
        {
            update1D(By, Az, cBB, cBA, {1, 0}, false, +1);
            update1D(Bx, Az, cBB, cBA, {0, 1}, false, -1);
        }
        else if constexpr(D == 3)
        {
            std::vector<Scalar>* A[3] = {&Ax, &Ay, &Az};
            std::vector<Scalar>* B[3] = {&Bx, &By, &Bz};

            for (int i = 0; i < 3; ++i)
            {
                int j = (i+1) % 3;
                int k = (i+2) % 3;
                update2D(*B[i], *A[j], *A[k], cBB, cBA,
                         basisVec<3>(j, 1), basisVec<3>(k, 1), false);
            }
        }
    }

    // One component of dA/dt = +/- rot(B)  (with one component of B = const) (Central FD)
    // sign: (sign of equation) * (sign of non-zero rot(B)-component)
    // shift: Skip first value in each line (A-field)
    void update1D(std::vector<Scalar>& A, std::vector<Scalar>& B,
                  std::vector<Scalar>& c1, std::vector<Scalar>& c2,
                  VecNi<D> dir, bool shift, int sign)
    {
        // 2D: Skip first line (B-field)
        auto n_shift = (1 - dir) * int(!shift);
        auto start = n_shift + dir * int(shift);

        forEachCell(RectNi<D>{start, N - 1 + start}, [&](VecNi<D> pos)
        {
            int i = to_idx(pos);
            int j1 = to_idx(pos - dir * int(shift)), j2 = to_idx(pos + dir * int(!shift));

            A[i] = A[i] * c1[i] + (B[j2] - B[j1]) * c2[i] * sign;
        });
    }

    // One component of dA/dt = rot(B) * sign  (Central FD)
    // (Example: dAz/dt = dBy/dx - dBx/dy => A = Az, B1 = Bx, B2 = By, sign = +1)
    // shift: Skip first row/col (A-field)
    void update2D(std::vector<Scalar>& A, std::vector<Scalar>& B1, std::vector<Scalar>& B2,
                  std::vector<Scalar>& c1, std::vector<Scalar>& c2,
                  VecNi<D> dir1, VecNi<D> dir2, bool shift)
    {
        int sign = 1 - 2*int(!shift);
        // dir1 = basisVec<3, int>((axis+1)%3), dir2 = basisVec<3, int>((axis+2)%3)

        // 3D: Skip first plane (B-field)
        auto n_shift = (1 - dir1 - dir2) * int(!shift);
        auto start = n_shift + (dir1 + dir2) * int(shift);

        forEachCell(RectNi<D>{start, N - 1 + start}, [&](VecNi<D> pos)
        {
            int i = to_idx(pos);
            int j1 = to_idx(pos - dir1 * int(shift)), j2 = to_idx(pos + dir1 * int(!shift));
            int k1 = to_idx(pos - dir2 * int(shift)), k2 = to_idx(pos + dir2 * int(!shift));

            A[i] = A[i] * c1[i] + ((B2[j2] - B2[j1]) - (B1[k2] - B1[k1])) * c2[i] * sign;
        });
    }

    // Sources + BCs
    void applySourcesA(Scalar q)
    {
        applyAbcs(abcsA);
        applyAdditiveSources(additiveSrcsA, Az, q);
        applyHardSources(hardSrcsA, Az, q);
    }

    // Sources + BCs
    void applySourcesB(Scalar q)
    {
        applyAbcs(abcsB);
        applyAdditiveSources(additiveSrcsB, By, q);
        applyHardSources(hardSrcsB, By, q);

        // Combined source for A and B
        for (auto& src : tfsfSrcs)
            applyTfsfSource(src, q);
    }

    void applyHardSources(std::vector<Source>& srcs, std::vector<Scalar>& F, Scalar q)
    {
        Scalar t = q;
        for (auto& src : srcs)
            F[src.cell] = src.f(t);
        Std::erase_if(srcs, [t](const Source& s) { return s.duration > 0 && t > s.duration; });
    }

    void applyAdditiveSources(std::vector<Source>& srcs, std::vector<Scalar>& F, Scalar q)
    {
        Scalar t = q;
        for (auto& src : srcs)
            F[src.cell] += src.f(t);
        Std::erase_if(srcs, [t](const Source& s) { return s.duration > 0 && t > s.duration; });
    }

    void applyAbc(ABC<D>& abc)
    {
        std::vector<Scalar>& F = *abc.field;

        // TODO: Calculate coeffs once, use for all nodes (assume same material)
        for (unsigned i = 0; i < abc.cells.size(); ++i)
        {
            int n0 = abc.cells[i];
            int n1 = to_idx(to_vec(n0) + abc.direction);

            if (abc.order == 0)
            {
                // Simple version (Assume Courant number = 1)
                F[n0] = std::exchange(abc.oldValues[i], F[n1]);
            }
            else if (abc.order == 1)
            {
                Scalar coeff = std::sqrt(cAB[n0] * cBA[n0]);
                coeff = (coeff - 1) / (coeff + 1);
                Scalar F1_old = std::exchange(abc.oldValues[i], F[n1]);
                F[n0] = F1_old + (F[n1] - F[n0]) * coeff;
            }
            else if (abc.order == 2)
            {
                int n[3] = {n0, n1, to_idx(to_vec(n0) + 2*abc.direction)};
                Scalar S = std::sqrt(cAB[n[0]] * cBA[n[0]]);
                Scalar factor = 1/S + 2 + S;
                Scalar coeff[3] = {-(1/S - 2 + S) / factor, -2*(S - 1/S) / factor, 4*(S + 1/S) / factor};
                Scalar* F_old1 = &abc.oldValues[6*i];
                Scalar* F_old2 = &abc.oldValues[6*i + 3];
                F[n[0]] = (F[n[2]] + F_old2[0]) * coeff[0]
                        + (F_old1[0] + F_old1[2] - F[n[1]] - F_old2[1]) * coeff[1]
                        + F_old1[1] * coeff[2] - F_old2[2];

                for (int j = 0; j < 3; ++j)
                    F_old2[j] = std::exchange(F_old1[j], F[n[j]]);
            }
        }
    }

    void applyAbcs(std::vector<ABC<D>>& abcs)
    {
        for (auto& abc : abcs)
            applyAbc(abc);
    }

    void applyTfsfSource(TFSF<D> src, Scalar q);


    // Field components parallel to the boundary
    std::array<std::vector<Scalar>*, 2>
        boundaryFields3D(Boundary b)
    {
        auto X = (int(b) % 2) ? std::array{&Bx, &By, &Bz}
                              : std::array{&Ax, &Ay, &Az};
        int t[2] = {(toAxis(b) + 1) % 3, (toAxis(b) + 2) % 3};

        return {X[t[0]], X[t[1]]};
    }

    // Field components parallel to the boundary
    // (Filters out unused fields for 2D and 1D)
    std::vector<std::vector<Scalar>*>
        boundaryFields(Boundary b)
    {
        std::vector<std::vector<Scalar>*> x;

        for (auto* p : boundaryFields3D(b))
            if (!p->empty())
                x.push_back(p);
        return x;

        // C++20:
        //auto v = std::views::filter(boundaryFields3D(b), [](auto* p){ return !p->empty(); });
        //return std::vector(v.begin(), v.end());
    }

    Scalar excite(Scalar t)
    {
        return fn::gauss(t);
    }

    Scalar excite(Scalar t, Scalar x_shift)  // Shifted wave f(t-x/c) for TF/SF
    {
        return excite(t - x_shift / S_c);  // Assumes empty space (c = c0)
    }

    void applyMaterials()
    {
        for (int i = 0; i < RectNi<D>{{}, N}.volume(); ++i)
        {
            // TODO: Interpolate materials
            UpdateCoeffs uc(mats[matIds[i]], refMat, dx, dt);
            cAA[i] = uc.cAA;
            cAB[i] = uc.cAB;
            cBB[i] = uc.cBB;
            cBA[i] = uc.cBA;
        }
    }

    void print()
    {
        ensureInitialState();

        if (timeStep % frameInterval) return;

        std::filesystem::create_directories("out/plot/2d");

        if constexpr(D == 1)
            print1D(RectNi<D>{{}, N}, Az);
        else if constexpr(D == 2)
            print2D(RectNi<D>{{}, N}, Az);
            //print1D(Rect2i{{0, N.y/2}, {N.x, N.y/2 + 1}}, Az);
        else
            //print3D(RectNi<D>{{}, N}, Az);
            print2D(Rect3i{{0, 0, N.z/2}, {N.x, N.y, N.z/2 + 1}}, Az);
    }

    // Assumes Rect is a 1D-slice (width 1) for D > 1
    void print1D(const RectNi<D>& rect, const std::vector<Scalar>& A)
    {
        std::ofstream out("out/fdtd1d." + std::to_string(frame++) + ".txt");

        forEachIndex(rect, [&](int i){ out << A[i] << "\n"; });
    }

    // Assumes Rect is a 2D-slice (thickness 1) for D = 3
    void print2D(const RectNi<D>& rect, const std::vector<Scalar>& A)
    {
        std::ofstream out("out/fdtd2d." + std::to_string(frame++) + ".txt");

        // Assume slice is in xy-plane
        Vec first = rect.min, last = rect.max - 1;
        for (int j = 0; j <= last[1] - first[1]; ++j)
        {
            for (int i = 0; i <= last[0] - first[0]; ++i)
                out << A[to_idx(first + resize<D>(Vec2i(i, j)))] << " ";
            out << "\n";
        }
        //for (...)
        //    forEachIndex(/*1D-slice*/, ...)
    }

    void print3D(const Rect3i& rect, const std::vector<Scalar>& A)
    {
    }
};


template <int D>
struct TFSF
{
    TFSF(RectNi<D> tfRegion, Scalar S_c)
      : bounds(tfRegion),
        sim1d(tfRegion.size()[0] + 20, S_c)
    {
        sim1d.addHardSource(0);
    }

    RectNi<D> bounds;

    FDTD<1> sim1d;
};


template <int D>
void FDTD<D>::applyTfsfSource(TFSF<D> src, Scalar q)
{
    // SF Correction

    // TF Update
}


} // namespace fdtd


using fdtd::FDTD;

} // namespace fidi


#endif // FIDITI_FDTD_HPP




/**
 * Notes:
 *
 * A grid of N cells is represented by Rect<D, int> = {0, N} with N+1 points 0,1,...,N
 * => Iterate from 0 to N-1 (forEachCell, forEachIndex)
 * => rect.volume() is the number of cells
 *
 *
 * A is assumed to be the electric E-field, B the magnetic H-field
 *
 * 3D: A-components are centered at cell edges, B at cell faces (dual edges)
 *
 * 2D/1D: TB(z) assumed (normal Az, transverse Bx, By) -> Later allow choosing TA(z)
 *
 * The EM update equations for E and H differ only in the sign:
 *  dD/dt = +rot(H)  vs  dB/dt = -rot(E)  (D = eps*E, B = mu*H)
 * Pulling the sign into the update coefficients (cAB or cBA) would make them symmetric?
 * => Could probably allow swapping E and H (A = H, B = E)
 *
 *
 *
 * To Do:
 *
 * Move sign of update-equations into coefficient cBA? -> Must add sign in ABC (cAB*cBA) ...
 *
 * Automate update equations (dA/dt = rot(B)):
 * - update2D(A[3], B[3], axis, cAA, cAB) -> Automatically calculates direction
 * - update1D(A[3], B[3], axisA, axisB, cAA, cAB) ?
 *
 * update1D -> update2D with zero placeholder for unused B-comp
 *
 * template <class Scalar>
 *
 * ABC: Exclude first (last) line on side parallel (orthogonal) to absorbed component (Corner in 2D)
 *
 * TFSF: Apply in 2 separate steps (A/B) later?
 *
 *
 *
 * Optimization:
 *
 * Store cell index offset for all grid directions (idxStride? N_offset?), avoid recomputing for each cell
 * (Use in to_idx(), and optimized advanceIndex(axis, d), or just add c + d*N_offset[i])
 * (advanceIndex: assert(idx_new >= 0 && idx_new < Rect{0,N}.volume()))
 *
 * Parallelize: Specialized FDTD::forEachCell() with OpenMP
 * Alternative: Do orthogonal 3D plane updates with 3 threads
**/
