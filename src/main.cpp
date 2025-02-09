#include <FiDiTi/FDTD.hpp>
#include <iostream>

using namespace fidi;

namespace
{
    Scalar toRad(Scalar deg)
    {
        return deg / 180 * constants::pi;
    }

    /*Scalar toDeg(Scalar rad)
    {
        return rad / constants::pi * 180;
    }*/

    template <int D>
    void setMaterials(FDTD<D>& sim)
    {
        fdtd::MaterialMap mats;
        mats.emplace(0, fdtd::Material{1,1});
        mats.emplace(1, fdtd::Material{1,1000000});
        sim.setMaterials(std::move(mats));
    }

    template <int D>
    void addReflectiveBoundary(FDTD<D>& sim, VecNi<D> N, std::vector<int>& matIds)
    {
        for (int i = 0; i < 2*D; ++i)
        {
            auto cells = sim.boundaryCells(fdtd::Boundary{i}, 2);
            for (auto c : cells)
                matIds[c] = 1;
        }
        sim.setCellMaterials(matIds);

        setMaterials(sim);
    }

    template <int D>
    void addScatterer(FDTD<D>& sim, VecNi<D> N, std::vector<int>& matIds)
    {
        forEachCell(Rect{N/5, N/3}, [&](VecNi<D> pos)
        {
            matIds[sim.to_idx(pos)] = 1;
        });
        sim.setCellMaterials(matIds);

        setMaterials(sim);
    }

    // 2 planes at an angle, parallel to z-axis -> Fill enclosed corners (at N.x) with material
    template <int D>
    void addOpening(FDTD<D>& sim, VecNi<D> N, std::vector<int>& matIds)
    {
        Scalar alpha = toRad(75);  // Normal angle, from x-axis (+/- at top/bottom)
        Scalar d = N.y * 0.36;     // Distance from corners

        using VecT = Vec<D, Scalar>;
        VecT n[2] = {resize<D>(Vec2<Scalar>{std::cos(alpha), std::sin(alpha)}),
                     resize<D>(Vec2<Scalar>{std::cos(-alpha), std::sin(-alpha)})};
        VecT c[2] = {VecT(N), VecT(project(N, 1))};
        Rect bounds{N*0 + 8, N-8};  // Place only in inner TFSF region

        forEachCell(Rect(N), [&](VecNi<D> pos)
        {
            VecT p{ pos };
            for (int i = 0; i < 2; ++i)
                if (dot(c[i] - p, n[i]) <= d  &&  bounds.contains(pos))
                    matIds[sim.to_idx(pos)] = 1;
        });
        sim.setCellMaterials(matIds);

        setMaterials(sim);
        static_assert(D > 1);
    }

    // Triangular obstacle, pointing in direction -x
    template <int D>
    void addWedge(FDTD<D>& sim, VecNi<D> N, std::vector<int>& matIds)
    {
        Scalar alpha = toRad(90+40);  // Normal angle of 2 faces (+/- at top/bottom)
        Scalar l = N.x * 0.12;        // Length

        using VecT = Vec<D, Scalar>;
        VecT c = VecT(N) / Scalar(2);  // Tip of the wedge
        VecT n[2] = {resize<D>(Vec2<Scalar>{std::cos(alpha), std::sin(alpha)}),
                     resize<D>(Vec2<Scalar>{std::cos(-alpha), std::sin(-alpha)})};

        forEachCell(Rect(N), [&](VecNi<D> pos)
        {
            VecT p{ pos };
            if (dot(p - c, n[0]) <= 0  &&  dot(p - c, n[1]) <= 0  &&  p[0] <= c[0] + l)
                matIds[sim.to_idx(pos)] = 1;
        });
        sim.setCellMaterials(matIds);

        setMaterials(sim);
        static_assert(D > 1);
    }
}


int main()
{
    constexpr int D = 2;

    VecNi<D> N = resize<D>(Vec3i(186, 120, 40));
    FDTD<D> sim(N);

    std::vector<int> matIds(Rect{N}.volume());

    // Set source type here:
    int s = 2;
    if(s==0) sim.addHardSource(N/2);                    // Point source
    if(s==1) sim.addPlaneSource(fdtd::Boundary::XMin);  // Array/plane of point sources
    if(s==2) sim.addTfsfSource(0, -1);                  // TF/SF plane source

    // Set boundary type here:
    int b = 1;
    if(b==0) addReflectiveBoundary<D>(sim, N, matIds);
    if(b==1) sim.addAbsorbingBoundary(2);

    //addScatterer<D>(sim, N, matIds);
    if constexpr(D > 1)
    {
        //addOpening(sim, N, matIds);
        //addWedge(sim, N, matIds);
    }

    par::set_num_threads(1);

    sim.run(152 + N[0]);
}
