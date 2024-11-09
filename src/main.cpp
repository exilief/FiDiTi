#include <FiDiTi/FDTD.hpp>
#include <iostream>

using namespace fidi;

namespace
{
    template <int D>
    void setMaterials(FDTD<D>& sim, VecNi<D> N)
    {
        std::vector<int> matIds(RectNi<D>{{}, N}.volume());
        for (int i = 0; i < 2*D; ++i)
        {
            auto cells = sim.boundaryCells(fdtd::Boundary{i}, 2);
            for (auto c : cells)
                matIds[c] = 1;
        }
        sim.setCellMaterials(std::move(matIds));

        fdtd::MaterialMap mats;
        mats.emplace(0, fdtd::Material{1,1});
        mats.emplace(1, fdtd::Material{1,1000000});
        sim.setMaterials(std::move(mats));
    }
}


int main()
{
    constexpr int D = 2;

    VecNi<D> N = resize<D>(Vec3i(52, 34, 8));
    FDTD<D> sim(N);
    sim.addHardSource(N/2);
    //sim.addPlaneSource(fdtd::Boundary::XMin);
    if(0)
        setMaterials<D>(sim, N);
    else
        sim.addAbsorbingBoundary(2);
    sim.run(152 + N[0]);
}
