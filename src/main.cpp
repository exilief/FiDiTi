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

    // Set source type here:
    int s = 0;
    if(s==0) sim.addHardSource(N/2);                    // Point source
    if(s==1) sim.addPlaneSource(fdtd::Boundary::XMin);  // Array/plane of point sources
    if(s==2) sim.addTfsfSource();                       // TF/SF plane source

    // Set boundary type here:
    int b = 1;
    if(b==0) setMaterials<D>(sim, N);  // Reflective boundary
    if(b==1) sim.addAbsorbingBoundary(2);

    sim.run(152 + N[0]);
}
