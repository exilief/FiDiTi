#pragma once

#include "FiDiTi/FDTD.hpp"

#include <lsToMesh.hpp>
#include <lsFromMesh.hpp>
#include <csDenseCellSet.hpp>

#include <cmath>

template <class Scalar, int D>
auto changeGridSpacing(viennals::SmartPointer<viennals::Domain<Scalar, D>> levelSet, Scalar multiplier)
{
    auto mesh = viennals::SmartPointer<viennals::Mesh<Scalar>>::New();
    viennals::ToSurfaceMesh<Scalar, D>(levelSet, mesh).apply();

    const auto& grid = levelSet->getGrid();
    auto min = grid.getMinBounds(), max = grid.getMaxBounds();
    auto bcs = grid.getBoundaryConditions();
    for (int i = 0; i < D; i++)
    {
        // Undo -1 change done in Grid constructor...
        if (bcs[i] == hrleGrid<D>::PERIODIC_BOUNDARY)
            min[i]++;

        min[i] = std::floor(min[i] / multiplier);
        max[i] = std::ceil(max[i] / multiplier);
    }
    hrleGrid<D> newGrid(&min[0], &max[0], grid.getGridDelta() * multiplier, &bcs[0]);

    auto newLevelSet = decltype(levelSet)::New(newGrid);
    viennals::FromSurfaceMesh<Scalar, D>(newLevelSet, mesh).apply();

    return newLevelSet;
}

template <class Scalar, int D>
auto changeGridSpacing(const std::vector<viennals::SmartPointer<viennals::Domain<Scalar, D>>>& levelSets, Scalar multiplier)
{
    std::decay_t<decltype(levelSets)> newLevelSets;
    for (const auto& ls : levelSets)
        newLevelSets.push_back(changeGridSpacing(ls, multiplier));
    return newLevelSets;
}

template <class Scalar, int D>
auto extractMaterials(viennacs::DenseCellSet<Scalar, D>& cellSet, fidi::Vec<D, int> N, fidi::Rect<D, int> innerBounds)
{
    const auto& csMatIds = *cellSet.getScalarData("Material");

    std::vector<int> csMatIdsInt(csMatIds.begin(), csMatIds.end());

    return fidi::resizeField(csMatIdsInt, innerBounds.size(), fidi::Rect{-innerBounds.min, N - innerBounds.min});
}

template <class Scalar, int D>
auto shrinkField(const std::vector<Scalar>& field, fidi::Vec<D, int> fieldSize, int width)
{
    fidi::Rect bounds(fieldSize - width);
    bounds.min += width;

    return fidi::resizeField(field, fieldSize, bounds);
}

template <class Scalar, int D>
void setSphereMaterial(std::vector<int>& field, fidi::VecNi<D> fieldSize, fidi::Rect<D, Scalar> bounds,
                       fidi::Vec<D, Scalar> center, Scalar radius, int material)
{
    fidi::RectNi<D> box(bounds);
    auto idxStride = indexStride(fieldSize);

    forEachCell(box, [&](fidi::VecNi<D> pos)
    {
        fidi::Vec<D, Scalar> p = pos;
        if (len_sq(p - center) <= radius*radius)
            field[dot(p, idxStride)] = material;
    });
}

template <class Scalar, int D>
auto print(const fidi::FDTD<D>& sim, viennacs::DenseCellSet<Scalar, D>& cellSet,
           fidi::Rect<D, int> innerBounds, int step)
{
    fidi::Vec<D, int> N = sim.gridSize();

    auto printField = [&] (const auto& field, std::string name)
    {
        if (field.empty()) return;
        auto* data = cellSet.addScalarData(name);
        *data = resizeField(field, N, innerBounds);
    };
    for (int i = 0; i < 3; ++i)
    {
        printField(sim.fieldA(i), "E_" + std::to_string(i+1));
        printField(sim.fieldB(i), "H_" + std::to_string(i+1));
    }
    cellSet.writeVTU("FDTD_cellSet_" + std::to_string(step) + ".vtu");
}

template <class Scalar, int D>
void runFDTD(viennacs::DenseCellSet<Scalar, D>& cellSet)
{
    using namespace fidi;

    auto numCells = cellSet.getNumberOfCells();
    auto dx = cellSet.getGridDelta();
    auto [min, max] = cellSet.getBoundingBox();

    Vec<D, int> csGridSize;
    for (int i = 0; i < D; ++i)
        csGridSize[i] = std::round((max[i] - min[i]) / dx);

    // Manually extend grid on all sides for absorbing boundary
    int extend = 5;
    Vec<D, int> N = csGridSize + 2*extend;
    Rect<D, int> innerBounds = {Vec<D, int>(extend), N - extend};

    FDTD<D> sim(N);

    fdtd::MaterialMap matMap;
    matMap.emplace(0, fdtd::Material{1,1});  // Vacuum
    matMap.emplace(1, fdtd::Material{1.5,1});  // {1,1000000}
    matMap.emplace(-1, fdtd::Material{1,1});

    sim.setMaterials(matMap);
    sim.setCellMaterials(extractMaterials(cellSet, N, innerBounds));

    //sim.addHardSource(N/2);
    sim.addTfsfSource(1);
    sim.addAbsorbingBoundary(2);

    int numSteps = 150;
    int frameInterval = 5;
    for (int q = 0; q < numSteps; ++q)
    {
        sim.step();

        if (!(q % frameInterval))
            print(sim, cellSet, innerBounds, q+1);
    }
    //sim.run(150);
}
