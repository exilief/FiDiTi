#include "FDTD.hpp"
#include "geometry.hpp"

#include <models/psSF6O2Etching.hpp>

#include <psProcess.hpp>
#include <psUtils.hpp>

namespace ps = viennaps;
namespace cs = viennacs;

int main(int argc, char *argv[]) {
  using Scalar = double;
  constexpr int D = 2;

  ps::Logger::setLogLevel(ps::LogLevel::INTERMEDIATE);
  omp_set_num_threads(4);

  // Parse the parameters
  ps::utils::Parameters params;
  if (argc > 1) {
    params.readConfigFile(argv[1]);
  } else {
    std::cout << "Usage: " << argv[0] << " <config file>" << std::endl;
    return 1;
  }

  bool withEtching = params.get("withEtching");

  // Geometry setup
  auto domain = makeGeometry<Scalar, D>(
      params.get("xExtent"), params.get("yExtent"),
      params.get("gridDelta"),
      params.get("bulkHeight"),
      params.get("numHoles"), params.get("holeWidth"), params.get("holeDepth"),
      params.get("passivationHeight"),
      params.get("antiReflectHeight1"), params.get("antiReflectHeight2"),
      params.get("maskHeight"),
      0. /*baseHeight*/,
      withEtching);

  if (withEtching)
  {
    // Save passivation & anti-reflect layers to re-add after etching
    auto [extraLayers, extraMaterials] = extractTopLevelSets(*domain, 3);

    // SF6O2 etching model
    auto model = ps::SmartPointer<ps::SF6O2Etching<Scalar, D>>::New(
        params.get("ionFlux"),
        params.get("etchantFlux"),
        params.get("oxygenFlux"),
        params.get("meanEnergy"),
        params.get("sigmaEnergy"),
        params.get("ionExponent") /*source power cosine distribution exponent*/,
        params.get("A_O") /*oxy sputter yield*/,
        params.get("etchStopDepth") /*max etch depth*/);

    // Process setup
    ps::Process<Scalar, D> process;
    process.setDomain(domain);
    process.setProcessModel(model);
    process.setMaxCoverageInitIterations(10);
    process.setNumberOfRaysPerPoint(params.get("raysPerPoint"));
    process.setProcessDuration(params.get("processTime"));


    domain->saveSurfaceMesh("initial.vtp");

    process.apply();

    domain->saveSurfaceMesh("final.vtp");


    // Remove mask
    domain->removeTopLevelSet();

    for (int i = 0; i < extraLayers.size(); ++i)
      domain->insertNextLevelSetAsMaterial(extraLayers[i], extraMaterials[i], false);
  }




  Scalar diodeHeight = params.get("bulkHeight") - params.get("holeDepth");
  Scalar solidHeight = params.get("bulkHeight") + params.get("passivationHeight") +
                       params.get("antiReflectHeight1") + params.get("antiReflectHeight2");
  Scalar height = solidHeight + params.get("airHeight") + params.get("gridDelta");

  auto levelSets = domain->getLevelSets();
  auto materialMap = domain->getMaterialMap();

  Scalar dxScale = 0.5;
  Scalar gridDelta = dxScale * params.get("gridDelta");

  if (dxScale != 1.0) levelSets = changeGridSpacing(levelSets, gridDelta);

  auto cellSet = cs::SmartPointer<cs::DenseCellSet<Scalar, D>>::New();
  cellSet->setCellSetPosition(/*isAboveSurface*/ true);
  cellSet->setCoverMaterial(int(ps::Material::Air));
  cellSet->fromLevelSets(levelSets, materialMap ? materialMap->getMaterialMap() : nullptr, height);
  cellSet->writeVTU("initial.vtu");

  //domain->generateCellSet(-5.0, ps::Material::Air, /*isAboveSurface*/ true)
  //auto& cellSet = domain->getCellSet();

  // Add lens
  using namespace fidi;
  Vec<D, int> csGridSize = getGridSize(*cellSet);
  Rect<D, Scalar> lensBounds(getBounds(*cellSet).size());
  lensBounds.min[D-1] = solidHeight;
  lensBounds.max[D-1] = lensBounds.min[D-1] + params.get("airHeight") + params.get("gridDelta");
  setSphereMaterial(*cellSet->getScalarData("Material"), csGridSize, lensBounds,
                    project(lensBounds.max / Scalar(2), D-1), height, gridDelta,
                    int(lensMaterial), int(ps::Material::Air));


  fdtd::MaterialMap matMap;
  matMap.emplace(int(ps::Material::Air),    fdtd::Material{1, 1});
  matMap.emplace(int(ps::Material::Si),     fdtd::Material{params.get("siliconPermittivity"), 1});
  matMap.emplace(int(ps::Material::SiN),    fdtd::Material{params.get("passivationPermittivity"), 1});
  matMap.emplace(int(lensMaterial),         fdtd::Material{params.get("lensPermittivity"), 1});
  matMap.emplace(int(antiReflectMaterial1), fdtd::Material{params.get("antiReflectPermittivity1"), 1});
  matMap.emplace(int(antiReflectMaterial2), fdtd::Material{params.get("antiReflectPermittivity2"), 1});


  runFDTD(*cellSet, std::move(matMap), params.get("fdtdSteps"), diodeHeight);
  std::cout << "FDTD done\n";

  cellSet->writeVTU("final.vtu");
}
