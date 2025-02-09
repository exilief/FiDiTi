#include "FDTD.hpp"
#include "geometry.hpp"

#include <models/psSF6O2Etching.hpp>

#include <psProcess.hpp>
#include <psUtils.hpp>

namespace ps = viennaps;
namespace cs = viennacs;

//const int lensMaterial = ps::Material::PolySi;

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

  bool withMask = false;

  // Geometry setup
  auto geometry = makeGeometry<Scalar, D>(
      params.get("xExtent"), params.get("yExtent"),
      params.get("gridDelta"),
      params.get("bulkHeight"),
      params.get("numHoles"),
      params.get("holeWidth"),
      params.get("holeDepth"),
      params.get("maskHeight"),
      0. /*baseHeight*/,
      withMask);

#if 0
  // SF6O2 etching
  auto model = ps::SmartPointer<ps::SF6O2Etching<Scalar, D>>::New(
      params.get("ionFlux"),
      params.get("etchantFlux"),
      params.get("oxygenFlux"),
      params.get("meanEnergy"),
      params.get("sigmaEnergy"),
      params.get("ionExponent") /*source power cosine distribution exponent*/,
      params.get("A_O") /*oxy sputter yield*/,
      params.get("etchStopDepth") /*max etch depth*/);

  // process setup
  ps::Process<Scalar, D> process;
  process.setDomain(geometry);
  process.setProcessModel(model);
  process.setMaxCoverageInitIterations(10);
  process.setNumberOfRaysPerPoint(params.get("raysPerPoint"));
  process.setProcessDuration(params.get("processTime"));
#endif

  geometry->saveSurfaceMesh("initial.vtp");

  //process.apply();

  geometry->saveSurfaceMesh("final.vtp");




  Scalar height = params.get("bulkHeight") + params.get("holeDepth") + params.get("airHeight") + params.get("gridDelta");

  auto levelSets = geometry->getLevelSets();
  auto materialMap = geometry->getMaterialMap();

  Scalar dxScale = 1.0;
  Scalar gridDelta = dxScale * params.get("gridDelta");

  if (dxScale != 1.0) levelSets = changeGridSpacing(levelSets, dxScale);

  auto cellSet = cs::SmartPointer<cs::DenseCellSet<Scalar, D>>::New();
  cellSet->setCellSetPosition(/*isAboveSurface*/ true);
  cellSet->setCoverMaterial(int(ps::Material::Air));
  cellSet->fromLevelSets(levelSets, materialMap ? materialMap->getMaterialMap() : nullptr, height);
  cellSet->writeVTU("initial.vtu");

  //geometry->generateCellSet(-5.0, ps::Material::Air, /*isAboveSurface*/ true)
  //auto& cellSet = geometry->getCellSet();

  // Add lens
  using namespace fidi;
  Vec<D, int> csGridSize = getGridSize(*cellSet);
  Rect<D, Scalar> lensBounds(getBounds(*cellSet).size());
  lensBounds.min[D-1] = params.get("bulkHeight") + params.get("holeDepth") + params.get("gridDelta");
  lensBounds.max[D-1] = lensBounds.min[D-1] + params.get("airHeight");
  setSphereMaterial(*cellSet->getScalarData("Material"), csGridSize, lensBounds,
                    project(lensBounds.max / Scalar(2), D-1), height, int(ps::Material::PolySi), gridDelta);


  fdtd::MaterialMap matMap;
  matMap.emplace(int(ps::Material::Air),     fdtd::Material{1, 1});
  matMap.emplace(int(ps::Material::Si),      fdtd::Material{params.get("siliconPermittivity"), 1});
  matMap.emplace(int(ps::Material::SiN),     fdtd::Material{params.get("passivationPermittivity"), 1});
  matMap.emplace(int(ps::Material::PolySi),  fdtd::Material{params.get("lensPermittivity"), 1});


  runFDTD(*cellSet, std::move(matMap));
  std::cout << "FDTD done\n";

  cellSet->writeVTU("final.vtu");
}
