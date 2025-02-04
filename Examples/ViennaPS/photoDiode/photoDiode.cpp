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

  // geometry setup
  auto geometry = makeGeometry<Scalar, D>(
      params.get("xExtent"), params.get("yExtent"),
      params.get("gridDelta"),
      params.get("bulkHeight"),
      params.get("numHoles"),
      params.get("holeWidth"),
      params.get("holeDepth"),
      params.get("maskHeight"),
      0. /*baseHeight*/,
      false /*create mask*/);
    //ps::Material::Si

#if 0
  // use pre-defined model SF6O2 etching model
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

  levelSets = changeGridSpacing(levelSets, 0.8);

  auto cellSet = cs::SmartPointer<cs::DenseCellSet<Scalar, D>>::New();
  cellSet->setCellSetPosition(/*isAboveSurface*/ true);
  cellSet->fromLevelSets(levelSets, materialMap ? materialMap->getMaterialMap() : nullptr, height);
  cellSet->writeVTU("initial.vtu");

  //geometry->generateCellSet(-5.0, ps::Material::Air, /*isAboveSurface*/ true)
  //auto& cellSet = geometry->getCellSet();

  // Add lens
  //setSphereMaterial(*cellSet->getScalarData("Material"), fieldSize, bounds, center, radius, material);

  runFDTD(*cellSet);
  std::cout << "FDTD done\n";

  cellSet->writeVTU("final.vtu");
}
