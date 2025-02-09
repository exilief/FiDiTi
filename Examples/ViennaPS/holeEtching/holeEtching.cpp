#include "FDTD.hpp"

#include <geometries/psMakeHole.hpp>
#include <models/psSF6O2Etching.hpp>

#include <psProcess.hpp>
#include <psUtils.hpp>

namespace ps = viennaps;
namespace cs = viennacs;

int main(int argc, char *argv[]) {
  using NumericType = double;
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
  auto geometry = ps::SmartPointer<ps::Domain<NumericType, D>>::New();
  ps::MakeHole<NumericType, D>(
      geometry, params.get("gridDelta") /* grid delta */,
      params.get("xExtent") /*x extent*/, params.get("yExtent") /*y extent*/,
      params.get("holeRadius") /*hole radius*/,
      params.get("maskHeight") /* mask height*/,
      params.get("taperAngle") /* tapering angle in degrees */,
      0 /* base height */, false /* periodic boundary */, true /*create mask*/,
      ps::Material::Si)
      .apply();

  // use pre-defined model SF6O2 etching model
  auto model = ps::SmartPointer<ps::SF6O2Etching<NumericType, D>>::New(
      params.get("ionFlux") /*ion flux*/,
      params.get("etchantFlux") /*etchant flux*/,
      params.get("oxygenFlux") /*oxygen flux*/,
      params.get("meanEnergy") /*mean energy*/,
      params.get("sigmaEnergy") /*energy sigma*/,
      params.get("ionExponent") /*source power cosine distribution exponent*/,
      params.get("A_O") /*oxy sputter yield*/,
      params.get("etchStopDepth") /*max etch depth*/);

  // process setup
  ps::Process<NumericType, D> process;
  process.setDomain(geometry);
  process.setProcessModel(model);
  process.setMaxCoverageInitIterations(10);
  process.setNumberOfRaysPerPoint(params.get("raysPerPoint"));
  process.setProcessDuration(params.get("processTime"));

  // print initial surface
  geometry->saveSurfaceMesh("initial.vtp");

  // run the process
  //process.apply();

  // print final surface
  geometry->saveSurfaceMesh("final.vtp");




  auto levelSets = geometry->getLevelSets();
  auto materialMap = geometry->getMaterialMap();

  levelSets = changeGridSpacing(levelSets, 0.5);

  // Create a cellSet
  auto cellSet = cs::SmartPointer<cs::DenseCellSet<NumericType, D>>::New();
  cellSet->setCellSetPosition(/*isAboveSurface*/ true);
  cellSet->fromLevelSets(levelSets, materialMap ? materialMap->getMaterialMap() : nullptr, params.get("maskHeight") + 5.);
  cellSet->writeVTU("initial.vtu");

  //geometry->generateCellSet(-5.0, ps::Material::Air, /*isAboveSurface*/ true)
  //auto& cellSet = geometry->getCellSet();

  runFDTD(*cellSet);
  std::cout << "FDTD done\n";

  cellSet->writeVTU("final.vtu");
}
