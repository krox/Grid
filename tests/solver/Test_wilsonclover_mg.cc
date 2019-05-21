/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/solver/Test_wilsonclover_mg.cc

    Copyright (C) 2015-2018

    Author: Daniel Richtmann <daniel.richtmann@ur.de>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
/*  END LEGAL */

#include <Grid/Grid.h>
#include <Test_multigrid_common.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;

// Enable control of nbasis from the compiler command line
// NOTE to self: Copy the value of CXXFLAGS from the makefile and call make as follows:
//   make CXXFLAGS="-DNBASIS=24 VALUE_OF_CXXFLAGS_IN_MAKEFILE" Test_wilsonclover_mg
#ifndef NBASIS
#define NBASIS 40
#endif

int main(int argc, char **argv) {

  Grid_init(&argc, &argv);

  GridCartesian *        FGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid);

  MGTestParams params;

  if(GridCmdOptionExists(argv, argv + argc, "--inputxml")) {
    std::string inputXml = GridCmdOptionPayload(argv, argv + argc, "--inputxml");
    assert(inputXml.length() != 0);

    XmlReader reader(inputXml);
    read(reader, "Params", params);
    std::cout << GridLogMessage << "Read in " << inputXml << std::endl;
  }

  // {
  //   XmlWriter writer("mg_params_template.xml");
  //   write(writer, "Params", params);
  //   std::cout << GridLogMessage << "Written mg_params_template.xml" << std::endl;
  // }

  checkParameterValidity(params);
  std::cout << params << std::endl;

  LevelInfo levelInfo(FGrid, params.mg);

  const int nbasis = NBASIS;
#if !defined(USE_TWOSPIN_COARSENING)
  static_assert((nbasis & 0x1) == 0, "");
#endif

  std::vector<int> fSeeds({1, 2, 3, 4});
  GridParallelRNG  fPRNG(FGrid);
  fPRNG.SeedFixedIntegers(fSeeds);

  LatticeFermion    src(FGrid);
  LatticeFermion    result(FGrid);
  LatticeGaugeField Umu(FGrid);

  if(params.test.sourceType == "ones")
    src = 1.;
  else if(params.test.sourceType == "random")
    random(fPRNG, src);
  else if(params.test.sourceType == "gaussian")
    gaussian(fPRNG, src);

  if(params.test.config != "foo") {
    FieldMetaData      header;
    IldgReader         _IldgReader;
    _IldgReader.open(params.test.config);
    _IldgReader.readConfiguration(Umu, header);
    _IldgReader.close();
  } else
    SU3::HotConfiguration(fPRNG, Umu);

  typename WilsonCloverFermionR::ImplParams implParams;
  WilsonAnisotropyCoefficients              anisParams;
  if(params.test.useAntiPeriodicBC)
    implParams.boundary_phases = {+1, +1, +1, -1};
  else
    implParams.boundary_phases = {+1, +1, +1, +1};

  WilsonCloverFermionR Dwc(Umu, *FGrid, *FrbGrid, params.test.mass, params.test.csw, params.test.csw, anisParams, implParams);

  MdagMLinearOperator<WilsonCloverFermionR, LatticeFermion> MdagMOpDwc(Dwc);

  TrivialPrecon<LatticeFermion> TrivialPrecon;

  std::cout << GridLogMessage << "**************************************************" << std::endl;
  std::cout << GridLogMessage << "Testing Multigrid for Wilson Clover" << std::endl;
  std::cout << GridLogMessage << "**************************************************" << std::endl;

#if defined (USE_TWOSPIN_COARSENING)
  auto MGPreconDwc = createMGInstance<vSpinColourVector,  vComplex, nbasis / 2, WilsonCloverFermionR>(params.mg, levelInfo, Dwc, Dwc);
#elif defined(USE_ONESPIN_COARSENING)
  auto MGPreconDwc = createMGInstance<vSpinColourVector,  vComplex, nbasis,     WilsonCloverFermionR>(params.mg, levelInfo, Dwc, Dwc);
#else // corresponds to original coarsening
  auto MGPreconDwc = createMGInstance<vSpinColourVector, vTComplex, nbasis,     WilsonCloverFermionR>(params.mg, levelInfo, Dwc, Dwc);
#endif

  bool  doRunChecks          = (GridCmdOptionExists(argv, argv + argc, "--runchecks")) ? true : false;
  RealD toleranceForMGChecks = (getPrecision<LatticeFermion>::value == 1) ? 1e-6 : 1e-13;

  MGPreconDwc->initialSetup();
  if(doRunChecks) MGPreconDwc->runChecks(toleranceForMGChecks);

  MGPreconDwc->iterativeSetup();
  if(doRunChecks) MGPreconDwc->runChecks(toleranceForMGChecks);

  auto outerSolverMaxIter = params.test.outerSolverMaxOuterIter * params.test.outerSolverMaxInnerIter;

  std::vector<std::unique_ptr<OperatorFunction<LatticeFermion>>> solversDwc;

  solversDwc.emplace_back(new ConjugateGradient<LatticeFermion>(params.test.outerSolverTol, 10000, false));
  solversDwc.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(params.test.outerSolverTol, 10000, TrivialPrecon, params.test.outerSolverMaxInnerIter, false));
  solversDwc.emplace_back(new FlexibleGeneralisedMinimalResidual<LatticeFermion>(params.test.outerSolverTol, outerSolverMaxIter, *MGPreconDwc, params.test.outerSolverMaxInnerIter, false));

  GridStopWatch solveTimer;
  for(auto const &solver : solversDwc) {
    solveTimer.Reset();
    std::cout << std::endl << "Starting with a new solver" << std::endl;
    result = zero;
    solveTimer.Start();
    (*solver)(MdagMOpDwc, src, result);
    solveTimer.Stop();
    std::cout << GridLogMessage << "Solver took: " << solveTimer.Elapsed() << std::endl;
  }

  MGPreconDwc->reportTimings();

  Grid_finalize();
}
