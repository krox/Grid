/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./benchmarks/Benchmark_blockorthogonalise.cc

    Copyright (C) 2015 - 2018

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
#include <Benchmark_helpers.h>

using namespace Grid;
using namespace Grid::QCD;
using namespace Grid::BenchmarkHelpers;

// Enable control of nbasis from the compiler command line
// NOTE to self: Copy the value of CXXFLAGS from the makefile and call make as follows:
//   make CXXFLAGS="-DNBASIS=24 VALUE_OF_CXXFLAGS_IN_MAKEFILE" Benchmark_blocksum
#ifndef NBASIS
#define NBASIS 32
#endif

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);

  std::vector<int> seeds({1, 2, 3, 4});

  // clang-format off
  const int        nBasis          = NBASIS; static_assert((nBasis & 0x1) == 0, "");
  const int        nb              = nBasis / 2;
  int              nIter           = readFromCommandLineInt(&argc, &argv, "--niter", 10);
  std::vector<int> blockSize       = readFromCommandLineIntVec(&argc, &argv, "--blocksize", std::vector<int>({4, 4, 4, 4}));
  bool             doPerfProfiling = readFromCommandLineToggle(&argc, &argv, "--perfprofiling");
  // clang-format on

  RealD mass = 0.5;

  std::vector<int> clatt = calcCoarseLattSize(GridDefaultLatt(), blockSize);

  typedef Aggregation<vSpinColourVector, vTComplex, nBasis> Aggregates;
  typedef Aggregates::CoarseVector                          CoarseVector;

  GridCartesian *        FGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridCartesian *        CGrid   = SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid);

  FGrid->show_decomposition();
  CGrid->show_decomposition();

  GridParallelRNG FPRNG(FGrid);
  GridParallelRNG CPRNG(CGrid);

  FPRNG.SeedFixedIntegers(seeds);
  CPRNG.SeedFixedIntegers(seeds);

  std::vector<LatticeFermion> BasisOriginal(nBasis, FGrid); for (auto & elem : BasisOriginal) random(FPRNG, elem);
  std::vector<LatticeFermion> BasisNew2Args(nBasis, FGrid); for (auto i = 0; i < BasisOriginal.size(); ++i) BasisNew2Args[i] = BasisOriginal[i];
  std::vector<LatticeFermion> BasisNew3Args(nBasis, FGrid); for (auto i = 0; i < BasisOriginal.size(); ++i) BasisNew3Args[i] = BasisOriginal[i];

  LatticeFermion FineXOriginal(FGrid);
  LatticeFermion FineYOriginal(FGrid);
  LatticeFermion FineXNew(FGrid);
  LatticeFermion FineYNew(FGrid);

  random(FPRNG, FineXOriginal);
  random(FPRNG, FineYOriginal);
  FineXNew = FineXOriginal;
  FineYNew = FineYOriginal;

  typedef Lattice<vTComplex> CoarseScalar;

  CoarseScalar InnerProdOriginal(CGrid);
  CoarseScalar InnerProdNew2Args(CGrid);
  CoarseScalar InnerProdNew3Args(CGrid);

  auto FSiteElems = getSiteElems<LatticeFermion>();
  auto CSiteElems = getSiteElems<CoarseScalar>();

  std::cout << FSiteElems << " " << CSiteElems << std::endl;

  auto FVolume = std::accumulate(FGrid->_fdimensions.begin(), FGrid->_fdimensions.end(), 1, std::multiplies<double>());
  auto CVolume = std::accumulate(CGrid->_fdimensions.begin(), CGrid->_fdimensions.end(), 1, std::multiplies<double>());

  // TODO: Calculate correct values here
  double flop = 1. * (1);
  double byte = 1. * (1);

  CoarseningLookUpTable lookUpTable(CGrid, FGrid);

  BenchmarkFunction(OriginalImpl::blockOrthogonalise, flop, byte, nIter, InnerProdOriginal, BasisOriginal);
  BenchmarkFunction(blockOrthogonalise,               flop, byte, nIter, InnerProdNew2Args, BasisNew2Args);
  BenchmarkFunction(blockOrthogonalise,               flop, byte, nIter, InnerProdNew3Args, BasisNew3Args, lookUpTable);

  for (auto i = 0; i < BasisOriginal.size(); ++i) printDeviationFromReference(BasisOriginal[i], BasisNew2Args[i]);
  for (auto i = 0; i < BasisOriginal.size(); ++i) printDeviationFromReference(BasisOriginal[i], BasisNew3Args[i]);

  BenchmarkFunction(OriginalImpl::blockInnerProduct, flop, byte, nIter, InnerProdOriginal, FineXOriginal, FineYOriginal);
  BenchmarkFunction(blockInnerProduct,               flop, byte, nIter, InnerProdNew3Args, FineXOriginal, FineYOriginal, lookUpTable);

  printDeviationFromReference(InnerProdOriginal, InnerProdNew3Args);

  BenchmarkFunction(OriginalImpl::blockNormalise, flop, byte, nIter, InnerProdOriginal, FineXOriginal);
  BenchmarkFunction(blockNormalise,              flop, byte, nIter, InnerProdNew3Args, FineXNew, lookUpTable);

  printDeviationFromReference(FineXOriginal, FineXNew);

  if (doPerfProfiling) {
    PerfProfileFunction(OriginalImpl::blockOrthogonalise, nIter, InnerProdOriginal, BasisOriginal);
    PerfProfileFunction(blockOrthogonalise,               nIter, InnerProdNew2Args, BasisNew2Args);
    PerfProfileFunction(blockOrthogonalise,               nIter, InnerProdNew3Args, BasisNew3Args, lookUpTable);

    PerfProfileFunction(OriginalImpl::blockInnerProduct, nIter, InnerProdOriginal, FineXOriginal, FineYOriginal);
    PerfProfileFunction(blockInnerProduct,               nIter, InnerProdNew3Args, FineXOriginal, FineYOriginal, lookUpTable);

    PerfProfileFunction(OriginalImpl::blockNormalise, nIter, InnerProdOriginal, FineXOriginal);
    PerfProfileFunction(blockNormalise,               nIter, InnerProdNew3Args, FineXNew, lookUpTable);
  }

  Grid_finalize();
}
