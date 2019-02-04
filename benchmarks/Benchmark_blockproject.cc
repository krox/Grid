/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./benchmarks/Benchmark_blockproject.cc

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
//   make CXXFLAGS="-DNBASIS=24 VALUE_OF_CXXFLAGS_IN_MAKEFILE" Benchmark_blockproject
#ifndef NBASIS
#define NBASIS 40
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

  LatticeGaugeField Umu(FGrid);
  SU3::HotConfiguration(FPRNG, Umu);

  WilsonFermionR                                      Dw(Umu, *FGrid, *FrbGrid, mass);
  MdagMLinearOperator<WilsonFermionR, LatticeFermion> MdagMOp(Dw);

  Aggregates Aggs(CGrid, FGrid, 0);
  Aggs.CreateSubspace(FPRNG, MdagMOp, nb);
  performChiralDoubling(Aggs.Subspace());

  LatticeFermion FineVec(FGrid);
  CoarseVector   CoarseVec(CGrid);
  CoarseVector   CoarseVecNew(CGrid);
  CoarseVector   CoarseVecOriginal(CGrid);
  CoarseVector   CoarseVecDiff(CGrid);
  CoarseVecNew      = zero;
  CoarseVec         = zero;
  CoarseVecOriginal = zero;

  random(FPRNG, FineVec);

  auto FSiteElems = getSiteElems<decltype(FineVec)>();
  auto CSiteElems = getSiteElems<decltype(CoarseVec)>();

  double FVolume = std::accumulate(FGrid->_fdimensions.begin(), FGrid->_fdimensions.end(), 1, std::multiplies<double>());
  double CVolume = std::accumulate(CGrid->_fdimensions.begin(), CGrid->_fdimensions.end(), 1, std::multiplies<double>());

  double flop = FVolume * (8 * FSiteElems) * nBasis;
  double byte = FVolume * (2 * 1 + 2 * FSiteElems) * nBasis * sizeof(Complex);

  CoarseningLookUpTable lookUpTable(CGrid, FGrid);

  BenchmarkFunction(OriginalImpl::blockProject, flop, byte, nIter, CoarseVecOriginal, FineVec, Aggs.Subspace());
  BenchmarkFunction(blockProject,               flop, byte, nIter, CoarseVec,         FineVec, Aggs.Subspace());
  BenchmarkFunction(blockProject,               flop, byte, nIter, CoarseVecNew,      FineVec, Aggs.Subspace(), lookUpTable);

  printDeviationFromReference(CoarseVecOriginal, CoarseVec);
  printDeviationFromReference(CoarseVecOriginal, CoarseVecNew);

  if (doPerfProfiling) {
    PerfProfileFunction(OriginalImpl::blockProject, nIter, CoarseVecOriginal, FineVec, Aggs.Subspace());
    PerfProfileFunction(blockProject,               nIter, CoarseVec,         FineVec, Aggs.Subspace());
    PerfProfileFunction(blockProject,               nIter, CoarseVecNew,      FineVec, Aggs.Subspace(), lookUpTable);
  }

  Grid_finalize();
}
