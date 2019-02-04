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

  std::vector<LatticeFermion> BasisOriginal(nBasis, FGrid); for (auto & elem : BasisOriginal) random(FPRNG, elem);
  std::vector<LatticeFermion> BasisNewFewerArgs(nBasis, FGrid); for (auto i = 0; i < BasisOriginal.size(); ++i) BasisNewFewerArgs[i] = BasisOriginal[i];
  std::vector<LatticeFermion> BasisNewMoreArgs(nBasis, FGrid); for (auto i = 0; i < BasisOriginal.size(); ++i) BasisNewMoreArgs[i] = BasisOriginal[i];

  LatticeFermion FineXOriginal(FGrid);
  LatticeFermion FineYOriginal(FGrid);
  LatticeFermion FineZOriginal(FGrid);
  LatticeFermion FineXNewFewerArgs(FGrid);
  LatticeFermion FineYNewFewerArgs(FGrid);
  LatticeFermion FineZNewFewerArgs(FGrid);
  LatticeFermion FineXNewMoreArgs(FGrid);
  LatticeFermion FineYNewMoreArgs(FGrid);
  LatticeFermion FineZNewMoreArgs(FGrid);

  random(FPRNG, FineXOriginal);
  random(FPRNG, FineYOriginal);
  random(FPRNG, FineZOriginal);
  FineXNewFewerArgs = FineXOriginal;
  FineYNewFewerArgs = FineYOriginal;
  FineZNewFewerArgs = FineZOriginal;
  FineXNewMoreArgs = FineXOriginal;
  FineYNewMoreArgs = FineYOriginal;
  FineZNewMoreArgs = FineZOriginal;

  typedef Lattice<vTComplex> CoarseScalar;

  CoarseScalar InnerProdOriginal(CGrid);
  CoarseScalar InnerProdNewFewerArgs(CGrid);
  CoarseScalar InnerProdNewMoreArgs(CGrid);

  auto FSiteVecElems = getSiteElems<LatticeFermion>();
  auto FSiteScalarElems = getSiteElems<decltype(innerProduct(FineXOriginal._odata[0], FineYOriginal._odata[0]))>();
  auto CSiteScalarElems = getSiteElems<CoarseScalar>();

  std::cout << FSiteVecElems << " " << CSiteScalarElems << std::endl;

  double FVolume = std::accumulate(FGrid->_fdimensions.begin(), FGrid->_fdimensions.end(), 1, std::multiplies<double>());
  double CVolume = std::accumulate(CGrid->_fdimensions.begin(), CGrid->_fdimensions.end(), 1, std::multiplies<double>());

  // These belong to blockInnerProduct
  double flopLocalInnerProduct = FVolume * (8 * FSiteVecElems - 2);
  double byteLocalInnerProduct = FVolume * (2 * FSiteVecElems + 1 * FSiteScalarElems) * sizeof(Complex);

  double flopBlockSum = FVolume * (2 * FSiteScalarElems);
  double byteBlockSum = FVolume * (2 * CSiteScalarElems + 1 * FSiteScalarElems) * sizeof(Complex);

  double flopCopy = CVolume * 0;
  double byteCopy = CVolume * (2 * CSiteScalarElems) * sizeof(Complex); // TODO: In terms of user intended byte, this should be zero

  double flopBlockInnerProduct = flopLocalInnerProduct + flopBlockSum + flopCopy;
  double byteBlockInnerProduct = byteLocalInnerProduct + byteBlockSum + byteCopy;

  // These belong to blockNormalise
  double flopPow = CVolume * (1); // TODO: I didn't know what to put for the 1
  double bytePow = CVolume * (2 * CSiteScalarElems) * sizeof(Complex);

  double flopBlockZAXPY = FVolume * (8 * FSiteVecElems);
  double byteBlockZAXPY = FVolume * (3 * FSiteVecElems + 1 * CSiteScalarElems) * sizeof(Complex); // TODO: Correct this number in the zaxpy benchmark file

  double flopBlockNormalise = flopBlockInnerProduct + flopPow + flopBlockZAXPY;
  double byteBlockNormalise = byteBlockInnerProduct + bytePow + byteBlockZAXPY;

  // These belong to blockOrthogonalise
  double flopMinus = CVolume * (6 * CSiteScalarElems) * sizeof(Complex); // TODO: This is actually a real multiplication
  double byteMinus = CVolume * (2 * CSiteScalarElems) * sizeof(Complex);

  double flopBlockOrthogonalise = flopBlockNormalise * nBasis + (flopBlockInnerProduct + flopMinus + flopBlockZAXPY) * nBasis * (nBasis - 1) / 2.;
  double byteBlockOrthogonalise = byteBlockNormalise * nBasis + (byteBlockInnerProduct + byteMinus + byteBlockZAXPY) * nBasis * (nBasis - 1) / 2.;

  CoarseningLookUpTable lookUpTable(CGrid, FGrid);

  BenchmarkFunction(OriginalImpl::blockOrthogonalise, flopBlockOrthogonalise, byteBlockOrthogonalise, nIter, InnerProdOriginal,     BasisOriginal);
  BenchmarkFunction(blockOrthogonalise,               flopBlockOrthogonalise, byteBlockOrthogonalise, nIter, InnerProdNewFewerArgs, BasisNewFewerArgs);
  BenchmarkFunction(blockOrthogonalise,               flopBlockOrthogonalise, byteBlockOrthogonalise, nIter, InnerProdNewMoreArgs,  BasisNewMoreArgs, lookUpTable);

  for (auto i = 0; i < BasisOriginal.size(); ++i) printDeviationFromReference(BasisOriginal[i], BasisNewFewerArgs[i]);
  for (auto i = 0; i < BasisOriginal.size(); ++i) printDeviationFromReference(BasisOriginal[i], BasisNewMoreArgs[i]);

  BenchmarkFunction(OriginalImpl::blockInnerProduct, flopBlockInnerProduct, byteBlockInnerProduct, nIter, InnerProdOriginal,     FineXOriginal,     FineYOriginal);
  BenchmarkFunction(blockInnerProduct,               flopBlockInnerProduct, byteBlockInnerProduct, nIter, InnerProdNewFewerArgs, FineXNewFewerArgs, FineYNewFewerArgs);
  BenchmarkFunction(blockInnerProduct,               flopBlockInnerProduct, byteBlockInnerProduct, nIter, InnerProdNewMoreArgs,  FineXNewMoreArgs,  FineYNewMoreArgs, lookUpTable);

  printDeviationFromReference(InnerProdOriginal, InnerProdNewFewerArgs);
  printDeviationFromReference(InnerProdOriginal, InnerProdNewMoreArgs);

  BenchmarkFunction(OriginalImpl::blockNormalise, flopBlockNormalise, byteBlockNormalise, nIter, InnerProdOriginal,     FineXOriginal);
  BenchmarkFunction(blockNormalise,               flopBlockNormalise, byteBlockNormalise, nIter, InnerProdNewFewerArgs, FineXNewFewerArgs);
  BenchmarkFunction(blockNormalise,               flopBlockNormalise, byteBlockNormalise, nIter, InnerProdNewMoreArgs,  FineXNewMoreArgs, lookUpTable);

  printDeviationFromReference(FineXOriginal, FineXNewFewerArgs);
  printDeviationFromReference(FineXOriginal, FineXNewMoreArgs);

  BenchmarkFunction(OriginalImpl::blockZAXPY, flopBlockZAXPY, byteBlockZAXPY, nIter, FineZOriginal,     InnerProdOriginal,     FineXOriginal,     FineYOriginal);
  BenchmarkFunction(blockZAXPY,               flopBlockZAXPY, byteBlockZAXPY, nIter, FineZNewFewerArgs, InnerProdNewFewerArgs, FineXNewFewerArgs, FineYNewFewerArgs);
  BenchmarkFunction(blockZAXPY,               flopBlockZAXPY, byteBlockZAXPY, nIter, FineZNewMoreArgs,  InnerProdNewMoreArgs,  FineXNewMoreArgs,  FineYNewMoreArgs, lookUpTable);

  printDeviationFromReference(FineZOriginal, FineZNewFewerArgs);

  if (doPerfProfiling) {
    PerfProfileFunction(OriginalImpl::blockOrthogonalise, nIter, InnerProdOriginal,     BasisOriginal);
    PerfProfileFunction(blockOrthogonalise,               nIter, InnerProdNewFewerArgs, BasisNewFewerArgs);
    PerfProfileFunction(blockOrthogonalise,               nIter, InnerProdNewMoreArgs,  BasisNewMoreArgs, lookUpTable);

    PerfProfileFunction(OriginalImpl::blockInnerProduct, nIter, InnerProdOriginal,     FineXOriginal,     FineYOriginal);
    PerfProfileFunction(blockInnerProduct,               nIter, InnerProdNewFewerArgs, FineXNewFewerArgs, FineYNewFewerArgs);
    PerfProfileFunction(blockInnerProduct,               nIter, InnerProdNewMoreArgs,  FineXNewMoreArgs,  FineYNewMoreArgs, lookUpTable);

    PerfProfileFunction(OriginalImpl::blockNormalise, nIter, InnerProdOriginal,     FineXOriginal);
    PerfProfileFunction(blockNormalise,               nIter, InnerProdNewFewerArgs, FineXNewFewerArgs);
    PerfProfileFunction(blockNormalise,               nIter, InnerProdNewMoreArgs,  FineXNewMoreArgs, lookUpTable);

    PerfProfileFunction(OriginalImpl::blockZAXPY, nIter, FineZOriginal,     InnerProdOriginal,     FineXOriginal,     FineYOriginal);
    PerfProfileFunction(blockZAXPY,               nIter, FineZNewFewerArgs, InnerProdNewFewerArgs, FineXNewFewerArgs, FineYNewFewerArgs);
    PerfProfileFunction(blockZAXPY,               nIter, FineZNewMoreArgs,  InnerProdNewMoreArgs,  FineXNewMoreArgs,  FineYNewMoreArgs, lookUpTable);
  }

  Grid_finalize();
}
