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

template<class vobj,class CComplex>
inline void blockInnerProductOriginal(Lattice<CComplex> &CoarseInner,
                              const Lattice<vobj> &fineX,
                              const Lattice<vobj> &fineY)
{
  // GridStopWatch totalTimer;
  // GridStopWatch preambleTimer;
  // GridStopWatch localInnerProductTimer;
  // GridStopWatch blockSumTimer;
  // GridStopWatch copyTimer;

  // totalTimer.Start();
  // preambleTimer.Start();
  typedef decltype(innerProduct(fineX._odata[0],fineY._odata[0])) dotp;

  GridBase *coarse(CoarseInner._grid);
  GridBase *fine  (fineX._grid);

  Lattice<dotp> fine_inner(fine); fine_inner.checkerboard = fineX.checkerboard;
  Lattice<dotp> coarse_inner(coarse);
  // preambleTimer.Stop();

  // Precision promotion?
  // localInnerProductTimer.Start();
  fine_inner = localInnerProduct(fineX,fineY);
  // localInnerProductTimer.Stop();
  // blockSumTimer.Start();
  blockSum(coarse_inner,fine_inner);
  // blockSumTimer.Stop();
  // copyTimer.Start();
  parallel_for(int ss=0;ss<coarse->oSites();ss++){
    CoarseInner._odata[ss] = coarse_inner._odata[ss];
  }
  // copyTimer.Stop();
  // totalTimer.Stop();


  // std::cout << "total             = " << totalTimer.Elapsed()             << " (" << 100. * totalTimer.useconds()/totalTimer.useconds()             << "% of total)" << std::endl;
  // std::cout << "preamble          = " << preambleTimer.Elapsed()          << " (" << 100. * preambleTimer.useconds()/totalTimer.useconds()          << "% of total)" << std::endl;
  // std::cout << "localInnerProduct = " << localInnerProductTimer.Elapsed() << " (" << 100. * localInnerProductTimer.useconds()/totalTimer.useconds() << "% of total)" << std::endl;
  // std::cout << "blockSum          = " << blockSumTimer.Elapsed()          << " (" << 100. * blockSumTimer.useconds()/totalTimer.useconds()          << "% of total)" << std::endl;
  // std::cout << "copy              = " << copyTimer.Elapsed()              << " (" << 100. * copyTimer.useconds()/totalTimer.useconds()              << "% of total)" << std::endl;
}

template<class vobj,class CComplex>
inline void blockNormaliseOriginal(Lattice<CComplex> &ip,Lattice<vobj> &fineX)
{
  // GridStopWatch totalTimer;
  // GridStopWatch preambleTimer;
  // GridStopWatch blockInnerProductTimer;
  // GridStopWatch powTimer;
  // GridStopWatch blockZAXPYTimer;

  // totalTimer.Start();
  // preambleTimer.Start();
  GridBase *coarse = ip._grid;
  Lattice<vobj> zz(fineX._grid); zz=zero; zz.checkerboard=fineX.checkerboard;
  // preambleTimer.Stop();

  // blockInnerProductTimer.Start();
  blockInnerProductOriginal(ip,fineX,fineX);
  // blockInnerProductTimer.Stop();

  // powTimer.Start();
  ip = pow(ip,-0.5);
  // powTimer.Stop();

  // blockZAXPYTimer.Start();
  blockZAXPY(fineX,ip,fineX,zz);
  // blockZAXPYTimer.Stop();
  // totalTimer.Stop();

  // std::cout << "total             = "   << totalTimer.Elapsed()             << " (" << 100. * totalTimer.useconds()/totalTimer.useconds()             << "% of total)" << std::endl;
  // std::cout << "preamble          = "   << preambleTimer.Elapsed()          << " (" << 100. * preambleTimer.useconds()/totalTimer.useconds()          << "% of total)" << std::endl;
  // std::cout << "blockInnerProduct = "   << blockInnerProductTimer.Elapsed() << " (" << 100. * blockInnerProductTimer.useconds()/totalTimer.useconds() << "% of total)" << std::endl;
  // std::cout << "pow               = "   << powTimer.Elapsed()               << " (" << 100. * powTimer.useconds()/totalTimer.useconds()               << "% of total)" << std::endl;
  // std::cout << "blockZAXPY        = " << blockZAXPYTimer.Elapsed()          << " (" << 100. * blockZAXPYTimer.useconds()/totalTimer.useconds()        << "% of total)" << std::endl;
}

template<class vobj,class CComplex>
inline void blockOrthogonaliseOriginal(Lattice<CComplex> &ip,std::vector<Lattice<vobj> > &Basis)
{
  // GridStopWatch totalTimer;
  // GridStopWatch preambleTimer;
  // GridStopWatch kernelTimer;
  // GridStopWatch blockInnerProductTimer;
  // GridStopWatch minusTimer;
  // GridStopWatch blockZAXPYTimer;
  // GridStopWatch blockNormaliseTimer;

  // totalTimer.Start();
  // preambleTimer.Start();
  GridBase *coarse = ip._grid;
  GridBase *fine   = Basis[0]._grid;

  int       nbasis = Basis.size() ;
  int  _ndimension = coarse->_ndimension;

  // checks
  subdivides(coarse,fine);
  for(int i=0;i<nbasis;i++){
    conformable(Basis[i]._grid,fine);
  }
  // preambleTimer.Stop();

  int iters = 0;
  // kernelTimer.Start();
  for(int v=0;v<nbasis;v++) {
    for(int u=0;u<v;u++) {
      //Inner product & remove component
      // blockInnerProductTimer.Start();
      blockInnerProductOriginal(ip,Basis[u],Basis[v]);
      // blockInnerProductTimer.Stop();
      // minusTimer.Start();
      ip = -ip;
      // minusTimer.Stop();
      // blockZAXPYTimer.Start();
      blockZAXPY<vobj,CComplex> (Basis[v],ip,Basis[u],Basis[v]);
      // blockZAXPYTimer.Stop();
      iters++;
    }
    // blockNormaliseTimer.Start();
    blockNormaliseOriginal(ip,Basis[v]);
    // blockNormaliseTimer.Stop();
  }
  // kernelTimer.Stop();
  // totalTimer.Stop();

  // std::cout << "total             = " << totalTimer.Elapsed()             << " (" << 100. * totalTimer.useconds()/totalTimer.useconds()              << "% of total)"  << std::endl;
  // std::cout << "preamble          = " << preambleTimer.Elapsed()          << " (" << 100. * preambleTimer.useconds()/totalTimer.useconds()           << "% of total)"  << std::endl;
  // std::cout << "kernel            = " << kernelTimer.Elapsed()            << " (" << 100. * kernelTimer.useconds()/totalTimer.useconds()             << "% of total)"  << std::endl;
  // std::cout << "blockInnerProduct = " << blockInnerProductTimer.Elapsed() << " (" << 100. * blockInnerProductTimer.useconds()/kernelTimer.useconds() << "% of kernel)" << std::endl;
  // std::cout << "minus             = " << minusTimer.Elapsed()             << " (" << 100. * minusTimer.useconds()/kernelTimer.useconds()             << "% of kernel)" << std::endl;
  // std::cout << "blockZAXPY        = " << blockZAXPYTimer.Elapsed()        << " (" << 100. * blockZAXPYTimer.useconds()/kernelTimer.useconds()        << "% of kernel)" << std::endl;
  // std::cout << "blockNormalise    = " << blockNormaliseTimer.Elapsed()    << " (" << 100. * blockNormaliseTimer.useconds()/kernelTimer.useconds()    << "% of kernel)" << std::endl;
  std::cout << "iters             = " << iters << std::endl;
}

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

  BenchmarkFunction(blockOrthogonaliseOriginal, flop, byte, nIter, InnerProdOriginal, BasisOriginal);
  BenchmarkFunction(blockOrthogonalise,         flop, byte, nIter, InnerProdNew2Args, BasisNew2Args);
  BenchmarkFunction(blockOrthogonalise,         flop, byte, nIter, InnerProdNew3Args, BasisNew3Args, lookUpTable);

  for (auto i = 0; i < BasisOriginal.size(); ++i) printDeviationFromReference(BasisOriginal[i], BasisNew2Args[i]);
  for (auto i = 0; i < BasisOriginal.size(); ++i) printDeviationFromReference(BasisOriginal[i], BasisNew3Args[i]);

  BenchmarkFunction(blockInnerProductOriginal, flop, byte, nIter, InnerProdOriginal, FineXOriginal, FineYOriginal);
  BenchmarkFunction(blockInnerProduct,         flop, byte, nIter, InnerProdNew3Args, FineXOriginal, FineYOriginal, lookUpTable);

  printDeviationFromReference(InnerProdOriginal, InnerProdNew3Args);

  BenchmarkFunction(blockNormaliseOriginal, flop, byte, nIter, InnerProdOriginal, FineXOriginal);
  BenchmarkFunction(blockNormalise,         flop, byte, nIter, InnerProdNew3Args, FineXNew, lookUpTable);

  printDeviationFromReference(FineXOriginal, FineXNew);

  if (doPerfProfiling) {
    PerfProfileFunction(blockOrthogonaliseOriginal, nIter, InnerProdOriginal, BasisOriginal);
    PerfProfileFunction(blockOrthogonalise,         nIter, InnerProdNew2Args, BasisNew2Args);
    PerfProfileFunction(blockOrthogonalise,         nIter, InnerProdNew3Args, BasisNew3Args, lookUpTable);

    PerfProfileFunction(blockInnerProductOriginal, nIter, InnerProdOriginal, FineXOriginal, FineYOriginal);
    PerfProfileFunction(blockInnerProduct,         nIter, InnerProdNew3Args, FineXOriginal, FineYOriginal, lookUpTable);

    PerfProfileFunction(blockNormaliseOriginal, nIter, InnerProdOriginal, FineXOriginal);
    PerfProfileFunction(blockNormalise,         nIter, InnerProdNew3Args, FineXNew, lookUpTable);
  }

  Grid_finalize();
}
