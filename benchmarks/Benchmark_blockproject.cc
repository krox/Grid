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
#define NBASIS 32
#endif

template<class vobj,class CComplex,int nbasis>
inline void blockProjectOriginal(Lattice<iVector<CComplex,nbasis > > &coarseData,
                                 const             Lattice<vobj>   &fineData,
                                 const std::vector<Lattice<vobj> > &Basis)
{
  GridBase * fine  = fineData._grid;
  GridBase * coarse= coarseData._grid;
  int  _ndimension = coarse->_ndimension;

  // checks
  assert( nbasis == Basis.size() );
  subdivides(coarse,fine);
  for(int i=0;i<nbasis;i++){
    conformable(Basis[i],fineData);
  }

  std::vector<int>  block_r      (_ndimension);

  for(int d=0 ; d<_ndimension;d++){
    block_r[d] = fine->_rdimensions[d] / coarse->_rdimensions[d];
    assert(block_r[d]*coarse->_rdimensions[d] == fine->_rdimensions[d]);
  }

  coarseData=zero;

  // Loop over coars parallel, and then loop over fine associated with coarse.
  parallel_for(int sf=0;sf<fine->oSites();sf++){

    int sc;
    std::vector<int> coor_c(_ndimension);
    std::vector<int> coor_f(_ndimension);
    Lexicographic::CoorFromIndex(coor_f,sf,fine->_rdimensions);
    for(int d=0;d<_ndimension;d++) coor_c[d]=coor_f[d]/block_r[d];
    Lexicographic::IndexFromCoor(coor_c,sc,coarse->_rdimensions);

PARALLEL_CRITICAL
    for(int i=0;i<nbasis;i++) {

      coarseData._odata[sc](i)=coarseData._odata[sc](i)
	+ innerProduct(Basis[i]._odata[sf],fineData._odata[sf]);

    }
  }
  return;
}

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);

  std::vector<int> seeds({1, 2, 3, 4});

  // clang-format off
  const int        nBasis    = NBASIS; static_assert((nBasis & 0x1) == 0, "");
  const int        nb        = nBasis / 2;
  int              nIter     = readFromCommandLineInt(&argc, &argv, "--niter", 10);
  std::vector<int> blockSize = readFromCommandLineIntVec(&argc, &argv, "--blocksize", std::vector<int>({4, 4, 4, 4}));
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
  performChiralDoubling(Aggs.subspace);

  LatticeFermion FineVec(FGrid);
  CoarseVector   CoarseVec(CGrid);
  CoarseVector   CoarseVecOriginal(CGrid);
  CoarseVector   CoarseVecDiff(CGrid);
  CoarseVec         = zero;
  CoarseVecOriginal = zero;

  random(FPRNG, FineVec);

  auto FSiteVecElems = Nc * Ns;
  auto CSiteVecElems = nBasis;

  auto FVolume = std::accumulate(FGrid->_fdimensions.begin(), FGrid->_fdimensions.end(), 1, std::multiplies<double>());
  auto CVolume = std::accumulate(CGrid->_fdimensions.begin(), CGrid->_fdimensions.end(), 1, std::multiplies<double>());

  double flop      = 8 * FSiteVecElems * nBasis * FVolume;
  double byte      = (2 * 1 + 2 * FSiteVecElems) * nBasis * FVolume * sizeof(Complex);

  BenchmarkFunction(blockProjectOriginal, flop, byte, nIter, CoarseVecOriginal, FineVec, Aggs.subspace);
  BenchmarkFunction(blockProject,         flop, byte, nIter, CoarseVec,         FineVec, Aggs.subspace);

  CoarseVecDiff = CoarseVecOriginal - CoarseVec;
  auto absDev   = norm2(CoarseVecDiff);
  auto relDev   = absDev / norm2(CoarseVecOriginal);
  std::cout << GridLogMessage << "absolute deviation = " << absDev << " | relative deviation = " << relDev << std::endl;

  Grid_finalize();
}
