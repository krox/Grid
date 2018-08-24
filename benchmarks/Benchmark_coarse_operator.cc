/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_coarse_operator.cc

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
#include <Benchmark_helpers.h>

using namespace Grid;
using namespace Grid::QCD;
using namespace Grid::BenchmarkHelpers;

// Enable control of nbasis from the compiler command line
// NOTE to self: Copy the value of CXXFLAGS from the makefile and call make as follows:
//   make CXXFLAGS="-DNBASIS=24 VALUE_OF_CXXFLAGS_IN_MAKEFILE" Benchmark_coarse_operator
#ifndef NBASIS
#define NBASIS 32
#endif

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);

  std::vector<int> seeds({1, 2, 3, 4});

  // clang-format off
  const int        nBasis    = NBASIS; static_assert((nBasis & 0x1) == 0, "");
  const int        nb        = nBasis / 2;
  int              nIter     = readFromCommandLineInt(&argc, &argv, "--niter", 10);
  // NB: blocksize not needed here, since we can use the --grid switch for the coarse grid directly
  // clang-format on

  typedef CoarsenedMatrix<vSpinColourVector, vTComplex, nBasis> CoarseDiracMatrix;
  typedef CoarseDiracMatrix::CoarseVector                       CoarseVector;

  GridCartesian *CGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());

  CGrid->show_decomposition();

  GridParallelRNG CPRNG(CGrid);
  CPRNG.SeedFixedIntegers(seeds);

  CoarseDiracMatrix CoarseMatrix(*CGrid);
  for(auto &elem : CoarseMatrix.A) random(CPRNG, elem);

  CoarseVector CoarseVecIn(CGrid);
  CoarseVector CoarseVecOut(CGrid);

  random(CPRNG, CoarseVecIn);
  random(CPRNG, CoarseVecOut);

  auto nStencil      = CoarseMatrix.geom.npoint;
  auto nAccum        = nStencil;
  auto CSiteVecElems = nBasis;

  auto CVolume = std::accumulate(CGrid->_fdimensions.begin(), CGrid->_fdimensions.end(), 1, std::multiplies<double>());

  double flopM      = (nStencil * (8 * CSiteVecElems * CSiteVecElems - 2 * CSiteVecElems) + nAccum * 2 * CSiteVecElems) * CVolume + 8 * CSiteVecElems * CVolume;
  double byteM      = (nStencil * (CSiteVecElems * CSiteVecElems + CSiteVecElems) + CSiteVecElems) * CVolume * sizeof(Complex) + CSiteVecElems * CVolume * sizeof(Complex);
  double footprintM = (2 * CSiteVecElems + nStencil * CSiteVecElems * CSiteVecElems) * CVolume * sizeof(Complex);

  double flopMdir      = (8 * CSiteVecElems * CSiteVecElems - 2 * CSiteVecElems) * CVolume;
  double byteMdir      = (CSiteVecElems * CSiteVecElems + 2 * CSiteVecElems) * CVolume * sizeof(Complex);
  double footprintMdir = (2 * CSiteVecElems + nStencil * CSiteVecElems * CSiteVecElems) * CVolume * sizeof(Complex);

  KernelPerf coarseM("coarseM", flopM, byteM, footprintM);
  KernelPerf coarseMdag("coarseMdag", flopM, byteM, footprintM); // TODO: with the current implementation of Mdag, this line is not correct
  KernelPerf coarseMdir("coarseMdir", flopMdir, byteMdir, footprintMdir);
  KernelPerf coarseMdiag("coarseMdiag", flopMdir, byteMdir, footprintMdir);

  {
    double start = usecond();
    for(int i = 0; i < nIter; ++i) CoarseMatrix.M(CoarseVecIn, CoarseVecOut);
    double stop = usecond();
    coarseM.reportPerformance(stop - start, nIter);
  }

  {
    double start = usecond();
    for(int i = 0; i < nIter; ++i) CoarseMatrix.Mdag(CoarseVecIn, CoarseVecOut);
    double stop = usecond();
    coarseMdag.reportPerformance(stop - start, nIter);
  }

  {
    double start = usecond();
    for(int i = 0; i < nIter; ++i) CoarseMatrix.Mdir(CoarseVecIn, CoarseVecOut, 2, 1);
    double stop = usecond();
    coarseMdir.reportPerformance(stop - start, nIter);
  }

  {
    double start = usecond();
    for(int i = 0; i < nIter; ++i) CoarseMatrix.Mdiag(CoarseVecIn, CoarseVecOut);
    double stop = usecond();
    coarseMdiag.reportPerformance(stop - start, nIter);
  }

  Grid_finalize();
}
