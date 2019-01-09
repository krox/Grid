/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./benchmarks/Benchmark_devel_new_coarsening.cc

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
//   make CXXFLAGS="-DNBASIS=24 VALUE_OF_CXXFLAGS_IN_MAKEFILE" Benchmark_devel_new_coarsening
#ifndef NBASIS
#define NBASIS 4
#endif

///////////////////////////////////////////////////////////////////////////////
//                       Definition of helper functions                      //
///////////////////////////////////////////////////////////////////////////////

template<class OldCoarsenedMatrix, class NewCoarsenedMatrix>
void
convertOldLayoutToNewLayout(typename OldCoarsenedMatrix::FermionField const &oldFormatIn, typename NewCoarsenedMatrix::FermionField &newFormatOut) {
  conformable(oldFormatIn._grid, newFormatOut._grid);
  GridBase * grid = oldFormatIn._grid;
  int Nbasis = OldCoarsenedMatrix::Nbasis;
  int Nb = Nbasis/2;
  parallel_for(int ss = 0; ss < grid->oSites(); ss++) {
    for(int n=0; n<Nbasis; n++) {
      int spinIndex = n / Nb;
      int colourIndex = n % Nb;
      newFormatOut._odata[ss]()(spinIndex)(colourIndex) = oldFormatIn._odata[ss](n)()()();
    }
  }
}

template<class OldCoarsenedMatrix, class NewCoarsenedMatrix>
void
convertOldLayoutToNewLayout(typename OldCoarsenedMatrix::LinkField const &oldFormatIn, typename NewCoarsenedMatrix::LinkField &newFormatOut) {
  conformable(oldFormatIn._grid, newFormatOut._grid);
  GridBase * grid = oldFormatIn._grid;
  int Nbasis = OldCoarsenedMatrix::Nbasis;
  int Nb = Nbasis/2;
  parallel_for(int ss = 0; ss < grid->oSites(); ss++) {
    for(int n1=0; n1<Nbasis; n1++) {
      for(int n2=0; n2<Nbasis; n2++) {
        int spinIndex1 = n1 / Nb;
        int spinIndex2 = n2 / Nb;
        int colourIndex1 = n1 % Nb;
        int colourIndex2 = n2 % Nb;
        newFormatOut._odata[ss]()(spinIndex1, spinIndex2)(colourIndex1, colourIndex2) = oldFormatIn._odata[ss](n1, n2)()()();
      }
    }
  }
}

template<class Field>
void printIndexRankInfo(std::string const &name) {
  typedef typename getVectorType<Field>::type vobj; // gives us the the type of the site object if it's a Lattice type

  std::cout << name << ": ColourN       = " << indexRank<ColourIndex,vobj>()  << std::endl;
  std::cout << name << ": ColourScalar  = " <<  isScalar<ColourIndex,vobj>()  << std::endl;
  std::cout << name << ": ColourVector  = " <<  isVector<ColourIndex,vobj>()  << std::endl;
  std::cout << name << ": ColourMatrix  = " <<  isMatrix<ColourIndex,vobj>()  << std::endl;

  std::cout << name << ": SpinN         = " << indexRank<SpinIndex,vobj>()    << std::endl;
  std::cout << name << ": SpinScalar    = " <<  isScalar<SpinIndex,vobj>()    << std::endl;
  std::cout << name << ": SpinVector    = " <<  isVector<SpinIndex,vobj>()    << std::endl;
  std::cout << name << ": SpinMatrix    = " <<  isMatrix<SpinIndex,vobj>()    << std::endl;

  std::cout << name << ": LorentzN      = " << indexRank<LorentzIndex,vobj>() << std::endl;
  std::cout << name << ": LorentzScalar = " <<  isScalar<LorentzIndex,vobj>() << std::endl;
  std::cout << name << ": LorentzVector = " <<  isVector<LorentzIndex,vobj>() << std::endl;
  std::cout << name << ": LorentzMatrix = " <<  isMatrix<LorentzIndex,vobj>() << std::endl;
}

int main(int argc, char **argv) {

  Grid_init(&argc, &argv);

  /////////////////////////////////////////////////////////////////////////////
  //                          Read from command line                         //
  /////////////////////////////////////////////////////////////////////////////

  const int nBasis           = NBASIS; static_assert((nBasis & 0x1) == 0, "");
  const int nB               = nBasis / 2;
  const int nIter            = readFromCommandLineInt(&argc, &argv, "--niter", 10);
  std::vector<int> blockSize = readFromCommandLineIntVec(&argc, &argv, "--blocksize", std::vector<int>({4, 4, 4, 4}));

  /////////////////////////////////////////////////////////////////////////////
  //                              General setup                              //
  /////////////////////////////////////////////////////////////////////////////

  std::vector<int> clatt = calcCoarseLattSize(GridDefaultLatt(), blockSize);

  GridCartesian *        FGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridCartesian *        CGrid   = SpaceTimeGrid::makeFourDimGrid(clatt, GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid);

  FGrid->show_decomposition();
  CGrid->show_decomposition();

  GridParallelRNG FPRNG(FGrid);
  GridParallelRNG CPRNG(CGrid);

  std::vector<int> seeds({1, 2, 3, 4});

  FPRNG.SeedFixedIntegers(seeds);
  CPRNG.SeedFixedIntegers(seeds);

  /////////////////////////////////////////////////////////////////////////////
  //                    Setup of Dirac Matrix and Operator                   //
  /////////////////////////////////////////////////////////////////////////////

  RealD mass=0.5;
  LatticeGaugeField Umu(FGrid); SU3::HotConfiguration(FPRNG,Umu);
  WilsonFermionR Dw(Umu,*FGrid,*FrbGrid,mass);
  MdagMLinearOperator<WilsonFermionR, LatticeFermion> LinOp(Dw);

  /////////////////////////////////////////////////////////////////////////////
  //                             Type definitions                            //
  /////////////////////////////////////////////////////////////////////////////

  typedef OriginalCoarseningPolicy<LatticeFermion, vComplex, nBasis> OldCoarseningPolicy;
  typedef TwoSpinCoarseningPolicy<LatticeFermion, vComplex, nB>      NewCoarseningPolicy;

  typedef CoarsenedMatrixUsingPolicies<OldCoarseningPolicy> OldCoarsenedMatrix;
  typedef CoarsenedMatrixUsingPolicies<NewCoarseningPolicy> NewCoarsenedMatrix;

  typedef OldCoarsenedMatrix::FermionField OldCoarseVector;
  typedef NewCoarsenedMatrix::FermionField NewCoarseVector;

  typedef AggregationUsingPolicies<OldCoarseningPolicy> OldAggregation;
  typedef AggregationUsingPolicies<NewCoarseningPolicy> NewAggregation;

  /////////////////////////////////////////////////////////////////////////////
  //                           Setup of Aggregation                          //
  /////////////////////////////////////////////////////////////////////////////

  OldAggregation oldAggregates(CGrid, FGrid, 0);
  NewAggregation newAggregates(CGrid, FGrid, 0);

  oldAggregates.CreateSubspaceRandom(FPRNG);
  for(int i=0; i<newAggregates._subspace.size(); i++) newAggregates._subspace[i] = oldAggregates._subspace[i];
  performChiralDoubling(oldAggregates._subspace);

  // /////////////////////////////////////////////////////////////////////////////
  // //         Artificial setup of CoarsenedMatrix for testing purposes        //
  // /////////////////////////////////////////////////////////////////////////////

  // OldCoarsenedMatrix oldCoarsenedMatrix(*CGrid);
  // NewCoarsenedMatrix newCoarsenedMatrix(*CGrid);

  // for(auto &elem : oldCoarsenedMatrix._Y) random(CPRNG, elem);
  // for(int i=0; i<oldCoarsenedMatrix._geom.npoint; i++) {
  //   convertOldLayoutToNewLayout<OldCoarsenedMatrix, NewCoarsenedMatrix>(oldCoarsenedMatrix._Y[i], newCoarsenedMatrix._Y[i]);
  // }

  /////////////////////////////////////////////////////////////////////////////
  //                     Correct setup of CoarsenedMatrix                    //
  /////////////////////////////////////////////////////////////////////////////

  OldCoarsenedMatrix oldCoarsenedMatrix(*CGrid);
  NewCoarsenedMatrix newCoarsenedMatrix(*CGrid);

  oldCoarsenedMatrix.CoarsenOperator(FGrid, LinOp, oldAggregates);
  newCoarsenedMatrix.CoarsenOperator(FGrid, LinOp, newAggregates);

  /////////////////////////////////////////////////////////////////////////////
  //                     Setup of vectors for Benchmarks                     //
  /////////////////////////////////////////////////////////////////////////////

  // NOTE: All "In" vectors set to random, all "Out" vectors set to zero

  OldCoarseVector oldCoarseVectorIn(CGrid);   random(CPRNG, oldCoarseVectorIn);
  NewCoarseVector newCoarseVectorIn(CGrid);   convertOldLayoutToNewLayout<OldCoarsenedMatrix, NewCoarsenedMatrix>(oldCoarseVectorIn, newCoarseVectorIn);
  OldCoarseVector oldCoarseVectorOut(CGrid);  oldCoarseVectorOut  = zero;
  NewCoarseVector newCoarseVectorOut(CGrid);  newCoarseVectorOut  = zero;
  NewCoarseVector tmpCoarseVector4Dev(CGrid); tmpCoarseVector4Dev = zero;

  LatticeFermion oldFineVectorIn(FGrid);  random(FPRNG, oldFineVectorIn);
  LatticeFermion newFineVectorIn(FGrid);  newFineVectorIn = oldFineVectorIn; // layout conversion not needed on fine grid
  LatticeFermion oldFineVectorOut(FGrid); oldFineVectorOut = zero;
  LatticeFermion newFineVectorOut(FGrid); newFineVectorOut = zero;

  /////////////////////////////////////////////////////////////////////////////
  //            Calculate performance figures for instrumentation            //
  /////////////////////////////////////////////////////////////////////////////

  double nStencil      = oldCoarsenedMatrix._geom.npoint;
  double nAccum        = nStencil;
  double FSiteElems    = Nc * Ns;
  double oldCSiteElems = OldCoarsenedMatrix::Nbasis;
  double newCSiteElems = NewCoarsenedMatrix::Ncs * NewCoarsenedMatrix::Nbasis;

  double FVolume = std::accumulate(FGrid->_fdimensions.begin(), FGrid->_fdimensions.end(), 1, std::multiplies<double>());
  double CVolume = std::accumulate(CGrid->_fdimensions.begin(), CGrid->_fdimensions.end(), 1, std::multiplies<double>());

  double flopOldM      = 1. * ((nStencil * (8 * oldCSiteElems * oldCSiteElems - 2 * oldCSiteElems) + nAccum * 2 * oldCSiteElems) + 8 * oldCSiteElems) * CVolume;
  double byteOldM      = 1. * ((nStencil * (oldCSiteElems * oldCSiteElems + oldCSiteElems) + oldCSiteElems) + oldCSiteElems) * CVolume * sizeof(Complex);
  double footprintOldM = 1. * (2 * oldCSiteElems + nStencil * oldCSiteElems * oldCSiteElems) * CVolume * sizeof(Complex);

  double flopNewM      = 1. * ((nStencil * (8 * newCSiteElems * newCSiteElems - 2 * newCSiteElems) + nAccum * 2 * newCSiteElems) + 8 * newCSiteElems) * CVolume;
  double byteNewM      = 1. * ((nStencil * (newCSiteElems * newCSiteElems + newCSiteElems) + newCSiteElems) + newCSiteElems) * CVolume * sizeof(Complex);
  double footprintNewM = 1. * (2 * newCSiteElems + nStencil * newCSiteElems * newCSiteElems) * CVolume * sizeof(Complex);

  double flopOldMdir      = 1. * (8 * oldCSiteElems * oldCSiteElems - 2 * oldCSiteElems) * CVolume;
  double byteOldMdir      = 1. * (oldCSiteElems * oldCSiteElems + 2 * oldCSiteElems) * CVolume * sizeof(Complex);
  double footprintOldMdir = 1. * (2 * oldCSiteElems + nStencil * oldCSiteElems * oldCSiteElems) * CVolume * sizeof(Complex);

  double flopNewMdir      = 1. * (8 * newCSiteElems * newCSiteElems - 2 * newCSiteElems) * CVolume;
  double byteNewMdir      = 1. * (newCSiteElems * newCSiteElems + 2 * newCSiteElems) * CVolume * sizeof(Complex);
  double footprintNewMdir = 1. * (2 * newCSiteElems + nStencil * newCSiteElems * newCSiteElems) * CVolume * sizeof(Complex);

  double flopOldMdiag      = flopOldMdir;
  double byteOldMdiag      = byteOldMdir;
  double footprintOldMdiag = footprintOldMdir;

  double flopNewMdiag      = flopNewMdir;
  double byteNewMdiag      = byteNewMdir;
  double footprintNewMdiag = footprintNewMdir;

  double flopOldProject      = 1. * (8 * FSiteElems) * nBasis * FVolume;
  double byteOldProject      = 1. * (2 * 1 + 2 * FSiteElems) * nBasis * FVolume * sizeof(Complex);
  double footprintOldProject = -1.;

  double flopNewProject      = 1. * (8 * FSiteElems) * nB * FVolume;
  double byteNewProject      = 1. * (2 * 1 + 2 * FSiteElems) * nB * FVolume * sizeof(Complex);
  double footprintNewProject = -1.;

  double flopOldPromote      = 1. * (8 * (nBasis - 1) + 6) * FSiteElems * FVolume;
  double byteOldPromote      = 1. * ((1 * 1 + 3 * FSiteElems) * (nBasis - 1) + (1 * 1 + 2 * FSiteElems) * 1) * FVolume * sizeof(Complex);
  double footprintOldPromote = 1. * (-1);

  double flopNewPromote      = 1. * (8 * (nB - 1) + 6) * FSiteElems * FVolume;
  double byteNewPromote      = 1. * ((1 * 1 + 3 * FSiteElems) * (nB - 1) + (1 * 1 + 2 * FSiteElems) * 1) * FVolume * sizeof(Complex);
  double footprintNewPromote = -1.;

  std::cout << GridLogMessage << "***************************************************************************" << std::endl;
  std::cout << GridLogMessage << "Performance figures and result comparison for benchmark coarse M" << std::endl;
  std::cout << GridLogMessage << "***************************************************************************" << std::endl;

  oldCoarseVectorOut = zero;
  newCoarseVectorOut = zero;
  BenchmarkFunction(oldCoarsenedMatrix.M, flopOldM, byteOldM, nIter, oldCoarseVectorIn, oldCoarseVectorOut);
  BenchmarkFunction(newCoarsenedMatrix.M, flopNewM, byteNewM, nIter, newCoarseVectorIn, newCoarseVectorOut);
  convertOldLayoutToNewLayout<OldCoarsenedMatrix, NewCoarsenedMatrix>(oldCoarseVectorOut, tmpCoarseVector4Dev);
  printDeviationFromReference(tmpCoarseVector4Dev, newCoarseVectorOut);

  std::cout << GridLogMessage << "***************************************************************************" << std::endl;
  std::cout << GridLogMessage << "Performance figures and result comparison for benchmark coarse Mdir" << std::endl;
  std::cout << GridLogMessage << "***************************************************************************" << std::endl;

  oldCoarseVectorOut = zero;
  newCoarseVectorOut = zero;
  BenchmarkFunction(oldCoarsenedMatrix.Mdir, flopOldMdir, byteOldMdir, nIter, oldCoarseVectorIn, oldCoarseVectorOut, 2, 1);
  BenchmarkFunction(newCoarsenedMatrix.Mdir, flopNewMdir, byteNewMdir, nIter, newCoarseVectorIn, newCoarseVectorOut, 2, 1);
  convertOldLayoutToNewLayout<OldCoarsenedMatrix, NewCoarsenedMatrix>(oldCoarseVectorOut, tmpCoarseVector4Dev);
  printDeviationFromReference(tmpCoarseVector4Dev, newCoarseVectorOut);

  std::cout << GridLogMessage << "***************************************************************************" << std::endl;
  std::cout << GridLogMessage << "Performance figures and result comparison for benchmark coarse Mdiag" << std::endl;
  std::cout << GridLogMessage << "***************************************************************************" << std::endl;

  oldCoarseVectorOut = zero;
  newCoarseVectorOut = zero;
  BenchmarkFunction(oldCoarsenedMatrix.Mdiag, flopOldMdiag, byteOldMdiag, nIter, oldCoarseVectorIn, oldCoarseVectorOut);
  BenchmarkFunction(newCoarsenedMatrix.Mdiag, flopNewMdiag, byteNewMdiag, nIter, newCoarseVectorIn, newCoarseVectorOut);
  convertOldLayoutToNewLayout<OldCoarsenedMatrix, NewCoarsenedMatrix>(oldCoarseVectorOut, tmpCoarseVector4Dev);
  printDeviationFromReference(tmpCoarseVector4Dev, newCoarseVectorOut);

  std::cout << GridLogMessage << "***************************************************************************" << std::endl;
  std::cout << GridLogMessage << "Performance figures and result comparison for benchmark projectToSubspace" << std::endl;
  std::cout << GridLogMessage << "***************************************************************************" << std::endl;

  oldCoarseVectorOut = zero;
  newCoarseVectorOut = zero;
  BenchmarkFunction(oldAggregates.ProjectToSubspace, flopOldProject, byteOldProject, nIter, oldCoarseVectorOut, oldFineVectorIn);
  BenchmarkFunction(newAggregates.ProjectToSubspace, flopNewProject, byteNewProject, nIter, newCoarseVectorOut, newFineVectorIn);
  convertOldLayoutToNewLayout<OldCoarsenedMatrix, NewCoarsenedMatrix>(oldCoarseVectorOut, tmpCoarseVector4Dev);
  printDeviationFromReference(tmpCoarseVector4Dev, newCoarseVectorOut);

  std::cout << GridLogMessage << "***************************************************************************" << std::endl;
  std::cout << GridLogMessage << "Performance figures and result comparison for benchmark promoteFromSubspace" << std::endl;
  std::cout << GridLogMessage << "***************************************************************************" << std::endl;

  oldFineVectorOut = zero;
  newFineVectorOut = zero;
  BenchmarkFunction(oldAggregates.PromoteFromSubspace, flopOldPromote, byteOldPromote, nIter, oldCoarseVectorIn, oldFineVectorOut);
  BenchmarkFunction(newAggregates.PromoteFromSubspace, flopNewPromote, byteNewPromote, nIter, newCoarseVectorIn, newFineVectorOut);
  printDeviationFromReference(oldFineVectorOut, newFineVectorOut);

  /////////////////////////////////////////////////////////////////////////////
  //                   Print info about types (for testing)                  //
  /////////////////////////////////////////////////////////////////////////////

  // printIndexRankInfo<LatticeFermion>("LatticeFermion");
  // printIndexRankInfo<LatticeHalfFermion>("LatticeHalfFermion");
  // printIndexRankInfo<LatticeGaugeField>("LatticeGaugeField");
  // printIndexRankInfo<LatticeDoubledGaugeField>("LatticeDoubledGaugeField");

  Grid_finalize();
}
