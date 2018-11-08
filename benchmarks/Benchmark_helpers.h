/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./benchmarks/Benchmark_helpers.h

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
#ifndef GRID_BENCHMARK_HELPERS_H
#define GRID_BENCHMARK_HELPERS_H

#include "Grid/util/Profiling.h"

namespace Grid {

namespace BenchmarkHelpers {

#define BenchmarkFunction(function, flop, byte, nIter, ...)     \
  do {                                                          \
    GridPerfMonitor perfMonitor(flop, byte);                    \
    perfMonitor.Start();                                        \
    for(int i = 0; i < nIter; ++i) {                            \
      __SSC_START;                                              \
      function(__VA_ARGS__);                                    \
      __SSC_STOP;                                               \
    }                                                           \
    perfMonitor.Stop(nIter);                                    \
    std::cout << GridLogPerformance << "Kernel "                \
              << std::setw(25) << std::right                    \
              << #function << ": " << perfMonitor << std::endl; \
  } while(0)

#define BenchmarkExpression(expression, flop, byte, nIter)        \
  do {                                                            \
    GridPerfMonitor perfMonitor(flop, byte);                      \
    perfMonitor.Start();                                          \
    for(int i = 0; i < nIter; ++i) {                              \
      __SSC_START;                                                \
      expression;                                                 \
      __SSC_STOP;                                                 \
    }                                                             \
    perfMonitor.Stop(nIter);                                      \
    std::cout << GridLogPerformance << "Kernel "                  \
              << std::setw(25) << std::right                      \
              << #expression << ": " << perfMonitor << std::endl; \
  } while(0)

#define PerfProfileFunction(function, nIter, ...)                                         \
  do {                                                                                    \
    std::string kernelName = "kernel."#function;                                          \
    std::string outputFile = kernelName + ".data";                                        \
    System::profile(kernelName, [&]() {                                                   \
      for(int i = 0; i < nIter; ++i) {                                                    \
        function(__VA_ARGS__);                                                            \
      }                                                                                   \
    });                                                                                   \
    std::cout << GridLogMessage << "Generated " << outputFile << std::endl;               \
    std::cout << GridLogMessage << "Use with: perf report -i " << outputFile <<std::endl; \
  } while(0)

#define PerfProfileExpression(expression, nIter)                                          \
  do {                                                                                    \
    std::string kernelName = "kernel."#expression;                                        \
    std::replace(kernelName.begin(), kernelName.end(), ' ', '_');                         \
    std::string outputFile = kernelName + ".data";                                        \
    System::profile(kernelName, [&]() {                                                   \
      for(int i = 0; i < nIter; ++i) {                                                    \
        expression;                                                                       \
      }                                                                                   \
    });                                                                                   \
    std::cout << GridLogMessage << "Generated " << outputFile << std::endl;               \
    std::cout << GridLogMessage << "Use with: perf report -i " << outputFile <<std::endl; \
  } while(0)

int readFromCommandLineInt(int *argc, char ***argv, const std::string &option, int defaultValue) {
  std::string arg;
  int         ret = defaultValue;
  if(GridCmdOptionExists(*argv, *argv + *argc, option)) {
    arg = GridCmdOptionPayload(*argv, *argv + *argc, option);
    GridCmdOptionInt(arg, ret);
  }
  return ret;
}

std::vector<int> readFromCommandLineIntVec(int *argc, char ***argv, const std::string &option, const std::vector<int> &defaultValues) {
  std::string      arg;
  std::vector<int> ret(defaultValues);
  if(GridCmdOptionExists(*argv, *argv + *argc, option)) {
    arg = GridCmdOptionPayload(*argv, *argv + *argc, option);
    GridCmdOptionIntVector(arg, ret);
  }
  return ret;
}

bool readFromCommandLineToggle(int *argc, char ***argv, const std::string &option) {
  std::string arg;
  bool ret = false;
  if(GridCmdOptionExists(*argv, *argv + *argc, option)) {
    ret = true;
  }
  return ret;
}

std::vector<int> calcCoarseLattSize(const std::vector<int> &fineLattSize, const std::vector<int> &blockSize) {
  std::vector<int> ret(fineLattSize);
  for(int d = 0; d < ret.size(); d++) {
    ret[d] /= blockSize[d];
  }
  return ret;
}

template<typename Field> void performChiralDoubling(std::vector<Field> &basisVectors) {
  assert(basisVectors.size() % 2 == 0);
  auto nb = basisVectors.size() / 2;

  for(int n = 0; n < nb; n++) {
    auto tmp1 = basisVectors[n];
    auto tmp2 = tmp1;
    G5C(tmp2, basisVectors[n]);
    axpby(basisVectors[n], 0.5, 0.5, tmp1, tmp2);
    axpby(basisVectors[n + nb], 0.5, -0.5, tmp1, tmp2);
    std::cout << GridLogMessage << "Chirally doubled vector " << n << ". "
              << "norm2(vec[" << n << "]) = " << norm2(basisVectors[n]) << ". "
              << "norm2(vec[" << n + nb << "]) = " << norm2(basisVectors[n + nb]) << std::endl;
  }
}

template<class Field>
void printDeviationFromReference(Field const &reference, Field const& result) {
  conformable(reference._grid, result._grid);
  Field diff(reference._grid);

  diff = reference - result;
  auto absDev   = norm2(diff);
  auto relDev   = absDev / norm2(reference);
  std::cout << GridLogMessage << "absolute deviation = " << absDev << " | relative deviation = " << relDev << std::endl;
}

template<class vobj>
int getSiteElems() {
  int ret = 1;

  typedef typename getVectorType<vobj>::type vobj_site; // gives us the the type of the site object if vobj is a Lattice type

  // TODO: This doesn't work for arbitry tensors. One would need to get the tensor recursion depth and then walk over these and multiply the values

  int _ColourN       = indexRank<ColourIndex,vobj_site>();
  int _ColourScalar  =  isScalar<ColourIndex,vobj_site>();
  int _ColourVector  =  isVector<ColourIndex,vobj_site>();
  int _ColourMatrix  =  isMatrix<ColourIndex,vobj_site>();

  int _SpinN         = indexRank<SpinIndex,vobj_site>();
  int _SpinScalar    =  isScalar<SpinIndex,vobj_site>();
  int _SpinVector    =  isVector<SpinIndex,vobj_site>();
  int _SpinMatrix    =  isMatrix<SpinIndex,vobj_site>();

  int _LorentzN      = indexRank<LorentzIndex,vobj_site>();
  int _LorentzScalar =  isScalar<LorentzIndex,vobj_site>();
  int _LorentzVector =  isVector<LorentzIndex,vobj_site>();
  int _LorentzMatrix =  isMatrix<LorentzIndex,vobj_site>();

  if (_ColourMatrix) ret *= _ColourN * _ColourN;
  else               ret *= _ColourN;

  if (_SpinMatrix) ret *= _SpinN * _SpinN;
  else             ret *= _SpinN;

  if (_LorentzMatrix) ret *= _LorentzN * _LorentzN;
  else                ret *= _LorentzN;

  return ret;
}

} // BenchmarkHelpers
} // Grid

#endif
