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

namespace Grid {

namespace BenchmarkHelpers {

#define BenchmarkFunction(function, flop, byte, nIter, ...)     \
  do {                                                          \
    GridPerfMonitor perfMonitor(flop, byte);                    \
    perfMonitor.Start();                                        \
    for(int i = 0; i < nIter; ++i) {                            \
      function(__VA_ARGS__);                                    \
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
      expression;                                                 \
    }                                                             \
    perfMonitor.Stop(nIter);                                      \
    std::cout << GridLogPerformance << "Kernel "                  \
              << std::setw(25) << std::right                      \
              << #expression << ": " << perfMonitor << std::endl; \
  } while(0)

class KernelPerf {
public:
  std::string name;
  double flop;
  double byte;
  double intensity;
  double footprint;

  KernelPerf(const std::string &_name, double _flop, double _byte, double _footprint)
    : name(_name)
    , flop(_flop)
    , byte(_byte)
    , intensity(flop/byte)
    , footprint(_footprint)
  {}

  void reportPerformance(double uSeconds, int nIter) {
    double time = uSeconds / nIter / 1e6;
    std::cout << GridLogMessage << std::scientific << "Kernel " << name << ":"
              << " Intensity[F/B] = " << intensity
              << " Time[s] = " << time
              << " Perf[F/s] = " << flop / time
              << " Traffic[B/s] = " << byte / time
              << " Footprint[B] = " << footprint << std::endl;
  }
};

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

} // BenchmarkHelpers
} // Grid

#endif
