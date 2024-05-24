/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#ifndef SWEC_HPP
#define SWEC_HPP 1

#include "core.hpp"
#include "platform.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "timeStepper.hpp"
#include "linAlg.hpp"

#define DSWEC LIBP_DIR"/solvers/SWEConvergence/"

using namespace libp;

class SWECSettings_t: public settings_t {
public:
  SWECSettings_t(comm_t _comm);
  void report();
  void parseFromFile(platformSettings_t& platformSettings,
                     meshSettings_t& meshSettings,
                     const std::string filename,
                     const std::string degree,
                     const std::string meshfile,
                     const std::string idx);
};

class SWEC_t: public solver_t {
public:
  mesh_t mesh;
  mesh_t meshPatch;

  int Nfields;
  int Ngrads;
  int cubature;
  int curvilinear;

  timeStepper_t timeStepper;

  ogs::halo_t fieldTraceHalo;
  ogs::halo_t gradTraceHalo;
  ogs::halo_t muTraceHalo;

  memory<dfloat> q;
  deviceMemory<dfloat> o_q;
  deviceMemory<dfloat> o_maxSpeed;

  memory<dfloat> gradq;
  deviceMemory<dfloat> o_gradq;

  memory<dfloat> mu;
  deviceMemory<dfloat> o_mu;

  deviceMemory<dfloat> o_Mq;

  kernel_t volumeKernel;
  kernel_t surfaceKernel;
  kernel_t cubatureVolumeKernel;
  kernel_t cubatureSurfaceKernel;

  kernel_t gradVolumeKernel;
  kernel_t gradSurfaceKernel;

  kernel_t initialConditionKernel;
  kernel_t maxWaveSpeedKernel;

  kernel_t viscosityKernel;
  kernel_t viscositySmoothKernel;

  SWEC_t() = default;
  SWEC_t(platform_t &_platform, mesh_t &_mesh,
              SWECSettings_t& _settings) {
    Setup(_platform, _mesh, _settings);
  }

  //setup
  void Setup(platform_t& _platform, mesh_t& _mesh,
             SWECSettings_t& _settings);

  void Run();

  void Report(dfloat time, int tstep);

  void PlotFields(memory<dfloat> Q, const std::string fileName);

  void rhsf(deviceMemory<dfloat>& o_q, deviceMemory<dfloat>& o_rhs, const dfloat time);

  dfloat MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T);
};

#endif
