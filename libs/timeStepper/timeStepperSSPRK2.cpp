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

#include "core.hpp"
#include "timeStepper.hpp"

namespace libp {

namespace TimeStepper {

ssprk2::ssprk2(dlong Nelements, dlong NhaloElements,
               int Np, int Nfields,
               platform_t& _platform, comm_t _comm):
  timeStepperBase_t(Nelements, NhaloElements, Np, Nfields,
                    _platform, _comm) {

  Nrk = 2;

  o_resq = platform.malloc<dfloat>(N);
  o_rhsq = platform.malloc<dfloat>(N);

  o_saveq = platform.malloc<dfloat>(N);
  o_q1 = platform.malloc<dfloat>(N);

  properties_t kernelInfo = platform.props(); //copy base occa properties from solver

  const int blocksize=256;


  kernelInfo["defines/" "p_blockSize"] = blocksize;

  updateKernel = platform.buildKernel(TIMESTEPPER_DIR "/okl/"
                                    "timeStepperSSPRK2.okl",
                                    "ssprk2Update",
                                    kernelInfo);

  dfloat _rka[2] = {0.0,1.0};
  dfloat _rkb[2] = {1,0};
  // added one more for advanced time step
  dfloat _rkc[2] = {0.0,1.0};

  rka.malloc(Nrk);
  rkb.malloc(Nrk);
  rkc.malloc(Nrk);
  rka.copyFrom(_rka);
  rkb.copyFrom(_rkb);
  rkc.copyFrom(_rkc);
}

void ssprk2::Run(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat start, dfloat end) {

  dfloat time = start;

  solver.Report(time,0);

  dfloat outputInterval=0.0;
  solver.settings.getSetting("OUTPUT INTERVAL", outputInterval);

  dfloat outputTime = time + outputInterval;

  int tstep=0;
  dfloat stepdt;
  while (time < end) {

    if (time<outputTime && time+dt>=outputTime) {

      //save current state
      o_saveq.copyFrom(o_q, N);

      stepdt = outputTime-time;

      //take small time step
      Step(solver, o_q, time, stepdt);

      //report state
      solver.Report(outputTime,tstep);

      //restore previous state
      o_q.copyFrom(o_saveq, N);

      outputTime += outputInterval;
    }

    //check for final timestep
    if (time+dt > end){
      stepdt = end-time;
    } else {
      stepdt = dt;
    }

    Step(solver, o_q, time, stepdt);
    time += stepdt;
    tstep++;
  }
}

void ssprk2::Step(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat time, dfloat _dt) {

  o_q1.copyFrom(o_q, N);
  for(int rk=0;rk<2;rk++){
  // Low storage explicit Runge Kutta (5 stages, 4th order)
    dfloat currentTime = time + (dfloat)rk*_dt;

    //evaluate ODE rhs = f(q,t)
    solver.rhsf(o_q, o_rhsq, currentTime);
    // update solution using Runge-Kutta
    updateKernel(N, _dt, rka[rk],o_rhsq, o_q, o_q1);
  }
}


/**************************************************/
/* PML version                                    */
/**************************************************/

/*
ssprk2_pml::ssprk2_pml(dlong _Nelements, dlong _NpmlElements, dlong _NhaloElements,
                      int _Np, int _Nfields, int _Npmlfields,
                      platform_t& _platform, comm_t _comm):
  ssprk2(_Nelements, _NhaloElements, _Np, _Nfields, _platform, _comm),
  Npml(_Npmlfields*_Np*_NpmlElements) {

  if (Npml) {
    memory<dfloat> pmlq(Npml,0.0);
    o_pmlq = platform.malloc<dfloat>(pmlq);

    o_respmlq = platform.malloc<dfloat>(Npml);
    o_rhspmlq = platform.malloc<dfloat>(Npml);
  }
}

void ssprk2_pml::Step(solver_t& solver, deviceMemory<dfloat> &o_q, dfloat time, dfloat _dt) {

  // Low storage explicit Runge Kutta (5 stages, 4th order)
  for(int rk=0;rk<Nrk;++rk){

    dfloat currentTime = time + rkc[rk]*_dt;

    //evaluate ODE rhs = f(q,t)
    solver.rhsf_pml(o_q, o_pmlq, o_rhsq, o_rhspmlq, currentTime);

    // update solution using Runge-Kutta
    updateKernel(N, _dt, rka[rk], rkb[rk],
                 o_rhsq, o_resq, o_q);

    if (Npml)
      updateKernel(Npml, _dt, rka[rk], rkb[rk],
                   o_rhspmlq, o_respmlq, o_pmlq);
  }
}
*/
} //namespace TimeStepper

} //namespace libp
