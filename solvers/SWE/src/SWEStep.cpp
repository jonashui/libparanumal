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

#include "SWE.hpp"

dfloat SWE_t::MaxWaveSpeed(){
  //wavespeed is constant 1 everywhere
  const dfloat vmax = 1.0;
  return vmax;
}

//evaluate ODE rhs = f(q,t)
void SWE_t::rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

  // extract q halo on DEVICE
  
  traceHalo.ExchangeStart(o_Q, 1);
  if (cubature) {
        cubatureVolumeKernel(mesh.Nelements,
                         mesh.o_cubvgeo,
                         mesh.o_cubvgeoCurv,
                         mesh.o_mapCurv,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubPDTs,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         T,
                         o_Q,
                         o_RHS);
  }
  else {
    volumeKernel(mesh.Nelements,
               mesh.o_vgeo,
               mesh.o_D,
               T,
               mesh.o_x,
               mesh.o_y,
               o_Q,
               o_RHS);
  }

  traceHalo.ExchangeFinish(o_Q, 1);

    if (cubature) {
            cubatureSurfaceKernel(mesh.Nelements,
                            mesh.o_cubsgeo,
                            mesh.o_cubsgeoCurv,
                            mesh.o_mapCurv,
                            mesh.o_vmapM,
                            mesh.o_vmapP,
                            mesh.o_EToB,
                            mesh.o_intInterp,
                            mesh.o_intLIFT,
                            mesh.o_intLIFTs,
                            mesh.o_intx,
                            mesh.o_inty,
                            mesh.o_intz,
                            T,
                            o_Q,
                            o_RHS);
    }
    else {
          surfaceKernel(mesh.Nelements,
                  mesh.o_sgeo,
                  mesh.o_LIFT,
                  mesh.o_vmapM,
                  mesh.o_vmapP,
                  mesh.o_EToB,
                  T,
                  mesh.o_x,
                  mesh.o_y,
                  mesh.o_z,
                  o_Q,
                  o_RHS);
  }
}
