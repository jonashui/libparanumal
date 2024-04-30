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

void SWE_t::Setup(platform_t& _platform, mesh_t& _mesh,
                        SWESettings_t& _settings){

  platform = _platform;
  mesh = _mesh;
  comm = _mesh.comm;
  settings = _settings;

  Nfields = (mesh.dim==3) ? 4:3;
  Ngrads = 6;

  dlong Nlocal = mesh.Nelements*mesh.Np*Nfields;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np*Nfields;

  dlong NlocalGrads = mesh.Nelements*mesh.Np*Ngrads;
  dlong NhaloGrads  = mesh.totalHaloPairs*mesh.Np*Ngrads;



  //Trigger JIT kernel builds
  ogs::InitializeKernels(platform, ogs::Dfloat, ogs::Add);

  cubature   = (settings.compareSetting("ADVECTION TYPE", "CUBATURE")) ? 1:0;

  if (cubature) {
    mesh.CubatureSetup();
    mesh.CubaturePhysicalNodes();
  }

  //setup linear algebra module
  platform.linAlg().InitKernels({"innerProd"});

  /*setup trace halo exchange */
  fieldTraceHalo = mesh.HaloTraceSetup(Nfields);
  gradTraceHalo  = mesh.HaloTraceSetup(Ngrads);

  //setup timeStepper
  if (settings.compareSetting("TIME INTEGRATOR","AB3")){
    timeStepper.Setup<TimeStepper::ab3>(mesh.Nelements,
                                        mesh.totalHaloPairs,
                                        mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","LSERK4")){
    timeStepper.Setup<TimeStepper::lserk4>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","DOPRI5")){
    timeStepper.Setup<TimeStepper::dopri5>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  } else if (settings.compareSetting("TIME INTEGRATOR","SSPRK2")){
    timeStepper.Setup<TimeStepper::ssprk2>(mesh.Nelements,
                                           mesh.totalHaloPairs,
                                           mesh.Np, Nfields, platform, comm);
  }

  // set penalty parameter
  dfloat Lambda2 = 0.5;
  dfloat p_grav = 9.81;

  // compute samples of q at interpolation nodes
  q.malloc(Nlocal+Nhalo);
  o_q = platform.malloc<dfloat>(q);

  gradq.malloc(NlocalGrads+NhaloGrads);
  o_gradq = platform.malloc<dfloat>(gradq);

  //storage for M*q during reporting
  o_Mq = platform.malloc<dfloat>(q);
  mesh.MassMatrixKernelSetup(Nfields); // mass matrix operator

  // OCCA build stuff
  properties_t kernelInfo = mesh.props; //copy base occa properties

  //add boundary data to kernel info
  std::string dataFileName;
  settings.getSetting("DATA FILE", dataFileName);
  kernelInfo["includes"] += dataFileName;

  kernelInfo["defines/" "p_Nfields"]= Nfields;
  kernelInfo["defines/" "p_Ngrads"]= Ngrads;

  const dfloat p_half = 1./2.;
  kernelInfo["defines/" "p_half"]= p_half;

  const dfloat p_quarter = 1./4.;
  kernelInfo["defines/" "p_quarter"]= p_quarter;

  kernelInfo["defines/" "p_grav"]= p_grav;

 
  

  int maxNodes = std::max(mesh.Np, (mesh.Nfp*mesh.Nfaces));
  kernelInfo["defines/" "p_maxNodes"]= maxNodes;

  int blockMax = 256;
  if (platform.device.mode() == "CUDA") blockMax = 512;

  int NblockV = std::max(1, blockMax/mesh.Np);
  kernelInfo["defines/" "p_NblockV"]= NblockV;

  int NblockS = std::max(1, blockMax/maxNodes);
  kernelInfo["defines/" "p_NblockS"]= NblockS;

  kernelInfo["defines/" "p_Lambda2"]= Lambda2;

if (cubature) {
    int cubMaxNodes = std::max(mesh.Np, (mesh.intNfp*mesh.Nfaces));
    kernelInfo["defines/" "p_cubMaxNodes"]= cubMaxNodes;
    int cubMaxNodes1 = std::max(mesh.Np, (mesh.intNfp));
    kernelInfo["defines/" "p_cubMaxNodes1"]= cubMaxNodes1;

    int cubNblockV = std::max(1, blockMax/mesh.cubNp);
    kernelInfo["defines/" "p_cubNblockV"]= cubNblockV;

    int cubNblockS = std::max(1, blockMax/cubMaxNodes);
    kernelInfo["defines/" "p_cubNblockS"]= cubNblockS;
  }


  
  //kernelInfo["defines/" "g"] = g;

  // set kernel name suffix
  std::string suffix;
  if(mesh.elementType==Mesh::TRIANGLES)
    suffix = "Tri2D";
  if(mesh.elementType==Mesh::QUADRILATERALS)
    suffix = "Quad2D";
  if(mesh.elementType==Mesh::TETRAHEDRA)
    suffix = "Tet3D";
  if(mesh.elementType==Mesh::HEXAHEDRA)
    suffix = "Hex3D";
  if(mesh.elementType==Mesh::CURVEDTRIANGLES)
    suffix = "Tri2DCurv";

  std::string oklFilePrefix = DSWE "/okl/";
  std::string oklFileSuffix = ".okl";

  std::string fileName, kernelName;

  // kernels from volume file
  fileName   = oklFilePrefix + "SWEGradVolume" + suffix + oklFileSuffix;
  kernelName = "SWEGradVolume" + suffix;

  gradVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                           kernelInfo);
  // kernels from surface file
  fileName   = oklFilePrefix + "SWEGradSurface" + suffix + oklFileSuffix;
  kernelName = "SWEGradSurface" + suffix;

  gradSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                           kernelInfo);


if (cubature) {
  // kernels from volume file
    fileName   = oklFilePrefix + "SWECubatureVolume" + suffix + oklFileSuffix;
    kernelName = "SWECubatureVolume" + suffix;

    cubatureVolumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
    // kernels from surface file
    fileName   = oklFilePrefix + "SWECubatureSurface" + suffix + oklFileSuffix;
    kernelName = "SWECubatureSurface" + suffix;

    cubatureSurfaceKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  }
  else {
    fileName   = oklFilePrefix + "SWEVolume" + suffix + oklFileSuffix;
    kernelName = "SWEVolume" + suffix;

    volumeKernel =  platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
    // kernels from surface file
    fileName   = oklFilePrefix + "SWESurface" + suffix + oklFileSuffix;
    kernelName = "SWESurface" + suffix;

    surfaceKernel = platform.buildKernel(fileName, kernelName,
                                         kernelInfo);
  }

  if (mesh.dim==2) {
    fileName   = oklFilePrefix + "SWEInitialCondition2D" + oklFileSuffix;
    kernelName = "SWEInitialCondition2D";
  } else {
    fileName   = oklFilePrefix + "SWEInitialCondition3D" + oklFileSuffix;
    kernelName = "SWEInitialCondition3D";
  }

  initialConditionKernel = platform.buildKernel(fileName, kernelName,
                                                  kernelInfo);
}
