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

#include "mesh.hpp"

namespace libp {

void mesh_t::CubatureSetupTri2DCurv(){

  /* Quadrature data */
  cubN = 3*(N+1); //cubature order
    CubatureNodesTri2D(cubN, cubNp, cubr, cubs, cubw);

  InterpolationMatrixTri2D(N, r, s, cubr, cubs, cubInterp);

  //cubature project cubProject = M^{-1} * cubInterp^T
  // Defined such that cubProject * cubW * cubInterp = Identity
  CubaturePmatrixTri2D(N, r, s, cubr, cubs, cubProject);

  //cubature derivates matrices, cubD: differentiate on cubature nodes
  // we dont use cubD on Tris/Tets  so skip computing
  

  // Instead, it's cheaper to:
  // make weak cubature derivatives cubPDT = cubProject * cubD^T 
  CubatureWeakDmatricesTri2D(N, r, s,
                             cubr, cubs,
                             cubPDT);
  cubPDrT = cubPDT + 0*cubNp*Np;
  cubPDsT = cubPDT + 1*cubNp*Np;

  // cubN+1 point Gauss-Legendre quadrature for surface integrals
  cubNq  = cubN+1;
  cubNfp = cubN+1;
  intNfp = cubN+1;
  JacobiGQ(0, 0, cubN, intr, intw);

  CubatureSurfaceMatricesTri2D(N, r, s, faceNodes,
                               intr, intw,
                               intInterp, intLIFT);

  // add compile time constants to kernels
  props["defines/" "p_cubNq"]= cubNq;
  props["defines/" "p_cubNp"]= cubNp;
  props["defines/" "p_intNfp"]= intNfp;
  props["defines/" "p_intNfpNfaces"]= intNfp*Nfaces;
  props["defines/" "p_cubNfp"]= cubNfp;

    // build transposes (we hold matrices as column major on device)
  memory<dfloat> cubProjectT(cubNp*Np);
  memory<dfloat> cubInterpT(cubNp*Np);
  linAlg_t::matrixTranspose(cubNp, Np, cubInterp, Np, cubInterpT, cubNp);
  linAlg_t::matrixTranspose(Np, cubNp, cubProject, cubNp, cubProjectT, Np);

  //pre-multiply cubProject by W on device
  for(int n=0;n<cubNp;++n){
    for(int m=0;m<Np;++m){
      cubProjectT[m+n*Np] *= cubw[n];
    }
  }

  memory<dfloat> cubPDTT(2*cubNp*Np);
  memory<dfloat> cubPDrTT = cubPDTT + 0*cubNp*Np;
  memory<dfloat> cubPDsTT = cubPDTT + 1*cubNp*Np;
  linAlg_t::matrixTranspose(Np, cubNp, cubPDrT, cubNp, cubPDrTT, Np);
  linAlg_t::matrixTranspose(Np, cubNp, cubPDsT, cubNp, cubPDsTT, Np);

  //pre-multiply cubPDT by W on device
  for(int n=0;n<cubNp;++n){
    for(int m=0;m<Np;++m){
      cubPDrTT[m+n*Np] *= cubw[n];
      cubPDsTT[m+n*Np] *= cubw[n];
    }
  }

  // build surface integration matrix transposes
  memory<dfloat> intLIFTT(Np*Nfaces*intNfp);
  memory<dfloat> intInterpT(Nfp*Nfaces*intNfp);
  linAlg_t::matrixTranspose(Np, Nfaces*intNfp, intLIFT, Nfaces*intNfp, intLIFTT, Np);
  linAlg_t::matrixTranspose(Nfaces*intNfp, Nfp, intInterp, Nfp, intInterpT, Nfaces*intNfp);

  o_cubInterp  = platform.malloc<dfloat>(Np*cubNp, cubInterpT);
  o_cubProject = platform.malloc<dfloat>(Np*cubNp, cubProjectT);

  o_cubPDT = platform.malloc<dfloat>(2*cubNp*Np, cubPDTT);

  o_intInterp = platform.malloc<dfloat>(Nfp*Nfaces*intNfp, intInterpT);
  o_intLIFT   = platform.malloc<dfloat>(Np*Nfaces*intNfp, intLIFTT);

  o_cubwJ = o_wJ;
  o_cubvgeo = o_vgeo;
  o_cubggeo = o_ggeo;
  o_cubsgeo = o_sgeo;


  cubwJCurv.malloc(Ncurv*cubNp);
  cubvgeoCurv.malloc(Ncurv*Nvgeo*cubNp);
  cubggeoCurv.malloc(Ncurv*Nggeo*cubNp);
  cubsgeoCurv.malloc(Ncurv*Nsgeo*cubNq*Nfaces);
  //temp arrays
  memory<dfloat> xre(Np);
  memory<dfloat> xse(Np);
  memory<dfloat> yre(Np);
  memory<dfloat> yse(Np);

  //geometric data for quadrature
  for(dlong e=0;e<Nelements;++e){ /* for each element */
    if (mapCurv[e] < 0) continue;
    dlong eC = mapCurv[e];
    for(int n=0;n<Np;++n){
        xre[n] = 0; xse[n] = 0;
        yre[n] = 0; yse[n] = 0;

        for(int m=0;m<Np;++m){
          int id = e*Np + m;
          xre[n] += Dr[n*Np+m]*x[id];
          xse[n] += Ds[n*Np+m]*x[id];
          yre[n] += Dr[n*Np+m]*y[id];
          yse[n] += Ds[n*Np+m]*y[id];
        }
    }

    //interpolate derivaties to cubature
      for(int n=0;n<cubNp;++n){
        dfloat xr = 0.0, xs = 0.0;
        dfloat yr = 0.0, ys = 0.0;
        for(int i=0;i<Np;++i){
          xr += cubInterp[n*Np+i]*xre[i];
          xs += cubInterp[n*Np+i]*xse[i];
          yr += cubInterp[n*Np+i]*yre[i];
          ys += cubInterp[n*Np+i]*yse[i];
        }

        /* compute geometric factors for affine coordinate transform*/
        dfloat J = xr*ys - xs*yr;

        LIBP_ABORT("Negative J found at element " << e, J<1e-8);
        //printf("J=%.15lf e=%d\n", J,e);

        dfloat rx =  ys/J;
        dfloat ry = -xs/J;
        dfloat sx = -yr/J;
        dfloat sy =  xr/J;


        /* store geometric factors */
        dlong base = Nvgeo*cubNp*eC + n;
        cubvgeoCurv[base + cubNp*RXID] = rx;
        cubvgeoCurv[base + cubNp*RYID] = ry;
        cubvgeoCurv[base + cubNp*SXID] = sx;
        cubvgeoCurv[base + cubNp*SYID] = sy;
        cubvgeoCurv[base + cubNp*JID]  = J; 

        /* store second order geometric factors */
        base = Nggeo*cubNp*eC + n;
        cubggeoCurv[base + cubNp*G00ID] = J*(rx*rx + ry*ry);
        cubggeoCurv[base + cubNp*G01ID] = J*(rx*sx + ry*sy);
        cubggeoCurv[base + cubNp*G11ID] = J*(sx*sx + sy*sy);

        cubwJCurv[cubNp*eC + n] = J;
    }
    //printf("_______________________\n");
    for(int f=0;f<Nfaces;++f){ // for each face
      for(int m=0;m<intNfp;++m){  // for each node on face

        //interpolate derivatives of physical coordinates
        dfloat xr = 0.0, xs = 0.0;
        dfloat yr = 0.0, ys = 0.0;
        for(int n=0;n<Nfp;++n){  // for each node on face
          /* volume index of face node */
          dlong idn = faceNodes[f*Nfp+n];
          xr += intInterp[n+m*Nfp+f*intNfp*Nfp]*xre[idn];
          xs += intInterp[n+m*Nfp+f*intNfp*Nfp]*xse[idn];
          yr += intInterp[n+m*Nfp+f*intNfp*Nfp]*yre[idn];
          ys += intInterp[n+m*Nfp+f*intNfp*Nfp]*yse[idn];
        }

        /* compute geometric factors for affine coordinate transform*/
        dfloat J = xr*ys - xs*yr;

        dfloat rx =  ys/J;
        dfloat ry = -xs/J;
        dfloat sx = -yr/J;
        dfloat sy =  xr/J;


        /* face f normal and length */
        dfloat nx=0.0, ny=0.0;
        switch(f){
        case 0: nx = -sx; ny = -sy; break;
        case 1: nx = rx+sx; ny = ry+sy; break;
        case 2: nx = -rx; ny = -ry; break;
        }
        dfloat  sJ = sqrt((nx)*(nx)+(ny)*(ny));
        nx /= sJ; ny /= sJ;
        sJ *= J;

        //printf("sJ=%lf, ",sJ);
        /* output index */
        dlong base = Nsgeo*(Nfaces*cubNfp*eC + cubNfp*f + m);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        cubsgeoCurv[base+NXID] = nx;
        cubsgeoCurv[base+NYID] = ny;
        cubsgeoCurv[base+SJID] = sJ;
        cubsgeoCurv[base+IJID] = 1./J;
      }
    }
    //printf("\n_______________________\n");
  }
  
  DmatrixTri2DCurv(N, r, s, cubr, cubs, cubD);
  invMMs.malloc(Np*Np*Ncurv);
  cubPDTs.malloc(2*Np*cubNp*Ncurv);
  intLIFTs.malloc(Nfaces*intNfp*Np*Ncurv);
  
  memory<dfloat> cubPDTTs(2*cubNp*Np*Ncurv);
  memory<dfloat> intLIFTTs(Np*Nfaces*intNfp*Ncurv);
  memory<dfloat> Je;
  memory<dfloat> invMMe;
  memory<dfloat> cubPDTe;
  memory<dfloat> intLIFTe;
  memory<dfloat> intLIFTTe;


  for(dlong e=0;e<Nelements;++e){
      if (mapCurv[e] < 0) continue;
      dlong eC = mapCurv[e];
      Je=cubvgeoCurv+Nvgeo*cubNp*eC+cubNp*JID;
      //printf("___________________________\n");
      //printf("e=%d\n",e);
      
      invMMe=invMMs+eC*Np*Np;
      invMassMatrixTri2DCurv(Np,cubNp,cubInterp,Je,cubw,invMMe);

      dlong dbase = 2*Np*cubNp*eC;
      cubPDTe=cubPDTs+dbase;

      
      CubatureWeakDmatricesTri2DCurv(N, r, s,
                             cubr, cubs, invMMe, cubD, cubPDTe);

      cubPDrT = cubPDTe + 0*cubNp*Np;
      cubPDsT = cubPDTe + 1*cubNp*Np;
  
      cubPDrTT = cubPDTTs+dbase + 0*cubNp*Np;
      cubPDsTT = cubPDTTs+dbase + 1*cubNp*Np;
      
      linAlg_t::matrixTranspose(Np, cubNp, cubPDrT, cubNp, cubPDrTT, Np);
      linAlg_t::matrixTranspose(Np, cubNp, cubPDsT, cubNp, cubPDsTT, Np);
      
      //pre-multiply cubPDT by W on device
      for(int n=0;n<cubNp;++n){
        for(int m=0;m<Np;++m){
          cubPDrTT[m+n*Np] *= cubw[n];
          cubPDsTT[m+n*Np] *= cubw[n];
        }
      }

      intLIFTe=intLIFTs+eC*Nfaces*intNfp*Np;
      CubatureSurfaceMatricesTri2DCurv(N, r, s, faceNodes,
                               intr, intw,
                               invMMe,intInterp,intLIFTe);
      intLIFTTe=intLIFTTs+eC*Nfaces*intNfp*Np;
      linAlg_t::matrixTranspose(Np, Nfaces*intNfp, intLIFTe, Nfaces*intNfp, intLIFTTe, Np);
  }
  
  o_cubPDTs = platform.malloc<dfloat>(2*cubNp*Np*Ncurv, cubPDTTs);
  o_intLIFTs = platform.malloc<dfloat>(Ncurv*Np*Nfaces*intNfp, intLIFTTs);

  o_cubwJCurv   = platform.malloc<dfloat>(Ncurv*cubNp, cubwJCurv);
  o_cubvgeoCurv = platform.malloc<dfloat>(Ncurv*Nvgeo*cubNp, cubvgeoCurv);
  o_cubggeoCurv = platform.malloc<dfloat>(Ncurv*Nggeo*cubNp, cubggeoCurv);
  o_cubsgeoCurv = platform.malloc<dfloat>(Ncurv*Nfaces*cubNq*Nsgeo, cubsgeoCurv);
}

} //namespace libp
