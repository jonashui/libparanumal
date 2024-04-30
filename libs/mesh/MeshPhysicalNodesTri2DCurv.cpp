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

void mesh_t::PhysicalNodesTri2DCurv(){

  x.malloc(Nelements*Np);
  y.malloc(Nelements*Np);

  mapCurv.malloc(Nelements);
  Ncurv=0;

  #pragma omp parallel for
  for(dlong e=0;e<Nelements;++e){ /* for each element */
    mapCurv[e]=-1;
    dlong id = e*Nverts+0;

    dfloat xe1 = EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = EX[id+1];
    dfloat xe3 = EX[id+2];

    dfloat ye1 = EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = EY[id+1];
    dfloat ye3 = EY[id+2];

    for(int n=0;n<Np;++n){ /* for each node */

      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = r[n];
      dfloat sn = s[n];

      /* physical coordinate of interpolation node */
      x[e*Np+n] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      y[e*Np+n] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
    }

    for(int f=0;f<Nfaces;++f){
      if(EToB[f+Nfaces*e] == 5) {
        mapCurv[e]=Ncurv;
        Ncurv +=1;
        dfloat x1=0.0; dfloat x2=0.0; dfloat y1=0.0; dfloat y2=0.0;
        memory<dfloat> vr;
        switch (f)
        {
        case 0:
          x1 = EX[id+0];
          x2 = EX[id+1];

          y1 = EY[id+0];
          y2 = EY[id+1];

          vr = r;
          break;
        case 1:
          x1 = EX[id+1];
          x2 = EX[id+2];

          y1 = EY[id+1];
          y2 = EY[id+2];

          vr = s;
          break;
        case 2:
          x1 = EX[id+0];
          x2 = EX[id+2];

          y1 = EY[id+0];
          y2 = EY[id+2];

          vr = s;
          break;
        }
      
        dfloat theta1 = atan2(y1-0.0, x1-0.0); 
        dfloat theta2 = atan2(y2-0.0, x2-0.0);

        //check to make sure they are in the same quadrant
        if ((theta2 > 0) && (theta1 < 0) &&  (fabs(theta1-theta2) > M_PI)) theta1 = theta1 + 2.0*M_PI;
        if ((theta1 > 0) && (theta2 < 0) &&  (fabs(theta1-theta2) > M_PI)) theta2 = theta2 + 2.0*M_PI;

        memory<dfloat> fdx(Nfp);
        memory<dfloat> fdy(Nfp);
        memory<dfloat> fr(Nfp);
        for (int i=0;i<Nfp;i++) {
          fr[i]= vr[faceNodes[f*Nfp+i]];
          dfloat theta = 0.5*theta1*(1.0-fr[i]) + 0.5*theta2*(1.0+fr[i]);
          
          fdx[i] = 0.0 + 1.0*cos(theta)-x[Np*e+faceNodes[f*Nfp+i]]; 
          fdy[i] = 0.0 + 1.0*sin(theta)-y[Np*e+faceNodes[f*Nfp+i]];
        }

        memory<dfloat> Vface, Vvol;
        Vandermonde1D(N, fr, Vface);
        Vandermonde1D(N, vr, Vvol);

        linAlg_t::matrixInverse(Nfp, Vface);

        memory<dfloat> vdx_temp(Nfp);
        memory<dfloat> vdy_temp(Nfp);
        memory<dfloat> vdx(Np);
        memory<dfloat> vdy(Np);
        for(int n=0;n<Nfp;++n){
            dfloat resx = 0.0; dfloat resy = 0.0;
            for(int i=0;i<Nfp;++i){
              resx += Vface[n*Nfp+i]*fdx[i];
              resy += Vface[n*Nfp+i]*fdy[i];
            }
            vdx_temp[n] = resx;
            vdy_temp[n] = resy;
        }
        for(int n=0;n<Np;++n){
            dfloat resx = 0.0; dfloat resy = 0.0;
            for(int i=0;i<Nfp;++i){
              resx += Vvol[n*Nfp+i]*vdx_temp[i];
              resy += Vvol[n*Nfp+i]*vdy_temp[i];
            }
            vdx[n] = resx;
            vdy[n] = resy;
        }

        for(int n=0;n<Np;++n) {
          if(fabs(1-vr[n]) > 1e-7) {
            dfloat blend=0.0;
            switch (f)
            {
            case 0:
              blend=-(r[n]+s[n])/(1.0-vr[n]);
              break;
            case 1:
              blend=(r[n]+1.0)/(1.0-vr[n]);
              break;
            case 2:
              blend=-(r[n]+s[n])/(1.0-vr[n]);
              break;
            }   
          x[e*Np+n] += blend*vdx[n];
          y[e*Np+n] += blend*vdy[n];
          }
        }
      }
    }
  }

  /*for(dlong e=0;e<Nelements;++e) {
    for(int n=0;n<Np;++n) {
      printf("x=%.15lf e=%d ", x[e*Np+n],e);
      printf("y=%.15lf e=%d\n", y[e*Np+n],e);
    }
    printf("_______________________\n");
  } */
  o_mapCurv = platform.malloc<dlong>(mapCurv);
  o_x = platform.malloc<dfloat>(x);
  o_y = platform.malloc<dfloat>(y);
}

} //namespace libp
