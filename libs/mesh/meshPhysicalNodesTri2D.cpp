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

void mesh_t::PhysicalNodesTri2D(){

  x.malloc(Nelements*Np);
  y.malloc(Nelements*Np);
  elementInfo2.malloc(Nelements);
  elementInfo3.malloc(Nelements);

  #pragma omp parallel for
  for(dlong e=0;e<Nelements;++e){ /* for each element */

    dlong id = e*Nverts+0;

    dfloat xe1 = EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = EX[id+1];
    dfloat xe3 = EX[id+2];

    dfloat ye1 = EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = EY[id+1];
    dfloat ye3 = EY[id+2];

    dfloat xavg=(xe1+xe2+xe3)/3.0;
    dfloat yavg=(ye1+ye2+ye3)/3.0;
    if( xavg <= 105) {
      elementInfo2[e]=1;
    } else {
      elementInfo2[e]=0;
    }

    if( xavg*xavg+yavg*yavg <= 2.5*2.5) {
      elementInfo3[e]=1;
    } else {
      elementInfo3[e]=0;
    }

    for(int n=0;n<Np;++n){ /* for each node */

      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = r[n];
      dfloat sn = s[n];

      /* physical coordinate of interpolation node */
      x[e*Np+n] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      y[e*Np+n] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
    }
  }

  o_x = platform.malloc<dfloat>(x);
  o_y = platform.malloc<dfloat>(y);
  o_elementInfo=platform.malloc<dlong>(elementInfo2);
  o_elementInfo2=platform.malloc<dlong>(elementInfo3);
}

} //namespace libp
