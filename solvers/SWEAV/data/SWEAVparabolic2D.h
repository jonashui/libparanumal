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

// Boundary conditions
/* wall 1, outflow 2 */
#define SWEAVDirichletConditions2D(bc, t, x, y, nx, ny, hM, qM, pM, hB, qB, pB) \
{                                             \
  if(bc==1){                                  \
    *(qB) = -qM+2.0*0.18;                     \
    *(pB) = -pM;                              \
  } else if(bc==2){                           \
    *(hB) = -hM+2.0*0.33;                     \
    *(pB) = -pM;                              \
  } else if(bc==3 || bc ==5){                 \
    *(hB) = hM;                               \
    *(qB) = qM - 2*nx*(nx*qM + ny*pM);        \
    *(pB) = pM - 2*ny*(nx*qM + ny*pM);        \
  }                                           \
}

#define SWEAVDerivativeConditions2D(bc, t, x, y, nx, ny, dhdxM, dhdyM, dudxM, dudyM, dvdxM, dvdyM, dhdxB, dhdyB, dudxB, dudyB, dvdxB, dvdyB) \
{                                       \
  if(bc==3 || bc == 5){                            \
    *dhdxB = dhdxM - 2*nx*(nx*dhdxM + ny*dhdyM);   \
    *dhdyB = dhdyM - 2*ny*(nx*dhdxM + ny*dhdyM);   \
    *dudxB = nx*nx*nx*nx*dudxM + 2*nx*nx*nx*ny*dvdxM + 2*nx*nx*ny*ny*dvdyM - 2*nx*ny*ny*ny*dudyM + ny*ny*ny*ny*dudxM;               \
    *dvdxB = -nx*nx*nx*nx*dvdxM + 2*ny*(dudxM-dvdyM)*nx*nx*nx + 2*ny*ny*(dudyM+dvdxM)*nx*nx + ny*ny*ny*ny*dvdxM;                       \
    *dudyB = nx*nx*nx*nx*dudyM + 2*ny*ny*(dudyM+dvdxM)*nx*nx - 2*ny*ny*ny*(dudxM-dvdyM)*nx - ny*ny*ny*ny*dudyM;                        \
    *dvdyB = nx*nx*nx*nx*dvdyM - 2*nx*nx*nx*ny*dvdxM + 2*nx*nx*ny*ny*dudxM + 2*nx*ny*ny*ny*dudyM + ny*ny*ny*ny*dvdyM;               \
  }                                                                                                                                  \
}
// Initial conditions

#define SWEAVInitialConditions2D(t, x, y, elementInfo, elementInfo2, h, q, p) \
{                                         \
if( (8.0<=x) && (x <= 12.0)) {          \
    *(h) =0.33-(0.2-0.05*(x-10)*(x-10));  \
} else {                                  \
    *(h) =0.33;                           \
}                                         \
    *(q)=0.18;                            \
    *(p)=0.0;                             \
}


#define SWEAVSourceTerms2D(xl, yl, t, h, q, p, s1, s2, s3) \
{                                       \
if( (8.0<=xl) && (xl <= 12.0)) {              \
    *s2=-p_grav*h*(-0.10*xl+1.0);        \
} else{                                 \
    *s2=0.0;                            \
}                                       \
  *s1=0.0;                              \
  *s3=0.0;                              \
}