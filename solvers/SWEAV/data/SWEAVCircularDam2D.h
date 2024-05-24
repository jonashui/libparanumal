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
    *(hB) = -hM+2.0*0.5;                      \
    *(qB) = -qM;                               \
    *(pB) = -pM;                              \
  } \
}


#define SWEAVDerivativeConditions2D(bc, t, x, y, nx, ny, dhdxM, dhdyM, dudxM, dudyM, dvdxM, dvdyM, dhdxB, dhdyB, dudxB, dudyB, dvdxB, dvdyB) \
{                                       \
}
// Initial conditions

/*
#define SWEAVInitialConditions2D(t, x, y, elementInfo, elementInfo2, h, q, p) \
{                                       \
    if(x*x+y*y <= 2.5*2.5) {            \
        *(h) = 2.5;                     \
    } else {                            \
        *(h) = 0.5;                     \
    }                                   \
    *(q)=0;                              \
    *(p)=0;                             \
}
*/
#define SWEAVInitialConditions2D(t, x, y, elementInfo, elementInfo2, h, q, p) \
{                                       \
    if(elementInfo2) {            \
        *(h) = 2.5;                     \
    } else {                            \
        *(h) = 0.5;                     \
    }                                   \
    *(q)=0;                              \
    *(p)=0;                             \
}


#define SWEAVSourceTerms2D(xl, yl, t, h, q, p, s1, s2, s3) \
{         \
  *s1=0.0; \
  *s2=0.0; \
  *s3=0.0; \
}