/*

The MIT License (MIT)

Copyright (c) 2.0017-2.002.02.0 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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
/* wall 1, outflow 2.0 */
#define SWECDirichletConditions2D(bc, t, x, y, nx, ny, hM, qM, pM, hB, qB, pB) \
{                                       \
  if(bc==5){                            \
    *(hB) = cos(-t + x)*cos(-t + y) + 2.0; \
    *(qB) = sin(x-t);                         \
    *(pB) = sin(y-t);                         \
  } else if(bc==1){     \
    *(hB) = cos(-t + x)*cos(-t + y) + 2.0; \
    *(qB) = sin(x-t);  \
    *(pB) = sin(y-t);  \
  }                                     \
}

/*
#define SWECDerivativeConditions2D(bc, t, x, y, nx, ny, dhdxM, dhdyM, dqdxM, dqdyM, dpdxM, dpdyM, dhdxB, dhdyB, dqdxB, dqdyB, dpdxB, dpdyB) \
{                                       \
  if(bc==1){                            \
    *dhdxB = sin(-x + t)*cos(-y + t);   \
    *dhdyB = sin(-y + t)*cos(-x + t);   \
    *dqdxB = cos(-x + t);               \
    *dqdyB = 0.0;                       \
    *dpdxB = 0.0;                       \
    *dpdyB = cos(-y + t);               \
  } else if(bc==5){                     \
    *dhdxB = sin(-x + t)*cos(-y + t);   \
    *dhdyB = sin(-y + t)*cos(-x + t);   \
    *dqdxB = cos(-x + t);               \
    *dqdyB = 0.0;                       \
    *dpdxB = 0.0;                       \
    *dpdyB = cos(-y + t);               \
  }                                     \
}
*/

/*
#define SWECDerivativeConditions2D(bc, t, x, y, nx, ny, dhdxM, dhdyM, dudxM, dudyM, dvdxM, dvdyM, dhdxB, dhdyB, dudxB, dudyB, dvdxB, dvdyB) \
{                                                                                                              \
  if(bc==1){                                                                                                   \
    *dhdxB = sin(-x + t)*cos(-y + t);                                                                          \
    *dhdyB = sin(-y + t)*cos(-x + t);                                                                          \
    *dudxB = (cos(-y + t) + 2.0*cos(-x + t))/((cos(-x + t)*cos(-y + t) + 2.0)*(cos(-x + t)*cos(-y + t) + 2.0));        \
    *dudyB = sin(-x + t)*sin(-y + t)*cos(-x + t)/((cos(-x + t)*cos(-y + t) + 2.0)*(cos(-x + t)*cos(-y + t) + 2.0));\
    *dvdxB = sin(-y + t)*sin(-x + t)*cos(-y + t)/((cos(-x + t)*cos(-y + t) + 2.0)*(cos(-x + t)*cos(-y + t) + 2.0));\
    *dvdyB = (cos(-x + t) + 2.0*cos(-y + t))/((cos(-x + t)*cos(-y + t) + 2.0)*(cos(-x + t)*cos(-y + t) + 2.0));        \
  } else if(bc==5){                                                                                            \
    *dhdxB = sin(-x + t)*cos(-y + t);                                                                          \
    *dhdyB = sin(-y + t)*cos(-x + t);                                                                          \
    *dudxB = (cos(-y + t) + 2.0*cos(-x + t))/((cos(-x + t)*cos(-y + t) + 2.0)*(cos(-x + t)*cos(-y + t) + 2.0));        \
    *dudyB = sin(-x + t)*sin(-y + t)*cos(-x + t)/((cos(-x + t)*cos(-y + t) + 2.0)*(cos(-x + t)*cos(-y + t) + 2.0));\
    *dvdxB = sin(-y + t)*sin(-x + t)*cos(-y + t)/((cos(-x + t)*cos(-y + t) + 2.0)*(cos(-x + t)*cos(-y + t) + 2.0));\
    *dvdyB = (cos(-x + t) + 2.0*cos(-y + t))/((cos(-x + t)*cos(-y + t) + 2.0)*(cos(-x + t)*cos(-y + t) + 2.0));        \
  }                                                                                                            \
} */



#define SWECDerivativeConditions2D(bc, t, x, y, nx, ny, dhdxM, dhdyM, dudxM, dudyM, dvdxM, dvdyM, dhdxB, dhdyB, dudxB, dudyB, dvdxB, dvdyB) \
{                                                                                               \
\
}

// Initial conditions

#define SWECInitialConditions2D(t, x, y,elementInfo, h, q, p) \
{                                       \
  *(h) = cos(-t + x)*cos(-t + y) + 2.0;      \
  *(q) = sin(x-t);                          \
  *(p) = sin(y-t);                          \
}

/*
#define SWECInitialConditions2D(t, x, y,elementInfo, h, q, p) \
{                                   \
    if(x < 0) {                   \
        *(h) = 10;                  \
    } else {                        \
        *(h) = 5;                   \
    }                               \
    *(q)=0;                         \
    *(p)=0;                          \
} */



#define SWECSourceTerms2D(xl, yl, t, s1, s2, s3) \
{                                       \
  *s1=-sin(t-xl)*cos(t-yl)-cos(t-xl)*sin(t-yl)+cos(t-xl)+cos(t-yl); \
  *s2=-cos(-xl+t)-sin(-xl+t)*sin(-xl+t)*sin(-xl+t)*cos(-yl+t)/ ((cos(-xl+t)*cos(-yl+t)+2)*(cos(-xl+t)*cos(-yl+t)+2) ) -2*sin(-xl + t)*cos(-xl + t)/(cos(-xl + t)*cos(-yl + t) + 2) +p_grav*(cos(-xl + t)*cos(-yl + t) + 2)*sin(-xl + t)*cos(-yl + t) -sin(-yl + t)*sin(-yl + t)*sin(-xl + t)*cos(-xl + t)/((cos(-xl + t)*cos(-yl + t) + 2)*(cos(-xl + t)*cos(-yl + t) + 2))-cos(-yl + t)*sin(-xl + t)/(cos(-xl + t)*cos(-yl + t) + 2); \
  *s3=-cos(-yl + t) - sin(-yl + t)*sin(-xl + t)*sin(-xl + t)*cos(-yl + t)/((cos(-xl + t)*cos(-yl + t) + 2)*(cos(-xl + t)*cos(-yl + t) + 2))- sin(-yl + t)*cos(-xl + t)/(cos(-xl + t)*cos(-yl + t) + 2)- sin(-yl + t)*sin(-yl + t)*sin(-yl + t)*cos(-xl + t)/((cos(-xl + t)*cos(-yl + t) + 2)*(cos(-xl + t)*cos(-yl + t) + 2))- 2*sin(-yl + t)*cos(-yl + t)/(cos(-xl + t)*cos(-yl + t) + 2)+ p_grav*(cos(-xl + t)*cos(-yl + t) + 2)*cos(-xl + t)*sin(-yl + t); \
}


/*
#define SWECSourceTerms2D(xl, yl, t, s1, s2, s3) \
{ \
  *s1 = -sin(-xl+t)*cos(-yl+t)-sin(-yl+t)*cos(-xl+t)+cos(-xl+t)+cos(-yl+t)+2.0*cos(-xl+t)*cos(-yl+t); \
  *s2 = -cos(-xl+t)-1.0/((cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0))*sin(-xl+t)*sin(-xl+t)*sin(-xl+t)*cos(-yl+t)-2.0/(cos(-xl+t)*cos(-yl+t)+2.0)*sin(-xl+t)*cos(-xl+t)+p_grav*(cos(-xl+t)*cos(-yl+t)+2.0)*sin(-xl+t)*cos(-yl+t)-1.0/((cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0))*sin(-xl+t)*sin(-yl+t)*sin(-yl+t)*cos(-xl+t)-1.0/(cos(-xl+t)*cos(-yl+t)+2.0)*sin(-xl+t)*cos(-yl+t)-sin(-xl+t)*cos(-yl+t)*(cos(-xl+t)/(cos(-xl+t)*cos(-yl+t)+2.0)+sin(-xl+t)*sin(-xl+t)/((cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0))*cos(-yl+t))-(cos(-xl+t)*cos(-yl+t)+2.0)*(sin(-xl+t)/(cos(-xl+t)*cos(-yl+t)+2.0)-3.0*cos(-xl+t)/((cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0))*sin(-xl+t)*cos(-yl+t)-2.0*sin(-xl+t)*sin(-xl+t)*sin(-xl+t)/((cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0))*cos(-yl+t)*cos(-yl+t))+1.0/((cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0))*sin(-xl+t)*sin(-yl+t)*sin(-yl+t)*cos(-xl+t)*cos(-xl+t)+1.0/(cos(-xl+t)*cos(-yl+t)+2.0)*sin(-xl+t)*cos(-yl+t)*cos(-xl+t);  \
  *s3 = -cos(-yl+t)-sin(-yl+t)*sin(-xl+t)*sin(-xl+t)*cos(-yl+t)/((cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0))-sin(-yl+t)*cos(-xl+t)/(cos(-xl+t)*cos(-yl+t)+2.0)-sin(-yl+t)*sin(-yl+t)*sin(-yl+t)*cos(-xl+t)/((cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0))-2.0*sin(-yl+t)*cos(-yl+t)/(cos(-xl+t)*cos(-yl+t)+2.0)+p_grav*(cos(-xl+t)*cos(-yl+t)+2.0)*sin(-yl+t)*cos(-xl+t)+1.0/((cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0))*sin(-yl+t)*sin(-xl+t)*sin(-xl+t)*cos(-yl+t)*cos(-yl+t)+1.0/(cos(-xl+t)*cos(-yl+t)+2.0)*sin(-yl+t)*cos(-xl+t)*cos(-yl+t)-sin(-yl+t)*cos(-xl+t)*(cos(-yl+t)/(cos(-xl+t)*cos(-yl+t)+2.0)+sin(-yl+t)*sin(-yl+t)/((cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0))*cos(-xl+t))-(cos(-xl+t)*cos(-yl+t)+2.0)*(sin(-yl+t)/(cos(-xl+t)*cos(-yl+t)+2.0)-3.0*cos(-yl+t)/((cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0))*sin(-yl+t)*cos(-xl+t)-2.0*sin(-yl+t)*sin(-yl+t)*sin(-yl+t)/((cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0)*(cos(-xl+t)*cos(-yl+t)+2.0))*cos(-xl+t)*cos(-xl+t));  \
} */



