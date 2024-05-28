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
void mesh_t::SetupVToETri2D(){

    Vcounts.malloc(3*(Nelements-totalRingElements));
    EToN.malloc(10*3*(Nelements-totalRingElements));

    if(Nelements < 1000) {
        VToE.malloc(3*Nelements*10);
        counts.calloc(3*Nelements);
    } else {
        VToE.malloc(Nelements*10);
        counts.calloc(Nelements);
    }

    #pragma omp parallel for
    for(dlong e=0;e<Nelements;++e) {
        dlong id = e*Nverts;

        for(int i=0;i<3;++i) {
            hlong v=EToV[id + i];
            VToE[v*10+counts[v]]=e;
            counts[v]++;
        }
    }

    #pragma omp parallel for
    for(dlong e=0;e<Nelements-totalRingElements;++e) {
        dlong id = e*Nverts;

        for(int i=0;i<3;++i) {
            hlong v=EToV[id + i];
            Vcounts[e*3+i]=counts[v];
            for(int j=0;j<counts[v];++j) {
                EToN[e*10*3+10*i+j]=VToE[v*10+j];
            }
        }
    }

    o_EToN=platform.malloc<dlong>(EToN);
    o_Vcounts=platform.malloc<dlong>(Vcounts);
}

} //namespace libp


/*
namespace libp {
void mesh_t::SetupVToETri2D(){

    Vcounts.malloc(3*(Nelements-totalRingElements));
    EToN.malloc(10*3*(Nelements-totalRingElements));

    #pragma omp parallel for
    for(dlong e=0;e<Nelements-totalRingElements;++e) {
        dlong id = e*Nverts;

        for(int i=0;i<3;++i) {
            Vcounts[e*3+i]=0;
            dfloat xe = EX[id+i];
            dfloat ye = EY[id+i];
            for(dlong c=0;c<Nelements;++c) {
                dlong id2 = c*Nverts;
                for(int j=0;j<3;++j) {
                dfloat xc = EX[id2+j];
                dfloat yc = EY[id2+j];
                if( (fabs(xe - xc) < 1e-10) && (fabs(ye - yc) < 1e-10) ) {
                    EToN[e*10*3+10*i+Vcounts[e*3+i]]=c;
                    Vcounts[e*3+i]++;
                }
                }
            }
        }
    }

    o_EToN=platform.malloc<dlong>(EToN);
    o_Vcounts=platform.malloc<dlong>(Vcounts);
}

} //namespace libp
*/



/*
#include "mesh.hpp"

namespace libp {
void mesh_t::SetupVToETri2D(){

    dlong totalcounts=0;
    Ncounts.malloc(Nelements-totalRingElements);
    Vcounts.malloc(3*(Nelements-totalRingElements));

    for(dlong e=0;e<Nelements-totalRingElements;++e) {
        dlong id = e*Nverts;
        Ncounts[e]=totalcounts;

        for(int i=0;i<3;++i) {
            Vcounts[e*3+i]=0;
            dfloat xe = EX[id+i];
            dfloat ye = EY[id+i];
            for(dlong c=0;c<Nelements;++c) {
                dlong id2 = c*Nverts;
                for(int j=0;j<3;++j) {
                dfloat xc = EX[id2+j];
                dfloat yc = EY[id2+j];
                if( (fabs(xe - xc) < 1e-10) && (fabs(ye - yc) < 1e-10) ) {
                    Vcounts[e*3+i]++;
                    totalcounts++;
                }
                }
            }
        }
    }
    printf("\n\ntotalcounts=%d\n\n",totalcounts);
    EToN.malloc(totalcounts);
    
    dlong counter=0;
    for(dlong e=0;e<Nelements-totalRingElements;++e) {
        dlong id = e*Nverts;
        for(int i=0;i<3;++i) {
            dfloat xe = EX[id+i];
            dfloat ye = EY[id+i];
            for(dlong c=0;c<Nelements;++c) {
                dlong id2 = c*Nverts;
                for(int j=0;j<3;++j) {
                dfloat xc = EX[id2+j];
                dfloat yc = EY[id2+j];
                if( (fabs(xe - xc) < 1e-10) && (fabs(ye - yc) < 1e-10) ) {
                    EToN[counter]=c;
                    counter++;
                }
                }
            }
        }
    }

    o_EToN=platform.malloc<dlong>(EToN);
    o_Ncounts=platform.malloc<dlong>(Ncounts);
    o_Vcounts=platform.malloc<dlong>(Vcounts);
}

} //namespace libp


*/