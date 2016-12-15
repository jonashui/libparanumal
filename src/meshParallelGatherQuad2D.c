#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"

// gathers to base nodes only
void meshParallelGatherQuad2D(mesh2D *mesh, dfloat *q,
			      dfloat *sendBuffer, dfloat *recvBuffer){

  int rank, size, Nbytes = sizeof(dfloat), tag = 111;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Status *recvStatuses = (MPI_Status*) calloc(size, sizeof(MPI_Status));
  MPI_Status *sendStatuses = (MPI_Status*) calloc(size, sizeof(MPI_Status));

  MPI_Request *recvRequests = (MPI_Request*) calloc(size, sizeof(MPI_Request));
  MPI_Request *sendRequests = (MPI_Request*) calloc(size, sizeof(MPI_Request));
  
  // grab node data to send
  for(iint n=0;n<mesh->gatherNsend;++n){
    // orders nodes as expected by base nodes on other processes
    iint nodeId = mesh->gatherSendIds[n]; 
    sendBuffer[n] = q[nodeId]; // only works for dfloats
  }

  // send non-base data out and recv other non-base data 
  int sendMessage = 0, recvMessage = 0;
  for(iint r=0;r<size;++r){
    if(r!=rank){
      if(mesh->gatherSendCounts[r]>0){
	MPI_Isend((char*)sendBuffer+Nbytes*mesh->gatherSendDispls[r],
		  Nbytes*mesh->gatherSendCounts[r],
		  MPI_CHAR, r, tag, MPI_COMM_WORLD, sendRequests+sendMessage);
	++sendMessage;
      }
      if(mesh->gatherRecvCounts[r]>0){
	MPI_Irecv((char*)recvBuffer+Nbytes*mesh->gatherRecvDispls[r],
		  Nbytes*mesh->gatherRecvCounts[r],
		  MPI_CHAR, r, tag, MPI_COMM_WORLD, recvRequests+recvMessage);
	++recvMessage;
      }
    }
  }

  // do gather operation for base nodes on this rank  (only local reduction)
  for(iint b=0;b<mesh->gatherNbaseRecv;++b){
    iint localStart = mesh->gatherLocalStarts[b];
    iint localEnd   = mesh->gatherLocalStarts[b+1];
    iint baseId     = mesh->gatherLocalBaseIds[b];
    dfloat qb = q[baseId];
    for(iint n=localStart;n<localEnd;++n){ // excludes self
      iint nodeId = mesh->gatherLocalSourceIds[n];
      qb += q[nodeId]; // only works for dfloats
    }
    q[baseId] = qb;
  }

  MPI_Waitall(recvMessage, recvRequests, recvStatuses);

  // complete gather operation for base nodes using incoming non-base node data
  for(iint b=0;b<mesh->gatherNbaseRecv;++b){
    iint sourceStart = mesh->gatherRecvStarts[b];
    iint sourceEnd   = mesh->gatherRecvStarts[b+1];
    iint baseId      = mesh->gatherRecvBaseIds[b];
    dfloat qb = q[baseId]; // need to find baseId
    for(iint n=sourceStart;n<sourceEnd;++n){ // excludes self
      iint nodeId = mesh->gatherRecvSourceIds[n];
      qb += recvBuffer[nodeId]; // only works for dfloats
    }
    q[baseId] = qb;
  }

  MPI_Waitall(sendMessage, sendRequests, sendStatuses);

  free(sendRequests); free(sendStatuses); free(recvRequests); free(recvStatuses);
}