#include "ellipticBenchmarkTri2D.h"

typedef struct {
  char x,y,z,w;
}char4;

void loadElementStiffnessMatricesTri2D(mesh_t *mesh, const char *options, int N){

  char fname[BUFSIZ];
  sprintf(fname, "sparseN%02d.dat", N);
  FILE *fp = fopen(fname, "r");

  char buf[BUFSIZ];
  if (fp == NULL) {
    printf("no file! %s\n", fname);
    exit(-1);
  }

  fgets(buf, BUFSIZ, fp); // read comment
  fgets(buf, BUFSIZ, fp);
  int maxNnzPerRow;
  sscanf(buf, "%d", &maxNnzPerRow);
  printf("maxNNz = %d Np = %d \n", maxNnzPerRow, mesh->Np);

  mesh->maxNnzPerRow = maxNnzPerRow;
  int paddedRowSize = 4*((mesh->maxNnzPerRow+3)/4);
  printf("maxNnzPerRow = %d, paddedNnzPerRow = %d\n", mesh->maxNnzPerRow, paddedRowSize);

  fgets(buf, BUFSIZ, fp); 
  mesh->Ind = (iint*) calloc(paddedRowSize*mesh->Np, sizeof(iint));                     
  for(int n=0;n<mesh->maxNnzPerRow*mesh->Np;++n){ 
    fscanf(fp, "%d ", mesh->Ind+n);                                             
  }   
  mesh->Srr =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  mesh->Srs =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  mesh->Sss =  (dfloat*) calloc(paddedRowSize*mesh->Np, sizeof(dfloat));
  fgets(buf, BUFSIZ, fp);
  for(int r=0;r<mesh->maxNnzPerRow;++r){
    for(int n=0;n<mesh->Np;++n){                                           
      int count = n + r*mesh->Np;
      int aa = fscanf(fp,  "%lf ", mesh->Srr+count); 
    }
  }
  fgets(buf, BUFSIZ, fp);
  for(int r=0;r<mesh->maxNnzPerRow;++r){
    for(int n=0;n<mesh->Np;++n){                                           
      int count = n + r*mesh->Np;
      int aa = fscanf(fp,  "%lf ", mesh->Srs+count); 
    }
  }
  fgets(buf, BUFSIZ, fp);
  for(int r=0;r<mesh->maxNnzPerRow;++r){
    for(int n=0;n<mesh->Np;++n){                                           
      int count = n + r*mesh->Np;
      int aa = fscanf(fp,  "%lf ", mesh->Sss+count); 
    }
  }
  fgets(buf, BUFSIZ, fp);
  
  /*char4 packing*/
  char4 * IndTchar = (char4*) calloc((paddedRowSize*mesh->Np)/4,sizeof(char4));
  for (iint m=0;m<paddedRowSize;m+=4) {  
    for (iint n=0;n<mesh->Np;n++) {
      char4 ind;
      ind.x = mesh->Ind[n + 0*mesh->Np + m*mesh->Np];
      ind.y = mesh->Ind[n + 1*mesh->Np + m*mesh->Np];
      ind.z = mesh->Ind[n + 2*mesh->Np + m*mesh->Np];
      ind.w = mesh->Ind[n + 3*mesh->Np + m*mesh->Np];
      
      IndTchar[n + mesh->Np*(m/4)] = ind;
    }
  }
  
  mesh->maxNnzPerRow = paddedRowSize;
  mesh->o_SrrT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Srr);
  mesh->o_SrsT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Srs);
  mesh->o_SssT = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Sss);
  mesh->o_Srr = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Srr);
  mesh->o_Srs = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Srs);
  mesh->o_Sss = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(dfloat), mesh->Sss);

  mesh->o_IndTchar = mesh->device.malloc(mesh->Np*mesh->maxNnzPerRow*sizeof(char), IndTchar);
  
  fclose(fp); 
}