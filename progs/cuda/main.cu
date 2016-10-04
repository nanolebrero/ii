//#include <cublas_v2.h> 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "avance.cu"
 
// Allocates a matrix with random float entries.
void ceroInit(float* data, int size)
{
   for (int i = 0; i < size; ++i)
   data[i] = 0.;
}
void gaussiana(float* data,  int nodos, int cuerdas)
{
       float cosa = 0;
     for(int i= 0; i< cuerdas; ++i){ 
     for(int j= 0; j < nodos; ++j){ 
      cosa = -0.02*((j-1*nodos/3)*(j-1*nodos/3)) ;     
     data[j+i*nodos] = 2.E-3*exp(cosa);
    // data[j] =( cosa-cosa*cosa/1024.);//^2/1024;
 //   printf ("data %e  \n",   data[j+i*nodos]) ;
   }
  }
}
void randommasa(float* data, int nodos,int cuerdas, float masa, float orden)
{
  for (int j=0; j< cuerdas; ++j) {
   for (int i = 0; i < nodos; ++i) {
 
   data[i+j*nodos] = masa - ((float) rand()/RAND_MAX)*orden*masa;
//    printf ("data %e  \n",   data[i]) ;
    }
  }
}
void friccionpuntas(float* data, int nodos,int cuerdas, float param, float orden)
{
   for (int j=0; j < cuerdas; ++j){
   for (int i = 0; i < nodos; ++i) {
   float factor = 0.2*exp(-orden*(i*i)) + exp(-orden*(i-nodos)*(i-nodos));
   data[i+j*nodos] += factor;
//    printf ("data %e  \n",   data[i]) ;
    }
  }
}
void random(float* data, int nodos,int cuerdas, int i)
{
        //  float cosa = 0;
     for(int j= 0; j < nodos; ++j){
       float cosa=data[j+i*nodos];
     data[j+i*nodos] = cosa*((float) rand()/RAND_MAX-0.5);
  }
}
 
/////////////////////////////////////////////////////////
// Program main
/////////////////////////////////////////////////////////
 
int
main(int argc, char** argv)
{

   // set seed for rand()
   //srand(2006);

 int npart=atoi(argv[1]);    // Cantidad de nodos en cada cuerda	
 int ncuerdas=atoi(argv[2]);  // cantidad de cuerdas
 int npasos=atoi(argv[3]);   // Cantidad de pasos
 float masatension=  atof(argv[4]); // Factor masa-tension
 float friccion=  atof(argv[5]);    // factor de fricción
 bool toco=false ;
   // 1. allocate host memory for matrices A and B
   unsigned int size_A = npart*ncuerdas;
   unsigned int mem_size_A = sizeof(float) * size_A;
   unsigned int buffer_salida = sizeof(float)*256*ncuerdas;
   unsigned int bcur = sizeof(bool)*ncuerdas;   

   float* X = (float*) malloc(mem_size_A);
   float* V = (float*) malloc(mem_size_A);
   float* F = (float*) malloc(mem_size_A);
   float* Fr = (float*) malloc(mem_size_A);
   float* M = (float*) malloc(mem_size_A);
   float* Fext = (float*) malloc(mem_size_A);
   float* salida=(float*) malloc(buffer_salida);
   bool* tococ=(bool*) malloc(bcur);
 
   // 8. allocate device memory
   float* d_X;
   float* d_V;
   float* d_F;
   float* d_Fr;
   float* d_M;
   float* d_Fext;
   float* d_salida;
   bool* d_tococ;
   cudaMalloc((void**) &d_X, mem_size_A);
   cudaMalloc((void**) &d_V, mem_size_A);
   cudaMalloc((void**) &d_F, mem_size_A);
   cudaMalloc((void**) &d_Fr, mem_size_A);
   cudaMalloc((void**) &d_M, mem_size_A);
   cudaMalloc((void**) &d_Fext, mem_size_A);
   cudaMalloc((void**) &d_salida, buffer_salida);
   cudaMalloc((void**) &d_tococ, bcur);

//    cublasHandle_t handle;
//    cublasCreate(&handle);

   // 2. initialize host memory
    ceroInit(X, size_A);
    ceroInit(V, size_A);
    ceroInit(F, size_A);
//    gaussiana(X, size_A);
    gaussiana(Fext, npart, ncuerdas);
    randommasa(M, npart, ncuerdas, masatension, 0.1);
    randommasa(Fr, npart, ncuerdas, friccion, 0.0);
    friccionpuntas(Fr, npart, ncuerdas, friccion, .10);
   

   // 9. copy host memory to device

   cudaMemcpy(d_X, X, mem_size_A, 
   cudaMemcpyHostToDevice);
   cudaMemcpy(d_V, V, mem_size_A, 
   cudaMemcpyHostToDevice);
   cudaMemcpy(d_F, F, mem_size_A, 
   cudaMemcpyHostToDevice);
   cudaMemcpy(d_Fr, Fr, mem_size_A, 
   cudaMemcpyHostToDevice);
   cudaMemcpy(d_M, M, mem_size_A, 
   cudaMemcpyHostToDevice);
   cudaMemcpy(d_Fext, Fext, mem_size_A, 
   cudaMemcpyHostToDevice);

 /*  for (uint nn=0; nn<npart; ++nn) {

    printf ("masa %e  \n", M[nn]) ;
    printf ("friccion %e  \n", Fr[nn]) ;

    } 
*/
    int sig=1;
   for (uint i = 0 ; i < npasos ; ++i) 
 { 
   
 
   // 5. perform the calculation
   // setup execution parameters
   dim3 threads(npart);
   dim3 grid(ncuerdas);
   // execute the kernel

  float random;
  for(int ii=0; ii < ncuerdas ; ++ii) {   

     tococ[ii] = 0;

      
 //    if(i > ii*450 + 10 && i < (ii + 1)*450 +20 ) {
     if(i == ii*400 + 10 ) { //&& i < (ii + 1)*150 +20 ) {
         tococ[ii]=1;
//      printf("toco la cuerda, %d %d", ii, i) ;
//      random(Fext, npart, ncuerdas,ii);
      }
    }
//   cudaMemcpy(d_Fext, Fext, mem_size_A, 
//   cudaMemcpyHostToDevice);

   cudaMemcpy(d_tococ, tococ, bcur, 
   cudaMemcpyHostToDevice);

/*    if(i % 100 == 0 ) {

       sig=-1*sig;
      }*/
       random=sig*((float) rand()/RAND_MAX-0.5);
   

   avance<<< grid, threads >>>(d_X, d_V, 
                                  d_F, d_Fext, d_salida, npart, d_M, d_Fr, d_tococ, random);
   

   // 11. copy result from device to host
   cudaMemcpy(salida, d_salida, buffer_salida, 
   cudaMemcpyDeviceToHost);
/*    for (uint j = 0 ; j< npart; ++j) {
    printf ("Xs %e  %d \n", X[j], j ) ; 
             }
*/
    for (uint j = 0 ; j< 256; ++j) {
       int sal =0 ;
       int sal2 =0 ;
         for (uint i=0; i< ncuerdas  ;++i) {
        if(((i+1) % 2) == 0){
         sal += 100000*salida[j+i*256];
          }
          else {
         sal2 += 100000*salida[j+i*256];
             }
                  }
    printf (" %d , %d  \n", sal , sal2) ; 
            }
   }
 
}

