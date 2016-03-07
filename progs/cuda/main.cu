
// Multiply two matrices A * B = C
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
      cosa = -0.01*((j-2*nodos/10)*(j-2*nodos/10)) ;     
     data[j+i*nodos] = 2.E-3*exp(cosa);
    // data[j] =( cosa-cosa*cosa/1024.);//^2/1024;
//    printf ("data %e  \n",   cosa) ;
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
   float factor = exp(-orden*(i*i)) + exp(-orden*(i-nodos)*(i-nodos));
   data[i+j*nodos] += factor;
//    printf ("data %e  \n",   data[i]) ;
    }
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

   float* X = (float*) malloc(mem_size_A);
   float* V = (float*) malloc(mem_size_A);
   float* F = (float*) malloc(mem_size_A);
   float* Fr = (float*) malloc(mem_size_A);
   float* M = (float*) malloc(mem_size_A);
   float* Fext = (float*) malloc(mem_size_A);
   float* salida=(float*) malloc(buffer_salida);

 
   // 8. allocate device memory
   float* d_X;
   float* d_V;
   float* d_F;
   float* d_Fr;
   float* d_M;
   float* d_Fext;
   float* d_salida;
   cudaMalloc((void**) &d_X, mem_size_A);
   cudaMalloc((void**) &d_V, mem_size_A);
   cudaMalloc((void**) &d_F, mem_size_A);
   cudaMalloc((void**) &d_Fr, mem_size_A);
   cudaMalloc((void**) &d_M, mem_size_A);
   cudaMalloc((void**) &d_Fext, mem_size_A);
   cudaMalloc((void**) &d_salida, buffer_salida);

//    cublasHandle_t handle;
//    cublasCreate(&handle);

   // 2. initialize host memory
    ceroInit(X, size_A);
    ceroInit(V, size_A);
    ceroInit(F, size_A);
//    gaussiana(X, size_A);
    gaussiana(Fext, npart, ncuerdas);
    randommasa(M, npart, ncuerdas, masatension, 0.1);
    randommasa(Fr, npart, ncuerdas, friccion, 0.1);
    friccionpuntas(Fr, npart, ncuerdas, friccion, 0.1);
   

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


   for (uint i = 0 ; i < npasos ; ++i) 
 { 
   
 
   // 5. perform the calculation
   // setup execution parameters
   dim3 threads(npart);
   dim3 grid(ncuerdas);
//     printf("pancho %e %e \n", masatension, friccion) ; 
   // execute the kernel
   if(10 < i && i < 12) {
   avance<<< grid, threads >>>(d_X, d_V, 
                                  d_F, d_Fext, d_salida, npart, d_M, d_Fr, 1);
 } else {
   avance<<< grid, threads >>>(d_X, d_V, 
                                 d_F, d_Fext, d_salida, npart, d_M, d_Fr, 0);

}
   


   // 11. copy result from device to host
   cudaMemcpy(salida, d_salida, buffer_salida, 
   cudaMemcpyDeviceToHost);
//   cudaMemcpy(X, d_X, mem_size_A, 
//   cudaMemcpyDeviceToHost);
/*    for (uint j = 0 ; j< npart; ++j) {
    printf ("Xs %e  %d \n", X[j], j ) ; 
             }
*/
    for (uint j = 0 ; j< 256; ++j) {
       int sal =0 ;
       int sal2 =0 ;
         for (uint i=0; i< ncuerdas/2 ;++i) {
         sal += 100000*salida[j+i*256];
         sal2 += 100000*salida[j+i*256*2];
                  }
    printf ("  %d , %d  \n", sal , sal2) ; 
            }
   }
 
}

