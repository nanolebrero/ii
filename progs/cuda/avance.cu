#define BLOCK_SIZE 1024 
#define ESCRIBE 1
__global__  void avance( float* X, float* V, float* F, const float* Fext, float* salida,const int nodos, const float* mt, const float* fr , bool toco )
{
      
    // Block index
    int num_cuer = blockIdx.x;
//     printf("Hola HOla %d \n", num_cuer) ;
      int nodo_en=threadIdx.x; 
    // Thread index
      int nodo = threadIdx.x+num_cuer*nodos;
      float vel=V[nodo] ;  
      float xr = X[nodo];
      float mtt=mt[nodo_en]/(0.2*num_cuer+1);
      float fext=Fext[nodo_en];
//       if(num_cuer==17)  xr=0;
//       if(num_cuer==0)     mtt=mt[nodo_en];
      float frr=1.-fr[nodo_en];
//     printf("frr %e \n", frr) ;

    for (int ii=0 ; ii < 256*ESCRIBE ; ++ii) { 
     if(nodo_en > 0 && nodo_en < nodos) {
     float fuerza = (X[nodo-1] + X[nodo+1] - 2.*xr);  
     if(toco && ii < 200) fuerza += fext ;

        vel += fuerza*mtt ;

        vel= vel*frr;

         xr = xr + vel ;

          X[nodo]=xr;
}
/*    if(num_cuer < 5) {
    if(nodo_en==nodos/2) {
     xr=0;

     }
     }*/
         __syncthreads();
    if(nodo_en == 0 || nodo_en==nodos)xr=0 ;
                             
 
          if (nodo_en == 7*nodos/12) {
           if( ii % ESCRIBE == 0) {
           salida[ii/ESCRIBE+num_cuer*256]=xr ;//sale/ESCRIBE;
//           printf("Hola HOla %e %d %d \n", salida[ii/ESCRIBE+num_cuer*256], num_cuer, nodos) ;
          }
          }
   }

           V[nodo] = vel ;

             

}

