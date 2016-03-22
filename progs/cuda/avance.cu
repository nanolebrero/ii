#define BLOCK_SIZE 1024 
#define ESCRIBE 1
__global__  void avance( float* X, float* V, float* F, const float* Fext, float* salida,const int nodos, const float* mt, const float* fr , bool* tococ )
{
      
    // Block index
    uint num_cuer = blockIdx.x;
    uint nodo_en=threadIdx.x; 
//     printf("num_cuer %d, \n", num_cuer);
    // Thread index
      uint nodo = threadIdx.x+num_cuer*nodos;
      float vel=V[nodo] ;  
      float xr = X[nodo];
      float mtt=mt[nodo_en]/(num_cuer+1.);
      float fext=Fext[nodo];
      bool toco=tococ[num_cuer];
      float frr=1.-fr[nodo_en];

        fext=fext*(0.8 + 0.2/mtt);
//        frr=frr*mtt;
//       printf("nodo,mttyfrr %d % , 
     for (uint ii=0 ; ii < 256*ESCRIBE ; ++ii) { 
        if(nodo_en > 0 && nodo_en < nodos) {
           float fuerza = (X[nodo-1] + X[nodo+1] - 2.*xr);  
           if(toco /*&& ii < 200*/) {

           fuerza += fext ;
//           printf("toco, %d %e \n ", nodo_en ,fext); 
          }
             
           vel += fuerza*mtt ;

           vel = vel*frr;

           xr = xr + vel ;

           X[nodo] = xr;

//         if(nodo_en == 1*nodos/4) printf(" bla %e \n",xr) ;
}
         __syncthreads();
    if(nodo_en == 0 || nodo_en==nodos)xr=0 ;
 
          if (nodo_en == 1*nodos/7) {
           if( ii % ESCRIBE == 0) {
           salida[ii/ESCRIBE+num_cuer*256]=xr ;//sale/ESCRIBE;
          }
          }
   }
/*         if (nodo_en == 1*nodos/4) {
          for(uint kk=0; kk < 256 ; ++kk) {

            printf("sale %d \n", salida[kk]) ;
           }
           }*/

           V[nodo] = vel ;

             

}

