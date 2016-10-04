#define BLOCK_SIZE 1024 
#define ESCRIBE 2
__global__  void avance( float* X, float* V, float* F, const float* Fext, float* salida,const int nodos, const float* mt, const float* fr , bool* tococ, const float random )
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
      float fext=Fext[nodo];//*random;
      bool toco=tococ[num_cuer];
      float frr=1.-fr[nodo_en];
      float fuerza_vieja=0. ;
       bool suelta ;
        fext=fext*(0.8 + 0.2/mtt);
//        frr=frr*mtt;
//       printf("nodo,mttyfrr %d % , 
     for (uint ii=0 ; ii < 256*ESCRIBE ; ++ii) { 
        if(nodo_en > 0 && nodo_en < nodos) {
               float fuerza ; 
          fuerza = (X[nodo-1] + X[nodo+1] - 2.*xr);
        // if(nodo_en > 2 && nodo_en < nodos-2)  fuerza = (X[nodo-3] + X[nodo+3] - 2.*xr);
        // if(nodo_en <= 2) fuerza = (0. + X[nodo+3] - 2.*xr); 
        // if(nodo_en >= nodos-2) fuerza = (0. + X[nodo-3] - 2.*xr); 

        //         fuerza +=-300*(pow((xr-X[nodo-1]),3) + pow((xr-X[nodo+1]),3)) ;
//           float fuerza = (X[nodo-1] + X[nodo+1] - 2.00001*xr); // PARECE ALGO 3D
           if(toco && ii < 2000) {

            //if(vel >= 0.02*random) fuerza += fext ;
             fuerza+= fext*10;//(-0.1*ii*ii+51.2*ii) ;
             // frr = 0.8;
         }
           
           vel += fuerza*mtt ;
           vel = vel*frr;
           xr = xr + vel ;
           X[nodo] = xr;
}
         __syncthreads();
    if(nodo_en == 0 || nodo_en==nodos) { 
             xr=0 ;
             X[nodo]=xr ;
          }
          if (nodo_en == 1*nodos/7) {
           if( ii % ESCRIBE == 0) {
           salida[ii/ESCRIBE+num_cuer*256]=xr ;//sale/ESCRIBE;
          }
          }
   }
           V[nodo] = vel ;
}

