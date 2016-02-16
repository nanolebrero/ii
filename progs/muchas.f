
           implicit real*4 (a-h,o-z)

            real*4 , dimension (:,:), allocatable :: x,y
            integer*4 listo
            real*4, dimension (:,:), allocatable :: fx,fy,vx,vy
            logical toco 
            real*4, dimension (:), allocatable::fact0,fact2,azar
            open(unit=56,file='geom.xyz')
            open(unit=23,file='salida')
            open(unit=55,file='ini')
           read(*,*) ncuerdas,npas,npart,nwrite
           allocate(fact0(ncuerdas),fact2(ncuerdas),azar(ncuerdas))

           do i=1,ncuerdas
             read(*,*) fact0(i),fact2(i)
           enddo

         allocate(x(npart,ncuerdas),y(npart,ncuerdas),fx(npart,ncuerdas)
     >    ,fy(npart,ncuerdas),vx(npart,ncuerdas) , vy(npart,ncuerdas))
           write(*,*) npart
            kwrite=0            
            x=0
            vy=0
            vx=0

              !!!!!!!!!!!! Precalculo de gaussiana de exitacion
            do i=1,npart
            fy(i,:) = 1.D-5*exp(-0.001*((i-npart/4)**2))
            enddo
              !!!!!!!!!!!!!!!
            do j=1,npas
            fx=0
            do kcu=1,ncuerdas
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculo fuerzas
            do i=2,npart-1
              fx(i,kcu)=x(i-1,kcu)+x(i+1,kcu)-2*x(i,kcu) !+ 0.0005*dx2
            enddo
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !!!!!!!!!!!!exitación
           toco=.false.
           if((j.gt.100.and.j.lt.400000).or.
     >    (j.gt.800000.and.j.lt.1500000).or.
     >  ( j.gt.1600000.and.j.lt.2500000).or.
     >   (j.gt.2900000.and.j.lt.3500000))
     > then
           toco=.true. 
          if(mod(j,250).eq.0) then !!!El valor aleatorio lo actualizo cada tanto para evitar ruido muy agudo
            azar(kcu)=rand()-0.5
c           fff=real(j)/1000.
c           fff=min(1.,fff)
c           azar(kcu) = azar(kcu)
           endif

          if(mod(j,250).lt.250) then   !!!! El valor aleatorio lo aplico 1/10 del tiempo
           do kk= 1,npart
           fx(kk,kcu)=fx(kk,kcu)+fy(kk,kcu)*azar(kcu)    !!!! Aplicando
           enddo
           endif
           endif
           if (toco) then   
           fact=fact0(kcu)
           else
            fact=fact0(kcu)*100    !!!!! Si dejo de tocar el apagamiento aumenta
            endif
           do i=1,npart
            vx(i,kcu)=vx(i,kcu)+fx(i,kcu)*fact2(kcu)  !!!! Actualizao velocidades por la fuerza
            vx(i,kcu) = vx(i,kcu)*(1.-fact)  !!!! Aplico fricción 
           enddo

            do i=1,npart
c            write(44,*) j, i, vx(i) 
            x(i,kcu)=x(i,kcu)+vx(i,kcu) !!!! Actualizo posiciones
            enddo
          
             x(1,kcu)=0.    !!!! Los extremos estan quietos
             x(npart,kcu)=0.

             enddo   ! Fin del loop x cuerda

             listo=0   !!!!Output en integer*4
             alisto=0  !!!! output en Real*4

             if (kwrite.eq.nwrite) then  !!! Cada cuanto escribe
              do kcu=1,ncuerdas
              do ii=npart/4-5,npart/4+5 !!! Acumulo sobre varios puentos para filtrar ultrasonido
              alisto=alisto + x(ii,kcu)
               enddo
               enddo

              listo=alisto*100000 !!!! Multiplico por un factor que de numeros razonables para el kwave
               kwrite=0
               write(23,*) listo 
             else
              kwrite=kwrite+1
            endif
           
           enddo
           end

