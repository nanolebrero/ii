
           implicit real*4 (a-h,o-z)

            real*4 , dimension (:), allocatable :: x,y
            integer*4 listo
            real*4, dimension (:), allocatable :: fx,fy,vx,vy
            logical toco 
            open(unit=56,file='geom.xyz')
            open(unit=23,file='salida')
            open(unit=55,file='ini')
           read(55,*) npas,npart,fact0,fact2,dxmax
           allocate(x(npart),y(npart),fx(npart),fy(npart)
     > ,vx(npart),vy(npart))
           write(*,*) npart
            kwrite=0            
            x=0
            vy=0
            vx=0

              !!!!!!!!!!!! Precalculo de gaussiana de exitacion
            do i=1,npart
            fy(i) = 1.D-5*exp(-0.001*
     >      ((i-npart/4)**2))
            enddo
              !!!!!!!!!!!!!!!
            do j=1,npas
            fx=0
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! calculo fuerzas
c            xcm=x(1)
c            fx(1) = x(2)-x(1)
            do i=2,npart-1
              fx(i)=x(i-1)+x(i+1)-2*x(i) !+ 0.0005*dx2
c              xcm=xcm + x(i)
            enddo
c            fx(npart)=x(npart-1)-x(npart)
c            xcm=xcm + x(npart)
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           !!!!!!!!!!!!exitación
           toco=.false.
           if((j.gt.100.and.j.lt.400000).
     >     or.(j.gt.800000.and.j.lt.1500000). 
     >     or.(j.gt.1600000.and.j.lt.2500000). 
     >     or.(j.gt.2900000.and.j.lt.3500000)) 
     >      then
           toco=.true. 
          if(mod(j,250).eq.0) then !!!El valor aleatorio lo actualizo cada tanto para evitar ruido muy agudo
            azar=rand()-0.5
           fff=real(j)/1000.
           fff=min(1.,fff)
           azar = fff*azar
           endif

          if(mod(j,250).lt.250) then   !!!! El valor aleatorio lo aplico 1/10 del tiempo
           do kk= 1,npart
           fx(kk)=fx(kk)+fy(kk)*azar    !!!! Aplicando
           enddo
           endif
           endif
           if (toco) then   
           fact=fact0
           else
            fact=fact0*100    !!!!! Si dejo de tocar el apagamiento aumenta
            endif
           do i=1,npart
            vx(i)=vx(i)+fx(i)*fact2  !!!! Actualizao velocidades por la fuerza
            vx(i) = vx(i)*(1.-fact)  !!!! Aplico fricción 
           enddo

            do i=1,npart
            x(i)=x(i)+vx(i)!-xcm !!!! Actualizo posiciones
            enddo

            x(1)=0.    !!!! Los extremos estan quietos
             x(npart)=0.

             listo=0   !!!!Output en integer*4
             alisto=0  !!!! output en Real*4
             if (kwrite.eq.1) then  !!! Cada cuanto escribe
             xcm=0
             do i=1,npart,10
              xcm=xcm + x(i)
             enddo
              xcm=10*xcm/npart

              do ii=npart/4-5,npart/4+5 !!! Acumulo sobre varios puentos para filtrar ultrasonido
              alisto=alisto + x(ii)!-10*xcm
               enddo
              listo=alisto*1000000 !!!! Multiplico por un factor que de numeros razonables para el kwave
               kwrite=0
               write(23,*) listo 
             else
              kwrite=kwrite+1
            endif
            
           enddo
           end

